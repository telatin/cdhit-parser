from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List
import csv
import os
import shutil
import subprocess
import tempfile

import click

from ._fasta import Sequence, read_fasta
from ._reader import Cluster, read_cdhit
from ._version import __version__

COMPRESSED_SUFFIXES = (".gz", ".bz2", ".xz", ".zst")


@dataclass(frozen=True)
class RenamedSequence:
    source_label: str
    renamed_name: str
    original: Sequence


@dataclass
class CategorySummary:
    label: str
    description: str
    output_path: Path
    clusters: int = 0
    output_records: int = 0
    source1_members: int = 0
    source2_members: int = 0


def has_cdhit_est() -> bool:
    return shutil.which("cd-hit-est") is not None


def input_prefix(path: Path) -> str:
    """
    Return a safe prefix derived from the input filename.
    """
    name = path.name
    for suffix in COMPRESSED_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    stem = Path(name).stem
    return stem.replace("___", "__")


def _record_name(prefix: str, seq_name: str) -> str:
    return f"{prefix}___{seq_name}"


def write_relabeled_fasta(input_path: Path, prefix: str, output_path: Path, append: bool) -> Dict[str, RenamedSequence]:
    """
    Relabel a FASTA file and append it to the concatenated input used for CD-HIT-EST.
    """
    mode = "a" if append else "w"
    renamed: Dict[str, RenamedSequence] = {}

    with output_path.open(mode, encoding="utf-8") as handle:
        for seq in read_fasta(input_path):
            renamed_name = _record_name(prefix, seq.name)
            if renamed_name in renamed:
                raise click.ClickException(
                    f"Duplicate sequence identifier after relabeling: {renamed_name}"
                )

            renamed[renamed_name] = RenamedSequence(prefix, renamed_name, seq)
            handle.write(f">{renamed_name}\n{seq.sequence}\n")

    return renamed


def classify_cluster(cluster: Cluster, prefix1: str, prefix2: str) -> tuple[str, int, int]:
    prefix1_marker = f"{prefix1}___"
    prefix2_marker = f"{prefix2}___"
    count1 = 0
    count2 = 0

    for member in cluster.sequences:
        if member.name.startswith(prefix1_marker):
            count1 += 1
        elif member.name.startswith(prefix2_marker):
            count2 += 1
        else:
            raise click.ClickException(
                f"Cluster member {member.name!r} does not belong to either input dataset."
            )

    if count1 == 1 and count2 == 1 and len(cluster) == 2:
        return "one_to_one_matches", count1, count2
    if count1 > 0 and count2 > 0:
        return "shared_multi", count1, count2
    if count1 > 0:
        return "only_file1", count1, count2
    return "only_file2", count1, count2


def sequences_for_cluster(
    cluster: Cluster,
    renamed_sequences: Dict[str, RenamedSequence],
    all_members: bool,
) -> List[Sequence]:
    if all_members:
        sequence_names = [member.name for member in cluster.sequences]
    else:
        if cluster.refname is None:
            raise click.ClickException(f"Cluster {cluster.name!r} does not define a representative sequence.")
        sequence_names = [cluster.refname]

    sequences: List[Sequence] = []
    for renamed_name in sequence_names:
        try:
            sequences.append(renamed_sequences[renamed_name].original)
        except KeyError as exc:
            raise click.ClickException(
                f"Could not restore original sequence for renamed identifier {renamed_name!r}."
            ) from exc
    return sequences


def write_fasta_output(path: Path, sequences: Iterable[Sequence]) -> int:
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        for sequence in sequences:
            handle.write(f"{sequence}\n")
            count += 1
    return count


def write_stats(path: Path, summaries: Iterable[CategorySummary]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "category",
                "description",
                "clusters",
                "output_records",
                "source1_members",
                "source2_members",
                "output_fasta",
            ]
        )
        for summary in summaries:
            writer.writerow(
                [
                    summary.label,
                    summary.description,
                    summary.clusters,
                    summary.output_records,
                    summary.source1_members,
                    summary.source2_members,
                    summary.output_path.name,
                ]
            )


def write_log(path: Path, lines: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.option("-1", "--fasta-one", "fasta_one", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path), help="First FASTA file.")
@click.option("-2", "--fasta-two", "fasta_two", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path), help="Second FASTA file.")
@click.option("-o", "--outdir", required=True, type=click.Path(file_okay=False, path_type=Path), help="Output directory.")
@click.option("-c", "--perc-identity", default=0.99, show_default=True, type=click.FloatRange(min=0.0, max=1.0, min_open=True), help="CD-HIT-EST sequence identity threshold.")
@click.option("--coverage", default=0.90, show_default=True, type=click.FloatRange(min=0.0, max=1.0, min_open=True), help="Coverage threshold applied to both -aS and -aL.")
@click.option("-d", "--description-length", default=0, show_default=True, type=int, help="Value passed to cd-hit-est -d.")
@click.option("--temp-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), help="Override the temporary directory base.")
@click.option("--keep-temp-dir", is_flag=True, default=False, help="Keep the temporary working directory.")
@click.option("--all-members", is_flag=True, default=False, help="Write all members of matching clusters instead of only the representative sequence.")
def compare_otus(
    fasta_one: Path,
    fasta_two: Path,
    outdir: Path,
    perc_identity: float,
    coverage: float,
    description_length: int,
    temp_dir: Path | None,
    keep_temp_dir: bool,
    all_members: bool,
) -> None:
    """
    Compare two OTU FASTA files with cd-hit-est.
    """
    if not has_cdhit_est():
        raise click.ClickException(
            "cd-hit-est is not available in PATH. Please install CD-HIT and try again."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    prefix1 = input_prefix(fasta_one)
    prefix2 = input_prefix(fasta_two)
    if prefix1 == prefix2:
        raise click.ClickException(
            "Input basenames resolve to the same prefix. Rename the inputs or adjust their filenames."
        )

    temp_base = temp_dir or Path(os.environ.get("TMP", "/tmp"))
    workdir = Path(tempfile.mkdtemp(prefix="cdhit_compare_otus_", dir=temp_base))
    concatenated_fasta = workdir / "concatenated.fasta"
    output_prefix = workdir / "cdhit_compare_otus"
    log_lines = [
        "cdhit-compare-otus run",
        f"fasta_one\t{fasta_one}",
        f"fasta_two\t{fasta_two}",
        f"outdir\t{outdir}",
        f"temp_dir\t{workdir}",
        f"prefix_one\t{prefix1}",
        f"prefix_two\t{prefix2}",
        f"perc_identity\t{perc_identity}",
        f"coverage\t{coverage}",
        f"description_length\t{description_length}",
        f"all_members\t{all_members}",
    ]

    stats_path = outdir / "compare_otus_stats.tsv"
    log_path = outdir / "compare_otus.log"

    summaries = {
        "one_to_one_matches": CategorySummary(
            label="one_to_one_matches",
            description="Exactly one member from each input file.",
            output_path=outdir / "one_to_one_matches.fasta",
        ),
        "shared_multi": CategorySummary(
            label="shared_multi",
            description="At least one member from both inputs and multiple members in at least one input.",
            output_path=outdir / "shared_multi.fasta",
        ),
        "only_file1": CategorySummary(
            label="only_file1",
            description=f"Clusters containing only sequences from {prefix1}.",
            output_path=outdir / f"only_{prefix1}.fasta",
        ),
        "only_file2": CategorySummary(
            label="only_file2",
            description=f"Clusters containing only sequences from {prefix2}.",
            output_path=outdir / f"only_{prefix2}.fasta",
        ),
    }

    output_sequences = {key: [] for key in summaries}

    try:
        renamed_sequences: Dict[str, RenamedSequence] = {}
        renamed_sequences.update(
            write_relabeled_fasta(fasta_one, prefix1, concatenated_fasta, append=False)
        )
        second_sequences = write_relabeled_fasta(
            fasta_two, prefix2, concatenated_fasta, append=True
        )
        overlap = set(renamed_sequences).intersection(second_sequences)
        if overlap:
            duplicated = ", ".join(sorted(overlap))
            raise click.ClickException(f"Duplicate renamed identifiers across inputs: {duplicated}")
        renamed_sequences.update(second_sequences)

        command = [
            "cd-hit-est",
            "-i",
            str(concatenated_fasta),
            "-o",
            str(output_prefix),
            "-c",
            str(perc_identity),
            "-aS",
            str(coverage),
            "-aL",
            str(coverage),
            "-d",
            str(description_length),
        ]
        log_lines.append(f"command\t{' '.join(command)}")

        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout.strip():
            log_lines.append("cd-hit-est stdout:")
            log_lines.extend(result.stdout.strip().splitlines())
        if result.stderr.strip():
            log_lines.append("cd-hit-est stderr:")
            log_lines.extend(result.stderr.strip().splitlines())

        for cluster in read_cdhit(f"{output_prefix}.clstr"):
            category, count1, count2 = classify_cluster(cluster, prefix1, prefix2)
            sequences = sequences_for_cluster(cluster, renamed_sequences, all_members=all_members)
            output_sequences[category].extend(sequences)
            summaries[category].clusters += 1
            summaries[category].output_records += len(sequences)
            summaries[category].source1_members += count1
            summaries[category].source2_members += count2

        for summary in summaries.values():
            written = write_fasta_output(summary.output_path, output_sequences[summary.label])
            if written != summary.output_records:
                raise click.ClickException(
                    f"Expected to write {summary.output_records} sequences to {summary.output_path.name}, wrote {written}."
                )
            log_lines.append(
                f"output\t{summary.label}\t{summary.output_path}\tclusters={summary.clusters}\tsequences={summary.output_records}"
            )

        write_stats(stats_path, summaries.values())
        log_lines.append(f"stats\t{stats_path}")
        log_lines.append(f"log\t{log_path}")
        click.echo(f"Results written to {outdir}")
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            log_lines.append("cd-hit-est stdout:")
            log_lines.extend(exc.stdout.strip().splitlines())
        if exc.stderr:
            log_lines.append("cd-hit-est stderr:")
            log_lines.extend(exc.stderr.strip().splitlines())
        raise click.ClickException("cd-hit-est failed; see compare_otus.log for details.") from exc
    finally:
        write_log(log_path, log_lines)
        if not keep_temp_dir:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            click.echo(f"Temporary directory kept at {workdir}")
