from pathlib import Path

from click.testing import CliRunner

from cdhit_reader._compare_otus import classify_cluster, compare_otus, input_prefix


class DummySequence:
    def __init__(self, name):
        self.name = name


class DummyCluster:
    def __init__(self, names, refname=None, name="Cluster 0"):
        self.sequences = [DummySequence(seq_name) for seq_name in names]
        self.refname = refname or names[0]
        self.name = name

    def __len__(self):
        return len(self.sequences)


def test_input_prefix_sanitizes_triple_underscore():
    assert input_prefix(Path("otu___table.fa.gz")) == "otu__table"


def test_classify_cluster_variants():
    assert classify_cluster(
        DummyCluster(["a___seq1", "b___seq2"]),
        "a",
        "b",
    ) == ("one_to_one_matches", 1, 1)
    assert classify_cluster(
        DummyCluster(["a___seq1", "a___seq2", "b___seq3"]),
        "a",
        "b",
    ) == ("shared_multi", 2, 1)
    assert classify_cluster(
        DummyCluster(["a___seq1"]),
        "a",
        "b",
    ) == ("only_file1", 1, 0)
    assert classify_cluster(
        DummyCluster(["b___seq1"]),
        "a",
        "b",
    ) == ("only_file2", 0, 1)


def test_compare_otus_cli_rejects_missing_cdhit_est(monkeypatch, tmp_path):
    runner = CliRunner()
    fasta_one = tmp_path / "one.fa"
    fasta_two = tmp_path / "two.fa"
    fasta_one.write_text(">a\nAAAA\n", encoding="utf-8")
    fasta_two.write_text(">b\nAAAA\n", encoding="utf-8")
    monkeypatch.setattr("cdhit_reader._compare_otus.has_cdhit_est", lambda: False)

    result = runner.invoke(
        compare_otus,
        [
            "-1",
            str(fasta_one),
            "-2",
            str(fasta_two),
            "-o",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert "cd-hit-est is not available" in result.output


def test_compare_otus_writes_outputs_with_restored_headers(monkeypatch, tmp_path):
    runner = CliRunner()
    fasta_one = tmp_path / "otu___table.fa"
    fasta_two = tmp_path / "db.fa"
    outdir = tmp_path / "out"
    tempdir = tmp_path / "tmp"
    tempdir.mkdir()
    fasta_one.write_text(">seqA comment one\nAAAA\n>seqOnly\nCCCC\n", encoding="utf-8")
    fasta_two.write_text(">seqB comment two\nAAAA\n>seqOther\nGGGG\n", encoding="utf-8")

    def fake_run(command, check, capture_output, text):
        concatenated = Path(command[command.index("-i") + 1]).read_text(encoding="utf-8")
        assert concatenated.count(">") == 4
        output_prefix = Path(command[command.index("-o") + 1])
        output_prefix.with_suffix(".clstr").write_text(
            ">Cluster 0\n"
            "0 4nt, >otu__table___seqA... *\n"
            "1 4nt, >db___seqB... at +/100.00%\n"
            ">Cluster 1\n"
            "0 4nt, >otu__table___seqOnly... *\n"
            ">Cluster 2\n"
            "0 4nt, >db___seqOther... *\n",
            encoding="utf-8",
        )
        return type("Completed", (), {"stdout": "", "stderr": ""})()

    monkeypatch.setattr("cdhit_reader._compare_otus.has_cdhit_est", lambda: True)
    monkeypatch.setattr("cdhit_reader._compare_otus.subprocess.run", fake_run)

    result = runner.invoke(
        compare_otus,
        [
            "-1",
            str(fasta_one),
            "-2",
            str(fasta_two),
            "-o",
            str(outdir),
            "--temp-dir",
            str(tempdir),
        ],
    )

    assert result.exit_code == 0
    assert "Results written to" in result.output
    assert ">seqA comment one\nAAAA\n" == (outdir / "one_to_one_matches.fasta").read_text(encoding="utf-8")
    assert (outdir / "shared_multi.fasta").read_text(encoding="utf-8") == ""
    assert ">seqOnly\nCCCC\n" == (outdir / "only_otu__table.fasta").read_text(encoding="utf-8")
    assert ">seqOther\nGGGG\n" == (outdir / "only_db.fasta").read_text(encoding="utf-8")

    stats = (outdir / "compare_otus_stats.tsv").read_text(encoding="utf-8")
    assert "one_to_one_matches" in stats
    assert "only_file1" in stats
    assert "only_file2" in stats
    log = (outdir / "compare_otus.log").read_text(encoding="utf-8")
    assert "prefix_one\totu__table" in log
    assert "prefix_two\tdb" in log
