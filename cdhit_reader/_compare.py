from email.policy import default
import sys
from statistics import mean

import click

from ._reader import read_cdhit
from ._version import __version__
from xopen import xopen
import os
import tempfile
import subprocess

def has_cdhit():
    cmd = ["cd-hit", "-h"]
    # Run cmd and save stdout as "cdout" variable, ignore exit status
    try:
        cdout = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.decode("utf-8")
        if "CD-HIT" in cdout:
            return True
    except Exception:
        return False
  
def read_fasta(path):
    name = None
    comment = ""
    seq = None
    with xopen(path) as fasta:
        for line in fasta:
            line = line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    yield name, comment, seq
                name = line[1:].split(" ")[0].split("\t")[0]
                comment = " ".join(line[1 + len(name):].split(" "))
                seq = ""
            else:
                seq += line
    yield name, comment, seq

def relabel_fasta(path, prefix, outpath):
    with xopen(outpath, "a") as out:
        for name, comment, seq in read_fasta(path):
            out.write(">{}{}{}{}\n{}\n".format(prefix, name, " " if len(comment) > 1 else "", comment, seq))

@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.argument("fasta1", type=click.Path(exists=True))
@click.argument("fasta2", type=click.Path(exists=True))
@click.option("--tag1",   help="Name of the first dataset")
@click.option("--tag2",  help="Name of the second dataset")
@click.option("--id", help="Identity threshold [default: 0.9]", default=0.95, type=float)
@click.option("--type", type=click.STRING, help="Type of the sequences (nucl or prot)")
@click.option("--tempdir",type=click.Path(exists=True), help="Temporary directory for intermediate files", default=tempfile.gettempdir())
@click.option("--verbose", default=False, is_flag=True, help="Show verbose information")
def compare(fasta1, fasta2, tag1: str, tag2: str, tempdir, type: str, id: float,verbose: bool):
    """
    Compare FASTA files

    \b
    Warning
    -------
    The commad line interface is in EXPERIMENTAL stage.
    """
    if not has_cdhit():
        click.echo("cd-hit is not installed. Please install it and try again.")
        sys.exit(1)
    # Generate temporary directory inside "tempdir"
    SEPARATOR = "_"
    tmp = tempfile.TemporaryDirectory(dir=tempdir)
    if verbose:
        print("Temporary directory: {}".format(tmp.name), file=sys.stderr)
    prefix1 = tag1 if tag1 else os.path.basename(fasta1).split("_")[0].split(".")[0] + SEPARATOR
    prefix2 = tag2 if tag2 else os.path.basename(fasta2).split("_")[0].split(".")[0] + SEPARATOR
    fasta_file = os.path.join(tmp.name, "seqs.fasta")
    clstr_file = os.path.join(tmp.name, "clusters.fasta")

    if type is None:
        type = "prot" if "faa" in fasta1 else "nucl"
        print("Type of sequences not specified, assuming {}".format(type), file=sys.stderr)

    if prefix1 == prefix2:
        print("Warning: prefixes are identical ({}): specify manual prefixes with --tag1 and --tag2".format(prefix1))
        exit(1)
    
    # Delete fasta_file if present:
    if os.path.exists(fasta_file):
        os.remove(fasta_file)
    relabel_fasta(fasta1, prefix1, fasta_file)
    relabel_fasta(fasta2, prefix2, fasta_file)
    if verbose:
        print("Relabeling {} to {} (prefix: {})".format(fasta1, fasta_file, prefix1), file=sys.stderr)
        print("Relabeling {} to {} (prefix: {})".format(fasta2, fasta_file, prefix2), file=sys.stderr)
        
    cmd = ["cd-hit" if type=="prot" else "cd-hit-est", "-i", fasta_file, "-o", clstr_file, "-c", str(id), "-d", "1000"]
    if verbose:
        print("Running {}".format(" ".join(cmd)), file=sys.stderr)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    stats = {prefix1: [], prefix2: [], "both": [], "multiple": [], "dupl": []}
    for cluster in read_cdhit(clstr_file + ".clstr"):

        if len(cluster) == 1:
            # Singleton
            pool = (cluster.refname).split(SEPARATOR)[0]
            seqname = cluster.refname[len(pool) + 1:]
            stats[pool+SEPARATOR].append(seqname)
        elif len(cluster) == 2:
            pool = (cluster.refname).split(SEPARATOR)[0]
            seqname = cluster.refname[len(pool) + 1:]
            # Pairwise comparison
            pair = []
            check = False
            for seq in cluster.sequences:
                sub_pool = (seq.name).split(SEPARATOR)[0]
                sub_seqname = (seq.name)[len(pool) + 1:]
                pair.append(sub_seqname)
                if sub_pool != pool:
                    check = True
            if check == True:    
                stats["both"].append(":".join(pair))
            else:
                stats["dupl"].append(":".join(pair))
        else:
            # Multiple sequences
            pass
    
    for key, list in stats.items():
        for i in list:
            print(key, i, sep="\t")


def show_hist(seq_lens):
    pass
