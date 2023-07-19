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

def split_cluster(cluster, tag1, tag2):
    pool1, pool2 = [], []
    for sequence in cluster.sequences:
        if sequence.name.startswith(tag1):
            pool1.append(sequence)
        elif sequence.name.startswith(tag2):
            pool2.append(sequence)
        else:
            raise ValueError("Sequence {} does not start with {} or {}".format(sequence.name, tag1, tag2))
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
    TAG1      = "1:::"
    TAG2      = "2:::"
    # tmp is a temporary directory not to be deleted
    tmp = tempfile.mkdtemp(dir=tempdir, prefix="cdhit_")

    #tmp = tempfile.TemporaryDirectory(dir=tempdir, prefix="cdhit_reader_", suffix=".tmp")
    if verbose:
        print("Temporary directory: {}".format(tmp), file=sys.stderr)

    prefix1 = tag1 if tag1 else os.path.basename(fasta1).split("_")[0].split(".")[0]
    prefix2 = tag2 if tag2 else os.path.basename(fasta2).split("_")[0].split(".")[0]

    fasta_file = os.path.join(tmp, "seqs.fasta")
    clstr_file = os.path.join(tmp, "clusters.fasta")

    if type is None:
        type = "prot" if "faa" in fasta1 else "nucl"
        print("Type of sequences not specified, assuming {}".format(type), file=sys.stderr)

    if prefix1 == prefix2:
        print("Warning: prefixes are identical ({}): specify manual prefixes with --tag1 and --tag2".format(prefix1))
        exit(1)
    
    # Delete fasta_file if present:
    if os.path.exists(fasta_file):
        os.remove(fasta_file)
    
    if verbose:
        print("Relabeling {} to {}".format(fasta1, fasta_file), file=sys.stderr)
    
    relabel_fasta(fasta1, TAG1, fasta_file)

    if verbose:
        print("Relabeling {} to {}".format(fasta2, fasta_file), file=sys.stderr)
    
    relabel_fasta(fasta2, TAG2, fasta_file)

    tags = {
        TAG1: prefix1,
        TAG2: prefix2
    } 
    if verbose:
        print("Relabeling {} to {} (prefix: {})".format(fasta1, fasta_file, prefix1), file=sys.stderr)
        print("Relabeling {} to {} (prefix: {})".format(fasta2, fasta_file, prefix2), file=sys.stderr)
        
    cmd = ["cd-hit" if type=="prot" else "cd-hit-est", "-i", fasta_file, "-o", clstr_file, "-c", str(id), "-d", "1000"]
    if verbose:
        print("Running {}".format(" ".join(cmd)), file=sys.stderr)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    stats = {TAG1: [], TAG2: [], "both": [], "multi": [], "dupl_" + TAG1: [],  "dupl_" + TAG2: []}
    for cluster in read_cdhit(clstr_file + ".clstr"):
        
        pool = (cluster.refname)[0:len(TAG1)]
        if pool != TAG1 and pool != TAG2:
            raise ValueError("Cluster {} does not start with {} or {}".format(cluster.refname, TAG1, TAG2))
        seqname = cluster.refname[len(pool) + 1:]
        if len(cluster) == 1:
            # Singleton
            stats[pool].append(seqname)
        elif len(cluster) == 2:
            # Pairwise comparison
            pair = []
            check = False
            for seq in cluster.sequences:
                sub_pool = (seq.name)[0:len(TAG1)]
                sub_seqname = (seq.name)[len(pool) + 1:]
                pair.append(sub_seqname)
                if sub_pool != pool:
                    check = True
            if check == True:    
                stats["both"].append(":".join(pair))
            else:
                stats["dupl_" + pool ].append(":".join(pair))
        else:
            
            stats["multi"].append(",".join(i.name for i in cluster.sequences))
    

    for key, list in stats.items():
        print("{} {}".format(key, len(list)), file=sys.stderr)
        for seqnames in list:
            key = tags[key] if key in tags else key
            key = "dupl_" + tags[key[-1*len(TAG1):]] if key[-1*len(TAG1):] in tags else key
            seqnames = seqnames.replace(TAG1, tags[TAG1] + "#").replace(TAG2, tags[TAG2] + "#")
            
            print(key, seqnames, sep="\t")

    #if not verbose:
    #    os.remove(tmp)

def show_hist(seq_lens):
    pass
