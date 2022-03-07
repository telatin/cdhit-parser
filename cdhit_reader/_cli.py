from email.policy import default
import sys
from statistics import mean

import click

from ._reader import read_cdhit
from ._version import __version__
#from ._writer import write_fasta


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.argument("clstr", type=click.Path(exists=True))
@click.option("--stats/--no-stats", default=True, help="Show sequence statistics.")
@click.option("--hist/--no-hist", default=False, help="Show histogram of sequence lengths.")
@click.option("--all", default=False, is_flag=True, help="Show all sequences in the cluster.")
def cli(clstr, stats: bool, hist: bool, all: bool):
    """
    Show information about CLSTR file

    \b
    Warning
    -------
    The commad line interface is in EXPERIMENTAL stage.  
    """

    nitems = 0
    nseqs = 0
    seq_lens = []
    for item in read_cdhit(clstr):
        seq_lens.append(len(item))
        nitems += 1
        nseqs  += len(item)
        if all:
            print(item)
            for s in item.sequences:
                print(f"   {s}")
    if stats:
        click.echo(f"Input file: {clstr}")
        click.echo(f"Number of clusters: {nitems}")
        click.echo(f"Total sequences: {nseqs}")
        
        msg = f"Cluster size: min {min(seq_lens)}, mean {mean(seq_lens):.2f}, max {max(seq_lens)}"
        click.echo(msg)

    if hist:
        show_hist(seq_lens)


def show_hist(seq_lens):
    pass