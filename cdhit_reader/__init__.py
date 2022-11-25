from ._cli import cli
from ._compare import compare
from ._reader import ParsingError, ClusterSequence, Cluster, ClstrReader, read_cdhit, SeqType, Strand
from ._fasta import Sequence, FastaReader, read_fasta
from ._testit import test
from ._version import __version__
#from ._writer import FASTAWriter, write_fasta

__all__ = [
    "ParsingError",
    "ClusterSequence",
    "Cluster",
    "ClstrReader",
    "read_cdhit",
    "FastaReader",
    "read_fasta",
    "SeqType",
    "Strand", 
    "Sequence",
    "FastaReader",
    "read_fasta",
    "__version__",
    "cli",
    "compare",
    "test",
]
