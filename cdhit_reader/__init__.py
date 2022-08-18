from ._cli import cli
from ._compare import compare
from ._reader import ParsingError, ClusterSequence, Cluster, ClstrReader, read_cdhit, SeqType, Strand
from ._testit import test
from ._version import __version__
#from ._writer import FASTAWriter, write_fasta

__all__ = [
    "ParsingError",
    "ClusterSequence",
    "Cluster",
    "ClstrReader",
    "read_cdhit",
    "SeqType",
    "Strand", # from reader
    "__version__",
    "cli",
    "compare",
    "test",
]
