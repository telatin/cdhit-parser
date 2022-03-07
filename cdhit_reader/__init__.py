from ._cli import cli
from ._reader import ParsingError, Cluster, ClstrReader, read_cdhit
from ._testit import test
from ._version import __version__
#from ._writer import FASTAWriter, write_fasta

__all__ = [
    "ParsingError", "Cluster", "ClstrReader", "read_cdhit", # from reader
    "__version__",
    "cli",
    "test",
]
