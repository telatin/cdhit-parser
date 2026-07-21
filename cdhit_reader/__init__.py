from importlib import import_module

from ._reader import ParsingError, ClusterSequence, Cluster, Clustering, ClstrReader, read_cdhit, SeqType, Strand
from ._fasta import FastaParsingError, Sequence, FastaReader, read_fasta
from ._version import __version__

_LAZY_EXPORTS = {
    "cli": ("._cli", "cli"),
    "compare": ("._compare", "compare"),
    "test": ("._testit", "test"),
}

__all__ = [
    "ParsingError",
    "FastaParsingError",
    "ClusterSequence",
    "Cluster",
    "Clustering",
    "ClstrReader",
    "read_cdhit",
    "FastaReader",
    "read_fasta",
    "SeqType",
    "Strand",
    "Sequence",
    "__version__",
    "cli",
    "compare",
    "test",
]


def __getattr__(name):
    if name in _LAZY_EXPORTS:
        module_name, attribute_name = _LAZY_EXPORTS[name]
        module = import_module(module_name, __name__)
        return getattr(module, attribute_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
