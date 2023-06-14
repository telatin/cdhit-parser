from __future__ import annotations
from pathlib import Path
from typing import IO, Iterator, List, Union
from enum import Enum
from more_itertools import peekable
from xopen import xopen
import re

__all__ = [
    "ParsingError",
    "ClusterSequence",
    "Cluster",
    "ClstrReader",
    "read_cdhit",
    "SeqType",
    "Strand",
    "FastaReader",
    "read_fasta",
]


class SeqType(Enum):
    """
    Sequence type.
    """

    PROTEIN = "Protein"
    NT = "DNA/RNA"
    NONE = "None"


class Strand(Enum):
    """
    Sequence type.
    """

    PLUS = "+"
    REVERSE = "-"
    NONE = "."


class ParsingError(Exception):
    """
    Parsing error.
    """

    def __init__(self, line_number: int):
        """
        Parameters
        ----------
        line_number
            Line number.
        """
        super().__init__(f"Line number {line_number}.")
        self._line_number = line_number

    @property
    def line_number(self) -> int:
        """
        Line number.

        Returns
        -------
        Line number.
        """
        return self._line_number


class ClusterSequence:
    """
    A single sequence of a cluster from line

    Attributes
    ----------
    line: str
    """

    def __init__(self, line: str):
        """
        Parameters
        ----------
        line
            Line.
        """
        self.line = line
        self.length = 0
        self.name = ""
        self.identity = 0
        self.is_ref = False
        self.seqtype = SeqType.NONE
        self.strand = Strand.NONE
        self.id = -1
        self.__parse()

    def __parse(self):
        """
        3       502nt, >IKXM6KN01CFAFI... at 1:502:1:503/+/97.81%
        """

        pattern = re.compile(
            r"(?P<id>\d+)\s+(?P<size>\d+)(?P<type>aa|nt), >(?P<name>.+?)\.\.\. (?P<attr>.+)"
        )
        attrpatt = re.compile(
              
            r"(?P<ref>\*|at) .*?(?P<strand>[+-]?)\/?(?P<percent>\d+\.?\d+)%"
        )

        match = pattern.search(self.line)

        self.seqtype = SeqType.PROTEIN if match["type"] == "aa" else SeqType.NT
        self.strand = Strand.NONE if self.seqtype == SeqType.PROTEIN else Strand.PLUS
        if match:
            self.name = match["name"]
            self.id = int(match["id"])
            self.length = int(match["size"])
            if match["attr"] == "*":
                self.is_ref = True
                self.identity = 100.0

            else:
                attrs = attrpatt.match(match["attr"])
                self.is_ref = False
                self.identity = float(attrs["percent"])
                self.strand = (
                    Strand.PLUS
                    if attrs["strand"] == "+"
                    else Strand.REVERSE
                    if attrs["strand"] == "-"
                    else Strand.NONE
                )

    def __repr__(self):
        """
        Returns a string representation of the ClusterSequence object.
        Args:
          self (ClusterSequence): The ClusterSequence object.
        Returns:
          str: A string representation of the ClusterSequence object.
        Examples:
          >>> cluster_sequence = ClusterSequence("seq1", "seq1", 0, 0.0, False, None, None)
          >>> cluster_sequence.__repr__()
          'ClusterSequence(id=seq1, name=seq1, length=0, identity=0.0, is_ref=False, seqtype=None, strand=None)'
        """
        return f"ClusterSequence(id={self.id}, name={self.name}, length={self.length}, identity={self.identity}, is_ref={self.is_ref}, seqtype={self.seqtype}, strand={self.strand})"


class Cluster:
    """
    Represents a cluster of sequences.
    Args:
      defline (str): The defline of the cluster.
      sequences (List[ClusterSequence]): The list of sequences in the cluster.
    Attributes:
      name (str): The defline of the cluster.
      sequences (List[ClusterSequence]): The list of sequences in the cluster.
      refname (str): The name of the reference sequence in the cluster.
    Examples:
      >>> cluster = Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])
      >>> cluster.name
      'Cluster_1'
      >>> cluster.sequences
      [ClusterSequence(id='seq1', name='seq1', length=0, identity=0.0, is_ref=False, seqtype=None, strand=None), ClusterSequence(id='seq2', name='seq2', length=0, identity=0.0, is_ref=False, seqtype=None, strand=None)]
      >>> cluster.refname
      None
    """

    def __init__(self, defline, sequences):
        """
        Initializes a Cluster object.
        Args:
          defline (str): The defline of the cluster.
          sequences (List[ClusterSequence]): The list of sequences in the cluster.
        Side Effects:
          Initializes the Cluster object.
        Examples:
          >>> cluster = Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])
          >>> cluster.name
          'Cluster_1'
          >>> cluster.sequences
          [ClusterSequence(id='seq1', name='seq1', length=0, identity=0.0, is_ref=False, seqtype=None, strand=None), ClusterSequence(id='seq2', name='seq2', length=0, identity=0.0, is_ref=False, seqtype=None, strand=None)]
        """
        self.name = defline
        self.sequences: List[ClusterSequence] = sequences
        self.refname = self._getref(sequences)

    def __repr__(self) -> str:
        """
        Returns a string representation of the Cluster object.
        Args:
          self (Cluster): The Cluster object.
        Returns:
          str: A string representation of the Cluster object.
        Examples:
          >>> cluster = Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])
          >>> cluster.__repr__()
          'Cluster(name=Cluster_1, len=2)'
        """
        return f"Cluster(name={self.name}, len={len(self.sequences)})"

    def _getref(self, sequences: List[ClusterSequence]) -> str:
        """
        Returns the name of the reference sequence in the cluster.
        Args:
          self (Cluster): The Cluster object.
          sequences (List[ClusterSequence]): The list of sequences in the cluster.
        Returns:
          str: The name of the reference sequence in the cluster.
        Examples:
          >>> cluster = Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])
          >>> cluster._getref([ClusterSequence("seq1"), ClusterSequence("seq2")])
          None
        """
        for seq in sequences:
            if seq.is_ref:
                return seq.name
        return None

    def __len__(self):
        """
        Returns the length of the cluster.
        Args:
          self (Cluster): The Cluster object.
        Returns:
          int: The length of the cluster.
        Examples:
          >>> cluster = Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])
          >>> cluster.__len__()
          2
        """
        return len(self.sequences)


class Clustering:
    """
    Represents a clustering of sequences.
    Args:
      name (str): The name of the clustering.
      clusters (List[Cluster]): The list of clusters in the clustering.
    Attributes:
      name (str): The name of the clustering.
      clusters (List[Cluster]): The list of clusters in the clustering.
      seqcluster (dict): A dictionary mapping sequence names to cluster names.
    Examples:
      >>> clustering = Clustering("clustering_1", [Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])])
      >>> clustering.name
      'clustering_1'
      >>> clustering.clusters
      [Cluster(name='Cluster_1', len=2)]
      >>> clustering.seqcluster
      {'seq1': 'Cluster_1', 'seq2': 'Cluster_1'}
    """

    def __init__(self, name, clusters):
        """
        Initializes a Clustering object.
        Args:
          name (str): The name of the clustering.
          clusters (List[Cluster]): The list of clusters in the clustering.
        Side Effects:
          Initializes the Clustering object.
        Examples:
          >>> clustering = Clustering("clustering_1", [Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])])
          >>> clustering.name
          'clustering_1'
          >>> clustering.clusters
          [Cluster(name='Cluster_1', len=2)]
        """
        self.name = name
        self.clusters = clusters
        self.seqcluster = self._todict()

    def __length__(self):
        """
        Returns the length of the clustering.
        Args:
          self (Clustering): The Clustering object.
        Returns:
          int: The length of the clustering.
        Examples:
          >>> clustering = Clustering("clustering_1", [Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])])
          >>> clustering.__length__()
          1
        """
        return len(self.clusters)

    def _todict(self):
        """
        Returns a dictionary mapping sequence names to cluster names.
        Args:
          self (Clustering): The Clustering object.
        Returns:
          dict: A dictionary mapping sequence names to cluster names.
        Examples:
          >>> clustering = Clustering("clustering_1", [Cluster("Cluster_1", [ClusterSequence("seq1"), ClusterSequence("seq2")])])
          >>> clustering._todict()
          {'seq1': 'Cluster_1', 'seq2': 'Cluster_1'}
        """
        for cluster in self.clusters:
            for seq in cluster.sequences:
                self.seqcluster[seq.name] = cluster.name


class FastaReader:
    """
    FASTA reader
    """

    def __init__(self, file: Union[str, Path, IO[str]]):
        """
        Parameters
        ----------
        file
            File path or IO stream.
        """
        if isinstance(file, str):
            file = Path(file)

        if isinstance(file, Path):
            file = xopen(file, "r")

        self._file = file
        self._seq = ""
        self._lines = peekable(line for line in file)
        self._line_number = 0

    def read_item(self) -> Cluster:
        """
        Get the next item.

        Returns
        -------
        Next item.
        """
        defline = self._next_defline()
        sequences = self._next_sequences()
        return Cluster(defline, sequences)

    def read_items(self) -> List[Cluster]:
        """
        Get the list of all items.

        Returns
        -------
        List of all items.
        """
        return list(self)

    def close(self):
        """
        Close the associated stream.
        """
        self._file.close()

    def _next_defline(self) -> str:
        """
        Returns the next defline in the FASTA file.
        Args:
          self (FastaReader): The FastaReader object.
        Returns:
          str: The next defline in the FASTA file.
        Raises:
          StopIteration: If the end of the file is reached.
        Examples:
          >>> fasta_reader = FastaReader(open("example.fasta"))
          >>> fasta_reader._next_defline()
          'seq1'
        """
        while True:
            line = next(self._lines)
            self._line_number += 1
            if line == "":
                raise StopIteration

            line = line.strip()
            if line.startswith(">"):
                return line[1:]
            if line != "":
                raise ParsingError(self._line_number)

    def _next_sequences(self) -> str:
        """
        Returns the next list of sequences in the FASTA file.
        Args:
          self (FastaReader): The FastaReader object.
        Returns:
          List[ClusterSequence]: The next list of sequences in the FASTA file.
        Raises:
          ParsingError: If the line is not a valid FASTA line.
        Examples:
          >>> fasta_reader = FastaReader(open("example.fasta"))
          >>> fasta_reader._next_sequences()
          [ClusterSequence(id='seq1', name='seq1', length=0, identity=0.0, is_ref=False, seqtype=None, strand=None), ClusterSequence(id='seq2', name='seq2', length=0, identity=0.0, is_ref=False, seqtype=None, strand=None)]
        """
        clusterSequences = []
        while True:
            line = next(self._lines)

            self._line_number += 1
            if line == "":
                raise ParsingError(self._line_number)

            line = line.strip()
            if not line.startswith(">"):
                clusterSequences.append(ClusterSequence(line.strip()))
                if self._sequence_continues():
                    continue
                return clusterSequences
            if line != "":
                raise ParsingError(self._line_number)

    def _sequence_continues(self):
        """
        Checks if the next line is a valid sequence line.
        Returns:
          bool: True if the next line is a valid sequence line, False otherwise.
        """
        try:
            next_line = self._lines.peek()
        except StopIteration:
            return False

        if next_line == "":
            return False
        next_line = next_line.strip()
        return len(next_line) > 0 and not next_line.startswith(">")

    def __iter__(self) -> Iterator[Cluster]:
        """
        Iterates over the FASTA file.
        Yields:
          Cluster: The next cluster in the FASTA file.
        """
        while True:
            try:
                yield self.read_item()
            except StopIteration:
                return

    def __enter__(self):
        """
        Enters the context manager.
        Returns:
          FastaReader: The FastaReader instance.
        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exits the context manager.
        Args:
          exception_type (type): The type of exception raised.
          exception_value (Exception): The exception raised.
          traceback (Traceback): The traceback of the exception.
        """
        del exception_type
        del exception_value
        del traceback
        self.close()


class ClstrReader:
    """
    CD-HIT (Clstr) reader.
    """

    def __init__(self, file: Union[str, Path, IO[str]]):
        """
        Parameters
        ----------
        file
            File path or IO stream.
        """
        if isinstance(file, str):
            file = Path(file)

        if isinstance(file, Path):
            file = xopen(file, "r")

        self._file = file
        self._clusterSequences = []
        self._lines = peekable(line for line in file)
        self._line_number = 0

    def read_item(self) -> Cluster:
        """
        Get the next item.

        Returns
        -------
        Next item.
        """
        defline = self._next_defline()
        sequences = self._next_sequences()
        return Cluster(defline, sequences)

    def read_items(self) -> List[Cluster]:
        """
        Get the list of all items.

        Returns
        -------
        List of all items.
        """
        return list(self)

    def close(self):
        """
        Close the associated stream.
        """
        self._file.close()

    def _next_defline(self) -> str:
        """
        Reads the next definition line from the Clstr file.
        Returns:
          str: The next definition line.
        Raises:
          StopIteration: If the end of the file is reached.
          ParsingError: If the line is not a valid definition line.
        """
        while True:
            line = next(self._lines)
            self._line_number += 1
            if line == "":
                raise StopIteration

            line = line.strip()
            if line.startswith(">"):
                return line[1:]
            if line != "":
                raise ParsingError(self._line_number)

    def _next_sequences(self) -> str:
        """
        Reads the next sequence lines from the Clstr file.
        Returns:
          List[ClusterSequence]: The next sequence lines.
        Raises:
          StopIteration: If the end of the file is reached.
          ParsingError: If the line is not a valid sequence line.
        """
        clusterSequences = []
        while True:
            line = next(self._lines)

            self._line_number += 1
            if line == "":
                raise ParsingError(self._line_number)

            line = line.strip()
            if not line.startswith(">"):
                clusterSequences.append(ClusterSequence(line.strip()))
                if self._sequence_continues():
                    continue
                return clusterSequences
            if line != "":
                raise ParsingError(self._line_number)

    def _sequence_continues(self):
        """
        Checks if the next line is a valid sequence line.
        Returns:
          bool: True if the next line is a valid sequence line, False otherwise.
        """
        try:
            next_line = self._lines.peek()
        except StopIteration:
            return False

        if next_line == "":
            return False
        next_line = next_line.strip()
        return len(next_line) > 0 and not next_line.startswith(">")

    def __iter__(self) -> Iterator[Cluster]:
        """
        Iterates over the Clstr file.
        Yields:
          Cluster: The next cluster in the Clstr file.
        """
        while True:
            try:
                yield self.read_item()
            except StopIteration:
                return

    def __enter__(self):
        """
        Enters the context manager.
        Returns:
          ClstrReader: The ClstrReader instance.
        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exits the context manager.
        Args:
          exception_type (type): The type of exception raised.
          exception_value (Exception): The exception raised.
          traceback (Traceback): The traceback of the exception.
        """
        del exception_type
        del exception_value
        del traceback
        self.close()


def read_cdhit(file: Union[str, Path, IO[str]]) -> ClstrReader:
    """
    Open a CD-HIT file for reading.

    Parameters
    ----------
    file
        File path or IO stream.

    Returns
    -------
    CD-HIT (Clstr) reader.
    """
    return ClstrReader(file)


def read_fasta(file: Union[str, Path, IO[str]]) -> FastaReader:
    """
    Open a FASTA file for reading.

    Parameters
    ----------
    file
        File path or IO stream.

    Returns
    -------
    FASTA reader.
    """
    return FastaReader(file)
