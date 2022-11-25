from __future__ import annotations
from pathlib import Path
from typing import IO, Iterator, List, Union
from enum import Enum
from more_itertools import peekable
from xopen import xopen
import re

__all__ = ["ParsingError", "ClusterSequence", "Cluster", "ClstrReader", "read_cdhit", "SeqType", "Strand", "FastaReader", "read_fasta"]

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
        pattern = re.compile(r'(?P<id>\d+)\s+(?P<size>\d+)(?P<type>aa|nt), >(?P<name>.+?)\.\.\. (?P<attr>.+)')
        attrpatt = re.compile(r'(?P<ref>\*|at) .*?(?P<strand>[+-]?)\/?(?P<percent>\d+\.\d+)%')
        match = pattern.search(self.line)

        self.seqtype = SeqType.PROTEIN if match["type"] == "aa" else SeqType.NT
        self.strand  = Strand.NONE if self.seqtype == SeqType.PROTEIN else Strand.PLUS
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
                self.strand = Strand.PLUS if attrs["strand"] == "+" else Strand.REVERSE if attrs["strand"] == "-" else Strand.NONE

    def __repr__(self):
        return f"ClusterSequence(id={self.id}, name={self.name}, length={self.length}, identity={self.identity}, is_ref={self.is_ref}, seqtype={self.seqtype}, strand={self.strand})"

 

class Cluster:
    def __init__(self, defline, sequences):
        self.name = defline
        self.sequences: List[ClusterSequence] = sequences
        self.refname = self._getref(sequences)
    
    def __repr__(self) -> str:
        return f"Cluster(name={self.name}, len={len(self.sequences)})"
    
    def _getref(self, sequences: List[ClusterSequence]) -> str:
        for seq in sequences:
            if seq.is_ref:
                return seq.name
        return None

    def __len__(self):
        return len(self.sequences)

class Clustering:
    def __init__(self, name, clusters):
        self.name = name
        self.clusters = clusters
        self.seqcluster = self._todict() 
    
    def __length__(self):
        return len(self.clusters)
    
    def _todict(self):
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
        clusterSequences = []
        while True:
            line = next(self._lines)
            
            self._line_number += 1
            if line == "":
                raise ParsingError(self._line_number)

            line = line.strip()
            if not line.startswith(">"):
                clusterSequences.append( ClusterSequence(line.strip()) )
                if self._sequence_continues():
                    continue
                return clusterSequences
            if line != "":
                raise ParsingError(self._line_number)
     

    def _sequence_continues(self):
        try:
            next_line = self._lines.peek()
        except StopIteration:
            return False

        if next_line == "":
            return False
        next_line = next_line.strip()
        return len(next_line) > 0 and not next_line.startswith(">")

    def __iter__(self) -> Iterator[Cluster]:
        while True:
            try:
                yield self.read_item()
            except StopIteration:
                return

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
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
        clusterSequences = []
        while True:
            line = next(self._lines)
            
            self._line_number += 1
            if line == "":
                raise ParsingError(self._line_number)

            line = line.strip()
            if not line.startswith(">"):
                clusterSequences.append( ClusterSequence(line.strip()) )
                if self._sequence_continues():
                    continue
                return clusterSequences
            if line != "":
                raise ParsingError(self._line_number)
     

    def _sequence_continues(self):
        try:
            next_line = self._lines.peek()
        except StopIteration:
            return False

        if next_line == "":
            return False
        next_line = next_line.strip()
        return len(next_line) > 0 and not next_line.startswith(">")

    def __iter__(self) -> Iterator[Cluster]:
        while True:
            try:
                yield self.read_item()
            except StopIteration:
                return

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
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