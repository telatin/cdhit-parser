from __future__ import annotations
from pathlib import Path
from typing import IO, Iterator, List, Union
from enum import Enum
from more_itertools import peekable
from xopen import xopen
import re

__all__ = ["Sequence", "FastaReader", "read_fasta"]
 
 

class Sequence:
    def __init__(self, name, sequence, comment=None, separator=" ", line_length=0):
        self.name = name
        self.sequence = sequence
        self.comment = comment
        self.separator = separator
        self.line_length = line_length
    
    def __repr__(self) -> str:
        comment_string = self.separator + self.comment if self.comment else ""
        sequence = self.sequence if self.line_length == 0 else self.wrap(self.sequence, self.line_length)
        return f">{self.name}{comment_string}\n{sequence}"
    

    def wrap(self, sequence, line_length):
        return "\n".join(sequence[i:i + line_length] for i in range(0, len(sequence), line_length))
    def __len__(self):
        return len(self.sequence)
 


class FastaReader:
    """
    FASTA reader
    """    
    def __init__(self, file: Union[str, Path, IO[str]], separator=" ", line_len=0):
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

        self.separator = separator
        self.line_len = line_len
        self._file = file
        self._seq = ""
        self._lines = peekable(line for line in file)
        self._line_number = 0

    def read_item(self) -> Sequence:
        """
        Get the next item.

        Returns
        -------
        Next item.
        """
        defline = self._next_defline()
        name = defline.split(maxsplit=1)[0]
        comment = defline.split(maxsplit=1)[1] if len(defline.split(maxsplit=1)) > 1 else None
        sequence = self._next_sequences()
        return Sequence(name, sequence, comment, separator=self.separator, line_length=self.line_len)

    def read_items(self) -> List[Sequence]:
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
                raise "Invalid FASTA file"

    def _next_sequences(self) -> str:
        seq = ""
        while True:
            line = next(self._lines).strip()
            
            self._line_number += 1
            if line == "":
                raise

            line = line.strip()
            if not line.startswith(">"):
                seq += line
                if self._sequence_continues():
                    continue
                return seq
            if line != "":
                raise 
     

    def _sequence_continues(self):
        try:
            next_line = self._lines.peek()
        except StopIteration:
            return False

        if next_line == "":
            return False
        next_line = next_line.strip()
        return len(next_line) > 0 and not next_line.startswith(">")

    def __iter__(self) -> Iterator[Sequence]:
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
            
def read_fasta(file: Union[str, Path, IO[str]], separator=" ", line_len = 0) -> FastaReader:
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
    return FastaReader(file, separator=separator, line_len=line_len)