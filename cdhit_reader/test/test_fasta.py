from io import StringIO

import pytest

from cdhit_reader._fasta import FastaParsingError, read_fasta


def test_fasta_reads_last_record_without_trailing_newline():
    records = list(read_fasta(StringIO(">seq1 comment\nACG\n>seq2\nTT")))

    assert [record.name for record in records] == ["seq1", "seq2"]
    assert [record.comment for record in records] == ["comment", None]
    assert [record.sequence for record in records] == ["ACG", "TT"]


def test_fasta_rejects_non_header_first_line():
    with pytest.raises(FastaParsingError) as excinfo:
        list(read_fasta(StringIO("seq1\nACG\n")))

    assert excinfo.value.line_number == 1


def test_fasta_rejects_header_without_sequence():
    with pytest.raises(FastaParsingError) as excinfo:
        list(read_fasta(StringIO(">seq1\n")))

    assert excinfo.value.line_number == 1


def test_fasta_rejects_header_without_sequence_before_next_record():
    with pytest.raises(FastaParsingError) as excinfo:
        list(read_fasta(StringIO(">seq1\n>seq2\nTT\n")))

    assert excinfo.value.line_number == 1
