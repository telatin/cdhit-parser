import os
from pathlib import Path

import pytest

from cdhit_reader import ParsingError, read_fasta, write_fasta


def _test_read_fasta_correct(filepath: Path):
    deflines = ["ID1", "ID2", "ID3", "ID4"]
    sequences = ["GAGUUA", "CAUAACAAATT", "AAGAA", "AAGAA"]

    f = read_fasta(filepath)
    item = f.read_item()
    assert item.defline == deflines[0]
    assert item.id == deflines[0]
    assert not item.has_desc
    assert item.sequence == sequences[0]

    item = f.read_item()
    assert item.defline == deflines[1]
    assert item.id == deflines[1]
    assert not item.has_desc
    assert item.sequence == sequences[1]

    item = f.read_item()
    assert item.defline == deflines[2]
    assert item.id == deflines[2]
    assert not item.has_desc
    assert item.sequence == sequences[2]

    item = f.read_item()
    assert item.defline == deflines[3]
    assert item.id == deflines[3]
    assert not item.has_desc
    assert item.sequence == sequences[3]

    with pytest.raises(StopIteration):
        f.read_item()

    f.close()

    f = read_fasta(filepath)
    for i, item in enumerate(f):
        assert item.defline == deflines[i]
        assert item.sequence == sequences[i]
    f.close()

    f = read_fasta(filepath)
    items = f.read_items()
    for i, defline in enumerate(deflines):
        assert items[i].defline == defline
        assert items[i].sequence == sequences[i]
    f.close()

    with read_fasta(filepath) as f:
        for i, item in enumerate(f):
            assert item.defline == deflines[i]
            assert item.sequence == sequences[i]

    with read_fasta(filepath) as f:
        f.close()


def test_read_fasta_correct1(correct1: Path):
    _test_read_fasta_correct(correct1)


def test_read_fasta_correct2(correct2: Path):
    _test_read_fasta_correct(correct2)


def test_read_fasta_protein(protein: Path, protein_gzip: Path):
    expected = """
QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE
KMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTS
VLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHP
FLFLIKHNPTNTIVYFGRYWSP
    """.replace(
        "\n", ""
    ).strip()

    with read_fasta(protein) as file:
        item = file.read_item()
        assert item.defline == "P01013 GENE X PROTEIN (OVALBUMIN-RELATED)"
        assert item.has_desc
        assert item.desc == "GENE X PROTEIN (OVALBUMIN-RELATED)"
        assert item.id == "P01013"
        assert item.sequence == expected

    assert len(list(read_fasta(protein))) == 1

    with read_fasta(protein_gzip) as file:
        item = file.read_item()
        assert item.defline == "P01013 GENE X PROTEIN (OVALBUMIN-RELATED)"
        assert item.has_desc
        assert item.desc == "GENE X PROTEIN (OVALBUMIN-RELATED)"
        assert item.id == "P01013"
        assert item.sequence == expected

    assert len(list(read_fasta(protein))) == 1


def test_read_fasta_damaged(damaged1, damaged2, damaged3):

    with read_fasta(damaged1) as f:
        with pytest.raises(ParsingError) as excinfo:
            f.read_item()
        e: ParsingError = excinfo.value
        assert e.line_number == 1

    with read_fasta(damaged2) as f:
        with pytest.raises(ParsingError) as excinfo:
            f.read_item()
        e: ParsingError = excinfo.value
        assert e.line_number == 2

    with read_fasta(damaged3) as f:
        f.read_item()
        with pytest.raises(ParsingError) as excinfo:
            f.read_item()
        e: ParsingError = excinfo.value
        assert e.line_number == 4


def test_write_fasta(tmp_path: Path):

    defline = ["defline1", "defline2 description"]
    sequence = ["ABCD", "ABCD" * 100]

    os.chdir(tmp_path)
    with write_fasta("output.faa") as writer:
        writer.write_item(defline[0], sequence[0])
        writer.write_item(defline[1], sequence[1])

    with read_fasta("output.faa") as reader:
        item = reader.read_item()
        assert item.defline == defline[0]
        assert item.sequence == sequence[0]

        item = reader.read_item()
        assert item.defline == defline[1]
        assert item.sequence == sequence[1]

    with write_fasta("output.faa.xz") as writer:
        writer.write_item(defline[0], sequence[0])
        writer.write_item(defline[1], sequence[1])

    with read_fasta("output.faa.xz") as reader:
        item = reader.read_item()
        assert item.defline == defline[0]
        assert item.sequence == sequence[0]

        item = reader.read_item()
        assert item.defline == defline[1]
        assert item.sequence == sequence[1]
