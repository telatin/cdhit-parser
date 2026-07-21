import gzip
import os
from pathlib import Path
from io import StringIO

import pytest
from cdhit_reader import ClusterSequence, ParsingError, read_cdhit, SeqType, Strand


def test_nt():
    input = "small_nt.clstr"

    cluster1names = ["seq1.A", "seq1.B", "seq1.C", "seq1.D"]
    cluster1strands = [Strand.PLUS, Strand.PLUS, Strand.PLUS, Strand.REVERSE]
    # List files in path
    path = Path(os.path.dirname(__file__))
    filePath = os.path.join(path, input)
    if not os.path.exists(filePath):
        pytest.skip("File not found: {}".format(filePath))

    n = 0
    for cluster in read_cdhit(filePath):
        # Cluster name is a progressive number from 0 to n-1
        assert cluster.name == f"Cluster {n}"
        assert len(cluster) > 0
        for item in cluster.sequences:
            assert item.seqtype == SeqType.NT
        assert _test_cluster_structure(cluster)
        if n == 0:
            # This is the first cluster, compare with the expected values
            assert len(cluster.sequences) == len(cluster1names)
            assert cluster.refname == cluster1names[0]
            assert [s.name for s in cluster.sequences] == cluster1names
            # Compare strands
            for s, st in zip(cluster.sequences, cluster1strands):
                assert s.strand == st
        n += 1


def test_aa():
    input = "small_aa.clstr"

    cluster2 = ["IBJJOHBJ_00001", "IBJJOHBJ_000F1"]
    # List files in path
    path = Path(os.path.dirname(__file__))
    filePath = os.path.join(path, input)
    if not os.path.exists(filePath):
        pytest.skip("File not found: {}".format(filePath))

    n = 0
    for cluster in read_cdhit(filePath):
        # Cluster name is a progressive number from 0 to n-1
        assert cluster.name == f"Cluster {n}"
        assert len(cluster) > 0
        for item in cluster.sequences:
            assert item.seqtype == SeqType.PROTEIN

        if n == 1:
            assert len(cluster.sequences) == len(cluster2)
            assert [s.name for s in cluster.sequences] == cluster2
            assert cluster.refname == cluster2[0]
        n += 1


def test_invalid_member_line_raises_parsing_error():
    with pytest.raises(ParsingError) as excinfo:
        list(read_cdhit(StringIO(">Cluster 0\n0 492nt, >seq1.A...\n")))

    assert excinfo.value.line_number == 2


def test_cluster_without_reference_has_none_refname():
    clusters = list(read_cdhit(StringIO(">Cluster 0\n0\t492nt, >seq1.A... at +/99.39%\n")))

    assert len(clusters) == 1
    assert clusters[0].refname is None


def test_reference_member_attributes():
    seq = ClusterSequence("0\t1016aa, >IBJJOHBJ_00007... *")

    assert seq.is_ref
    assert seq.identity == 100.0
    assert seq.id == 0
    assert seq.length == 1016


def test_member_with_coordinate_attributes():
    # Attribute format produced with alignment coordinates
    seq = ClusterSequence("3\t502nt, >IKXM6KN01CFAFI... at 1:502:1:503/+/97.81%")

    assert not seq.is_ref
    assert seq.id == 3
    assert seq.length == 502
    assert seq.name == "IKXM6KN01CFAFI"
    assert seq.identity == 97.81
    assert seq.strand == Strand.PLUS


def test_member_with_integer_percent():
    seq = ClusterSequence("0\t100nt, >x... at +/100%")

    assert not seq.is_ref
    assert seq.identity == 100.0
    assert seq.strand == Strand.PLUS


def test_aa_member_has_no_strand():
    seq = ClusterSequence("1\t366aa, >IBJJOHBJ_000F1... at 98.91%")

    assert seq.seqtype == SeqType.PROTEIN
    assert seq.strand == Strand.NONE
    assert seq.identity == 98.91


def test_blank_lines_around_clusters_are_tolerated():
    content = "\n>Cluster 0\n0\t492nt, >seq1.A... *\n\n>Cluster 1\n0\t100nt, >s2... *\n\n"
    clusters = list(read_cdhit(StringIO(content)))

    assert [cluster.name for cluster in clusters] == ["Cluster 0", "Cluster 1"]
    assert [cluster.sequences[0].name for cluster in clusters] == ["seq1.A", "s2"]


def test_crlf_line_endings():
    content = ">Cluster 0\r\n0\t492nt, >seq1.A... *\r\n1\t492nt, >seq1.B... at +/99.39%\r\n"
    clusters = list(read_cdhit(StringIO(content)))

    assert len(clusters) == 1
    assert clusters[0].sequences[1].identity == 99.39


def test_gzipped_clstr(tmp_path):
    gz_path = tmp_path / "tiny.clstr.gz"
    with gzip.open(gz_path, "wt") as handle:
        handle.write(">Cluster 0\n0\t492nt, >seq1.A... *\n")

    clusters = list(read_cdhit(gz_path))
    assert clusters[0].sequences[0].name == "seq1.A"


def test_merge_clstr_file_parses():
    # Real-world output of clstr_merge.pl, with dotted sequence names
    input_path = Path(__file__).resolve().parents[2] / "data" / "merge-clust.clstr"
    clusters = list(read_cdhit(input_path))

    assert len(clusters) == 8
    assert sum(len(cluster) for cluster in clusters) == 18
    assert clusters[0].refname == "ciao.IBJJOHBJ_00007"
    for cluster in clusters:
        for seq in cluster.sequences:
            assert seq.seqtype == SeqType.PROTEIN


def _test_cluster_structure(cluster):
    # name
    assert cluster.name is not None
    # refname
    assert cluster.refname is not None
    assert cluster.refname in [s.name for s in cluster.sequences]
    # sequences
    assert len(cluster.sequences) > 0
    assert [_test_seq(s) for s in cluster.sequences]
    return True


def _test_seq(s):
    assert s.id is not None
    assert int(s.id) >= 0
    assert s.identity >= 0
    assert s.identity <= 100
    assert s.length > 0
    assert int(s.length) > 0

    assert s.seqtype is not None
    assert s.seqtype == SeqType.NT or s.seqtype == SeqType.PROTEIN

    assert s.strand is not None
    assert s.strand == Strand.PLUS or s.strand == Strand.REVERSE
    return True
