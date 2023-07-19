import os
from pathlib import Path

import pytest
import sys
from cdhit_reader import ParsingError, read_cdhit, SeqType, Strand

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
        assert _test_cluster_structure(cluster) == True
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
        
def _test_cluster_structure(cluster):
    # name
    assert cluster.name is not None
    # refname
    assert cluster.refname is not None
    assert cluster.refname in [s.name for s in cluster.sequences]
    # sequences
    assert len(cluster.sequences) > 0
    assert [_test_seq(s)  == True for s in cluster.sequences]
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