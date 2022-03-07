import os
from pathlib import Path

import pytest
from cdhit_reader import ParsingError, read_cdhit, SeqType

def test_nt():
    input = "small_nt.clstr"
    
    cluster1 = ["seq1.A", "seq1.B", "seq1.C", "seq1.D"]
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

        if n == 0:
            assert len(cluster.sequences) == len(cluster1)
            assert [s.name for s in cluster.sequences] == cluster1
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
        n += 1
        