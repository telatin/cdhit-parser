from pathlib import Path

from cdhit_reader import read_cdhit, read_fasta


def test_cl_py_example_matches_expected_cluster_counts():
    input_clst = Path(__file__).resolve().parents[2] / "data" / "aa.clstr"

    total = 0
    members = 0
    rendered = []
    for cluster in read_cdhit(input_clst):
        total += 1
        members += len(cluster)
        rendered.append(repr(cluster))

    assert total == 7
    assert members == 10
    assert rendered[0] == "Cluster(name=Cluster 0, len=1)"
    assert rendered[-1] == "Cluster(name=Cluster 6, len=1)"


def test_fa_py_example_matches_expected_sequence_counts():
    input_fasta = Path(__file__).resolve().parents[2] / "data" / "test.fa.gz"

    total = 0
    comments = 0
    bp = 0
    summaries = []
    for sequence in read_fasta(input_fasta, line_len=60):
        summaries.append((sequence.name, sequence.comment))
        total += 1
        bp += len(sequence)
        if sequence.comment:
            comments += 1

    assert total == 3
    assert comments == 2
    assert bp == 29
    assert summaries == [
        ("name", "comment"),
        ("name", "tabbed\tcomment multi"),
        ("nocomments", None),
    ]
