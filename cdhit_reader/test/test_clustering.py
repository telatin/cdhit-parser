from io import StringIO
from pathlib import Path

from cdhit_reader import Clustering

DATA = Path(__file__).resolve().parents[2] / "data" / "aa.clstr"


def test_from_file_builds_reverse_index():
    clustering = Clustering.from_file(DATA)

    assert clustering.name == "aa"
    assert len(clustering) == 7
    assert len(clustering.seqcluster) == 10
    assert clustering.seqcluster["IBJJOHBJ_00007"] == "Cluster 0"
    assert clustering.seqcluster["IBJJOHBJ_00001"] == "Cluster 1"
    assert clustering.seqcluster["CBJJOHBJ_000C6"] == "Cluster 2"


def test_constructor_with_explicit_name():
    clusters = Clustering.from_file(DATA).clusters
    clustering = Clustering("my_run", clusters)

    assert clustering.name == "my_run"
    assert len(clustering) == 7


def test_from_file_name_override():
    clustering = Clustering.from_file(DATA, name="custom")

    assert clustering.name == "custom"


def test_from_stream_has_empty_default_name():
    clustering = Clustering.from_file(StringIO(">Cluster 0\n0\t492nt, >seq1.A... *\n"))

    assert clustering.name == ""
    assert clustering.seqcluster == {"seq1.A": "Cluster 0"}


def test_iter_yields_all_clusters():
    clustering = Clustering.from_file(DATA)
    names = [cluster.name for cluster in clustering]

    assert names == [f"Cluster {n}" for n in range(7)]


def test_every_sequence_is_indexed():
    clustering = Clustering.from_file(DATA)
    for cluster in clustering:
        for seq in cluster.sequences:
            assert clustering.seqcluster[seq.name] == cluster.name


def test_empty_clustering():
    clustering = Clustering("empty", [])

    assert len(clustering) == 0
    assert clustering.seqcluster == {}
    assert list(clustering) == []


def test_repr():
    clustering = Clustering.from_file(DATA)

    assert repr(clustering) == "Clustering(name=aa, clusters=7, sequences=10)"
