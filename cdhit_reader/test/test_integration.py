from pathlib import Path

from cdhit_reader import __version__, read_cdhit


def test_public_api_smoke_on_aa_fixture(capsys):
    input_path = Path(__file__).resolve().parents[2] / "data" / "aa.clstr"

    clusters = read_cdhit(input_path).read_items()
    assert len(clusters) == 7
    assert sum(len(cluster) for cluster in clusters) == 10

    clust_count = 0
    for cluster in read_cdhit(input_path):
        clust_count += 1
        print(f"\n[{clust_count}] {cluster.name} refSequence={cluster.refname} size={len(cluster)}")
        for member_count, member in enumerate(cluster.sequences, start=1):
            print(
                f"{'REF' if member.is_ref else 'MEM'}\t{member_count}: "
                f"{member.name} ({member.length}) %ID={member.identity}%"
            )

    print(f"\nCD-HIT_reader version: {__version__}")

    captured = capsys.readouterr().out
    assert "[1] Cluster 0 refSequence=IBJJOHBJ_00007 size=1" in captured
    assert "[7] Cluster 6 refSequence=IBJJOHBJ_00004 size=1" in captured
    assert "REF\t1: IBJJOHBJ_00001 (371) %ID=100.0%" in captured
    assert "MEM\t2: IBJJOHBJ_000F1 (366) %ID=98.91%" in captured
    assert f"CD-HIT_reader version: {__version__}" in captured
