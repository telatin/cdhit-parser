from cdhit_reader._compare import TAG1, TAG2, has_cdhit, strip_tag


def test_strip_tag_keeps_full_sequence_name():
    assert strip_tag(TAG1 + "IBJJOHBJ_00007") == "IBJJOHBJ_00007"
    assert strip_tag(TAG2 + "seq.with.dots") == "seq.with.dots"


def test_has_cdhit_returns_bool():
    assert isinstance(has_cdhit(), bool)
