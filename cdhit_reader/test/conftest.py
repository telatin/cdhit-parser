import pytest


@pytest.fixture
def correct1(tmp_path):
    return _write_file(tmp_path, "correct1.faa")


@pytest.fixture
def correct2(tmp_path):
    return _write_file(tmp_path, "correct2.faa")


@pytest.fixture
def damaged1(tmp_path):
    return _write_file(tmp_path, "damaged1.faa")


@pytest.fixture
def damaged2(tmp_path):
    return _write_file(tmp_path, "damaged2.faa")


@pytest.fixture
def damaged3(tmp_path):
    return _write_file(tmp_path, "damaged3.faa")


@pytest.fixture
def protein(tmp_path):
    return _write_file(tmp_path, "protein.faa")


@pytest.fixture
def protein_gzip(tmp_path):
    return _write_file(tmp_path, "protein.faa.gz")


def _write_file(path, filename):
    import importlib_resources as pkg_resources

    import cdhit_reader

    content = pkg_resources.read_binary(cdhit_reader.test, filename)

    with open(path / filename, "wb") as f:
        f.write(content)

    return path / filename
