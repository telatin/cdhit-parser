__all__ = ["test"]


def test(verbose=True) -> int:
    """
    Run tests to verify this package's integrity.

    Parameters
    ----------
    verbose
        ``True`` to show diagnostic. Defaults to ``True``.

    Returns
    -------
    Exit code: ``0`` for success.
    """

    args = ["--doctest-modules", "-k", "not test_testit"]
    if not verbose:
        args += ["--quiet"]

    args += ["--pyargs", __name__.split(".")[0]]
    return __import__("pytest").main(args)
