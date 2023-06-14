import re
from os.path import join

from setuptools import find_packages


def get():
    """
    Gets the version of the package.
    Args:
      None
    Returns:
      str: The version of the package.
    Examples:
      >>> get()
      '1.0.0'
    """
    pkgnames = find_packages()
    if len(pkgnames) == 0:
        return "unknown"
    pkgname = pkgnames[0]
    content = open(join(pkgname, "_version.py")).read()
    c = re.compile(r"__version__ *= *('[^']+'|\"[^\"]+\")")
    m = c.search(content)
    if m is None:
        return "unknown"
    return m.groups()[0][1:-1]
