"""The GREAT 3D suite for integrative multi-omics 3D epigenomics."""

# Add imports here
import time

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

__version__ = get_versions()['version']
del get_versions


def timing(func):
    """Wrapper to time functions.py.

    Works as a decorator, adapted from https://stackoverflow.com/questions/5478351/python-time-measure-function
    """
    def wrap(*args):
        time1 = time.time()
        ret = func(*args)
        time2 = time.time()
        print(f"{func.__name__} function took {(time2 - time1) * 1000.0:.3f} ms")
        return ret
    return wrap
