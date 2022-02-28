"""A set of tools and pipelines for molecular dynamics simulations of protein systems."""

# Add imports here
from .disorder import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
