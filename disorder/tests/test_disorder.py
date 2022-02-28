"""
Unit and regression test for the disorder package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import disorder


def test_disorder_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "disorder" in sys.modules
