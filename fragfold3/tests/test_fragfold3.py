"""
Unit and regression test for the fragfold3 package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import fragfold3


def test_fragfold3_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "fragfold3" in sys.modules
