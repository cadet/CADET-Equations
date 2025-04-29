import os
import pytest

@pytest.fixture
def test_dir():
    """Returns the path to the test data directory."""
    return os.path.dirname(__file__)