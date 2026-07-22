"""
@author: jmbr
"""

import os

import pytest
import streamlit as st


@pytest.fixture
def test_dir():
    """Returns the path to the test data directory."""
    return os.path.dirname(__file__)


@pytest.fixture(autouse=True)
def skip_logo(monkeypatch):
    monkeypatch.setattr(st, "logo", lambda *args, **kwargs: None)
