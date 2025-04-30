# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import streamlit as st
import os
import pytest


@pytest.fixture
def test_dir():
    """Returns the path to the test data directory."""
    return os.path.dirname(__file__)


@pytest.fixture(autouse=True)
def skip_logo(monkeypatch):
    monkeypatch.setattr(st, "logo", lambda *args, **kwargs: None)