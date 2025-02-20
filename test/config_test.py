# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

from streamlit.testing.v1 import AppTest


def test_streamlit_app():
    at = AppTest.from_file("../app.py")
    at.run()
    assert not at.exception