# -*- coding: utf-8 -*-
"""
SPDX-License-Identifier: GPL-3.0-only
Copyright (c) 2026 Jan Michael Breuer
See LICENSE file for details.


Unit tests for rerender_variables() from Equation-Generator.py.
Covers CADET format, Legacy format, and error handling.
"""

import pytest
from importlib import import_module

# Equation-Generator.py has a hyphenated name, so we use importlib
_eg = import_module("Equation-Generator")
rerender_variables = _eg.rerender_variables


# %% CADET format

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("input_str, expected_frag", [
    (r"c^{\b}",  r"\mathrm{b}"),
    (r"c^{\p}",  r"\mathrm{p}"),
    (r"c^{\s}",  r"\mathrm{s}"),
    (r"c^{\l}",  r"\mathrm{\ell}"),
])
def test_rerender_CADET_superscript_substitution(input_str, expected_frag):
    """CADET format should wrap short LaTeX phase commands in \\mathrm."""
    result = rerender_variables(input_str, "CADET")
    assert expected_frag in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_rerender_CADET_preserves_longer_commands():
    """CADET format should not modify longer LaTeX commands like \\partial."""
    assert rerender_variables(r"\partial", "CADET") == r"\partial"


# %% Legacy format

@pytest.mark.ci
@pytest.mark.unit_test
def test_rerender_Legacy_bulk_conc():
    """Legacy format should strip the superscript from bulk concentration."""
    assert rerender_variables(r"c^{\b}", "Legacy") == "c"


@pytest.mark.ci
@pytest.mark.unit_test
def test_rerender_Legacy_particle_conc():
    """Legacy format should move the particle phase to a subscript."""
    assert rerender_variables(r"c^{\p}_{i}", "Legacy") == r"c_{\mathrm{p},i}"


@pytest.mark.ci
@pytest.mark.unit_test
def test_rerender_Legacy_solid_conc():
    """Legacy format should replace solid concentration c^s with q."""
    result = rerender_variables(r"c^{\s}", "Legacy")
    assert "q" in result


# %% Error handling

@pytest.mark.ci
@pytest.mark.unit_test
def test_rerender_invalid_format_raises():
    """An unsupported variable format should raise ValueError."""
    with pytest.raises(ValueError, match="not supported"):
        rerender_variables(r"c^{\b}", "InvalidFormat")
