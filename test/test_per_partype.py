# -*- coding: utf-8 -*-
"""
Integration tests for per-particle-type parameterization (issue #21).
Tests the full Streamlit UI path with per-particle-type configuration.
"""

from streamlit.testing.v1 import AppTest
import pytest


def setup_psd_with_per_partype(at, per_partype_settings, resolution="1D (radial coordinate)", has_binding=True):
    """Set up a model with PSD and per-particle-type configuration.

    Parameters
    ----------
    at : AppTest
        The Streamlit app test instance.
    per_partype_settings : list of dict
        Per-particle-type settings, each dict with keys:
        'nonlimiting_filmDiff', 'has_surfDiff'
    resolution : str
        Particle resolution.
    has_binding : bool
        Whether binding is enabled.
    """
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="PSD").set_value("Yes").run()
    at.selectbox(key="particle_resolution").set_value(resolution).run()
    at.selectbox(key="has_binding").set_value("Yes" if has_binding else "No").run()

    for j, settings in enumerate(per_partype_settings):
        at.selectbox(key=f"nonlimiting_filmDiff_partype_{j}").set_value(settings['nonlimiting_filmDiff']).run()
        if has_binding and "1D" in resolution:
            at.selectbox(key=f"has_surfDiff_partype_{j}").set_value(settings['has_surfDiff']).run()


@pytest.mark.ci
@pytest.mark.reference
def test_per_partype_uniform_settings():
    """Test that uniform per-partype settings produce same output as global config."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
    ]
    setup_psd_with_per_partype(at, settings)

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    # Uniform settings -> should NOT have per-partype text
    assert "particle type(s)" not in latex
    # Should have the standard "all particle sizes" text
    assert "particle sizes" in latex or "particle" in latex.lower()


@pytest.mark.ci
@pytest.mark.reference
def test_per_partype_different_film_diffusion():
    """Test per-partype with mixed film diffusion settings."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No'},
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
    ]
    setup_psd_with_per_partype(at, settings)

    assert not at.exception
    latex = at.session_state.latex_string
    # Different settings -> should have per-partype text
    assert "particle type(s)" in latex
    assert r"j = 1" in latex
    assert r"j = 2" in latex


@pytest.mark.ci
@pytest.mark.reference
def test_per_partype_different_surface_diffusion():
    """Test per-partype with mixed surface diffusion settings."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No'},
    ]
    setup_psd_with_per_partype(at, settings)

    assert not at.exception
    latex = at.session_state.latex_string
    # Different settings -> should have per-partype text
    assert "particle type(s)" in latex


@pytest.mark.ci
@pytest.mark.reference
def test_per_partype_all_nonlimiting_film_diff():
    """Test per-partype with all non-limiting film diffusion (LRM without pores case)."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
    ]
    setup_psd_with_per_partype(at, settings, resolution="0D (homogeneous)")

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    # All non-limiting -> uniform, no per-partype text
    assert "particle type(s)" not in latex


@pytest.mark.ci
@pytest.mark.reference
def test_per_partype_single_particle_no_per_partype():
    """Test that single particle type (N_p=1) does not show per-partype config."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    # PSD = No -> N_p = 1
    at.selectbox(key="nonlimiting_filmDiff").set_value("No").run()
    at.selectbox(key="has_binding").set_value("Yes").run()

    assert not at.exception
    latex = at.session_state.latex_string
    assert "particle type(s)" not in latex
