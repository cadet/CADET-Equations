# -*- coding: utf-8 -*-
"""Tests for per-particle-type parameterization of film diffusion and surface diffusion."""

from streamlit.testing.v1 import AppTest
import pytest


def setup_grm_psd(at, per_partype_settings):
    """Set up a GRM model with PSD and per-particle-type configuration.

    Parameters
    ----------
    at : AppTest
        The Streamlit app test instance.
    per_partype_settings : list of dict
        Per-particle-type settings, each dict with keys:
        'nonlimiting_filmDiff', 'has_surfDiff'
    """
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="PSD").set_value("Yes").run()
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    at.selectbox(key="has_binding").set_value("Yes").run()

    for jj, settings in enumerate(per_partype_settings):
        at.selectbox(key=f"nonlimiting_filmDiff_partype_{jj}").set_value(settings['nonlimiting_filmDiff']).run()
        if 'has_surfDiff' in settings:
            at.selectbox(key=f"has_surfDiff_partype_{jj}").set_value(settings['has_surfDiff']).run()


@pytest.mark.ci
def test_psd_uniform_settings():
    """PSD with uniform per-partype settings should produce consistent equations."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_psd(at, [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    assert r"k^{\mathrm{f}}" in latex


@pytest.mark.ci
def test_psd_different_film_diffusion():
    """PSD with different film diffusion settings per particle type."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_psd(at, [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No'},
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    # Should have particle type labels since types differ
    assert "particle type" in latex.lower() or r"j = " in latex


@pytest.mark.ci
def test_psd_different_surface_diffusion():
    """PSD with different surface diffusion settings per particle type."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_psd(at, [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex


@pytest.mark.ci
def test_psd_all_nonlimiting():
    """PSD where all particle types have non-limiting film diffusion."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_psd(at, [
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    # Non-limiting means no film diffusion coefficient
    assert r"k^{\mathrm{f}}" not in latex


@pytest.mark.ci
def test_single_particle_unchanged():
    """Single particle type should still use global settings."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="PSD").set_value("No").run()
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    at.selectbox(key="nonlimiting_filmDiff").set_value("No").run()
    at.selectbox(key="has_binding").set_value("Yes").run()
    at.selectbox(key="has_surfDiff").set_value("Yes").run()

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    assert r"k^{\mathrm{f}}" in latex
