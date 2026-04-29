# -*- coding: utf-8 -*-
"""Tests for per-particle-type parameterization of film diffusion and surface diffusion."""

from streamlit.testing.v1 import AppTest
import pytest


def setup_grm_ptd(at, per_partype_settings):
    """Set up a GRM model with PTD (dev mode) and per-particle-type configuration.

    Parameters
    ----------
    at : AppTest
        The Streamlit app test instance.
    per_partype_settings : list of dict
        Per-particle-type settings, each dict with keys:
        'nonlimiting_filmDiff', 'has_surfDiff'
    """
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="dev_mode").set_value("On").run()
    n_p = len(per_partype_settings)
    at.number_input(key=r"N^\mathrm{p}").set_value(n_p).run()

    for jj in range(n_p):
        at.selectbox(key=f"parType_{jj+1}_resolution").set_value("1D (radial coordinate)").run()
        at.selectbox(key=f"parType_{jj+1}_nonlimiting_filmDiff").set_value(per_partype_settings[jj]['nonlimiting_filmDiff']).run()

    at.selectbox(key="has_binding").set_value("Yes").run()

    for jj, settings in enumerate(per_partype_settings):
        if 'has_surfDiff' in settings:
            at.selectbox(key=f"parType_{jj+1}_has_surfDiff").set_value(settings['has_surfDiff']).run()


@pytest.mark.ci
def test_ptd_uniform_settings():
    """PTD with uniform per-partype settings should produce consistent equations."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_ptd(at, [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    assert r"k^{\mathrm{f}}" in latex


@pytest.mark.ci
def test_ptd_different_film_diffusion():
    """PTD with different film diffusion settings per particle type."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_ptd(at, [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No'},
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    # Should have particle type labels since types differ
    assert "particle type" in latex.lower() or r"j = " in latex


@pytest.mark.ci
def test_ptd_different_surface_diffusion():
    """PTD with different surface diffusion settings per particle type."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_ptd(at, [
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes'},
        {'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex


@pytest.mark.ci
def test_ptd_all_nonlimiting():
    """PTD where all particle types have non-limiting film diffusion."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    setup_grm_ptd(at, [
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
        {'nonlimiting_filmDiff': 'Yes', 'has_surfDiff': 'No'},
    ])

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    # Non-limiting means no film diffusion coefficient
    assert r"k^{\mathrm{f}}" not in latex


@pytest.mark.ci
def test_psd_shared_config():
    """PSD in advanced mode should use shared config for all particle types."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="PSD").set_value("Particle size distribution").run()
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    at.selectbox(key="particle_nonlimiting_filmDiff").set_value("No").run()
    at.selectbox(key="has_binding").set_value("Yes").run()
    at.selectbox(key="particle_has_surfDiff").set_value("Yes").run()

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    assert r"k^{\mathrm{f}}" in latex


@pytest.mark.ci
def test_single_particle_unchanged():
    """Single particle type should still use global settings."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="PSD").set_value("Yes").run()
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    at.selectbox(key="particle_nonlimiting_filmDiff").set_value("No").run()
    at.selectbox(key="has_binding").set_value("Yes").run()
    at.selectbox(key="particle_has_surfDiff").set_value("Yes").run()

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    assert r"k^{\mathrm{f}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_format_partype_set_single():
    """format_partype_set with a single index."""
    from src.model_column import Column
    assert Column.format_partype_set([1]) == r"$j = 1$"


@pytest.mark.ci
@pytest.mark.unit_test
def test_format_partype_set_multiple():
    """format_partype_set with multiple indices."""
    from src.model_column import Column
    assert Column.format_partype_set([1, 3]) == r"$j \in \{1, 3\}$"
