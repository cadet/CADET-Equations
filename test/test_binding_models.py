# -*- coding: utf-8 -*-
"""
Integration tests for explicit binding models (issue #15).
"""

from streamlit.testing.v1 import AppTest
import pytest


def setup_app_with_binding(binding_model="Arbitrary", particle_resolution="1D (radial coordinate)", req_binding="Kinetic"):
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    assert not at.exception

    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="particle_resolution").set_value(particle_resolution).run()
    at.selectbox(key="has_binding").set_value("Yes").run()
    assert not at.exception

    at.selectbox(key="req_binding").set_value(req_binding).run()
    at.selectbox(key="binding_model").set_value(binding_model).run()
    assert not at.exception

    return at


@pytest.mark.ci
def test_arbitrary_binding_default():
    """Default binding model should produce f^bind."""
    at = setup_app_with_binding("Arbitrary")
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{bind}}" in latex


@pytest.mark.ci
def test_linear_binding_1D():
    """Linear binding model should produce k_a and k_d terms."""
    at = setup_app_with_binding("Linear")
    latex = at.session_state.latex_string
    assert r"k^{\mathrm{a}}" in latex
    assert r"k^{\mathrm{d}}" in latex
    assert r"f^{\mathrm{bind}}" not in latex


@pytest.mark.ci
def test_langmuir_binding_1D():
    """Langmuir binding model should produce q_max and summation."""
    at = setup_app_with_binding("Langmuir")
    latex = at.session_state.latex_string
    assert r"q^{\mathrm{max}}" in latex
    assert r"k^{\mathrm{a}}" in latex
    assert r"f^{\mathrm{bind}}" not in latex


@pytest.mark.ci
def test_sma_binding_1D():
    """SMA binding model should produce nu, sigma, Lambda, and auxiliary equations."""
    at = setup_app_with_binding("SMA")
    latex = at.session_state.latex_string
    assert r"\nu" in latex
    assert r"\Lambda" in latex
    assert r"\bar{q}" in latex
    assert r"q^{\mathrm{ref}}" in latex
    assert r"c^{\mathrm{ref}}" in latex
    assert r"f^{\mathrm{bind}}" not in latex


@pytest.mark.ci
def test_linear_binding_0D():
    """Linear binding should work with 0D (homogeneous) particles."""
    at = setup_app_with_binding("Linear", particle_resolution="0D (homogeneous)")
    latex = at.session_state.latex_string
    assert r"k^{\mathrm{a}}" in latex
    assert r"k^{\mathrm{d}}" in latex


@pytest.mark.ci
def test_langmuir_binding_0D():
    """Langmuir binding should work with 0D (homogeneous) particles."""
    at = setup_app_with_binding("Langmuir", particle_resolution="0D (homogeneous)")
    latex = at.session_state.latex_string
    assert r"q^{\mathrm{max}}" in latex


@pytest.mark.ci
def test_sma_binding_rapid_equilibrium():
    """SMA binding with rapid equilibrium should still produce SMA terms."""
    at = setup_app_with_binding("SMA", req_binding="Rapid-equilibrium")
    latex = at.session_state.latex_string
    assert r"\bar{q}" in latex
    assert r"\Lambda" in latex


@pytest.mark.ci
def test_binding_model_not_shown_without_advanced_mode():
    """Binding model selector should not appear without advanced mode."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    selectbox_keys = [box.key for box in at.selectbox]
    assert "binding_model" not in selectbox_keys


@pytest.mark.ci
def test_binding_model_not_shown_without_binding():
    """Binding model selector should not appear when binding is disabled."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="has_binding").set_value("No").run()
    assert not at.exception
    selectbox_keys = [box.key for box in at.selectbox]
    assert "binding_model" not in selectbox_keys
