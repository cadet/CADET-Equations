# -*- coding: utf-8 -*-
"""
Tests for dynamic reaction terms
"""

from streamlit.testing.v1 import AppTest
import pytest


def setup_app_with_advanced_mode(add_particles="Yes", particle_resolution="1D (radial coordinate)", has_binding="Yes"):
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    assert not at.exception

    if add_particles == "Yes":
        at.selectbox(key="add_particles").set_value("Yes").run()
        at.selectbox(key="particle_resolution").set_value(particle_resolution).run()
        at.selectbox(key="has_binding").set_value(has_binding).run()

    return at


@pytest.mark.ci
def test_bulk_reaction_no_particles():
    at = setup_app_with_advanced_mode(add_particles="No")
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex


@pytest.mark.ci
def test_bulk_reaction_with_particles():
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex


@pytest.mark.ci
def test_particle_liquid_reaction_1D():
    at = setup_app_with_advanced_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{p}}" in latex


@pytest.mark.ci
def test_particle_solid_reaction_1D():
    at = setup_app_with_advanced_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{s}}" in latex


@pytest.mark.ci
def test_particle_liquid_reaction_0D():
    at = setup_app_with_advanced_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{p}}" in latex


@pytest.mark.ci
def test_particle_solid_reaction_0D():
    at = setup_app_with_advanced_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{s}}" in latex


@pytest.mark.ci
def test_all_reactions_enabled():
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{b}}" in latex
    assert r"f^{\mathrm{react},\mathrm{p}}" in latex
    assert r"f^{\mathrm{react},\mathrm{s}}" in latex


@pytest.mark.ci
def test_no_reactions_by_default():
    at = setup_app_with_advanced_mode()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" not in latex


@pytest.mark.ci
def test_reactions_only_in_advanced_mode():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    selectbox_keys = [box.key for box in at.selectbox]
    assert "has_reaction_bulk" not in selectbox_keys


@pytest.mark.ci
def test_bulk_reaction_0D_tank():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="column_resolution").set_value("0D (Homogeneous Tank)").run()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex
