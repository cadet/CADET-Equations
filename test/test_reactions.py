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


# === Rapid-equilibrium reaction unit tests ===

@pytest.mark.ci
def test_req_reaction_equation_functions():
    """Unit tests for rapid-equilibrium reaction equation snippet functions."""
    from src import equations as eq

    # Bulk constraint
    result = eq.req_reaction_bulk_constraint()
    assert r"g^{\mathrm{react,eq},\b}" in result

    # Particle liquid constraint - single and multi particle
    result_single = eq.req_reaction_particle_liquid_constraint(singleParticle=True)
    assert r"g^{\mathrm{react,eq},\p}_{k}" in result_single
    result_multi = eq.req_reaction_particle_liquid_constraint(singleParticle=False)
    assert r"g^{\mathrm{react,eq},\p}_{j,k}" in result_multi

    # Particle solid constraint - single and multi particle
    result_single = eq.req_reaction_particle_solid_constraint(singleParticle=True)
    assert r"g^{\mathrm{react,eq},\s}_{k}" in result_single
    result_multi = eq.req_reaction_particle_solid_constraint(singleParticle=False)
    assert r"g^{\mathrm{react,eq},\s}_{j,k}" in result_multi

    # Conserved moiety equations
    result = eq.conserved_moiety_equation_bulk()
    assert r"M^{\b}" in result

    result_single = eq.conserved_moiety_equation_particle_liquid(singleParticle=True)
    assert r"M^{\p}_{l,i}" in result_single
    result_multi = eq.conserved_moiety_equation_particle_liquid(singleParticle=False)
    assert r"M^{\p}_{j,l,i}" in result_multi

    result_single = eq.conserved_moiety_equation_particle_solid(singleParticle=True)
    assert r"M^{\s}_{l,i}" in result_single
    result_multi = eq.conserved_moiety_equation_particle_solid(singleParticle=False)
    assert r"M^{\s}_{j,l,i}" in result_multi


# === Rapid-equilibrium reaction integration tests ===

@pytest.mark.ci
def test_req_reaction_bulk():
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    # Should have equilibrium constraint, not kinetic reaction term
    assert r"g^{\mathrm{react,eq}" in latex
    assert r"M^{\mathrm{b}}" in latex
    assert r"f^{\mathrm{react},\mathrm{b}}" not in latex


@pytest.mark.ci
def test_req_reaction_bulk_no_particles():
    at = setup_app_with_advanced_mode(add_particles="No")
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{b}}" in latex
    assert r"M^{\mathrm{b}}" in latex


@pytest.mark.ci
def test_req_reaction_particle_liquid_1D():
    at = setup_app_with_advanced_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{p}}" in latex
    assert r"M^{\mathrm{p}}" in latex
    assert r"f^{\mathrm{react},\mathrm{p}}" not in latex


@pytest.mark.ci
def test_req_reaction_particle_solid_1D():
    at = setup_app_with_advanced_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_solid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{s}}" in latex
    assert r"M^{\mathrm{s}}" in latex
    assert r"f^{\mathrm{react},\mathrm{s}}" not in latex


@pytest.mark.ci
def test_req_reaction_particle_liquid_0D():
    at = setup_app_with_advanced_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{p}}" in latex
    assert r"M^{\mathrm{p}}" in latex


@pytest.mark.ci
def test_req_reaction_particle_solid_0D():
    at = setup_app_with_advanced_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_solid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{s}}" in latex
    assert r"M^{\mathrm{s}}" in latex


@pytest.mark.ci
def test_req_reaction_all():
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_solid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{b}}" in latex
    assert r"g^{\mathrm{react,eq},\mathrm{p}}" in latex
    assert r"g^{\mathrm{react,eq},\mathrm{s}}" in latex
    assert r"f^{\mathrm{react}" not in latex


@pytest.mark.ci
def test_req_reaction_kinetics_mode_not_shown_without_reaction():
    at = setup_app_with_advanced_mode()
    assert not at.exception
    selectbox_keys = [box.key for box in at.selectbox]
    assert "req_reaction_bulk" not in selectbox_keys
    assert "req_reaction_particle_liquid" not in selectbox_keys
    assert "req_reaction_particle_solid" not in selectbox_keys


@pytest.mark.ci
def test_req_reaction_assumptions():
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    at.toggle(key="model_assumptions").set_value(True).run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert "conserved moieties" in latex


@pytest.mark.ci
def test_mixed_kinetic_and_req_reactions():
    """Test that kinetic and rapid-equilibrium reactions can coexist in different phases."""
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    # bulk stays kinetic (default)
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    # Bulk should have kinetic reaction term
    assert r"f^{\mathrm{react}" in latex
    # Particle liquid should have equilibrium constraint
    assert r"g^{\mathrm{react,eq},\mathrm{p}}" in latex


@pytest.mark.ci
def test_req_reaction_bulk_0D_tank():
    """Test rapid-equilibrium bulk reaction in 0D tank model."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="column_resolution").set_value("0D (Homogeneous Tank)").run()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{b}}" in latex
    assert r"M^{\mathrm{b}}" in latex
    assert r"f^{\mathrm{react},\mathrm{b}}" not in latex


@pytest.mark.ci
def test_req_reaction_with_psd():
    """Test rapid-equilibrium reactions with particle size distribution (N_p > 1)."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="PSD").set_value("Yes").run()
    at.selectbox(key="has_binding").set_value("Yes").run()
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_solid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{p}}" in latex
    assert r"g^{\mathrm{react,eq},\mathrm{s}}" in latex
    assert r"M^{\mathrm{p}}" in latex
    assert r"M^{\mathrm{s}}" in latex


@pytest.mark.ci
def test_req_reaction_assumptions_multiple_phases():
    """Test assumptions text lists all phases with rapid-equilibrium reactions."""
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    at.toggle(key="model_assumptions").set_value(True).run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert "bulk liquid" in latex
    assert "particle liquid" in latex
    assert "conserved moieties" in latex


@pytest.mark.ci
def test_req_reaction_symbol_table():
    """Test that symbol table contains rapid-equilibrium reaction symbols."""
    at = setup_app_with_advanced_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    at.toggle(key="sym_table").set_value(True).run()
    assert not at.exception
