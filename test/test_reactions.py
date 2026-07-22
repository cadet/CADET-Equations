"""
Tests for dynamic reaction terms
"""

import pytest
from streamlit.testing.v1 import AppTest


def setup_app_with_dev_mode(add_particles="Yes", particle_resolution="1D (radial coordinate)", has_binding="Yes"):
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    at.button(key="dev_mode_button").click().run()
    assert not at.exception

    if add_particles == "Yes":
        # dev_mode uses number_input for N_p
        at.number_input(key=r"N^\mathrm{p}").set_value(1).run()
        at.selectbox(key="particle_resolution").set_value(particle_resolution).run()
        at.selectbox(key="has_binding").set_value(has_binding).run()

    return at


@pytest.mark.ci
@pytest.mark.unit_test
def test_bulk_reaction_no_particles():
    at = setup_app_with_dev_mode(add_particles="No")
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_bulk_reaction_with_particles():
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_liquid_reaction_1D():
    at = setup_app_with_dev_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{p}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_solid_reaction_1D():
    at = setup_app_with_dev_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{s}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_liquid_reaction_0D():
    at = setup_app_with_dev_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{p}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_solid_reaction_0D():
    at = setup_app_with_dev_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{s}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_all_reactions_enabled():
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react},\mathrm{b}}" in latex
    assert r"f^{\mathrm{react},\mathrm{p}}" in latex
    assert r"f^{\mathrm{react},\mathrm{s}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_no_reactions_by_default():
    at = setup_app_with_dev_mode()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" not in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reactions_only_in_dev_mode():
    """Rapid-equilibrium and particle reaction selectboxes should not appear without dev_mode."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    selectbox_keys = [box.key for box in at.selectbox]
    assert "req_reaction_bulk" not in selectbox_keys
    assert "has_reaction_particle_liquid" not in selectbox_keys


@pytest.mark.ci
@pytest.mark.unit_test
def test_bulk_reaction_0D_tank():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.button(key="dev_mode_button").click().run()
    at.selectbox(key="column_type").set_value("Mixed tank").run()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex


# === Rapid-equilibrium reaction unit tests ===


@pytest.mark.ci
@pytest.mark.unit_test
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
@pytest.mark.unit_test
def test_req_reaction_bulk():
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    # Should have equilibrium constraint, not kinetic reaction term
    assert r"g^{\mathrm{react,eq}" in latex
    assert r"M^{\mathrm{b}}" in latex
    assert r"f^{\mathrm{react},\mathrm{b}}" not in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_bulk_no_particles():
    at = setup_app_with_dev_mode(add_particles="No")
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{b}}" in latex
    assert r"M^{\mathrm{b}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_particle_liquid_1D():
    at = setup_app_with_dev_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{p}}" in latex
    assert r"M^{\mathrm{p}}" in latex
    assert r"f^{\mathrm{react},\mathrm{p}}" not in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_particle_solid_1D():
    at = setup_app_with_dev_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_solid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{s}}" in latex
    assert r"M^{\mathrm{s}}" in latex
    assert r"f^{\mathrm{react},\mathrm{s}}" not in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_particle_liquid_0D():
    at = setup_app_with_dev_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{p}}" in latex
    assert r"M^{\mathrm{p}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_particle_solid_0D():
    at = setup_app_with_dev_mode(particle_resolution="0D (homogeneous)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_solid").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{s}}" in latex
    assert r"M^{\mathrm{s}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_all():
    at = setup_app_with_dev_mode()
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
@pytest.mark.unit_test
def test_req_reaction_kinetics_mode_not_shown_without_reaction():
    at = setup_app_with_dev_mode()
    assert not at.exception
    selectbox_keys = [box.key for box in at.selectbox]
    assert "req_reaction_bulk" not in selectbox_keys
    assert "req_reaction_particle_liquid" not in selectbox_keys
    assert "req_reaction_particle_solid" not in selectbox_keys


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_assumptions():
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    at.toggle(key="model_assumptions").set_value(True).run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert "conserved moieties" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_mixed_kinetic_and_req_reactions():
    """Test that kinetic and rapid-equilibrium reactions can coexist in different phases."""
    at = setup_app_with_dev_mode()
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
@pytest.mark.unit_test
def test_req_reaction_bulk_0D_tank():
    """Test rapid-equilibrium bulk reaction in 0D tank model."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.button(key="dev_mode_button").click().run()
    at.selectbox(key="column_type").set_value("Mixed tank").run()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"g^{\mathrm{react,eq},\mathrm{b}}" in latex
    assert r"M^{\mathrm{b}}" in latex
    assert r"f^{\mathrm{react},\mathrm{b}}" not in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_req_reaction_with_psd():
    """Test rapid-equilibrium reactions with particle size distribution (N_p > 1)."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.button(key="dev_mode_button").click().run()
    # dev_mode uses number_input for N_p
    at.number_input(key=r"N^\mathrm{p}").set_value(2).run()
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
@pytest.mark.unit_test
def test_req_reaction_assumptions_multiple_phases():
    """Test assumptions text lists all phases with rapid-equilibrium reactions."""
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="req_reaction_bulk").set_value("Rapid-equilibrium").run()
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_liquid").set_value("Rapid-equilibrium").run()
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="req_reaction_particle_solid").set_value("Rapid-equilibrium").run()
    at.toggle(key="model_assumptions").set_value(True).run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert "bulk liquid" in latex
    assert "particle liquid" in latex
    assert "particle solid" in latex
    assert "conserved moieties" in latex


# === Explicit reaction model unit tests ===


@pytest.mark.ci
@pytest.mark.unit_test
def test_mass_action_law_equation_functions():
    """Unit tests for Mass Action Law equation snippet functions."""
    from src import equations as eq

    for phase in ["bulk", "particle_liquid", "particle_solid"]:
        result = eq.mass_action_law_definition(phase)
        assert result is not None
        main_eq, flux_eq = result
        assert r"S^{" in main_eq
        assert r"\varphi^{" in main_eq
        assert r"k^{\mathrm{fwd}" in flux_eq
        assert r"k^{\mathrm{bwd}" in flux_eq
        assert r"e^{\mathrm{fwd}" in flux_eq
        assert r"e^{\mathrm{bwd}" in flux_eq


@pytest.mark.ci
@pytest.mark.unit_test
def test_michaelis_menten_equation_functions():
    """Unit tests for Michaelis-Menten equation snippet functions."""
    from src import equations as eq

    for phase in ["bulk", "particle_liquid", "particle_solid"]:
        result = eq.michaelis_menten_definition(phase)
        assert result is not None
        main_eq, flux_eq = result
        assert r"S^{" in main_eq
        assert r"\nu^{" in main_eq
        assert r"v^{\mathrm{max}" in flux_eq
        assert r"K^{\mathrm{M}" in flux_eq


@pytest.mark.ci
@pytest.mark.unit_test
def test_reaction_model_definition_dispatcher():
    """Test that reaction_model_definition dispatches correctly."""
    from src import equations as eq

    assert eq.reaction_model_definition("Mass Action Law", "bulk") is not None
    assert eq.reaction_model_definition("Michaelis Menten", "bulk") is not None
    assert eq.reaction_model_definition("Arbitrary", "bulk") is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_reaction_model_assumptions():
    """Test reaction model assumption strings."""
    from src import equations as eq

    assert eq.reaction_model_assumptions("Arbitrary") == ["No custom assumptions"]

    mal_asmpts = eq.reaction_model_assumptions("Mass Action Law")
    assert any("mass action" in a for a in mal_asmpts)

    mm_asmpts = eq.reaction_model_assumptions("Michaelis Menten")
    assert any("Michaelis-Menten" in a for a in mm_asmpts)


# === Explicit reaction model integration tests ===


@pytest.mark.ci
@pytest.mark.unit_test
def test_mass_action_law_bulk():
    """Test Mass Action Law reaction model with bulk reactions."""
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Mass Action Law").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex
    assert r"S^{\mathrm{b}}" in latex
    assert r"\varphi^{\mathrm{b}}" in latex
    assert r"k^{\mathrm{fwd}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_michaelis_menten_bulk():
    """Test Michaelis-Menten reaction model with bulk reactions."""
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Michaelis Menten").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex
    assert r"S^{\mathrm{b}}" in latex
    assert r"\nu^{\mathrm{b}}" in latex
    assert r"v^{\mathrm{max}" in latex
    assert r"K^{\mathrm{M}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_mass_action_law_particle_liquid():
    """Test Mass Action Law with particle liquid reactions."""
    at = setup_app_with_dev_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Mass Action Law").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"S^{\mathrm{p}}" in latex
    assert r"\varphi^{\mathrm{p}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_michaelis_menten_particle_solid():
    """Test Michaelis-Menten with particle solid reactions."""
    at = setup_app_with_dev_mode(particle_resolution="1D (radial coordinate)")
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Michaelis Menten").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"S^{\mathrm{s}}" in latex
    assert r"\nu^{\mathrm{s}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_mass_action_law_all_phases():
    """Test Mass Action Law with all reaction phases enabled."""
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="has_reaction_particle_liquid").set_value("Yes").run()
    at.selectbox(key="has_reaction_particle_solid").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Mass Action Law").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"S^{\mathrm{b}}" in latex
    assert r"S^{\mathrm{p}}" in latex
    assert r"S^{\mathrm{s}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_mass_action_law_assumptions():
    """Test that Mass Action Law assumptions appear when enabled."""
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Mass Action Law").run()
    at.toggle(key="model_assumptions").set_value(True).run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert "mass action" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_arbitrary_reaction_no_definition_section():
    """Test that Arbitrary model does not show a reaction definition section."""
    at = setup_app_with_dev_mode()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"f^{\mathrm{react}" in latex
    assert r"\varphi^{" not in latex
    assert r"\nu^{" not in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_mass_action_law_no_particles():
    """Test Mass Action Law with only bulk reaction and no particles."""
    at = setup_app_with_dev_mode(add_particles="No")
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Mass Action Law").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"S^{\mathrm{b}}" in latex
    assert r"\varphi^{\mathrm{b}}" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_vars_mass_action_law_liquid():
    """Direct Particle unit test: Mass Action Law vars for particle liquid reaction."""
    from src.model_particle import Particle

    p = Particle(
        geometry="Sphere",
        has_core=False,
        var_format="CADET",
        resolution="1D",
        has_binding=True,
        req_binding=False,
        has_mult_bnd_states=False,
        has_surfDiff=True,
        nonlimiting_filmDiff=False,
        interstitial_volume_resolution="1D",
        column_type="Axial",
        single_partype=True,
        PTD=False,
        binding_model="Linear",
        has_reaction_liquid=True,
        has_reaction_solid=False,
        req_reaction_liquid=False,
        req_reaction_solid=False,
        reaction_model="Mass Action Law",
    )
    symbols = [v["Symbol"] for v in p.vars_and_params]
    assert any(r"S^{\mathrm{p}}" in s for s in symbols)
    assert any(r"k^{\mathrm{fwd},\mathrm{p}}" in s for s in symbols)
    assert any(r"k^{\mathrm{bwd},\mathrm{p}}" in s for s in symbols)
    assert any(r"e^{\mathrm{fwd},\mathrm{p}}" in s for s in symbols)
    assert any(r"e^{\mathrm{bwd},\mathrm{p}}" in s for s in symbols)


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_vars_michaelis_menten_liquid():
    """Direct Particle unit test: Michaelis-Menten vars for particle liquid reaction."""
    from src.model_particle import Particle

    p = Particle(
        geometry="Sphere",
        has_core=False,
        var_format="CADET",
        resolution="1D",
        has_binding=True,
        req_binding=False,
        has_mult_bnd_states=False,
        has_surfDiff=True,
        nonlimiting_filmDiff=False,
        interstitial_volume_resolution="1D",
        column_type="Axial",
        single_partype=True,
        PTD=False,
        binding_model="Linear",
        has_reaction_liquid=True,
        has_reaction_solid=False,
        req_reaction_liquid=False,
        req_reaction_solid=False,
        reaction_model="Michaelis Menten",
    )
    symbols = [v["Symbol"] for v in p.vars_and_params]
    assert any(r"S^{\mathrm{p}}" in s for s in symbols)
    assert any(r"v^{\mathrm{max},\mathrm{p}}" in s for s in symbols)
    assert any(r"K^{\mathrm{M},\mathrm{p}}" in s for s in symbols)
    assert any(r"N^{\mathrm{sub},\mathrm{p}}" in s for s in symbols)


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_vars_mass_action_law_solid():
    """Direct Particle unit test: Mass Action Law vars for particle solid reaction."""
    from src.model_particle import Particle

    p = Particle(
        geometry="Sphere",
        has_core=False,
        var_format="CADET",
        resolution="1D",
        has_binding=True,
        req_binding=False,
        has_mult_bnd_states=False,
        has_surfDiff=True,
        nonlimiting_filmDiff=False,
        interstitial_volume_resolution="1D",
        column_type="Axial",
        single_partype=True,
        PTD=False,
        binding_model="Linear",
        has_reaction_liquid=False,
        has_reaction_solid=True,
        req_reaction_liquid=False,
        req_reaction_solid=False,
        reaction_model="Mass Action Law",
    )
    symbols = [v["Symbol"] for v in p.vars_and_params]
    assert any(r"S^{\mathrm{s}}" in s for s in symbols)
    assert any(r"k^{\mathrm{fwd},\mathrm{s}}" in s for s in symbols)
    assert any(r"k^{\mathrm{bwd},\mathrm{s}}" in s for s in symbols)


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_vars_michaelis_menten_solid():
    """Direct Particle unit test: Michaelis-Menten vars for particle solid reaction."""
    from src.model_particle import Particle

    p = Particle(
        geometry="Sphere",
        has_core=False,
        var_format="CADET",
        resolution="1D",
        has_binding=True,
        req_binding=False,
        has_mult_bnd_states=False,
        has_surfDiff=True,
        nonlimiting_filmDiff=False,
        interstitial_volume_resolution="1D",
        column_type="Axial",
        single_partype=True,
        PTD=False,
        binding_model="Linear",
        has_reaction_liquid=False,
        has_reaction_solid=True,
        req_reaction_liquid=False,
        req_reaction_solid=False,
        reaction_model="Michaelis Menten",
    )
    symbols = [v["Symbol"] for v in p.vars_and_params]
    assert any(r"S^{\mathrm{s}}" in s for s in symbols)
    assert any(r"v^{\mathrm{max},\mathrm{s}}" in s for s in symbols)
    assert any(r"K^{\mathrm{M},\mathrm{s}}" in s for s in symbols)
    assert any(r"N^{\mathrm{sub},\mathrm{s}}" in s for s in symbols)


@pytest.mark.ci
@pytest.mark.unit_test
def test_michaelis_menten_0D_tank():
    """Test Michaelis-Menten in 0D tank model with bulk reaction."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.button(key="dev_mode_button").click().run()
    at.selectbox(key="column_type").set_value("Mixed tank").run()
    at.selectbox(key="has_reaction_bulk").set_value("Yes").run()
    at.selectbox(key="reaction_model").set_value("Michaelis Menten").run()
    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\nu^{\mathrm{b}}" in latex
    assert r"K^{\mathrm{M}" in latex
