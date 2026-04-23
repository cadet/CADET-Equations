# -*- coding: utf-8 -*-
"""
Integration tests for per-component parameterization (issue #7).
Tests the full Streamlit UI path with per-component configuration.
"""

from streamlit.testing.v1 import AppTest
import pytest


def setup_grm_with_per_component(at, n_comp, per_comp_settings):
    """Set up a GRM model with per-component configuration.

    Parameters
    ----------
    at : AppTest
        The Streamlit app test instance.
    n_comp : int
        Number of components.
    per_comp_settings : list of dict
        Per-component settings, each dict with keys:
        'req_binding', 'nonlimiting_filmDiff', 'has_surfDiff', 'has_mult_bnd_states'
    """
    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="N_c_choice").set_value(n_comp).run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    at.selectbox(key="has_binding").set_value("Yes").run()

    for i, settings in enumerate(per_comp_settings):
        at.selectbox(key=f"req_binding_comp_{i}").set_value(settings['req_binding']).run()
        at.selectbox(key=f"nonlimiting_filmDiff_comp_{i}").set_value(settings['nonlimiting_filmDiff']).run()
        at.selectbox(key=f"has_surfDiff_comp_{i}").set_value(settings['has_surfDiff']).run()
        at.selectbox(key=f"has_mult_bnd_states_comp_{i}").set_value(settings['has_mult_bnd_states']).run()


@pytest.mark.ci
@pytest.mark.reference
def test_per_component_uniform_settings():
    """Test that uniform per-component settings produce a single equation block."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'req_binding': 'Kinetic', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes', 'has_mult_bnd_states': 'No'},
        {'req_binding': 'Kinetic', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes', 'has_mult_bnd_states': 'No'},
    ]
    setup_grm_with_per_component(at, 2, settings)

    assert not at.exception
    latex = at.session_state.latex_string
    assert r"\begin{align}" in latex
    # Uniform settings -> should NOT have "For component(s)" text
    assert "For component(s)" not in latex


@pytest.mark.ci
@pytest.mark.reference
def test_per_component_different_binding_kinetics():
    """Test per-component with mixed kinetic / rapid-equilibrium binding."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'req_binding': 'Kinetic', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No', 'has_mult_bnd_states': 'No'},
        {'req_binding': 'Rapid-equilibrium', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No', 'has_mult_bnd_states': 'No'},
    ]
    setup_grm_with_per_component(at, 2, settings)

    assert not at.exception
    latex = at.session_state.latex_string
    # Different settings -> should have per-component text
    assert "For component(s)" in latex
    assert r"i = 1" in latex
    assert r"i = 2" in latex


@pytest.mark.ci
@pytest.mark.reference
def test_per_component_three_components_two_groups():
    """Test 3 components grouped into 2 groups by shared settings."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'req_binding': 'Kinetic', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No', 'has_mult_bnd_states': 'No'},
        {'req_binding': 'Rapid-equilibrium', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No', 'has_mult_bnd_states': 'No'},
        {'req_binding': 'Kinetic', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No', 'has_mult_bnd_states': 'No'},
    ]
    setup_grm_with_per_component(at, 3, settings)

    assert not at.exception
    latex = at.session_state.latex_string
    assert "For component(s)" in latex
    # Components 1 and 3 share kinetic binding
    assert r"i \in \{1, 3\}" in latex
    # Component 2 has rapid-equilibrium
    assert r"i = 2" in latex


@pytest.mark.ci
@pytest.mark.reference
def test_per_component_different_surface_diffusion():
    """Test per-component with mixed surface diffusion settings."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    settings = [
        {'req_binding': 'Kinetic', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'Yes', 'has_mult_bnd_states': 'No'},
        {'req_binding': 'Kinetic', 'nonlimiting_filmDiff': 'No', 'has_surfDiff': 'No', 'has_mult_bnd_states': 'No'},
    ]
    setup_grm_with_per_component(at, 2, settings)

    assert not at.exception
    latex = at.session_state.latex_string
    assert "For component(s)" in latex


@pytest.mark.ci
@pytest.mark.reference
def test_per_component_0d_particle_resolution():
    """Test per-component with 0D (homogeneous) particle resolution."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="N_c_choice").set_value(2).run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="particle_resolution").set_value("0D (homogeneous)").run()
    at.selectbox(key="has_binding").set_value("Yes").run()

    at.selectbox(key=f"req_binding_comp_0").set_value("Kinetic").run()
    at.selectbox(key=f"nonlimiting_filmDiff_comp_0").set_value("No").run()
    at.selectbox(key=f"has_mult_bnd_states_comp_0").set_value("No").run()

    at.selectbox(key=f"req_binding_comp_1").set_value("Rapid-equilibrium").run()
    at.selectbox(key=f"nonlimiting_filmDiff_comp_1").set_value("No").run()
    at.selectbox(key=f"has_mult_bnd_states_comp_1").set_value("No").run()

    assert not at.exception
    latex = at.session_state.latex_string
    assert "For component(s)" in latex


@pytest.mark.ci
@pytest.mark.reference
def test_arbitrary_n_c_no_per_component():
    """Test that N_c_choice = Arbitrary does not enable per-component config."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()
    at.selectbox(key="N_c_choice").set_value("Arbitrary").run()
    at.selectbox(key="add_particles").set_value("Yes").run()
    at.selectbox(key="has_binding").set_value("Yes").run()

    assert not at.exception
    latex = at.session_state.latex_string
    # Should not have per-component text
    assert "For component(s)" not in latex
