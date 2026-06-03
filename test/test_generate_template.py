# -*- coding: utf-8 -*-
"""
Tests for src/generate_template.py.
"""

from streamlit.testing.v1 import AppTest
import pytest


def setup_app(**overrides):
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    for key, value in overrides.items():
        at.selectbox(key=key).set_value(value).run()
        assert not at.exception
    return at


def get_script(at):
    return at.session_state.template_script


@pytest.mark.ci
class TestGenerateTemplateIntegration:

    def test_template_in_session_state(self):
        at = setup_app()
        assert "template_script" in at.session_state

    def test_default_grm_template(self):
        at = setup_app(
            add_particles="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
        )
        script = get_script(at)
        assert "def get_unit_operation(" in script
        assert "'COLUMN_MODEL_1D'" in script
        assert "'has_pore_diffusion'" in script
        assert "'has_film_diffusion'" in script
        assert "return unit" in script
        _assert_valid_python(script)

    def test_no_particles(self):
        at = setup_app(add_particles="No")
        script = get_script(at)
        assert "'npartype'" not in script
        assert "'par_geom'" not in script
        assert "'col_porosity'] = 1.0" in script
        _assert_valid_python(script)

    def test_lrmp_template(self):
        at = setup_app(
            add_particles="Yes",
            has_binding="Yes",
            particle_resolution="0D (homogeneous)",
            particle_nonlimiting_filmDiff="No",
        )
        script = get_script(at)
        assert "'has_film_diffusion'] = 1" in script
        assert "'has_pore_diffusion'] = 0" in script
        _assert_valid_python(script)

    def test_lrm_without_pores(self):
        at = setup_app(
            add_particles="Yes",
            has_binding="Yes",
            particle_resolution="0D (homogeneous)",
            particle_nonlimiting_filmDiff="Yes",
        )
        script = get_script(at)
        assert "'has_film_diffusion'] = 0" in script
        assert "'total_porosity'" in script
        _assert_valid_python(script)

    def test_linear_binding(self):
        at = setup_app(
            advanced_mode="On",
            PSD="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
            binding_model="Linear",
        )
        script = get_script(at)
        assert "'LINEAR'" in script
        assert "'lin_ka'" in script
        assert "'lin_kd'" in script
        _assert_valid_python(script)

    def test_langmuir_binding(self):
        at = setup_app(
            advanced_mode="On",
            PSD="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
            binding_model="Langmuir",
        )
        script = get_script(at)
        assert "'MULTI_COMPONENT_LANGMUIR'" in script
        assert "'mcl_ka'" in script
        assert "'mcl_qmax'" in script
        _assert_valid_python(script)

    def test_sma_binding(self):
        at = setup_app(
            advanced_mode="On",
            PSD="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
            binding_model="SMA",
        )
        script = get_script(at)
        assert "'STERIC_MASS_ACTION'" in script
        assert "'sma_ka'" in script
        assert "'sma_lambda'" in script
        assert "'sma_nu'" in script
        assert "ncomp=4" in script
        _assert_valid_python(script)

    def test_no_binding(self):
        at = setup_app(
            add_particles="Yes",
            has_binding="No",
        )
        script = get_script(at)
        assert "'adsorption_model'" not in script
        _assert_valid_python(script)

    def test_no_axial_dispersion(self):
        at = setup_app(
            add_particles="No",
            has_axial_dispersion="No",
        )
        script = get_script(at)
        assert "'col_dispersion'] = 0.0" in script
        _assert_valid_python(script)

    def test_surface_diffusion(self):
        at = setup_app(
            add_particles="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
            particle_has_surfDiff="Yes",
        )
        script = get_script(at)
        assert "'has_surface_diffusion'] = 1" in script
        assert "'surface_diffusion'" in script
        _assert_valid_python(script)

    def test_cstr_template(self):
        at = setup_app(column_type="Mixed tank")
        assert not at.exception
        script = get_script(at)
        assert "'CSTR'" in script
        assert "'init_volume'" in script
        assert "'col_length'" not in script
        assert "'velocity'" not in script
        _assert_valid_python(script)

    def test_radial_template(self):
        at = setup_app(column_type="Radial flow cylinder")
        script = get_script(at)
        assert "'RADIAL_COLUMN_MODEL_1D'" in script
        assert "'col_radius_inner'" in script
        assert "'col_radius_outer'" in script
        assert "'velocity_coeff'" in script
        _assert_valid_python(script)

    def test_frustum_template(self):
        at = setup_app(column_type="Frustum")
        script = get_script(at)
        assert "'FRUSTUM_COLUMN_MODEL_1D'" in script
        assert "'col_length'" in script
        assert "'col_radius_inlet'" in script
        assert "'col_radius_outlet'" in script
        _assert_valid_python(script)

    def test_rapid_equilibrium_binding(self):
        at = setup_app(
            advanced_mode="On",
            PSD="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
            binding_model="Linear",
            req_binding="Rapid-equilibrium",
        )
        script = get_script(at)
        assert "'is_kinetic': 0" in script
        _assert_valid_python(script)

    def test_generated_script_runs(self):
        at = setup_app(
            add_particles="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
        )
        script = get_script(at)
        namespace = {}
        exec(script, namespace)
        result = namespace['get_unit_operation']()
        assert isinstance(result, dict)
        assert result['UNIT_TYPE'] == 'COLUMN_MODEL_1D'
        assert 'particle_type_000' in result
        assert result['ncomp'] == 1

    def test_generated_script_runs_with_ncomp(self):
        at = setup_app(
            add_particles="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
        )
        script = get_script(at)
        namespace = {}
        exec(script, namespace)
        result = namespace['get_unit_operation'](ncomp=3)
        assert result['ncomp'] == 3
        assert len(result['init_c']) == 3
        par = result['particle_type_000']
        assert len(par['film_diffusion']) == 3
        assert len(par['pore_diffusion']) == 3
        assert len(par['nbound']) == 3

    def test_sma_script_runs(self):
        at = setup_app(
            advanced_mode="On",
            PSD="Yes",
            has_binding="Yes",
            particle_resolution="1D (radial coordinate)",
            binding_model="SMA",
        )
        script = get_script(at)
        namespace = {}
        exec(script, namespace)
        result = namespace['get_unit_operation']()
        assert result['ncomp'] == 4
        par = result['particle_type_000']
        assert par['adsorption_model'] == 'STERIC_MASS_ACTION'
        assert len(par['adsorption']['sma_ka']) == 4
        assert par['adsorption']['sma_ka'][0] == 0.0


@pytest.mark.ci
@pytest.mark.unit_test
class TestGenerateCrystallizationTemplateIntegration:

    @staticmethod
    def _setup_cry_app(**overrides):
        at = AppTest.from_file("../Equation-Generator.py")
        at.run()
        assert not at.exception
        at.button(key="dev_mode_button").click().run()
        assert not at.exception
        at.button(key="model_type_crystallization_button").click().run()
        assert not at.exception
        for key, value in overrides.items():
            at.selectbox(key=key).set_value(value).run()
            assert not at.exception
        return at

    def test_dpfr_without_axial_dispersion(self):
        at = self._setup_cry_app(
            cry_column_type="DPFR",
            cry_has_axial_dispersion="No",
        )
        script = get_script(at)
        assert "'LUMPED_RATE_MODEL_WITHOUT_PORES'" in script
        assert "'col_dispersion'] = 0.0" in script
        _assert_valid_python(script)

    def test_cstr_template(self):
        at = self._setup_cry_app(cry_column_type="CSTR")
        script = get_script(at)
        assert "'CSTR'" in script
        assert "'reaction_model'] = 'CRYSTALLIZATION'" in script
        assert "cry['cry_mode']" in script
        _assert_valid_python(script)


def _assert_valid_python(script):
    try:
        compile(script, "<template>", "exec")
    except SyntaxError as e:
        pytest.fail(f"Generated script has syntax error: {e}\n\nScript:\n{script}")
