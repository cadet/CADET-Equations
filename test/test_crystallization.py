# -*- coding: utf-8 -*-
"""
Unit tests for crystallization / population balance model equations and UI.
"""

import pytest
from streamlit.testing.v1 import AppTest

from src import equations as eq


# %% CRY_MODE helpers

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("mode, expected", [
    (1, True), (2, False), (3, True), (4, False),
    (5, True), (6, False), (7, True),
])
def test_cry_has_primary_formation(mode, expected):
    assert eq.cry_has_primary_formation(mode) == expected


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("mode, expected", [
    (1, False), (2, True), (3, True), (4, False),
    (5, False), (6, True), (7, True),
])
def test_cry_has_aggregation(mode, expected):
    assert eq.cry_has_aggregation(mode) == expected


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("mode, expected", [
    (1, False), (2, False), (3, False), (4, True),
    (5, True), (6, True), (7, True),
])
def test_cry_has_fragmentation(mode, expected):
    assert eq.cry_has_fragmentation(mode) == expected


# %% Constitutive relations

@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_supersaturation():
    result = eq.cry_supersaturation()
    assert r"c_{\mathrm{eq}}" in result
    assert r"s =" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("size_dep, expected_frag, not_expected", [
    (True, r"\gamma", []),
    (False, r"k_g", [r"\gamma"]),
])
def test_cry_growth_rate(size_dep, expected_frag, not_expected):
    result = eq.cry_growth_rate(size_dep)
    assert expected_frag in result
    for frag in not_expected:
        assert frag not in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_primary_nucleation():
    result = eq.cry_primary_nucleation()
    assert r"k_p" in result
    assert r"s^u" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_secondary_nucleation():
    result = eq.cry_secondary_nucleation()
    assert r"k_b" in result
    assert r"M^k" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_total_nucleation():
    result = eq.cry_total_nucleation()
    assert r"B_p" in result
    assert r"B_s" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_suspension_density():
    result = eq.cry_suspension_density()
    assert r"k_v" in result
    assert r"\rho" in result
    assert r"\int" in result


# %% PBE equations

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_primary, has_gd, has_agg, has_frag", [
    (True, False, False, False),
    (True, True, False, False),
    (False, False, True, False),
    (False, False, False, True),
    (True, True, True, True),
])
def test_cry_pbe_cstr(has_primary, has_gd, has_agg, has_frag):
    result = eq.cry_pbe_cstr(has_primary, has_gd, has_agg, has_frag)
    assert r"\begin{align}" in result
    if has_primary:
        assert r"v_G" in result
    if has_gd:
        assert r"D_g" in result
    if has_agg:
        assert r"B_{\mathrm{agg}}" in result
    if has_frag:
        assert r"B_{\mathrm{frag}}" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_primary, has_ax_disp, has_gd, has_agg, has_frag", [
    (True, True, False, False, False),
    (True, False, True, False, False),
    (False, True, False, True, False),
    (False, True, False, False, True),
    (True, True, True, True, True),
])
def test_cry_pbe_dpfr(has_primary, has_ax_disp, has_gd, has_agg, has_frag):
    result = eq.cry_pbe_dpfr(has_primary, has_ax_disp, has_gd, has_agg, has_frag)
    assert r"\begin{align}" in result
    assert r"v_{\mathrm{ax}}" in result
    if has_ax_disp:
        assert r"D_{\mathrm{ax}}" in result
    if has_primary:
        assert r"v_G" in result
    if has_agg:
        assert r"B_{\mathrm{agg}}" in result
    if has_frag:
        assert r"B_{\mathrm{frag}}" in result


# %% Mass balance

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_primary", [True, False])
def test_cry_mass_balance_cstr(has_primary):
    result = eq.cry_mass_balance_cstr(has_primary)
    assert r"F_{\mathrm{in}}" in result
    if has_primary:
        assert r"\rho k_v" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_primary, has_ax_disp", [
    (True, True), (True, False), (False, True),
])
def test_cry_mass_balance_dpfr(has_primary, has_ax_disp):
    result = eq.cry_mass_balance_dpfr(has_primary, has_ax_disp)
    assert r"v_{\mathrm{ax}}" in result
    if has_ax_disp:
        assert r"D_{\mathrm{ax}}" in result
    if has_primary:
        assert r"\rho k_v" in result


# %% Boundary conditions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_primary, has_gd", [
    (True, True), (True, False),
])
def test_cry_pbe_bc_internal(has_primary, has_gd):
    result = eq.cry_pbe_bc_internal(has_primary, has_gd)
    assert r"\begin{align}" in result
    if has_gd:
        assert r"D_g" in result
    else:
        assert r"v_G(x_c)" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_pbe_bc_internal_no_primary():
    result = eq.cry_pbe_bc_internal(False, False)
    assert result == ""


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_ax_disp, expected_frag", [
    (True, r"D_{\mathrm{ax}}"),
    (False, r"n_{\mathrm{in},x}"),
])
def test_cry_pbe_bc_external_dpfr(has_ax_disp, expected_frag):
    result = eq.cry_pbe_bc_external_dpfr(has_ax_disp)
    assert expected_frag in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_ax_disp, expected_frag", [
    (True, r"D_{\mathrm{ax}}"),
    (False, r"c_{\mathrm{in}}"),
])
def test_cry_solute_bc_dpfr(has_ax_disp, expected_frag):
    result = eq.cry_solute_bc_dpfr(has_ax_disp)
    assert expected_frag in result


# %% Aggregation

@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_aggregation_birth_death():
    result = eq.cry_aggregation_birth_death()
    assert r"B_{\mathrm{agg}}" in result
    assert r"D_{\mathrm{agg}}" in result
    assert r"\beta" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("kernel_idx", [0, 1, 2, 3, 4])
def test_cry_aggregation_kernel(kernel_idx):
    result = eq.cry_aggregation_kernel(kernel_idx)
    assert r"\beta_0" in result
    assert r"\begin{align}" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_aggregation_kernel_unknown_falls_back():
    result = eq.cry_aggregation_kernel(99)
    assert r"\beta_0" in result


# %% Fragmentation

@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_fragmentation_birth_death():
    result = eq.cry_fragmentation_birth_death()
    assert r"B_{\mathrm{frag}}" in result
    assert r"D_{\mathrm{frag}}" in result
    assert r"S(" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_selection_function():
    result = eq.cry_selection_function()
    assert r"S_0" in result
    assert r"\alpha" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_cry_breakage_function():
    result = eq.cry_breakage_function()
    assert r"\gamma" in result
    assert r"\lambda" in result


# %% Assumptions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("col_type, cry_mode", [
    ("CSTR", 1), ("DPFR", 1), ("CSTR", 2),
    ("CSTR", 4), ("CSTR", 7),
])
def test_cry_assumptions(col_type, cry_mode):
    result = eq.cry_assumptions(col_type, cry_mode)
    assert len(result) >= 3
    if col_type == "CSTR":
        assert any("perfectly mixed" in a for a in result)
    else:
        assert any("1D axial" in a for a in result)
    if eq.cry_has_primary_formation(cry_mode):
        assert any("nucleation" in a for a in result)
    if eq.cry_has_aggregation(cry_mode):
        assert any("aggregation" in a for a in result)
    if eq.cry_has_fragmentation(cry_mode):
        assert any("fragmentation" in a for a in result)


# %% Streamlit UI smoke tests

@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_cstr_primary_formation():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.selectbox(key="model_type").set_value("Crystallization").run()
    assert not at.exception
    assert "Crystallization" in at.session_state.latex_string


@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_dpfr():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    at.selectbox(key="cry_column_type").set_value("DPFR").run()
    assert not at.exception
    assert "DPFR" in at.session_state.latex_string


@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_aggregation_mode():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    at.selectbox(key="cry_mode").set_value("Aggregation").run()
    assert not at.exception
    assert "aggregation" in at.session_state.latex_string.lower()


@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_fragmentation_mode():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    at.selectbox(key="cry_mode").set_value("Fragmentation").run()
    assert not at.exception
    assert "fragmentation" in at.session_state.latex_string.lower()


@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_all_mechanisms():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    at.selectbox(key="cry_mode").set_value(
        "Primary particle formation + Aggregation + Fragmentation").run()
    assert not at.exception
    latex = at.session_state.latex_string.lower()
    assert "aggregation" in latex
    assert "fragmentation" in latex


@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_secondary_nucleation():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    at.selectbox(key="cry_has_secondary_nucleation").set_value("Yes").run()
    assert not at.exception
    assert "k_b" in at.session_state.latex_string


@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_size_dependent_growth():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    at.selectbox(key="cry_size_dependent_growth").set_value("Yes").run()
    assert not at.exception
    assert r"\gamma" in at.session_state.latex_string


@pytest.mark.ci
@pytest.mark.unit_test
def test_crystallization_growth_dispersion():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    at.selectbox(key="cry_has_growth_dispersion").set_value("Yes").run()
    assert not at.exception
    assert "D_g" in at.session_state.latex_string


@pytest.mark.ci
@pytest.mark.unit_test
def test_switching_back_to_chromatography():
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    at.selectbox(key="model_type").set_value("Crystallization").run()
    assert not at.exception
    at.selectbox(key="model_type").set_value("Chromatography").run()
    assert not at.exception
    assert "Crystallization" not in at.session_state.latex_string
