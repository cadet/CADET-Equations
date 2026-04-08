# -*- coding: utf-8 -*-
"""
SPDX-License-Identifier: GPL-3.0-only
Copyright (c) 2026 Jan Michael Breuer
See LICENSE file for details.

Unit tests for src/equations.py.
Each function is tested in isolation with relevant parameter combinations.
"""

import pytest
from unittest.mock import MagicMock
from src import equations as eq


# %% Helpers

def _make_particle(**overrides):
    """Create a minimal Particle-like mock object.

    Avoids importing the Streamlit-dependent Particle dataclass from
    Equation-Generator.py while providing the attributes that equation
    functions expect.
    """
    defaults = dict(
        geometry="Sphere",
        resolution="1D",
        has_core=False,
        has_binding=True,
        req_binding=False,
        has_mult_bnd_states=False,
        has_surfDiff=False,
        nonlimiting_filmDiff=False,
        surface_volume_ratio=3,
        interstitial_volume_resolution="1D",
        single_partype=True,
        PTD=False,
    )
    defaults.update(overrides)
    return MagicMock(**defaults)


# %% HRM_asmpt

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("column_resolution, N_p, has_binding, expected_keyword", [
    ("0D (Homogeneous Tank)", 0, False, "tank"),
    ("1D (axial coordinate)", 0, False, "column"),
    ("1D (axial coordinate)", 1, True, "spherical"),
    ("1D (axial coordinate)", 0, False, "no solid phase"),
])
def test_HRM_asmpt_device_and_particle_text(column_resolution, N_p, has_binding, expected_keyword):
    """Verify that HRM assumptions adapt wording to device type and particle presence."""
    result = eq.HRM_asmpt(N_p=N_p, nonlimiting_filmDiff=False,
                          has_binding=has_binding, has_surfDiff=False,
                          column_resolution=column_resolution)
    combined = " ".join(result)
    assert expected_keyword in combined


@pytest.mark.ci
@pytest.mark.unit_test
def test_HRM_asmpt_binding_toggle():
    """Binding-related assumptions should only appear when has_binding is True."""
    with_bnd = eq.HRM_asmpt(N_p=1, nonlimiting_filmDiff=False,
                             has_binding=True, has_surfDiff=False,
                             column_resolution="1D (axial coordinate)")
    without_bnd = eq.HRM_asmpt(N_p=1, nonlimiting_filmDiff=False,
                                has_binding=False, has_surfDiff=False,
                                column_resolution="1D (axial coordinate)")
    assert "binding model" in " ".join(with_bnd)
    assert "binding model" not in " ".join(without_bnd)


@pytest.mark.ci
@pytest.mark.unit_test
def test_HRM_asmpt_filters_empty_strings():
    """Empty assumption strings resulting from disabled options should be removed."""
    result = eq.HRM_asmpt(N_p=0, nonlimiting_filmDiff=False,
                          has_binding=False, has_surfDiff=False,
                          column_resolution="1D (axial coordinate)")
    assert "" not in result


# %% Interstitial volume continuum assumptions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("resolution, expected_keyword", [
    ("3D", "continuum"),
    ("2D", "radially symmetric"),
    ("1D", "radially symmetric and homogeneous"),
    ("0D", "spatially homogeneous"),
])
def test_int_vol_continuum_asmpt_per_resolution(resolution, expected_keyword):
    """Each resolution should produce its characteristic assumption text."""
    result = eq.int_vol_continuum_asmpt(resolution, N_p=1, nonlimiting_filmDiff=False)
    combined = " ".join(result)
    assert expected_keyword in combined


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("N_p, expected_keyword", [
    (0, "tank"),
    (1, "interstitial volume"),
])
def test_int_vol_0D_volume_naming(N_p, expected_keyword):
    """0D assumption text refers to 'tank' without particles and 'interstitial volume' with."""
    result = eq.int_vol_0DContinuum_asmpt(N_p=N_p, nonlimiting_filmDiff=False)
    combined = " ".join(result)
    assert expected_keyword in combined


@pytest.mark.ci
@pytest.mark.unit_test
def test_int_vol_3D_multiple_particle_types():
    """Multiple particle types (N_p > 1) should mention representative radii."""
    result = eq.int_vol_3DContinuum_asmpt(N_p=2, nonlimiting_filmDiff=False)
    combined = " ".join(result)
    assert "N^{\\mathrm{p}}" in combined


# %% Particle assumptions

@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_asmpts_base():
    """Base particle assumptions should reference the outer particle surface."""
    result = eq.particle_asmpts()
    assert len(result) == 1
    assert "outer particle surface" in result[0]


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_surfDiff, expected_keyword", [
    (True, "surface diffusion"),
    (False, "no surface diffusion"),
])
def test_particle_1D_surfDiff_toggle(has_surfDiff, expected_keyword):
    """1D particle assumptions should reflect whether surface diffusion is active."""
    result = eq.particle_1D_asmpt(has_surfDiff=has_surfDiff)
    combined = " ".join(result)
    assert expected_keyword in combined


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_0D_infinite_diffusion():
    """0D particles should state that diffusion is infinitely fast."""
    result = eq.particle_0D_asmpt()
    assert "infinitely fast" in " ".join(result)


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("resolution, has_surfDiff", [
    ("0D", False),
    ("1D", True),
    ("1D", False),
])
def test_particle_asmpt_router(resolution, has_surfDiff):
    """The router function should delegate to the correct resolution-specific function."""
    result = eq.particle_asmpt(resolution, has_surfDiff=has_surfDiff)
    if resolution == "0D":
        assert result == eq.particle_0D_asmpt()
    else:
        assert result == eq.particle_1D_asmpt(has_surfDiff=has_surfDiff)


# %% Equation term generators (bulk transport)

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("term_func, latex_fragment", [
    (eq.bulk_time_derivative, r"\frac{\partial c^{\b}_i}{\partial t}"),
    (eq.solid_time_derivative, r"\frac{\partial c^{\s}_i}{\partial t}"),
    (eq.axial_convection, r"\partial z"),
    (eq.axial_dispersion, r"D^{\mathrm{ax}}"),
    (eq.radial_dispersion, r"D^{\mathrm{rad}}"),
    (eq.angular_dispersion, r"D^{\mathrm{ang}}"),
])
def test_transport_term_without_eps(term_func, latex_fragment):
    """Transport terms without porosity should contain their characteristic LaTeX fragment."""
    assert latex_fragment in term_func()


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("term_func", [
    eq.bulk_time_derivative,
    eq.solid_time_derivative,
    eq.axial_convection,
    eq.axial_dispersion,
    eq.radial_dispersion,
    eq.angular_dispersion,
])
def test_transport_term_with_eps(term_func):
    """When porosity is provided, all transport terms should include it."""
    result = term_func(eps=r"\varepsilon")
    assert r"\varepsilon" in result


# %% Film diffusion term

@pytest.mark.ci
@pytest.mark.unit_test
def test_int_filmDiff_single_particle():
    """Single-particle film diffusion should not contain a summation."""
    particle = _make_particle(surface_volume_ratio=3, resolution="1D")
    result = eq.int_filmDiff_term(particle, 1, 1, singleParticle=True,
                                   nonLimitingFilmDiff=False, hasSurfDiff=False)
    assert r"k^{\mathrm{f}}" in result
    assert r"\sum" not in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_int_filmDiff_multiple_particles():
    """Multiple-particle film diffusion should contain a summation."""
    particle = _make_particle(surface_volume_ratio=3, resolution="1D")
    result = eq.int_filmDiff_term(particle, 1, 3, singleParticle=False,
                                   nonLimitingFilmDiff=False, hasSurfDiff=False)
    assert r"\sum" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_int_filmDiff_nonlimiting_0D_returns_empty():
    """Non-limiting film diffusion with 0D resolution should produce an empty string."""
    particle = _make_particle(resolution="0D")
    result = eq.int_filmDiff_term(particle, 1, 1, singleParticle=True,
                                   nonLimitingFilmDiff=True, hasSurfDiff=False)
    assert result == ""


@pytest.mark.ci
@pytest.mark.unit_test
def test_int_filmDiff_nonlimiting_1D_substitutes_BC():
    """Non-limiting film diffusion in 1D substitutes the BC into the term, removing k_f."""
    particle = _make_particle(surface_volume_ratio=3, resolution="1D")
    result = eq.int_filmDiff_term(particle, 1, 1, singleParticle=False,
                                   nonLimitingFilmDiff=True, hasSurfDiff=False)
    assert r"D^{\mathrm{p}}" in result
    assert r"k^{\mathrm{f}}" not in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_int_filmDiff_nonlimiting_with_surfDiff():
    """Non-limiting film diffusion with surface diffusion should include D^s."""
    particle = _make_particle(surface_volume_ratio=3, resolution="1D")
    result = eq.int_filmDiff_term(particle, 1, 1, singleParticle=False,
                                   nonLimitingFilmDiff=True, hasSurfDiff=True)
    assert r"D^{\mathrm{s}}" in result


# %% Interstitial volume boundary conditions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("resolution, hasAxDisp, expected, not_expected", [
    ("1D", True,  [r"D^{\mathrm{ax}}", r"z=0", r"z=L"], []),
    ("1D", False, [r"c_{\mathrm{in}"],                    [r"z=L"]),
    ("2D", True,  [r"D^{\mathrm{rad}}", r"R^{\mathrm{c}}"], []),
    ("3D", True,  [r"D^{\mathrm{ang}}", r"2\pi"],         []),
])
def test_int_vol_BC(resolution, hasAxDisp, expected, not_expected):
    """Boundary conditions should include/exclude terms based on resolution and dispersion."""
    result = eq.int_vol_BC(resolution, hasAxialDispersion=hasAxDisp)
    for frag in expected:
        assert frag in result
    for frag in not_expected:
        assert frag not in result


# %% Interstitial volume initial conditions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("resolution, includeParLiquid, has_bulk, has_par", [
    ("1D", False, True,  False),
    ("1D", True,  True,  True),
    ("2D", False, True,  False),
    ("3D", True,  True,  True),
])
def test_int_vol_initial(resolution, includeParLiquid, has_bulk, has_par):
    """Initial conditions should include bulk and optionally particle liquid concentrations."""
    result = eq.int_vol_initial(resolution, includeParLiquid=includeParLiquid)
    assert (r"c^{\b}" in result) == has_bulk
    assert (r"c^{\p}" in result) == has_par


@pytest.mark.ci
@pytest.mark.unit_test
def test_int_vol_initial_2D_domain():
    """2D initial conditions domain should reference column radius."""
    result = eq.int_vol_initial("2D", includeParLiquid=False)
    assert r"R^{\mathrm{c}}" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_int_vol_initial_3D_domain():
    """3D initial conditions domain should include angular component."""
    result = eq.int_vol_initial("3D", includeParLiquid=False)
    assert r"2\pi" in result


# %% Domain functions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("resolution, with_time, expected, not_expected", [
    ("0D", True,  [r"T^\mathrm{end}"],              [r"L"]),
    ("1D", True,  [r"T^\mathrm{end}", r"(0, L)"],   []),
    ("2D", True,  [r"R^\mathrm{c}"],                []),
    ("3D", True,  [r"2\pi"],                         []),
    ("1D", False, [r"(0, L)"],                       [r"T^\mathrm{end}"]),
])
def test_int_vol_domain(resolution, with_time, expected, not_expected):
    """Interstitial volume domain should grow with resolution dimension."""
    result = eq.int_vol_domain(resolution, with_time_domain=with_time)
    for frag in expected:
        assert frag in result
    for frag in not_expected:
        assert frag not in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("hasCore, with_par_index, expected, not_expected", [
    (False, False, [r"(0, R^{\mathrm{p}})"],  [r"R^{\mathrm{pc}}", r"_{j}"]),
    (True,  False, [r"R^{\mathrm{pc}}"],       []),
    (False, True,  [r"_{j}"],                  []),
])
def test_particle_domain(hasCore, with_par_index, expected, not_expected):
    """Particle domain should reflect core presence and particle index."""
    result = eq.particle_domain("1D", hasCore=hasCore, with_par_index=with_par_index)
    for frag in expected:
        assert frag in result
    for frag in not_expected:
        assert frag not in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("col_res, par_res, hasCore, with_time, expected, not_expected", [
    ("1D", "1D", False, True,  [r"(0, L)", r"R^{\mathrm{p}}"],        []),
    ("0D", "0D", False, True,  [r"T^\mathrm{end}"],                   [r"R^{\mathrm{p}}"]),
    ("2D", "1D", False, True,  [r"R^\mathrm{c}"],                     []),
    ("1D", "1D", True,  True,  [r"R^{\mathrm{pc}}"],                  []),
    ("1D", "1D", False, False, [],                                     [r"T^\mathrm{end}"]),
])
def test_full_particle_conc_domain(col_res, par_res, hasCore, with_time, expected, not_expected):
    """Full particle concentration domain should compose column and particle domains correctly."""
    result = eq.full_particle_conc_domain(col_res, par_res, hasCore=hasCore, with_time_domain=with_time)
    for frag in expected:
        assert frag in result
    for frag in not_expected:
        assert frag not in result


# %% Particle transport

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("resolution, geometry, has_binding, expected_frag", [
    ("0D", "Sphere", True,  r"\begin{align}"),
    ("0D", "Sphere", False, r"\begin{align}"),
    ("1D", "Sphere", True,  r"r^2"),
])
def test_particle_transport_resolution_and_geometry(resolution, geometry, has_binding, expected_frag):
    """Particle transport should select the correct equation form for given resolution and geometry."""
    particle = _make_particle(resolution=resolution, geometry=geometry)
    result = eq.particle_transport(
        particle, singleParticle=True, nonlimiting_filmDiff=False,
        has_surfDiff=False, has_binding=has_binding, req_binding=False,
        has_mult_bnd_states=False)
    assert expected_frag in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_single_removes_j_index():
    """Single-particle transport equations should not contain the particle-type index j."""
    particle = _make_particle(resolution="1D", geometry="Sphere")
    result = eq.particle_transport(
        particle, singleParticle=True, nonlimiting_filmDiff=False,
        has_surfDiff=False, has_binding=True, req_binding=False,
        has_mult_bnd_states=False)
    assert "j,i" not in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_multiple_keeps_j_index():
    """Multi-particle-type transport equations should contain the particle-type index j."""
    particle = _make_particle(resolution="1D", geometry="Sphere")
    result = eq.particle_transport(
        particle, singleParticle=False, nonlimiting_filmDiff=False,
        has_surfDiff=False, has_binding=True, req_binding=False,
        has_mult_bnd_states=False)
    assert "j,i" in result or "j}" in result


# %% Particle transport radial (geometry-specific)

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("geometry, expected_frag", [
    ("Sphere",   r"r^2"),
    ("Cylinder", r"\frac{1}{r}"),
    ("Slab",     r"\frac{\partial }{\partial r}"),
])
def test_particle_transport_radial_geometry(geometry, expected_frag):
    """Radial transport should use the correct differential operator for each geometry."""
    result = eq.particle_transport_radial(
        geometry, has_surfDiff=False, has_binding=True,
        req_binding=False, has_mult_bnd_states=False)
    assert expected_frag in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_radial_surfDiff_toggle():
    """Surface diffusion coefficient should appear only when has_surfDiff is True."""
    with_sd = eq.particle_transport_radial(
        "Sphere", has_surfDiff=True, has_binding=True,
        req_binding=False, has_mult_bnd_states=False)
    without_sd = eq.particle_transport_radial(
        "Sphere", has_surfDiff=False, has_binding=True,
        req_binding=False, has_mult_bnd_states=False)
    assert r"D_{j,i}^{\s}" in with_sd
    assert r"D_{j,i}^{\s}" not in without_sd


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_radial_req_binding():
    """Rapid-equilibrium binding should set the solid phase LHS to zero."""
    result = eq.particle_transport_radial(
        "Sphere", has_surfDiff=False, has_binding=True,
        req_binding=True, has_mult_bnd_states=False)
    # The solid equation starts with "0 &="
    assert r"0 " in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_radial_no_binding():
    """Without binding, the solid phase concentration should not appear."""
    result = eq.particle_transport_radial(
        "Sphere", has_surfDiff=False, has_binding=False,
        req_binding=False, has_mult_bnd_states=False)
    assert r"c^{\s}" not in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_radial_mult_bnd_states():
    """Multiple bound states should introduce a sum over bound state index k."""
    result = eq.particle_transport_radial(
        "Sphere", has_surfDiff=False, has_binding=True,
        req_binding=False, has_mult_bnd_states=True)
    assert r"N^{\mathrm{b}}" in result


# %% Particle transport homogeneous

@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_homogeneous_liquid_binding_toggle():
    """Homogeneous liquid transport should include binding term only when binding is active."""
    with_bnd = eq.particle_transport_homogeneous_liquid(
        has_binding=True, req_binding=False, has_mult_bnd_states=False)
    without_bnd = eq.particle_transport_homogeneous_liquid(
        has_binding=False, req_binding=False, has_mult_bnd_states=False)
    assert r"f^{\mathrm{bind}}" in with_bnd
    assert r"f^{\mathrm{bind}}" not in without_bnd


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_homogeneous_solid_req_binding():
    """Rapid-equilibrium solid transport should have zero on the LHS."""
    result = eq.particle_transport_homogeneous_solid(
        req_binding=True, has_mult_bnd_states=False)
    assert "0" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_homogeneous_combined():
    """Combined homogeneous transport should wrap both equations in align environment."""
    result = eq.particle_transport_homogeneous(
        has_binding=True, req_binding=False, has_mult_bnd_states=False)
    assert r"\begin{align}" in result
    assert r"\end{align}" in result


# %% Particle boundary conditions

@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_boundary_0D_returns_empty():
    """0D particles have no spatial boundary, so boundary conditions should be empty."""
    particle = _make_particle(resolution="0D")
    result = eq.particle_boundary(
        particle, singleParticle=True, nonlimiting_filmDiff=False,
        has_surfDiff=False, has_binding=True, req_binding=False,
        has_mult_bnd_states=False)
    assert result == ""


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_core, expected_frag", [
    (True,  r"R^{\mathrm{pc}}"),
    (False, r"|_{r=0}"),
])
def test_particle_boundary_core_toggle(has_core, expected_frag):
    """Inner boundary location should depend on particle core presence."""
    particle = _make_particle(resolution="1D", has_core=has_core)
    result = eq.particle_boundary(
        particle, singleParticle=True, nonlimiting_filmDiff=False,
        has_surfDiff=False, has_binding=True, req_binding=False,
        has_mult_bnd_states=False)
    assert expected_frag in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_boundary_nonlimiting():
    """Non-limiting film diffusion should set particle concentration equal to bulk at boundary."""
    particle = _make_particle(resolution="1D", has_core=False)
    result = eq.particle_boundary(
        particle, singleParticle=True, nonlimiting_filmDiff=True,
        has_surfDiff=False, has_binding=True, req_binding=False,
        has_mult_bnd_states=False)
    assert r"c^{\b}_i" in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_boundary_single_particle_removes_j():
    """Single-particle boundary conditions should not contain the particle type index."""
    particle = _make_particle(resolution="1D", has_core=False)
    result = eq.particle_boundary(
        particle, singleParticle=True, nonlimiting_filmDiff=False,
        has_surfDiff=False, has_binding=True, req_binding=False,
        has_mult_bnd_states=False)
    assert "j," not in result


# %% Particle initial conditions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("includeParLiquid, has_p, has_s", [
    (True,  True,  True),
    (False, False, True),
])
def test_particle_initial_liquid_toggle(includeParLiquid, has_p, has_s):
    """Particle initial conditions should include liquid phase only when requested."""
    domain = r"$(0, R^{\mathrm{p}})$"
    result = eq.particle_initial(domain, singleParticle=True, includeParLiquid=includeParLiquid)
    assert (r"c^{\p}" in result) == has_p
    assert (r"c^{\s}" in result) == has_s


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_initial_single_particle_removes_comma_j():
    """Single-particle initial conditions should not contain ',j' subscript patterns."""
    domain = r"$(0, R^{\mathrm{p}})$"
    result = eq.particle_initial(domain, singleParticle=True, includeParLiquid=True)
    assert ",j" not in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_initial_multiple_particles_keeps_j():
    """Multi-particle-type initial conditions should retain the particle-type index j."""
    domain = r"$(0, R^{\mathrm{p}}_{j})$"
    result = eq.particle_initial(domain, singleParticle=False, includeParLiquid=True)
    assert "j," in result


# %% Reaction terms

@pytest.mark.ci
@pytest.mark.unit_test
def test_bulk_reaction_term():
    """Bulk reaction term should contain the reaction function for bulk phase."""
    result = eq.bulk_reaction_term()
    assert r"f^{\mathrm{react},\b}" in result
    assert r"\vec{c}^{\b}" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("singleParticle, expected_idx", [
    (True,  r"_{i}"),
    (False, r"_{j,i}"),
])
def test_particle_liquid_reaction_term(singleParticle, expected_idx):
    """Particle liquid reaction term should use correct index based on particle count."""
    result = eq.particle_liquid_reaction_term(singleParticle)
    assert r"f^{\mathrm{react},\p}" in result
    assert expected_idx in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("singleParticle, expected_idx", [
    (True,  r"_{i}"),
    (False, r"_{j,i}"),
])
def test_particle_solid_reaction_term(singleParticle, expected_idx):
    """Particle solid reaction term should use correct index based on particle count."""
    result = eq.particle_solid_reaction_term(singleParticle)
    assert r"f^{\mathrm{react},\s}" in result
    assert expected_idx in result


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_homogeneous_liquid_reaction():
    """Homogeneous liquid transport should include reaction term when enabled."""
    with_react = eq.particle_transport_homogeneous_liquid(
        has_binding=True, req_binding=False, has_mult_bnd_states=False, has_reaction_liquid=True)
    without_react = eq.particle_transport_homogeneous_liquid(
        has_binding=True, req_binding=False, has_mult_bnd_states=False, has_reaction_liquid=False)
    assert r"f^{\mathrm{react},\p}" in with_react
    assert r"f^{\mathrm{react},\p}" not in without_react


@pytest.mark.ci
@pytest.mark.unit_test
def test_particle_transport_homogeneous_solid_reaction():
    """Homogeneous solid transport should include reaction term when enabled."""
    with_react = eq.particle_transport_homogeneous_solid(
        req_binding=False, has_mult_bnd_states=False, has_reaction_solid=True)
    without_react = eq.particle_transport_homogeneous_solid(
        req_binding=False, has_mult_bnd_states=False, has_reaction_solid=False)
    assert r"f^{\mathrm{react},\s}" in with_react
    assert r"f^{\mathrm{react},\s}" not in without_react


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("geometry", ["Sphere", "Cylinder", "Slab"])
def test_particle_transport_radial_reaction(geometry):
    """Radial transport should include reaction terms when enabled for all geometries."""
    with_react = eq.particle_transport_radial(
        geometry, has_surfDiff=False, has_binding=True,
        req_binding=False, has_mult_bnd_states=False,
        has_reaction_liquid=True, has_reaction_solid=True)
    without_react = eq.particle_transport_radial(
        geometry, has_surfDiff=False, has_binding=True,
        req_binding=False, has_mult_bnd_states=False,
        has_reaction_liquid=False, has_reaction_solid=False)
    assert r"f^{\mathrm{react},\p}" in with_react
    assert r"f^{\mathrm{react},\s}" in with_react
    assert r"f^{\mathrm{react}" not in without_react


# %% Binding model terms

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("PTD, expected_idx", [
    (True, "j,i"),
    (False, "i"),
])
def test_binding_term_arbitrary(PTD, expected_idx):
    """Arbitrary binding term should use f^bind with correct indices."""
    result = eq.binding_term_arbitrary(PTD=PTD)
    assert r"f^{\mathrm{bind}}" in result
    assert expected_idx in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("PTD, expected_idx", [
    (True, "j,i"),
    (False, "i"),
])
def test_binding_term_linear(PTD, expected_idx):
    """Linear binding term should contain k_a, k_d, c^p, and c^s."""
    result = eq.binding_term_linear(PTD=PTD)
    assert r"k^{\mathrm{a}}" in result
    assert r"k^{\mathrm{d}}" in result
    assert r"c^{\p}" in result
    assert r"c^{\s}" in result
    assert expected_idx in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("PTD, expected_idx", [
    (True, "j,i"),
    (False, "i"),
])
def test_binding_term_langmuir(PTD, expected_idx):
    """Langmuir binding term should contain k_a, k_d, q_max, and summation."""
    result = eq.binding_term_langmuir(PTD=PTD)
    assert r"k^{\mathrm{a}}" in result
    assert r"k^{\mathrm{d}}" in result
    assert r"q^{\mathrm{max}}" in result
    assert r"\sum" in result
    assert expected_idx in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("PTD, expected_idx", [
    (True, "j,i"),
    (False, "i"),
])
def test_binding_term_sma(PTD, expected_idx):
    """SMA binding term should contain k_a, k_d, nu, bar_q, and reference concentrations."""
    result = eq.binding_term_sma(PTD=PTD)
    assert r"k^{\mathrm{a}}" in result
    assert r"k^{\mathrm{d}}" in result
    assert r"\nu" in result
    assert r"\bar{q}" in result
    assert r"q^{\mathrm{ref}}" in result
    assert r"c^{\mathrm{ref}}" in result
    assert expected_idx in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("PTD", [True, False])
def test_sma_free_binding_sites(PTD):
    """SMA free binding sites should contain Lambda, nu, sigma, and summation."""
    result = eq.sma_free_binding_sites(PTD=PTD)
    assert r"\bar{q}" in result
    assert r"\Lambda" in result
    assert r"\nu" in result
    assert r"\sigma" in result
    assert r"\sum" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("PTD", [True, False])
def test_sma_electroneutrality(PTD):
    """SMA electroneutrality should contain Lambda, nu, and summation."""
    result = eq.sma_electroneutrality(PTD=PTD)
    assert r"\Lambda" in result
    assert r"\nu" in result
    assert r"\sum" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("model", eq.BINDING_MODELS)
def test_get_binding_term(model):
    """get_binding_term should return a non-empty string for all known models."""
    result = eq.get_binding_term(model, PTD=True)
    assert isinstance(result, str)
    assert len(result) > 0


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_binding_term_unknown_model_falls_back():
    """Unknown binding model should fall back to arbitrary."""
    result = eq.get_binding_term("UnknownModel", PTD=True)
    assert r"f^{\mathrm{bind}}" in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("binding_model, expected_fragment", [
    ("Arbitrary", r"f^{\mathrm{bind}}"),
    ("Linear", r"k^{\mathrm{a}}"),
    ("Langmuir", r"q^{\mathrm{max}}"),
    ("SMA", r"\bar{q}"),
])
def test_particle_transport_radial_binding_model(binding_model, expected_fragment):
    """Radial transport should use correct binding term based on binding_model."""
    result = eq.particle_transport_radial(
        "Sphere", has_surfDiff=False, has_binding=True,
        req_binding=False, has_mult_bnd_states=False,
        binding_model=binding_model)
    assert expected_fragment in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("binding_model, expected_fragment", [
    ("Arbitrary", r"f^{\mathrm{bind}}"),
    ("Linear", r"k^{\mathrm{a}}"),
    ("Langmuir", r"q^{\mathrm{max}}"),
    ("SMA", r"\bar{q}"),
])
def test_particle_transport_homogeneous_binding_model(binding_model, expected_fragment):
    """Homogeneous transport should use correct binding term based on binding_model."""
    result = eq.particle_transport_homogeneous(
        has_binding=True, req_binding=False, has_mult_bnd_states=False,
        binding_model=binding_model)
    assert expected_fragment in result


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("binding_model, expected_fragment", [
    ("Arbitrary", r"f^{\mathrm{bind}}"),
    ("Linear", r"k^{\mathrm{d}}"),
])
def test_particle_transport_full_binding_model(binding_model, expected_fragment):
    """particle_transport should pass binding_model through to sub-functions."""
    particle = _make_particle(resolution="1D", geometry="Sphere")
    result = eq.particle_transport(
        particle, singleParticle=True, nonlimiting_filmDiff=False,
        has_surfDiff=False, has_binding=True, req_binding=False,
        has_mult_bnd_states=False, PTD=False, binding_model=binding_model)
    assert expected_fragment in result


# %% Binding model assumptions

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("binding_model, expected_keyword", [
    ("Arbitrary", "No custom assumptions"),
    ("Linear", "independently"),
    ("Langmuir", "finite number of binding sites"),
    ("SMA", "electroneutrality"),
])
def test_binding_model_assumptions(binding_model, expected_keyword):
    """Each binding model should return its characteristic assumptions."""
    result = eq.binding_model_assumptions(binding_model)
    assert result is not None
    combined = " ".join(result)
    assert expected_keyword in combined


@pytest.mark.ci
@pytest.mark.unit_test
def test_binding_model_assumptions_unknown_returns_none():
    """Unknown binding model should return None."""
    result = eq.binding_model_assumptions("UnknownModel")
    assert result is None


# %% Surface diffusion for Cylinder and Slab geometries

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("geometry, expected_frag", [
    ("Cylinder", r"\frac{1}{r}"),
    ("Slab", r"\frac{\partial }{\partial r}"),
])
def test_particle_transport_radial_surfDiff_cylinder_slab(geometry, expected_frag):
    """Cylinder and Slab geometries with surface diffusion should include D^s term."""
    result = eq.particle_transport_radial(
        geometry, has_surfDiff=True, has_binding=True,
        req_binding=False, has_mult_bnd_states=False)
    assert r"D_{j,i}^{\s}" in result
    assert expected_frag in result
