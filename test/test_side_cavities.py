"""
Tests for the side-cavity particle model (Model PC of the modeling paper).

The pore phase is partitioned into a main pore network (compartment index 0)
and N^sc side-cavity types (m = 1, ..., N^sc). Only the main pore network is
radially resolved and only it touches the particle surface; the cavities are
reached through cavity exchange.
"""

import pytest

from src import equations as eq
from src.model_particle import Particle


def _particle(resolution="1D", n_cavities=2, models=("Arbitrary", "Arbitrary"), reqs=(False, False), **kw):
    defaults = dict(
        geometry="Sphere",
        has_core=False,
        var_format="CADET",
        resolution=resolution,
        has_binding=True,
        req_binding=False,
        has_mult_bnd_states=False,
        has_surfDiff=False,
        nonlimiting_filmDiff=False,
        interstitial_volume_resolution="1D",
        has_side_cavities=True,
        n_side_cavities=n_cavities,
        side_cavity_binding_models=models,
        side_cavity_req_binding=reqs,
    )
    defaults.update(kw)
    return Particle(**defaults)


def _transport(particle, **kw):
    opts = dict(
        singleParticle=True,
        nonlimiting_filmDiff=False,
        has_surfDiff=particle.has_surfDiff,
        has_binding=True,
        req_binding=False,
        has_mult_bnd_states=False,
    )
    opts.update(kw)
    return eq.particle_transport(particle, **opts)


class TestSideCavityTransport:
    def test_main_pore_is_radially_resolved(self):
        result = _transport(_particle("1D"))
        assert r"\frac{\partial c^{\p}_{0,i}}{\partial t}" in result
        assert r"D^{\p}_{i} \frac{\partial c^{\p}_{0,i}}{\partial r}" in result

    def test_homogeneous_main_pore_has_film_source_instead_of_flux(self):
        """Without radial resolution the surface flux becomes a source term."""
        result = _transport(_particle("0D"))
        assert r"d^{\mathrm{mp}} \varepsilon^{\mathrm{mp}} \frac{\partial c^{\p}_{0,i}}{\partial t}" in result
        assert r"k^{\mathrm{f}}_{i} \left( c^{\b}_{i} - c^{\p}_{0,i} \right)" in result
        assert r"\frac{\partial }{\partial r}" not in result

    @pytest.mark.parametrize("resolution", ["1D", "0D"])
    def test_cavity_exchange_sink_in_main_pore(self, resolution):
        result = _transport(_particle(resolution))
        assert r"\sum_{m=1}^{N^{\mathrm{sc}}}" in result
        assert r"k^{\mathrm{sc}}_{m,i} \left( c^{\p}_{0,i} - c^{\p}_{m,i} \right)" in result
        assert (
            r"\frac{d^{\mathrm{sc}}_{m} \varepsilon^{\mathrm{sc}}_{m}}{d^{\mathrm{mp}} \varepsilon^{\mathrm{mp}}}"
            in result
        )

    @pytest.mark.parametrize("resolution", ["1D", "0D"])
    def test_one_equation_pair_per_cavity_type(self, resolution):
        result = _transport(_particle(resolution, n_cavities=3, models=("Arbitrary",) * 3, reqs=(False,) * 3))
        for m in (1, 2, 3):
            assert rf"\frac{{\partial c^{{\p}}_{{{m},i}}}}{{\partial t}}" in result
            assert rf"k^{{\mathrm{{sc}}}}_{{{m},i}}" in result

    def test_cavities_use_main_pore_porosity_only_in_main_pore(self):
        result = _transport(_particle("1D"))
        assert r"\varepsilon^{\mathrm{mp}}" in result
        assert r"\varepsilon^{\mathrm{sc}}_{1}" in result
        # the single-compartment particle porosity must not appear any more
        assert r"\varepsilon^{\mathrm{p}}" not in result

    def test_per_cavity_binding_models(self):
        """Remark (PCc): each cavity type may carry its own binding model."""
        result = _transport(_particle("1D", models=("Linear", "Langmuir")))
        # cavity 1 linear: no capacity term; cavity 2 Langmuir: capacity term
        assert r"q^{\mathrm{max}}_{2,i}" in result
        assert r"q^{\mathrm{max}}_{1,i}" not in result

    def test_langmuir_sum_index_does_not_clash_with_cavity_index(self):
        result = _transport(_particle("1D", models=("Langmuir", "Langmuir")))
        # the competition sum runs over n, the cavity sum over m
        assert r"\sum_{n=0}^{N^{\mathrm{c}} - 1}" in result
        assert r"\sum_{m=0}" not in result

    def test_rapid_equilibrium_cavity_gives_algebraic_constraint(self):
        result = _transport(_particle("1D", models=("Arbitrary", "Arbitrary"), reqs=(False, True)))
        assert r"0 &= f^{\mathrm{bind}}_{2,i}" in result
        # and the conserved moiety appears on the left-hand side
        assert (
            r"\frac{1 - \varepsilon^{\mathrm{sc}}_{2}}{\varepsilon^{\mathrm{sc}}_{2}} \frac{\partial c^{\s}_{2,i}}"
            in result
        )


class TestSideCavityBoundary:
    def _bc(self, particle, **kw):
        opts = dict(
            singleParticle=True,
            nonlimiting_filmDiff=False,
            has_surfDiff=particle.has_surfDiff,
            has_binding=True,
            req_binding=False,
            has_mult_bnd_states=False,
        )
        opts.update(kw)
        return eq.particle_boundary(particle, **opts)

    def test_surface_flux_carries_main_pore_volume_fraction(self):
        result = self._bc(_particle("1D"))
        assert r"d^{\mathrm{mp}} \varepsilon^{\mathrm{mp}}" in result
        assert r"k^{\mathrm{f}}_{i} \left( c^{\b}_{i}" in result

    def test_only_main_pore_has_boundary_conditions(self):
        result = self._bc(_particle("1D", n_cavities=2))
        assert r"c^{\p}_{1,i}" not in result
        assert r"c^{\p}_{2,i}" not in result

    def test_homogeneous_particle_has_no_boundary_conditions(self):
        assert self._bc(_particle("0D")) == ""

    def test_core_shifts_inner_boundary(self):
        result = self._bc(_particle("1D", has_core=True))
        assert r"_{r=R^{\mathrm{c}}}" in result


class TestSideCavityCoupling:
    def test_bulk_couples_to_main_pore_network(self):
        particle = _particle("1D")
        term = eq.int_filmDiff_term(particle, 1, 1, True, False, False)
        assert r"d^{\mathrm{mp}}" in term
        assert r"c^{\p}_{0,i}" in term

    def test_bulk_coupling_unchanged_without_cavities(self):
        particle = _particle("1D", n_cavities=0, models=(), reqs=(), has_side_cavities=False)
        term = eq.int_filmDiff_term(particle, 1, 1, True, False, False)
        assert r"d^{\mathrm{mp}}" not in term
        assert r"c^{\p}_{i}" in term

    def test_multiple_particle_types_carry_both_fractions(self):
        particle = _particle("1D")
        term = eq.int_filmDiff_term(particle, 1, r"N^{\mathrm{p}}", False, False, False)
        assert r"d_{j} d^{\mathrm{mp}}_{j}" in term
        assert r"c^{\p}_{j,0,i}" in term


class TestSideCavityParameters:
    @pytest.mark.parametrize(
        "symbol",
        [
            r"N^{\mathrm{sc}}",
            r"\varepsilon^{\mathrm{mp}}",
            r"\varepsilon^{\mathrm{sc}}_{m}",
            r"d^{\mathrm{mp}}",
            r"d^{\mathrm{sc}}_{m}",
            r"k^{\mathrm{sc}}_{m,i}",
        ],
    )
    def test_parameter_is_listed(self, symbol):
        symbols = [entry["Symbol"] for entry in _particle("1D").vars_and_params]
        assert symbol in symbols

    def test_cavity_exchange_is_a_first_order_rate(self):
        """Remark (PCa): k^sc is volumetric, so 1/s rather than m/s."""
        entry = next(e for e in _particle("1D").vars_and_params if e["Symbol"] == r"k^{\mathrm{sc}}_{m,i}")
        assert entry["Unit"] == r"\frac{1}{s}"

    def test_parameters_absent_without_cavities(self):
        particle = _particle("1D", n_cavities=0, models=(), reqs=(), has_side_cavities=False)
        symbols = [entry["Symbol"] for entry in particle.vars_and_params]
        assert r"k^{\mathrm{sc}}_{m,i}" not in symbols
