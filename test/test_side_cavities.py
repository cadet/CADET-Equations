"""
Tests for the side-cavity particle model (Model PC of the modeling paper).

The pore phase is partitioned into a main pore network (compartment index 0)
and N^sc side-cavity types (m = 1, ..., N^sc). Only the main pore network is
radially resolved and only it touches the particle surface; the cavities are
reached through cavity exchange.
"""

from typing import ClassVar

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
        assert r"\frac{\partial c^{\p}_{0,i}}{\partial t}" in result
        assert (
            r"\frac{3}{R^{\mathrm{p}} d^{\mathrm{mp}} \varepsilon^{\mathrm{mp}}} "
            r"k^{\mathrm{f}}_{i} \left( c^{\b}_{i} - c^{\p}_{0,i} \right)" in result
        )
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

    def test_isotherm_arguments_are_compartment_local(self):
        """A compartment's isotherm only sees the concentrations of that compartment."""
        result = _transport(_particle("1D"))
        for m in (0, 1, 2):
            args = rf"\left( \vec{{c}}^{{\p}}_{{{m}}}, \vec{{c}}^{{\s}}_{{{m}}} \right)"
            assert rf"f^{{\mathrm{{bind}}}}_{{{m},i}}{args}" in result

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


class TestSideCavityWithoutBinding:
    """Without binding there is no solid phase anywhere in the pore space."""

    @pytest.mark.parametrize("resolution", ["1D", "0D"])
    def test_no_solid_equations(self, resolution):
        particle = _particle(resolution, has_binding=False)
        result = _transport(particle, has_binding=False)
        assert r"c^{\s}" not in result
        assert r"f^{\mathrm{bind}}" not in result

    @pytest.mark.parametrize("resolution", ["1D", "0D"])
    def test_cavities_still_exchange(self, resolution):
        particle = _particle(resolution, has_binding=False)
        result = _transport(particle, has_binding=False)
        for m in (1, 2):
            assert rf"k^{{\mathrm{{sc}}}}_{{{m},i}}" in result
            assert rf"\frac{{\partial c^{{\p}}_{{{m},i}}}}{{\partial t}}" in result

    def test_no_porosity_weighting_without_a_solid_phase(self):
        result = _transport(_particle("1D", has_binding=False), has_binding=False)
        # the (1 - eps)/eps prefactor only multiplies a binding term
        assert r"\frac{1 - \varepsilon^{\mathrm{mp}}}{\varepsilon^{\mathrm{mp}}}" not in result
        assert r"\frac{1 - \varepsilon^{\mathrm{sc}}_{1}}{\varepsilon^{\mathrm{sc}}_{1}}" not in result

    def test_solid_phase_is_listed_only_with_binding(self):
        without = [e["Symbol"] for e in _particle("1D", has_binding=False).vars_and_params]
        with_ = [e["Symbol"] for e in _particle("1D").vars_and_params]
        assert r"c^{\mathrm{s}}_{i}" not in without
        assert r"c^{\mathrm{s}}_{i}" in with_


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
        """k^f is per total external particle surface, so only the index changes."""
        particle = _particle("1D")
        term = eq.int_filmDiff_term(particle, 1, 1, True, False, False)
        assert r"c^{\p}_{0,i}" in term
        assert r"\frac{3}{R^{\mathrm{p}}}" in term
        assert r"d^{\mathrm{mp}}" not in term

    def test_bulk_coupling_unchanged_without_cavities(self):
        particle = _particle("1D", n_cavities=0, models=(), reqs=(), has_side_cavities=False)
        term = eq.int_filmDiff_term(particle, 1, 1, True, False, False)
        assert r"d^{\mathrm{mp}}" not in term
        assert r"c^{\p}_{i}" in term

    def test_multiple_particle_types_carry_only_the_type_fraction(self):
        particle = _particle("1D")
        term = eq.int_filmDiff_term(particle, 1, r"N^{\mathrm{p}}", False, False, False)
        assert r"3d_{j}" in term
        assert r"d^{\mathrm{mp}}_{j}" not in term
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


class TestSideCavityAppConfiguration:
    """End-to-end checks on configurations that the app must handle gracefully."""

    @staticmethod
    def _app(**session_state):
        from streamlit.testing.v1 import AppTest

        at = AppTest.from_file("../Equation-Generator.py")
        at.run()
        assert not at.exception, at.exception
        # mirror the app's own configuration upload, which writes straight
        # into the session state rather than driving the widgets
        for key, value in session_state.items():
            at.session_state[key] = value
        at.run()
        assert not at.exception, at.exception
        return at

    BASE: ClassVar[dict] = {
        "dev_mode": True,
        "column_type": "Axial flow cylinder",
        "column_resolution": "1D (axial coordinate)",
        "particle_resolution": "0D (homogeneous)",
        r"N^\mathrm{p}": 1,
        "has_binding": "No",
        "particle_has_side_cavities": "Yes",
        "particle_n_side_cavities": 2,
        "particle_nonlimiting_filmDiff": "Kinetic",
    }

    def test_homogeneous_particle_without_binding_keeps_its_cavities(self):
        at = self._app(**self.BASE)
        latex = at.session_state.latex_string
        assert r"k^{\mathrm{sc}}_{1,i}" in latex
        assert r"k^{\mathrm{sc}}_{2,i}" in latex
        assert r"f^{\mathrm{bind}}" not in latex
        assert r"c^{\mathrm{s}}" not in latex

    def test_rapid_equilibrium_film_diffusion_falls_back_with_a_warning(self):
        at = self._app(**{**self.BASE, "particle_nonlimiting_filmDiff": "Rapid-equilibrium"})
        assert any("radially homogeneous side-cavity particle" in warning.value for warning in at.warning)
        latex = at.session_state.latex_string
        # the particle survives and keeps the kinetic film diffusion coupling
        assert r"k^{\mathrm{f}}_{i}" in latex
        assert r"k^{\mathrm{sc}}_{1,i}" in latex
        # and the bulk is not also collapsed into the without-pores notation
        assert r"c^{\mathrm{\ell}}" not in latex


class TestSideCavityRapidEquilibriumFilmDiffusion:
    """Rapid equilibrium across the film is orthogonal to the cavities in 1D."""

    def test_radial_main_pore_takes_the_bulk_concentration(self):
        particle = _particle("1D", nonlimiting_filmDiff=True)
        result = eq.particle_boundary(
            particle,
            singleParticle=True,
            nonlimiting_filmDiff=True,
            has_surfDiff=False,
            has_binding=True,
            req_binding=False,
            has_mult_bnd_states=False,
        )
        assert r"\left. c^{\p}_{0,i} \right|_{r = R^{\mathrm{p}}} &= c^{\b}_{i}" in result
        assert r"k^{\mathrm{f}}" not in result

    def test_cavities_keep_their_own_exchange(self):
        """The cavities never touch the surface, so they are unaffected."""
        particle = _particle("1D", nonlimiting_filmDiff=True)
        result = eq.particle_transport(
            particle,
            singleParticle=True,
            nonlimiting_filmDiff=True,
            has_surfDiff=False,
            has_binding=True,
            req_binding=False,
            has_mult_bnd_states=False,
        )
        assert r"k^{\mathrm{sc}}_{1,i}" in result
        assert r"k^{\mathrm{sc}}_{2,i}" in result

    def test_bulk_substitutes_the_main_pore_flux(self):
        """d^mp sits in the prefactor, so only the main pore porosity is substituted."""
        particle = _particle("1D", nonlimiting_filmDiff=True)
        term = eq.int_filmDiff_term(particle, 1, 1, True, True, False)
        assert r"\frac{3}{R^{\mathrm{p}}}" in term
        assert r"d^{\mathrm{mp}} \varepsilon^{\mathrm{mp}} D^{\mathrm{p}}_{i}" in term
        assert r"c^{\p}_{0,i}" in term
        assert r"k^{\mathrm{f}}" not in term


class TestSideCavityModelName:
    @staticmethod
    def _name(**session_state):
        from streamlit.testing.v1 import AppTest

        at = AppTest.from_file("../Equation-Generator.py")
        at.run()
        for key, value in session_state.items():
            at.session_state[key] = value
        at.run()
        assert not at.exception, at.exception
        latex = at.session_state.latex_string
        return latex.split(r"\section*{")[1].split("}")[0]

    BASE: ClassVar[dict] = {
        "dev_mode": True,
        "column_type": "Axial flow cylinder",
        "column_resolution": "1D (axial coordinate)",
        r"N^\mathrm{p}": 1,
        "has_binding": "Yes",
        "particle_has_side_cavities": "Yes",
        "particle_n_side_cavities": 2,
        "particle_nonlimiting_filmDiff": "Kinetic",
    }

    def test_radial_particle(self):
        name = self._name(**self.BASE, particle_resolution="1D (radial coordinate)")
        assert name == "General Rate Model with side cavities"

    def test_homogeneous_particle_reads_as_a_list(self):
        """The name already carries a "with", so the second clause uses "and"."""
        name = self._name(**self.BASE, particle_resolution="0D (homogeneous)")
        assert name == "Lumped Rate Model with Pores and side cavities"

    def test_absent_without_cavities(self):
        name = self._name(
            **{**self.BASE, "particle_has_side_cavities": "No"},
            particle_resolution="1D (radial coordinate)",
        )
        assert "side cavities" not in name


class TestSideCavitySolverAvailability:
    """The side-cavity model has no implementation in any of the solvers."""

    BASE: ClassVar[dict] = {
        "dev_mode": True,
        "column_type": "Axial flow cylinder",
        "column_resolution": "1D (axial coordinate)",
        r"N^\mathrm{p}": 1,
        "has_binding": "Yes",
        "particle_resolution": "1D (radial coordinate)",
        "particle_nonlimiting_filmDiff": "Kinetic",
        "particle_n_side_cavities": 2,
    }

    @staticmethod
    def _badges(**session_state):
        from streamlit.testing.v1 import AppTest

        at = AppTest.from_file("../Equation-Generator.py")
        at.run()
        for key, value in session_state.items():
            at.session_state[key] = value
        at.run()
        assert not at.exception, at.exception
        blob = " ".join(str(block.value) for block in at.markdown if "CADET-" in str(block.value))
        return {
            tool: "not supported" in blob[blob.index(tool) :][:400]
            for tool in ("CADET-Core", "CADET-Process", "CADET-Semi-Analytic")
        }

    def test_no_solver_supports_side_cavities(self):
        badges = self._badges(**self.BASE, particle_has_side_cavities="Yes")
        assert all(badges.values()), badges

    def test_the_same_model_without_cavities_is_supported(self):
        badges = self._badges(**self.BASE, particle_has_side_cavities="No")
        assert not any(badges.values()), badges


class TestSideCavityWidgetLabel:
    def test_label_states_what_the_model_enables(self):
        from streamlit.testing.v1 import AppTest

        at = AppTest.from_file("../Equation-Generator.py")
        at.run()
        at.session_state["dev_mode"] = True
        at.session_state[r"N^\mathrm{p}"] = 1
        at.run()
        label = next(box.label for box in at.selectbox if box.key == "particle_has_side_cavities")
        assert (
            label
            == "Add side-cavities (enables simultaneous component-specific pore accessibility and competitive binding)"
        )
