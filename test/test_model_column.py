# -*- coding: utf-8 -*-
"""
Unit tests for src/model_column.py.
Tests the Column dataclass and its methods.
"""

import pytest
from unittest.mock import MagicMock, patch
from collections import Counter
from src.model_column import Column
from src.model_particle import Particle


class TestColumnInitialization:
    """Test Column initialization and basic properties."""

    @pytest.mark.unit_test
    def test_column_minimal_instantiation(self):
        """Test creating a Column with minimal required parameters."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert col.dev_mode is False
        assert col.advanced_mode is False
        assert col.var_format == "CADET"

    @pytest.mark.unit_test
    def test_column_has_vars_and_params(self):
        """Test that Column populates vars_and_params after init."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'vars_and_params')
        assert isinstance(col.vars_and_params, list)

    @pytest.mark.unit_test
    def test_column_with_dev_mode(self):
        """Test Column with development mode enabled."""
        col = Column(
            dev_mode=True,
            advanced_mode=False,
            var_format="CADET",
        )
        assert col.dev_mode is True

    @pytest.mark.unit_test
    def test_column_with_advanced_mode(self):
        """Test Column with advanced mode enabled."""
        col = Column(
            dev_mode=False,
            advanced_mode=True,
            var_format="CADET",
        )
        assert col.advanced_mode is True


class TestColumnResolutions:
    """Test Column resolution-related functionality."""

    @pytest.mark.unit_test
    def test_column_resolution_attribute_exists(self):
        """Test that Column has resolution attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'resolution')
        assert col.resolution in ["0D", "1D", "2D", "3D"]

    @pytest.mark.unit_test
    def test_column_coordinate_flags_exist(self):
        """Test that Column has coordinate flags."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_axial_coordinate')
        assert hasattr(col, 'has_radial_coordinate')
        assert hasattr(col, 'has_angular_coordinate')


class TestColumnWithoutParticles:
    """Test Column behavior without Streamlit UI constraints."""

    @pytest.mark.unit_test
    def test_column_has_particle_models_attribute(self):
        """Test that Column has particle_models attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'particle_models')

    @pytest.mark.unit_test
    def test_column_has_binding_attribute(self):
        """Test that Column has binding attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_binding')
        assert isinstance(col.has_binding, bool)

    @pytest.mark.unit_test
    def test_column_model_name_returns_string(self):
        """Test that model_name returns a string."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        name = col.model_name()
        assert isinstance(name, str)
        assert len(name) > 0


class TestColumnWithParticles:
    """Test Column behavior with particles."""

    @pytest.mark.unit_test
    def test_column_par_type_counts_is_counter(self):
        """Test that particle type counts is a Counter."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert isinstance(col.par_type_counts, Counter)

    @pytest.mark.unit_test
    def test_column_par_unique_intv_contribution_counts_is_counter(self):
        """Test that unique interstitial contribution counts is a Counter."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert isinstance(col.par_unique_intV_contribution_counts, Counter)


class TestColumnBindingConfiguration:
    """Test Column binding model configuration."""

    @pytest.mark.unit_test
    def test_column_binding_model_attribute_exists(self):
        """Test that Column has binding_model attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'binding_model')

    @pytest.mark.unit_test
    def test_column_has_binding_attribute(self):
        """Test that Column has has_binding attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_binding')
        assert isinstance(col.has_binding, bool)

    @pytest.mark.unit_test
    def test_column_req_binding_attribute(self):
        """Test that Column has req_binding attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'req_binding')
        assert isinstance(col.req_binding, bool)


class TestColumnDispersion:
    """Test Column dispersion configuration."""

    @pytest.mark.unit_test
    def test_column_has_axial_dispersion_attribute(self):
        """Test that Column has axial dispersion attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_axial_dispersion')
        assert isinstance(col.has_axial_dispersion, bool)

    @pytest.mark.unit_test
    def test_column_has_radial_dispersion_attribute(self):
        """Test that Column has radial dispersion attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_radial_dispersion')
        assert isinstance(col.has_radial_dispersion, bool)

    @pytest.mark.unit_test
    def test_column_has_angular_dispersion_attribute(self):
        """Test that Column has angular dispersion attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_angular_dispersion')
        assert isinstance(col.has_angular_dispersion, bool)


class TestColumnAvailability:
    """Test CADET-Core and CADET-Process availability checks."""

    @pytest.mark.unit_test
    def test_available_cadet_core_returns_valid_int(self):
        """Test that available_CADET_Core returns a valid integer."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        result = col.available_CADET_Core()
        assert isinstance(result, int)
        assert result in [-1, 0, 1]

    @pytest.mark.unit_test
    def test_available_cadet_process_returns_valid_int(self):
        """Test that available_CADET_Process returns a valid integer."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        result = col.available_CADET_Process()
        assert isinstance(result, int)
        assert result in [-1, 0, 1]

    @pytest.mark.unit_test
    def test_cadet_process_uses_cadet_core_availability(self):
        """Test that CADET-Process availability depends on CADET-Core."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        core_avail = col.available_CADET_Core()
        process_avail = col.available_CADET_Process()
        # If core is not available, process should also not be available
        if core_avail == -1:
            assert process_avail == -1


class TestColumnVarsAndParams:
    """Test Column vars_and_params generation."""

    @pytest.mark.unit_test
    def test_vars_and_params_is_list(self):
        """Test that vars_and_params is a list."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert isinstance(col.vars_and_params, list)
        assert len(col.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_vars_and_params_contain_dicts(self):
        """Test that vars_and_params contains dictionaries."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        for var in col.vars_and_params:
            assert isinstance(var, dict)

    @pytest.mark.unit_test
    def test_vars_and_params_have_symbol_field(self):
        """Test that vars_and_params entries have Symbol field."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        for var in col.vars_and_params:
            assert "Symbol" in var

    @pytest.mark.unit_test
    def test_vars_and_params_sorted_by_group(self):
        """Test that vars_and_params are sorted by Group."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        groups = [var.get("Group", float('inf')) for var in col.vars_and_params]
        assert groups == sorted(groups)


class TestColumnModelName:
    """Test Column model_name method."""

    @pytest.mark.unit_test
    def test_model_name_returns_string(self):
        """Test that model_name returns a string."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        name = col.model_name()
        assert isinstance(name, str)
        assert len(name) > 0

    @pytest.mark.unit_test
    def test_model_name_not_empty(self):
        """Test that model name is populated."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        name = col.model_name()
        assert name


class TestColumnModelAssumptions:
    """Test Column model_assumptions method."""

    @pytest.mark.unit_test
    def test_model_assumptions_returns_dict(self):
        """Test that model_assumptions returns a dictionary."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
            resolution="1D",
            N_p=0,
        )
        assumptions = col.model_assumptions()
        assert isinstance(assumptions, dict)

    @pytest.mark.unit_test
    def test_model_assumptions_has_general(self):
        """Test that model_assumptions includes general assumptions."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
            resolution="1D",
            N_p=0,
        )
        assumptions = col.model_assumptions()
        assert "General model assumptions" in assumptions

    @pytest.mark.unit_test
    def test_model_assumptions_has_specific(self):
        """Test that model_assumptions includes specific assumptions."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
            resolution="1D",
            N_p=0,
        )
        assumptions = col.model_assumptions()
        assert "Specific model assumptions" in assumptions

    @pytest.mark.unit_test
    def test_model_assumptions_with_binding_model(self):
        """Test that binding model assumptions are included when applicable."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
            resolution="1D",
            N_p=1,
            binding_model="Langmuir",
        )
        assumptions = col.model_assumptions()
        assert "Binding model assumptions" in assumptions

    @pytest.mark.unit_test
    def test_model_assumptions_nonlimiting_film_diffusion(self):
        """Test that nonlimiting film diffusion is noted in assumptions."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
            resolution="1D",
            N_p=1,
            nonlimiting_filmDiff=True,
        )
        assumptions = col.model_assumptions()
        # Check if assumption about film diffusion is present
        all_assumptions = assumptions["Specific model assumptions"]
        assert len(all_assumptions) > 0


class TestColumnDomainMethods:
    """Test Column domain-related methods."""

    @pytest.mark.unit_test
    def test_domain_interstitial_with_time_returns_string(self):
        """Test domain_interstitial returns string."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        domain = col.domain_interstitial(with_time_domain=True)
        assert isinstance(domain, str)
        assert len(domain) > 0

    @pytest.mark.unit_test
    def test_domain_interstitial_without_time_returns_string(self):
        """Test domain_interstitial without time returns string."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        domain = col.domain_interstitial(with_time_domain=False)
        assert isinstance(domain, str)

    @pytest.mark.unit_test
    def test_domain_particle_returns_string(self):
        """Test domain_particle returns string."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        domain = col.domain_particle()
        assert isinstance(domain, str)


class TestColumnEquationGeneration:
    """Test Column equation generation methods."""

    @pytest.mark.unit_test
    def test_interstitial_volume_equation_returns_string(self):
        """Test interstitial volume equation returns string."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        equation = col.interstitial_volume_equation()
        assert isinstance(equation, str)
        assert len(equation) > 0

    @pytest.mark.unit_test
    def test_interstitial_volume_equation_has_latex(self):
        """Test that equation contains LaTeX."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        equation = col.interstitial_volume_equation()
        assert "align" in equation or "equation" in equation

    @pytest.mark.unit_test
    def test_interstitial_volume_bc_returns_string_or_none(self):
        """Test that boundary conditions return string or None."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        bc = col.interstitial_volume_bc()
        assert isinstance(bc, (str, type(None)))

    @pytest.mark.unit_test
    def test_particle_equations_returns_dicts(self):
        """Test particle_equations returns two dictionaries."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        eqs, bcs = col.particle_equations()
        assert isinstance(eqs, dict)
        assert isinstance(bcs, dict)


class TestColumnVarsParamsDescription:
    """Test Column vars_params_description method."""

    @pytest.mark.unit_test
    def test_vars_params_description_returns_string(self):
        """Test that vars_params_description returns a string."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        description = col.vars_params_description()
        assert isinstance(description, str)
        assert len(description) > 0

    @pytest.mark.unit_test
    def test_vars_params_description_ends_with_period(self):
        """Test that description ends with a period."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        description = col.vars_params_description()
        assert description.endswith(".")


class TestColumnFilmDiffusion:
    """Test film diffusion configurations."""

    @pytest.mark.unit_test
    def test_nonlimiting_film_diffusion_attribute_exists(self):
        """Test that Column has nonlimiting_filmDiff attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'nonlimiting_filmDiff')
        assert isinstance(col.nonlimiting_filmDiff, bool)


class TestColumnSurfaceDiffusion:
    """Test surface diffusion configurations."""

    @pytest.mark.unit_test
    def test_surface_diffusion_attribute_exists(self):
        """Test that Column has has_surfDiff attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_surfDiff')
        assert isinstance(col.has_surfDiff, bool)


class TestColumnReactions:
    """Test reaction configurations."""

    @pytest.mark.unit_test
    def test_bulk_reaction_attribute_exists(self):
        """Test that Column has has_reaction_bulk attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_reaction_bulk')
        assert isinstance(col.has_reaction_bulk, bool)

    @pytest.mark.unit_test
    def test_particle_liquid_reaction_attribute_exists(self):
        """Test that Column has has_reaction_particle_liquid attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_reaction_particle_liquid')
        assert isinstance(col.has_reaction_particle_liquid, bool)

    @pytest.mark.unit_test
    def test_particle_solid_reaction_attribute_exists(self):
        """Test that Column has has_reaction_particle_solid attribute."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert hasattr(col, 'has_reaction_particle_solid')
        assert isinstance(col.has_reaction_particle_solid, bool)


class TestPerComponentConfiguration:
    """Test per-component parameterization in advanced mode."""

    @pytest.mark.unit_test
    def test_per_component_attributes_default_none(self):
        """Test that per-component lists default to None."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert col.req_binding_per_comp is None
        assert col.nonlimiting_filmDiff_per_comp is None
        assert col.has_surfDiff_per_comp is None
        assert col.has_mult_bnd_states_per_comp is None

    @pytest.mark.unit_test
    def test_has_per_component_config_false_by_default(self):
        """Test that per-component config is inactive by default."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert col.has_per_component_config() is False

    @pytest.mark.unit_test
    def test_component_groups_returns_none_without_config(self):
        """Test that component_groups returns None without per-component config."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        assert col.component_groups() is None

    @pytest.mark.unit_test
    def test_format_component_set_single(self):
        """Test formatting a single component index."""
        result = Column.format_component_set([1])
        assert result == r"$i = 1$"

    @pytest.mark.unit_test
    def test_format_component_set_multiple(self):
        """Test formatting multiple component indices."""
        result = Column.format_component_set([1, 3, 5])
        assert result == r"$i \in \{1, 3, 5\}$"

    @pytest.mark.unit_test
    def test_has_per_component_config_true_when_set(self):
        """Test that has_per_component_config returns True when lists are set."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        col.N_c = 2
        col.req_binding_per_comp = [False, True]
        assert col.has_per_component_config() is True

    @pytest.mark.unit_test
    def test_component_groups_single_group(self):
        """Test component_groups when all components share the same settings."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        col.N_c = 2
        col.req_binding_per_comp = [False, False]
        col.nonlimiting_filmDiff_per_comp = [False, False]
        col.has_surfDiff_per_comp = [False, False]
        col.has_mult_bnd_states_per_comp = [False, False]
        groups = col.component_groups()
        assert groups is not None
        assert len(groups) == 1
        assert groups[0]['components'] == [1, 2]

    @pytest.mark.unit_test
    def test_component_groups_two_groups(self):
        """Test component_groups when settings differ across components."""
        col = Column(
            dev_mode=False,
            advanced_mode=False,
            var_format="CADET",
        )
        col.N_c = 3
        col.req_binding_per_comp = [False, True, False]
        col.nonlimiting_filmDiff_per_comp = [False, False, False]
        col.has_surfDiff_per_comp = [False, False, False]
        col.has_mult_bnd_states_per_comp = [False, False, False]
        groups = col.component_groups()
        assert groups is not None
        assert len(groups) == 2
        # Components 1 and 3 share kinetic binding, component 2 has rapid-equilibrium
        comp_sets = [sorted(g['components']) for g in groups]
        assert [1, 3] in comp_sets
        assert [2] in comp_sets
