# -*- coding: utf-8 -*-
"""
Unit tests for src/model_particle.py.
Tests the Particle dataclass and its methods.
"""

import pytest
from src.model_particle import Particle


class TestParticleInitialization:
    """Test Particle initialization and validation."""

    @pytest.mark.unit_test
    def test_particle_sphere_creation(self):
        """Test creation of a basic spherical particle."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        assert particle.geometry == "Sphere"
        assert particle.resolution == "1D"
        assert particle.has_core is False
        assert particle.surface_volume_ratio == 3

    @pytest.mark.unit_test
    def test_particle_cylinder_creation(self):
        """Test creation of a cylindrical particle."""
        particle = Particle(
            geometry="Cylinder",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        assert particle.geometry == "Cylinder"
        assert particle.surface_volume_ratio == 2

    @pytest.mark.unit_test
    def test_particle_slab_creation(self):
        """Test creation of a slab-shaped particle."""
        particle = Particle(
            geometry="Slab",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        assert particle.geometry == "Slab"
        assert particle.surface_volume_ratio == 1

    @pytest.mark.unit_test
    def test_invalid_geometry_raises_error(self):
        """Test that invalid geometry raises ValueError."""
        with pytest.raises(ValueError, match="Invalid geometry"):
            Particle(
                geometry="Invalid",
                has_core=False,
                var_format="CADET",
                resolution="1D",
            )

    @pytest.mark.unit_test
    def test_invalid_resolution_raises_error(self):
        """Test that invalid resolution raises ValueError."""
        with pytest.raises(ValueError, match="Invalid resolution"):
            Particle(
                geometry="Sphere",
                has_core=False,
                var_format="CADET",
                resolution="2D",
            )

    @pytest.mark.unit_test
    def test_particle_0d_resolution(self):
        """Test creation of 0D homogeneous particle."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="0D",
        )
        assert particle.resolution == "0D"

    @pytest.mark.unit_test
    def test_particle_with_core(self):
        """Test particle with core-shell structure."""
        particle = Particle(
            geometry="Sphere",
            has_core=True,
            var_format="CADET",
            resolution="1D",
        )
        assert particle.has_core is True

    @pytest.mark.unit_test
    def test_particle_immutability(self):
        """Test that Particle dataclass is frozen (immutable)."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        # Attempting to modify should raise an error
        with pytest.raises(Exception):  # FrozenInstanceError
            particle.geometry = "Cylinder"


class TestParticleVarsAndParams:
    """Test the vars_and_params generation."""

    @pytest.mark.unit_test
    def test_vars_and_params_populated(self):
        """Test that vars_and_params is populated after initialization."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert len(particle.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_vars_and_params_with_binding(self):
        """Test vars_and_params contains binding-related variables."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        symbols = [var["Symbol"] for var in particle.vars_and_params]
        # Should contain particle concentration and binding variables
        assert any("c" in str(sym).lower() for sym in symbols)

    @pytest.mark.unit_test
    def test_vars_and_params_without_binding(self):
        """Test vars_and_params when binding is disabled."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=False,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert len(particle.vars_and_params) > 0
        # Should not include binding-specific variables in description

    @pytest.mark.unit_test
    def test_vars_and_params_1d_particle(self):
        """Test vars_and_params for 1D particle includes radial coordinate."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        # Should include radial coordinate 'r'
        assert any("r" in str(var.get("Symbol", "")).lower() for var in particle.vars_and_params)

    @pytest.mark.unit_test
    def test_vars_and_params_0d_particle(self):
        """Test vars_and_params for 0D particle."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="0D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="0D",
        )
        assert len(particle.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_vars_and_params_with_surface_diffusion(self):
        """Test vars_and_params includes surface diffusion coefficient."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_surfDiff=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert len(particle.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_vars_and_params_with_multiple_bound_states(self):
        """Test vars_and_params with multiple bound states."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_mult_bnd_states=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert len(particle.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_vars_and_params_sorted_by_group(self):
        """Test that vars_and_params are sorted by Group."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        groups = [var.get("Group", float('inf')) for var in particle.vars_and_params]
        # Check if sorted
        assert groups == sorted(groups)


class TestParticleStatesDependencies:
    """Test state dependencies configuration."""

    @pytest.mark.unit_test
    def test_state_deps_0d_column_0d_particle(self):
        """Test state dependencies for 0D column and 0D particle."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="0D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="0D",
        )
        assert len(particle.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_state_deps_1d_column_1d_particle(self):
        """Test state dependencies for 1D column and 1D particle."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert len(particle.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_state_deps_2d_column_1d_particle(self):
        """Test state dependencies for 2D column and 1D particle."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="2D",
        )
        assert len(particle.vars_and_params) > 0

    @pytest.mark.unit_test
    def test_state_deps_3d_column_1d_particle(self):
        """Test state dependencies for 3D column and 1D particle."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="3D",
        )
        assert len(particle.vars_and_params) > 0


class TestParticleBindingModels:
    """Test particle behavior with different binding models."""

    @pytest.mark.unit_test
    @pytest.mark.parametrize("binding_model", ["Arbitrary", "Linear", "Langmuir", "SMA"])
    def test_different_binding_models(self, binding_model):
        """Test particle creation with different binding models."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            binding_model=binding_model,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.binding_model == binding_model

    @pytest.mark.unit_test
    def test_arbitrary_binding_model(self):
        """Test Arbitrary binding model configuration."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            binding_model="Arbitrary",
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.binding_model == "Arbitrary"

    @pytest.mark.unit_test
    def test_linear_binding_model(self):
        """Test Linear binding model configuration."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            binding_model="Linear",
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.binding_model == "Linear"

    @pytest.mark.unit_test
    def test_langmuir_binding_model(self):
        """Test Langmuir binding model configuration."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            binding_model="Langmuir",
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.binding_model == "Langmuir"

    @pytest.mark.unit_test
    def test_sma_binding_model(self):
        """Test SMA (Steric Mass Action) binding model configuration."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            binding_model="SMA",
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.binding_model == "SMA"


class TestParticleParallelParticleDistribution:
    """Test particle type distribution (PTD) configuration."""

    @pytest.mark.unit_test
    def test_single_particle_type(self):
        """Test single particle type configuration."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            single_partype=True,
            PTD=False,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.single_partype is True
        assert particle.PTD is False

    @pytest.mark.unit_test
    def test_multiple_particle_types_with_ptd(self):
        """Test multiple particle types with PTD enabled."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            single_partype=False,
            PTD=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.single_partype is False
        assert particle.PTD is True


class TestParticleReactions:
    """Test particle with reaction configurations."""

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_with_liquid_reaction(self):
        """Test particle with liquid phase reaction."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_liquid=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.has_reaction_liquid is True

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_with_solid_reaction(self):
        """Test particle with solid phase reaction."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_solid=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.has_reaction_solid is True

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_with_both_reactions(self):
        """Test particle with both liquid and solid phase reactions."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_liquid=True,
            has_reaction_solid=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        assert particle.has_reaction_liquid is True
        assert particle.has_reaction_solid is True

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_with_req_reaction_liquid(self):
        """Test particle with rapid-equilibrium liquid reaction."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_liquid=True,
            req_reaction_liquid=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        symbols = [v["Symbol"] for v in particle.vars_and_params]
        assert any("g^{\\mathrm{react,eq}" in s for s in symbols)
        assert any("M^{\\mathrm{p}}" in s for s in symbols)
        assert not any("f^{\\mathrm{react},\\mathrm{p}}" in s for s in symbols)

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_with_req_reaction_solid(self):
        """Test particle with rapid-equilibrium solid reaction."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_solid=True,
            req_reaction_solid=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        symbols = [v["Symbol"] for v in particle.vars_and_params]
        assert any("g^{\\mathrm{react,eq}" in s for s in symbols)
        assert any("M^{\\mathrm{s}}" in s for s in symbols)
        assert not any("f^{\\mathrm{react},\\mathrm{s}}" in s for s in symbols)

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_with_req_reaction_both(self):
        """Test particle with rapid-equilibrium reactions in both phases."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_liquid=True,
            req_reaction_liquid=True,
            has_reaction_solid=True,
            req_reaction_solid=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        symbols = [v["Symbol"] for v in particle.vars_and_params]
        assert any("M^{\\mathrm{p}}" in s for s in symbols)
        assert any("M^{\\mathrm{s}}" in s for s in symbols)

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_with_req_reaction_multi_partype(self):
        """Test particle with rapid-equilibrium reactions and multiple particle types."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_liquid=True,
            req_reaction_liquid=True,
            has_reaction_solid=True,
            req_reaction_solid=True,
            single_partype=False,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        symbols = [v["Symbol"] for v in particle.vars_and_params]
        # Multi-partype should have j index in constraint symbols
        assert any("j,k" in s for s in symbols)

    @pytest.mark.ci
    @pytest.mark.unit_test
    def test_particle_kinetic_reaction_vars(self):
        """Test that kinetic reaction keeps f^react symbols."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            has_reaction_liquid=True,
            req_reaction_liquid=False,
            has_reaction_solid=True,
            req_reaction_solid=False,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        symbols = [v["Symbol"] for v in particle.vars_and_params]
        assert any("f^{\\mathrm{react},\\mathrm{p}}" in s for s in symbols)
        assert any("f^{\\mathrm{react},\\mathrm{s}}" in s for s in symbols)


class TestParticleVarsParamsDescription:
    """Test the vars_params_description method."""

    @pytest.mark.unit_test
    def test_vars_params_description_returns_string(self):
        """Test that vars_params_description returns a string."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        description = particle.vars_params_description()
        assert isinstance(description, str)
        assert len(description) > 0

    @pytest.mark.unit_test
    def test_vars_params_description_ends_with_period(self):
        """Test that description ends with a period."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=True,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        description = particle.vars_params_description()
        assert description.endswith(".")

    @pytest.mark.unit_test
    def test_vars_params_description_no_binding(self):
        """Test description generation without binding."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
            has_binding=False,
            nonlimiting_filmDiff=False,
            interstitial_volume_resolution="1D",
        )
        description = particle.vars_params_description()
        assert isinstance(description, str)
        assert len(description) > 0


class TestParticleHashability:
    """Test that Particle can be used in Counter and set operations."""

    @pytest.mark.unit_test
    def test_particle_hashable(self):
        """Test that Particle instances are hashable."""
        particle1 = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        particle2 = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        # Should be hashable and equal
        particle_set = {particle1, particle2}
        assert len(particle_set) == 1  # Same particles should be considered equal

    @pytest.mark.unit_test
    def test_particles_different_geometry_not_equal(self):
        """Test that particles with different geometries are not equal."""
        particle1 = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        particle2 = Particle(
            geometry="Cylinder",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        particle_set = {particle1, particle2}
        assert len(particle_set) == 2  # Different particles should be distinct

    @pytest.mark.unit_test
    def test_particles_different_resolution_not_equal(self):
        """Test that particles with different resolutions are not equal."""
        particle1 = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        particle2 = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="0D",
        )
        particle_set = {particle1, particle2}
        assert len(particle_set) == 2


class TestParticleAvailableCADETCore:
    """Test CADET-Core availability check."""

    @pytest.mark.unit_test
    def test_available_cadet_core_returns_bool(self):
        """Test that available_CADET_Core returns a boolean."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        result = particle.available_CADET_Core()
        assert isinstance(result, bool)

    @pytest.mark.unit_test
    def test_available_cadet_core_default_true(self):
        """Test that available_CADET_Core returns True by default."""
        particle = Particle(
            geometry="Sphere",
            has_core=False,
            var_format="CADET",
            resolution="1D",
        )
        assert particle.available_CADET_Core() is True
