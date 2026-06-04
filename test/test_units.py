# -*- coding: utf-8 -*-
"""
Unit tests for the unit system module.
"""

import pytest
from src.units import (
    get_unit, get_conversion_factor, format_conversion_factor,
    UNIT_SYSTEMS, AVAILABLE_SYSTEMS, CONVERSION_FACTORS,
)


@pytest.mark.ci
@pytest.mark.unit_test
class TestGetUnit:

    def test_si_time(self):
        assert get_unit("time", "SI") == r"s"

    def test_si_length(self):
        assert get_unit("length", "SI") == r"m"

    def test_si_concentration_molar(self):
        assert get_unit("concentration_molar", "SI") == r"\frac{mol}{m^3}"

    def test_si_dimensionless(self):
        assert get_unit("dimensionless", "SI") == r"-"

    def test_si_diffusion(self):
        assert get_unit("diffusion", "SI") == r"\frac{m^2}{s}"

    def test_si_velocity(self):
        assert get_unit("velocity", "SI") == r"\frac{m}{s}"

    def test_si_volumetric_flow(self):
        assert get_unit("volumetric_flow", "SI") == r"\frac{m^3}{s}"

    def test_si_volume(self):
        assert get_unit("volume", "SI") == r"m^3"

    def test_si_concentration_mass(self):
        assert get_unit("concentration_mass", "SI") == r"\frac{kg}{m^3}"

    def test_si_reaction_rate_molar(self):
        assert get_unit("reaction_rate_molar", "SI") == r"\frac{mol}{m^3 \cdot s}"

    def test_si_rate_first_order(self):
        assert get_unit("rate_first_order", "SI") == r"\frac{1}{s}"

    def test_si_rate_second_order(self):
        assert get_unit("rate_second_order", "SI") == r"\frac{m^3}{mol \cdot s}"

    def test_si_rate_nth_order(self):
        assert get_unit("rate_nth_order", "SI") == r"\frac{1}{s} \cdot \left(\frac{m^3}{mol}\right)^{n}"

    def test_si_nucleation_rate(self):
        assert get_unit("nucleation_rate", "SI") == r"\frac{1}{m^3 \cdot s}"

    def test_si_number_density(self):
        assert get_unit("number_density", "SI") == r"\frac{1}{m \cdot m^3}"

    def test_si_inverse_length(self):
        assert get_unit("inverse_length", "SI") == r"\frac{1}{m}"

    def test_default_system_is_si(self):
        assert get_unit("time") == get_unit("time", "SI")

    def test_unknown_system_raises(self):
        with pytest.raises(ValueError, match="Unknown unit system"):
            get_unit("time", "INVALID")

    def test_unknown_key_raises(self):
        with pytest.raises(KeyError, match="Unknown unit key"):
            get_unit("nonexistent_unit", "SI")


@pytest.mark.ci
@pytest.mark.unit_test
class TestCGSSystem:

    def test_cgs_time(self):
        assert get_unit("time", "CGS") == r"s"

    def test_cgs_length(self):
        assert get_unit("length", "CGS") == r"cm"

    def test_cgs_volume(self):
        assert get_unit("volume", "CGS") == r"cm^3"

    def test_cgs_concentration_molar(self):
        assert get_unit("concentration_molar", "CGS") == r"\frac{mol}{cm^3}"

    def test_cgs_concentration_mass(self):
        assert get_unit("concentration_mass", "CGS") == r"\frac{g}{cm^3}"

    def test_cgs_velocity(self):
        assert get_unit("velocity", "CGS") == r"\frac{cm}{s}"

    def test_cgs_diffusion(self):
        assert get_unit("diffusion", "CGS") == r"\frac{cm^2}{s}"


@pytest.mark.ci
@pytest.mark.unit_test
class TestPracticalSystem:

    def test_practical_time(self):
        assert get_unit("time", "Practical") == r"min"

    def test_practical_length(self):
        assert get_unit("length", "Practical") == r"cm"

    def test_practical_volume(self):
        assert get_unit("volume", "Practical") == r"mL"

    def test_practical_concentration_molar(self):
        assert get_unit("concentration_molar", "Practical") == r"mM"

    def test_practical_concentration_mass(self):
        assert get_unit("concentration_mass", "Practical") == r"\frac{g}{L}"

    def test_practical_velocity(self):
        assert get_unit("velocity", "Practical") == r"\frac{cm}{min}"

    def test_practical_volumetric_flow(self):
        assert get_unit("volumetric_flow", "Practical") == r"\frac{mL}{min}"

    def test_practical_reaction_rate(self):
        assert get_unit("reaction_rate_molar", "Practical") == r"\frac{mM}{min}"


@pytest.mark.ci
@pytest.mark.unit_test
class TestConversionFactors:

    def test_si_to_si_is_unity(self):
        for key in UNIT_SYSTEMS["SI"]:
            assert get_conversion_factor(key, "SI") == 1.0

    def test_cgs_length(self):
        assert get_conversion_factor("length", "CGS") == pytest.approx(1e2)

    def test_cgs_volume(self):
        assert get_conversion_factor("volume", "CGS") == pytest.approx(1e6)

    def test_cgs_concentration_molar(self):
        assert get_conversion_factor("concentration_molar", "CGS") == pytest.approx(1e-6)

    def test_cgs_concentration_mass(self):
        assert get_conversion_factor("concentration_mass", "CGS") == pytest.approx(1e-3)

    def test_cgs_diffusion(self):
        assert get_conversion_factor("diffusion", "CGS") == pytest.approx(1e4)

    def test_cgs_rate_first_order_unchanged(self):
        assert get_conversion_factor("rate_first_order", "CGS") == pytest.approx(1.0)

    def test_cgs_rate_nth_order_is_none(self):
        assert get_conversion_factor("rate_nth_order", "CGS") is None

    def test_practical_time(self):
        assert get_conversion_factor("time", "Practical") == pytest.approx(1.0 / 60)

    def test_practical_concentration_molar_unity(self):
        assert get_conversion_factor("concentration_molar", "Practical") == pytest.approx(1.0)

    def test_practical_concentration_mass_unity(self):
        assert get_conversion_factor("concentration_mass", "Practical") == pytest.approx(1.0)

    def test_practical_velocity(self):
        assert get_conversion_factor("velocity", "Practical") == pytest.approx(6e3)

    def test_practical_rate_first_order(self):
        assert get_conversion_factor("rate_first_order", "Practical") == pytest.approx(60.0)

    def test_practical_volumetric_flow(self):
        assert get_conversion_factor("volumetric_flow", "Practical") == pytest.approx(6e7)

    def test_unknown_system_raises(self):
        with pytest.raises(ValueError, match="Unknown unit system"):
            get_conversion_factor("time", "INVALID")

    def test_unknown_key_raises(self):
        with pytest.raises(KeyError, match="Unknown unit key"):
            get_conversion_factor("nonexistent", "CGS")

    def test_all_systems_have_factors_for_all_keys(self):
        si_keys = set(UNIT_SYSTEMS["SI"].keys())
        for name in CONVERSION_FACTORS:
            assert set(CONVERSION_FACTORS[name].keys()) == si_keys, (
                f"System '{name}' is missing conversion factors"
            )


@pytest.mark.ci
@pytest.mark.unit_test
class TestFormatConversionFactor:

    def test_unity(self):
        assert format_conversion_factor(1.0) == "1"

    def test_integer(self):
        assert format_conversion_factor(60.0) == "60"

    def test_scientific(self):
        result = format_conversion_factor(1e-6)
        assert "1e-06" in result or "1E-06" in result or "1e-6" in result

    def test_none(self):
        result = format_conversion_factor(None)
        assert "n" in result


@pytest.mark.ci
@pytest.mark.unit_test
class TestRegistry:

    def test_si_in_available(self):
        assert "SI" in AVAILABLE_SYSTEMS

    def test_cgs_in_available(self):
        assert "CGS" in AVAILABLE_SYSTEMS

    def test_practical_in_available(self):
        assert "Practical" in AVAILABLE_SYSTEMS

    def test_all_systems_have_same_keys(self):
        si_keys = set(UNIT_SYSTEMS["SI"].keys())
        for name, system in UNIT_SYSTEMS.items():
            assert set(system.keys()) == si_keys, (
                f"System '{name}' has different keys than SI"
            )
