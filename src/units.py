# -*- coding: utf-8 -*-
"""
Unit system definitions and conversion utilities.

Each unit system maps a canonical unit key to a LaTeX string.
Conversion factors are stored relative to SI:
    value_in_target = value_in_SI * CONVERSION_FACTORS[target][key]
"""

from typing import Dict, Optional


UNIT_SYSTEMS: Dict[str, Dict[str, str]] = {
    "SI": {
        "dimensionless":          r"-",
        "time":                   r"s",
        "length":                 r"m",
        "volume":                 r"m^3",
        "velocity":               r"\frac{m}{s}",
        "velocity_coeff":         r"\frac{m^2}{s}",
        "diffusion":              r"\frac{m^2}{s}",
        "volumetric_flow":        r"\frac{m^3}{s}",
        "concentration_molar":    r"\frac{mol}{m^3}",
        "concentration_molar_vec":r"[\frac{mol}{m^3}]",
        "concentration_mass":     r"\frac{kg}{m^3}",
        "reaction_rate_molar":    r"\frac{mol}{m^3 \cdot s}",
        "rate_first_order":       r"\frac{1}{s}",
        "rate_second_order":      r"\frac{m^3}{mol \cdot s}",
        "rate_nth_order":         r"\frac{1}{s} \cdot \left(\frac{m^3}{mol}\right)^{n}",
        "nucleation_rate":        r"\frac{1}{m^3 \cdot s}",
        "number_density":         r"\frac{1}{m \cdot m^3}",
        "inverse_length":         r"\frac{1}{m}",
    },
    "CGS": {
        "dimensionless":          r"-",
        "time":                   r"s",
        "length":                 r"cm",
        "volume":                 r"cm^3",
        "velocity":               r"\frac{cm}{s}",
        "velocity_coeff":         r"\frac{cm^2}{s}",
        "diffusion":              r"\frac{cm^2}{s}",
        "volumetric_flow":        r"\frac{cm^3}{s}",
        "concentration_molar":    r"\frac{mol}{cm^3}",
        "concentration_molar_vec":r"[\frac{mol}{cm^3}]",
        "concentration_mass":     r"\frac{g}{cm^3}",
        "reaction_rate_molar":    r"\frac{mol}{cm^3 \cdot s}",
        "rate_first_order":       r"\frac{1}{s}",
        "rate_second_order":      r"\frac{cm^3}{mol \cdot s}",
        "rate_nth_order":         r"\frac{1}{s} \cdot \left(\frac{cm^3}{mol}\right)^{n}",
        "nucleation_rate":        r"\frac{1}{cm^3 \cdot s}",
        "number_density":         r"\frac{1}{cm \cdot cm^3}",
        "inverse_length":         r"\frac{1}{cm}",
    },
    "Practical": {
        "dimensionless":          r"-",
        "time":                   r"min",
        "length":                 r"cm",
        "volume":                 r"mL",
        "velocity":               r"\frac{cm}{min}",
        "velocity_coeff":         r"\frac{cm^2}{min}",
        "diffusion":              r"\frac{cm^2}{min}",
        "volumetric_flow":        r"\frac{mL}{min}",
        "concentration_molar":    r"mM",
        "concentration_molar_vec":r"[mM]",
        "concentration_mass":     r"\frac{g}{L}",
        "reaction_rate_molar":    r"\frac{mM}{min}",
        "rate_first_order":       r"\frac{1}{min}",
        "rate_second_order":      r"\frac{1}{mM \cdot min}",
        "rate_nth_order":         r"\frac{1}{min} \cdot \left(\frac{1}{mM}\right)^{n}",
        "nucleation_rate":        r"\frac{1}{mL \cdot min}",
        "number_density":         r"\frac{1}{cm \cdot mL}",
        "inverse_length":         r"\frac{1}{cm}",
    },
}

AVAILABLE_SYSTEMS = list(UNIT_SYSTEMS.keys())

# Conversion factors from SI to each system.
# value_in_target = value_in_SI * factor
# None means the factor depends on a model parameter (e.g. reaction order n).
CONVERSION_FACTORS: Dict[str, Dict[str, Optional[float]]] = {
    "SI": {key: 1.0 for key in UNIT_SYSTEMS["SI"]},
    "CGS": {
        "dimensionless":          1.0,
        "time":                   1.0,            # 1 s = 1 s
        "length":                 1e2,             # 1 m = 100 cm
        "volume":                 1e6,             # 1 m³ = 10⁶ cm³
        "velocity":               1e2,             # 1 m/s = 100 cm/s
        "velocity_coeff":         1e4,             # 1 m²/s = 10⁴ cm²/s
        "diffusion":              1e4,             # 1 m²/s = 10⁴ cm²/s
        "volumetric_flow":        1e6,             # 1 m³/s = 10⁶ cm³/s
        "concentration_molar":    1e-6,            # 1 mol/m³ = 10⁻⁶ mol/cm³
        "concentration_molar_vec":1e-6,
        "concentration_mass":     1e-3,            # 1 kg/m³ = 10⁻³ g/cm³
        "reaction_rate_molar":    1e-6,            # 1 mol/(m³·s) = 10⁻⁶ mol/(cm³·s)
        "rate_first_order":       1.0,             # 1/s = 1/s
        "rate_second_order":      1e6,             # 1 m³/(mol·s) = 10⁶ cm³/(mol·s)
        "rate_nth_order":         None,            # depends on reaction order n
        "nucleation_rate":        1e-6,            # 1/(m³·s) = 10⁻⁶/(cm³·s)
        "number_density":         1e-8,            # 1/(m·m³) = 10⁻⁸/(cm·cm³)
        "inverse_length":         1e-2,            # 1/m = 10⁻²/cm
    },
    "Practical": {
        "dimensionless":          1.0,
        "time":                   1.0 / 60,        # 1 s = 1/60 min
        "length":                 1e2,              # 1 m = 100 cm
        "volume":                 1e6,              # 1 m³ = 10⁶ mL
        "velocity":               6e3,              # 1 m/s = 6000 cm/min
        "velocity_coeff":         6e5,              # 1 m²/s = 6×10⁵ cm²/min
        "diffusion":              6e5,              # 1 m²/s = 6×10⁵ cm²/min
        "volumetric_flow":        6e7,              # 1 m³/s = 6×10⁷ mL/min
        "concentration_molar":    1.0,              # 1 mol/m³ = 1 mM
        "concentration_molar_vec":1.0,
        "concentration_mass":     1.0,              # 1 kg/m³ = 1 g/L
        "reaction_rate_molar":    60.0,             # 1 mol/(m³·s) = 60 mM/min
        "rate_first_order":       60.0,             # 1/s = 60/min
        "rate_second_order":      60.0,             # 1 m³/(mol·s) = 60/(mM·min)
        "rate_nth_order":         None,             # depends on reaction order n
        "nucleation_rate":        6e-5,             # 1/(m³·s) = 6×10⁻⁵/(mL·min)
        "number_density":         1e-8,             # 1/(m·m³) = 10⁻⁸/(cm·mL)
        "inverse_length":         1e-2,             # 1/m = 10⁻²/cm
    },
}


def get_unit(unit_key: str, system: str = "SI") -> str:
    """Return the LaTeX string for *unit_key* in the given *system*."""
    if system not in UNIT_SYSTEMS:
        raise ValueError(
            f"Unknown unit system '{system}'. "
            f"Available systems: {AVAILABLE_SYSTEMS}"
        )
    units = UNIT_SYSTEMS[system]
    if unit_key not in units:
        raise KeyError(
            f"Unknown unit key '{unit_key}' in system '{system}'. "
            f"Available keys: {list(units.keys())}"
        )
    return units[unit_key]


def get_conversion_factor(unit_key: str, target: str) -> Optional[float]:
    """Return the factor to convert a value from SI to *target* system.

    Returns None when the factor depends on a model parameter (e.g.
    reaction order *n*).
    """
    if target not in CONVERSION_FACTORS:
        raise ValueError(
            f"Unknown unit system '{target}'. "
            f"Available systems: {AVAILABLE_SYSTEMS}"
        )
    factors = CONVERSION_FACTORS[target]
    if unit_key not in factors:
        raise KeyError(
            f"Unknown unit key '{unit_key}' in system '{target}'."
        )
    return factors[unit_key]


def format_conversion_factor(factor: Optional[float]) -> str:
    """Format a conversion factor as a short LaTeX string for display."""
    if factor is None:
        return r"\text{depends on } n"
    if factor == 1.0:
        return r"1"
    if factor == int(factor) and abs(factor) < 1e4:
        return str(int(factor))
    return f"{factor:.4g}"


__all__ = [
    "UNIT_SYSTEMS",
    "AVAILABLE_SYSTEMS",
    "CONVERSION_FACTORS",
    "get_unit",
    "get_conversion_factor",
    "format_conversion_factor",
]
