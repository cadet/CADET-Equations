# -*- coding: utf-8 -*-
"""
This script implements the `Particle` dataclass and related logic.
"""

from dataclasses import dataclass, field
from collections import Counter
from typing import List, Literal, Optional
import re
import streamlit as st

from src import equations as eq
from src.utils import format_variables

# immutable & hashable dataclass -> unique particle type counter
@dataclass(frozen=True)
class Particle:
    """Minimal description of a particle type used by the generator.

    Attributes are intentionally concise; see field names for intent.
    This dataclass is immutable and hashable so it can be used with
    collections. Counter to count distinct particle types.
    """

    geometry: Literal["Sphere", "Cylinder", "Slab"]
    has_core: bool
    var_format: Literal
    resolution: Literal["1D", "0D"]
    has_binding: bool =None # and thus solid phase
    req_binding: bool = None
    has_mult_bnd_states: bool = None
    has_surfDiff: bool = None
    nonlimiting_filmDiff: bool = None
    surface_volume_ratio: float = None
    interstitial_volume_resolution: str = None
    single_partype: bool = True
    PTD: bool = False
    binding_model: str = "Arbitrary"
    has_reaction_liquid: bool = False
    has_reaction_solid: bool = False
    # volume fraction ?
    # binding -> is_kinetic, nBound
    vars_and_params: List[dict] = field(default_factory=list, init=False, compare=False, hash=False)

    def __post_init__(self):
        valid_geometries = {"Sphere", "Cylinder", "Slab"}
        valid_resolutions = {"1D", "0D"}

        if self.geometry not in valid_geometries:
            raise ValueError(
                f"Invalid geometry: {self.geometry}. Must be one of {valid_geometries}.")
        if self.resolution not in valid_resolutions:
            raise ValueError(
                f"Invalid resolution: {self.resolution}. Must be one of {valid_resolutions}.")

        if self.geometry == "Sphere":
            object.__setattr__(self, 'surface_volume_ratio', 3)
        elif self.geometry == "Cylinder":
            object.__setattr__(self, 'surface_volume_ratio', 2)
        elif self.geometry == "Slab":
            object.__setattr__(self, 'surface_volume_ratio', 1)

        vars_and_params_ = []
        
        state_deps = r"t"
        if self.interstitial_volume_resolution == "1D":
            state_deps = r"t, z"
        if self.interstitial_volume_resolution == "2D":
            state_deps = r"t, z, \rho"
        if self.interstitial_volume_resolution == "3D":
            state_deps = r"t, z, \rho, \phi"
        if self.resolution == "1D":
            state_deps += r", r"

        state_deps += "; i" if self.single_partype else "; j, i"

        if not (self.nonlimiting_filmDiff and self.resolution == "0D"):
            symbol_name_ = r"c^{\p}_{i}" if self.single_partype else r"c^{\p}_{j,i}"
            vars_and_params_.append({"Group" : 1, "Symbol": symbol_name_, "Description": r"particle liquid concentration", "Unit": r"\frac{mol}{m^3}", "Dependence" : state_deps, "Domain" : eq.full_particle_conc_domain(column_resolution=self.interstitial_volume_resolution, particle_resolution=self.resolution, hasCore=self.has_core, with_par_index=False, with_time_domain=True)})
            
        if self.resolution == "1D":
            vars_and_params_.append({"Group" : 0, "Symbol": r"r", "Description": r"radial particle coordinate", "Unit": r"m", "Dependence": r"\text{independent variable}", "Property": r""})
            symbol_name_ = r"D^\mathrm{p}_{i}" if self.single_partype else r"D^\mathrm{p}_{j,i}"
            vars_and_params_.append({"Group" : 6.1, "Symbol": symbol_name_, "Description": r"particle diffusion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": r"i" if self.single_partype else r"j,i", "Property": r"> 0"})

            if not self.nonlimiting_filmDiff:
                vars_and_params_.append({"Group" : 4, "Symbol": r"\varepsilon^{\mathrm{p}}" if self.single_partype else r"\varepsilon^{\mathrm{p}}_{j}", "Description": r"particle porosity", "Unit": r"-", "Dependence": r"\text{constant}" if self.single_partype else r"j", "Property": r"\in (0, 1)"})

        if self.has_binding:
            symbol_name_ = r"c^{\s}_{i}" if self.single_partype else r"c^{\s}_{j,i}"
            vars_and_params_.append({"Group" : 1, "Symbol": symbol_name_, "Description": r"particle solid concentration", "Unit": r"\frac{mol}{m^3}", "Dependence" : state_deps, "Domain" : eq.full_particle_conc_domain(column_resolution=self.interstitial_volume_resolution, particle_resolution=self.resolution, hasCore=self.has_core, with_par_index=False, with_time_domain=True)})

            if self.binding_model == "Arbitrary":
                symbol_name_ = r"f^\mathrm{bind}_{j,i}" if self.PTD else r"f^\mathrm{bind}_{i}"
                dep_ = r"\vec{c}^\mathrm{p}, \vec{c}^\mathrm{s}; j, i" if self.PTD else r"\vec{c}^\mathrm{p}, \vec{c}^\mathrm{s}; i"
                vars_and_params_.append({"Group" : 10, "Symbol": symbol_name_, "Description": r"adsorption isotherm function", "Unit": r"\frac{1}{s}", "Dependence": dep_})
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"\vec{c}^\mathrm{p}", "Description": r"particle liquid components vector", "Unit": r"[\frac{mol}{m^3}]", "Dependence": state_deps})
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"\vec{c}^\mathrm{s}", "Description": r"particle solid components vector", "Unit": r"[\frac{mol}{m^3}]", "Dependence": state_deps})
            else:
                idx_ = "i" if self.single_partype else "j,i"
                vars_and_params_.append({"Group" : 10, "Symbol": r"k^{\mathrm{a}}_{" + idx_ + r"}", "Description": r"adsorption rate constant", "Unit": r"\frac{m^3}{mol \cdot s}", "Dependence": idx_})
                vars_and_params_.append({"Group" : 10, "Symbol": r"k^{\mathrm{d}}_{" + idx_ + r"}", "Description": r"desorption rate constant", "Unit": r"\frac{1}{s}", "Dependence": idx_})

            if self.binding_model == "Langmuir":
                idx_ = "i" if self.single_partype else "j,i"
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"q^{\mathrm{max}}_{" + idx_ + r"}", "Description": r"maximum binding capacity", "Unit": r"\frac{mol}{m^3}", "Dependence": idx_})

            if self.binding_model == "SMA":
                idx_ = "i" if self.single_partype else "j,i"
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"\nu_{" + idx_ + r"}", "Description": r"characteristic charge", "Unit": r"-", "Dependence": idx_})
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"\sigma_{" + idx_ + r"}", "Description": r"steric factor", "Unit": r"-", "Dependence": idx_})
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"\Lambda", "Description": r"ionic capacity (binding site concentration)", "Unit": r"\frac{mol}{m^3}", "Dependence": r"\text{constant}"})
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"q^{\mathrm{ref}}", "Description": r"reference solid phase concentration", "Unit": r"\frac{mol}{m^3}", "Dependence": r"\text{constant}"})
                vars_and_params_.append({"Group" : 10.1, "Symbol": r"c^{\mathrm{ref}}", "Description": r"reference liquid phase concentration", "Unit": r"\frac{mol}{m^3}", "Dependence": r"\text{constant}"})
            
            if self.has_surfDiff:
                symbol_name_ = r"D^\mathrm{s}_{i}" if self.single_partype else r"D^\mathrm{s}_{j,i}"
                vars_and_params_.append({"Group" : 6.1, "Symbol": symbol_name_, "Description": r"surface diffusion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": r"i" if self.single_partype else r"j,i", "Property": r"\geq 0"})
            
            if self.has_mult_bnd_states:
                symbol_name_ = r"N^{\mathrm{b}}_{i}" if self.single_partype else r"N^{\mathrm{b}}_{j,i}"
                vars_and_params_.append({"Group" : 11, "Symbol": symbol_name_, "Description": r"number of bound states", "Unit": r"-", "Dependence": r"i" if self.single_partype else r"j,i"})
            
            if self.has_core:
                self.vars_and_params.append({"Group" : 0.1, "Symbol": r"R^\mathrm{pc}", "Description": r"particle core radius", "Unit": r"-", "Dependence": r"-", "Property": r"\in (0, R^\mathrm{p})"})

        if self.has_reaction_liquid:
            symbol_name_ = r"f^{\mathrm{react},\p}_{i}" if self.single_partype else r"f^{\mathrm{react},\p}_{j,i}"
            dep_ = r"\vec{c}^{\p}, \vec{c}^{\s}; i" if self.single_partype else r"\vec{c}^{\p}, \vec{c}^{\s}; j, i"
            vars_and_params_.append({"Group" : 10.2, "Symbol": symbol_name_, "Description": r"particle liquid phase reaction function", "Unit": r"\frac{mol}{m^3 \cdot s}", "Dependence": dep_})

        if self.has_reaction_solid:
            symbol_name_ = r"f^{\mathrm{react},\s}_{i}" if self.single_partype else r"f^{\mathrm{react},\s}_{j,i}"
            dep_ = r"\vec{c}^{\p}, \vec{c}^{\s}; i" if self.single_partype else r"\vec{c}^{\p}, \vec{c}^{\s}; j, i"
            vars_and_params_.append({"Group" : 10.3, "Symbol": symbol_name_, "Description": r"particle solid phase reaction function", "Unit": r"\frac{mol}{m^3 \cdot s}", "Dependence": dep_})

        if not self.single_partype:
            vars_and_params_.append({"Group" : -0.1, "Symbol": r"j", "Description": r"particle type index", "Unit": r"-", "Dependence": r"-", "Property": r""})
            vars_and_params_.append({"Group" : 1.9, "Symbol": r"d_j", "Description": r"particle type volume fraction", "Unit": r"-", "Dependence": r"j", "Property": r""})

        for var_ in vars_and_params_:
            var_["Symbol"] = format_variables(var_["Symbol"], self.var_format)
        
        vars_and_params_ = sorted(vars_and_params_, key=lambda x: x['Group'])

        object.__setattr__(self, 'vars_and_params', vars_and_params_)

    
    def available_CADET_Core(self):
        return True

    def vars_params_description(self):

        description_ = ""

        idx = 1
        num_VP = len(self.vars_and_params)

        for thing in self.vars_and_params:

            if thing.get("Group", -1) < 0:
                num_VP -= 1
                continue

            if not idx == 1:
                description_ += ", " if idx < num_VP else ", and "
                
            description_ += r"$" + thing["Symbol"]

            if not thing.get("Domain", "-") == "-":
                description_ += r"\colon " + re.sub(r"\$", "", thing["Domain"]) + r" \mapsto \mathbb{R}"

            description_ += thing.get("Property", "")

            description_ += r"$"
            
            description_ += " is the " + thing["Description"]

            idx += 1
                
        return description_ + "."