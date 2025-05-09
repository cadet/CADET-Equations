# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import streamlit as st
import re
from collections import Counter
from dataclasses import dataclass, field
import json
import subprocess
import tempfile
import pandas as pd

from typing import List
from typing import Literal
from collections import OrderedDict

from src import equations as eq
from src import load_CADET_h5
from src import ui_config


def rerender_variables(input_str:str, var_format:int, mul_pars:bool=False):
    if var_format == "CADET":
        input_str = re.sub(r"\\l(?![a-zA-Z])", r"\\mathrm{\\ell}", input_str)
        input_str = re.sub(r"\\b(?![a-zA-Z])", r"\\mathrm{b}", input_str)
        input_str = re.sub(r"\\p(?![a-zA-Z])", r"\\mathrm{p}", input_str)
        input_str = re.sub(r"\\s(?![a-zA-Z])", r"\\mathrm{s}", input_str)
    elif var_format == "Legacy":
        input_str = re.sub(r"c\^\{\\l\}(?![a-zA-Z])", r"c", input_str)
        input_str = re.sub(r"c\^\{\\b\}(?![a-zA-Z])", r"c", input_str)

        input_str = re.sub(
            r"V\^\{\\l\}(?![a-zA-Z])", r"V^{\\mathrm{\\ell}}", input_str)
        input_str = re.sub(
            r"V\^\{\\b\}(?![a-zA-Z])", r"V^{\\mathrm{\\ell}}", input_str)

        input_str = re.sub(
            r"c\^\{\\p\}_\{i\}(?![a-zA-Z])", r"c_{\\mathrm{p},i}", input_str)
        input_str = re.sub(
            r"c\^\{\\p\}_i(?![a-zA-Z])", r"c_{\\mathrm{p},i}", input_str)
        input_str = re.sub(
            r"c\^\{\\p\}_\{j,i\}(?![a-zA-Z])", r"c_{\\mathrm{p},j,i}", input_str)
        input_str = re.sub(
            r"c\^\{\\p\}(?![a-zA-Z])", r"c_{\\mathrm{p}}", input_str)
        input_str = re.sub(
            r"c\}\^\{\\p\}(?![a-zA-Z])", r"c}_{\\mathrm{p}}", input_str)
        input_str = re.sub(r"c\}\^\{\\l\}(?![a-zA-Z])", r"c}", input_str)

        input_str = re.sub(r"c\^\{\\s\}(?![a-zA-Z])", r"q", input_str)
        input_str = re.sub(r"c\}\^\{\\s\}(?![a-zA-Z])", r"q}", input_str)

        input_str = re.sub(r"\\p(?![a-zA-Z])", r"\\mathrm{p}", input_str)
        input_str = re.sub(r"\\s(?![a-zA-Z])", r"\\mathrm{s}", input_str)
    else:
        raise ValueError(f"Format {var_format} is not supported.")

    return input_str


# immutable & hashable dataclass -> unique particle type counter
@dataclass(frozen=True)
class Particle:

    geometry: Literal["Sphere", "Cylinder", "Slab"]
    has_core: bool
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
    # volume fraction ?
    # binding -> is_kinetic, nBound
    vars_and_params = []

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

        if self.has_binding:
            symbol_name_ = r"c^{\s}_{i}" if self.single_partype else r"c^{\s}_{j,i}"
            vars_and_params_.append({"Group" : 1, "Symbol": symbol_name_, "Description": r"particle solid concentration", "Unit": r"\frac{mol}{m^3}", "Dependence" : state_deps, "Domain" : eq.full_particle_conc_domain(column_resolution=self.interstitial_volume_resolution, particle_resolution=self.resolution, hasCore=self.has_core, with_par_index=False, with_time_domain=True)})
            symbol_name_ = r"f^\mathrm{bind}_{j,i}" if self.PTD else r"f^\mathrm{bind}_{i}"
            dep_ = r"\vec{c}^\mathrm{p}, \vec{c}^\mathrm{s}; j, i" if self.PTD else r"\vec{c}^\mathrm{p}, \vec{c}^\mathrm{s}; i"
            vars_and_params_.append({"Group" : 10, "Symbol": symbol_name_, "Description": r"adsorption isotherm function", "Unit": r"\frac{1}{s}", "Dependence": dep_})
            vars_and_params_.append({"Group" : 10.1, "Symbol": r"\vec{c}^\mathrm{p}", "Description": r"particle liquid components vector", "Unit": r"[\frac{mol}{m^3}]", "Dependence": state_deps})
            vars_and_params_.append({"Group" : 10.1, "Symbol": r"\vec{c}^\mathrm{s}", "Description": r"particle solid components vector", "Unit": r"[\frac{mol}{m^3}]", "Dependence": state_deps})
            
            if not (self.nonlimiting_filmDiff and self.resolution == "0D"):
                vars_and_params_.append({"Group" : 4, "Symbol": r"\varepsilon^{\mathrm{p}}" if self.single_partype else r"\varepsilon^{\mathrm{p}}_{j}", "Description": r"particle porosity", "Unit": r"-", "Dependence": r"\text{constant}" if self.single_partype else r"j", "Property": r"\in (0, 1)"})
            
            if self.has_surfDiff:
                symbol_name_ = r"D^\mathrm{s}_{i}" if self.single_partype else r"D^\mathrm{s}_{j,i}"
                vars_and_params_.append({"Group" : 6.1, "Symbol": symbol_name_, "Description": r"surface diffusion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": r"i" if self.single_partype else r"j,i", "Property": r"\geq 0"})
            
            if self.has_mult_bnd_states:
                symbol_name_ = r"N^{\mathrm{b}}_{i}" if self.single_partype else r"N^{\mathrm{b}}_{j,i}"
                vars_and_params_.append({"Group" : 11, "Symbol": symbol_name_, "Description": r"number of bound states", "Unit": r"-", "Dependence": r"i" if self.single_partype else r"j,i"})
            
            if self.has_core:
                self.vars_and_params.append({"Group" : 0.1, "Symbol": r"R^\mathrm{pc}", "Description": r"particle core radius", "Unit": r"-", "Dependence": r"-", "Property": r"\in (0, R^\mathrm{p})"})

        if not self.single_partype:
            vars_and_params_.append({"Group" : -0.1, "Symbol": r"j", "Description": r"particle type index", "Unit": r"-", "Dependence": r"-", "Property": r""})
            vars_and_params_.append({"Group" : 1.9, "Symbol": r"d_j", "Description": r"particle type volume fraction", "Unit": r"-", "Dependence": r"particle type", "Property": r""})

        for var_ in vars_and_params_:
            var_["Symbol"] = rerender_variables(var_["Symbol"], var_format_)
        
        vars_and_params_ = sorted(vars_and_params_, key=lambda x: x['Group'])

        object.__setattr__(self, 'vars_and_params', vars_and_params_)

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

@dataclass
class Column:

    dev_mode: bool  # dev mode including untested, unstable and wip features
    advanced_mode: bool  # Ask for detailed parameter inputs
    resolution: Literal = ""
    N_c: int = -1
    N_p: int = -1

    has_axial_coordinate: bool = False
    has_radial_coordinate: bool = False
    has_angular_coordinate: bool = False
    has_axial_dispersion: bool = False
    has_radial_dispersion: bool = False
    has_angular_dispersion: bool = False

    # List[bool] = None # TODO make this possible per component and particle type in advanced mode
    nonlimiting_filmDiff: bool = False

    particle_models: List[Particle] = None
    # counts per unique particle type (geometry, has_core, resolution)
    par_type_counts: Counter[Particle] = field(default_factory=Counter)
    # puts particle types together that have a similar contribution to the interstitial volume equation and counts them
    par_unique_intV_contribution_counts: Counter[Particle] = field(
        default_factory=Counter)
    has_surfDiff: bool = False
    has_binding: bool = True  # and thus solid phase
    req_binding: bool = False
    has_mult_bnd_states: bool = False
    PTD: bool = False
    # nCompReq: int
    # nCompKin: int
    # noPoresButReqBinding: bool

    has_filter: bool = False

    vars_and_params = List[dict]

    def __post_init__(self):

        # Configure interstitial column transport
        st.sidebar.write("Configure interstitial volume model")
        col1, col2 = st.columns(2)

        with col1:
            self.resolution = re.search(r'\dD', st.sidebar.selectbox("Column resolution", [
                                        "1D (axial coordinate)", "0D (Homogeneous Tank)", "2D (axial and radial coordinate)", "3D (axial, radial and angular coordinate)"], key=r"column_resolution")).group()

        valid_resolutions = {"3D", "2D", "1D", "0D"}

        if self.resolution not in valid_resolutions:
            raise ValueError(
                f"Invalid resolution: {self.resolution}. Must be one of {valid_resolutions}.")
        if int(re.search("\\d", self.resolution).group()) == 0:
            self.has_filter = st.sidebar.selectbox(
                    "Add flow rate filter", ["No", "Yes"], key=r"has_filter") == "Yes"
        if int(re.search("\\d", self.resolution).group()) > 0:
            self.has_axial_coordinate = True
            self.has_axial_dispersion = True
        if int(re.search("\\d", self.resolution).group()) > 1:
            self.has_radial_coordinate = True
            self.has_radial_dispersion = True
        if int(re.search("\\d", self.resolution).group()) > 2:
            self.has_angular_coordinate = True
            self.has_angular_dispersion = True

        self.N_c = st.sidebar.number_input(
            "Number of components", key=r"N_c", min_value=1, step=1) if dev_mode_ else -1

        if self.has_axial_coordinate:
            with col2:
                self.has_axial_dispersion = st.sidebar.selectbox(
                    "Add axial Dispersion", ["No", "Yes"], key=r"has_axial_dispersion") == "Yes"

        if self.dev_mode:
            if self.has_radial_coordinate:
                self.has_radial_dispersion = st.sidebar.selectbox(
                    "Add radial Dispersion", ["Yes", "No"], key=r"has_radial_dispersion") == "Yes"
            if self.has_angular_coordinate:
                self.has_angular_dispersion = st.sidebar.selectbox(
                    "Add angular Dispersion", ["Yes", "No"], key=r"has_angular_dispersion") == "Yes"

        # Configure particle model
        st.sidebar.write("Configure particle model")

        self.N_p = int(st.sidebar.selectbox("Add particles", ["No", "Yes"], key=r"add_particles") == "Yes")

        if self.N_p > 0:
            if dev_mode_:
                self.N_p = st.sidebar.number_input("Number of particle types", key=r"N^\mathrm{p}", min_value=0, step=1)
            elif advanced_mode_:
                self.N_p = 1 + int(st.sidebar.selectbox("Particle size distribution", ["No", "Yes"], key=r"PSD") == "Yes")
                if self.N_p > 1:
                    self.PTD = st.sidebar.selectbox("Particle type distribution (binding per j)", ["No", "Yes"], key=r"PTD") == "Yes"
                
        if self.N_p > 0:

            self.particle_models = []

            for j in range(self.N_p):

                if self.N_p > 1 and dev_mode_:
                    st.sidebar.write(f"Configure particle type {j + 1}")

                col1, col2 = st.columns(2)

                with col1:

                    if dev_mode_:  # multiple particle types
                        resolution = re.search(r'\dD', st.sidebar.selectbox(f"Select spatial resolution of particle type {j + 1}", [
                                               "1D (radial coordinate)", "0D (homogeneous)"], key=f"parType_{j+1}_resolution")).group()
                        has_core = st.sidebar.selectbox(f"Choose if particle type {j + 1} is a core-shell particle (i.e. " + r"$R^\mathrm{pc} > 0$)", [
                                                       "No core-shell", "Has core-shell"], key=f"parType_{j+1}_has_core") == "Has core-shell"
                        geometry = st.sidebar.selectbox(f"Select geometry of particle type {j + 1}", [
                                                        "Sphere", "Cylinder", "Slab"], key=f"parType_{j+1}__geometry")
                    elif j == 0:
                        resolution = re.search(r'\dD', st.sidebar.selectbox(f"Select spatial resolution of particles", [
                                               "1D (radial coordinate)", "0D (homogeneous)"], key=r"particle_resolution")).group()
                        has_core = st.sidebar.selectbox(f"Add impenetrable core-shell to particles (i.e. " + r"$R^\mathrm{pc} > 0$)", [
                                                       "No", "Yes"], key=f"particle_has_core") == "Yes" if (resolution == "1D" and advanced_mode_) else False
                        geometry = "Sphere"

                with col2:
                    if j == 0:  # todo make this configurable for every particle type
                        self.nonlimiting_filmDiff = st.sidebar.selectbox(
                            "Non-limiting film diffusion", ["No", "Yes"], key=r"nonlimiting_filmDiff") == "Yes"

                if j == 0:  # todo make this configurable for every particle type
                    # Configure binding model
                    st.sidebar.write("Configure binding model")

                    self.has_binding = st.sidebar.selectbox(
                        "Add binding", ["No", "Yes"], key=r"has_binding") == "Yes"

                    if self.has_binding:
                        
                        self.req_binding = st.sidebar.selectbox("Binding kinetics mode", [
                                                                "Kinetic", "Rapid-equilibrium"], key=r"req_binding") == "Rapid-equilibrium"
                        self.has_mult_bnd_states = st.sidebar.selectbox("Add multiple bound states", [
                                                                        "No", "Yes"], key=r"has_mult_bnd_states") == "Yes" if advanced_mode_ else False
                        self.has_surfDiff = st.sidebar.selectbox("Add surface diffusion", [
                                                                 "No", "Yes"], key=r"has_surfDiff") == "Yes" if resolution == "1D" else False

                self.particle_models.append(
                    Particle(
                        geometry=geometry,
                        resolution=resolution,
                        has_core=has_core,
                        has_binding=self.has_binding,
                        req_binding=self.req_binding,
                        has_mult_bnd_states=self.has_mult_bnd_states,
                        has_surfDiff=self.has_surfDiff,
                        nonlimiting_filmDiff=self.nonlimiting_filmDiff,
                        interstitial_volume_resolution=self.resolution,
                        single_partype=(self.N_p == 1),
                        PTD=self.PTD
                    )
                )

            # We need to count and thus sort particle_models:
            #  Particle types need individual particle equations, which is why we count them
            #  Only specific differences lead to changes in the interstitial volume equations: geometry if kinetic film diffusion. else geometry + resolution
            #  Sorting by (geometry, resolution) ensures we are using the same indices across the interstitial volume and particle eqquations.
            self.particle_models = sorted(self.particle_models, key=lambda particle: (
                particle.geometry, particle.resolution))
            self.par_type_counts = Counter(self.particle_models)
            if self.nonlimiting_filmDiff:
                self.par_unique_intV_contribution_counts = Counter(
                    (particle.geometry, particle.resolution) for particle in self.particle_models)
            else:
                self.par_unique_intV_contribution_counts = Counter(
                    particle.geometry for particle in self.particle_models)

        else:
            self.has_binding = False

        self.vars_and_params()

    def vars_and_params(self):

        without_pores_ = self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D"

        state_deps = r"t"
        param_deps = r"\text{constant}"
        if self.resolution == "1D":
            state_deps = r"t, z"
        if self.resolution == "2D":
            state_deps = r"t, z, \rho"
            param_deps = r"\rho"
        if self.resolution == "3D":
            state_deps = r"t, z, \rho, \phi"
            param_deps = r"z, \rho, \phi"

        state_deps += "; i"
        param_deps_comp = r"i" if param_deps == r"\text{constant}" else param_deps + "; i"

        self.vars_and_params = [
            {"Group" : 0, "Symbol": r"t", "Description": r"time coordinate", "Unit": r"s", "Dependence": r"\text{independent variable}", "Property": r"\in (0, T^{\mathrm{end}})"},
            {"Group" : 1, "Symbol": r"c^{\b}_i", "Description": r"bulk liquid concentration", "Unit": r"\frac{mol}{m^3}", "Dependence" : state_deps, "Domain": eq.int_vol_domain(self.resolution)},
            {"Group" : -1, "Symbol": r"T^{\mathrm{end}}", "Description": r"process end time", "Unit": r"s", "Dependence": r"\text{constant}", "Property": r" > 0"},
            {"Group" : -0.2, "Symbol": r"N^\mathrm{c}", "Description": r"number of components", "Unit": r"-", "Dependence": r"-", "Property": r"\in \mathbb{N}"},
            {"Group" : -0.1, "Symbol": r"i", "Description": r"component index", "Unit": r"s", "Dependence": r"-", "Property": r"-"},
            ]

        if self.has_axial_dispersion:
            if not without_pores_:
                self.vars_and_params.append({"Group" : 6, "Symbol": r"D^\mathrm{ax}_i", "Description": r"axial dispersion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})
            else:
                self.vars_and_params.append({"Group" : 6, "Symbol": r"\tilde{D}^\mathrm{ax}_i", "Description": r"apparent axial dispersion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})
        if self.has_radial_dispersion:
            self.vars_and_params.append({"Group" : 6, "Symbol": r"D^\mathrm{rad}_i", "Description": r"radial dispersion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})
        if self.has_angular_dispersion:
            self.vars_and_params.append({"Group" : 6, "Symbol": r"D^\mathrm{ang}_i", "Description": r"angular dispersion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})

        if self.resolution == "0D":
            self.vars_and_params.append({"Group" : 2, "Symbol": r"Q^\mathrm{in}", "Description": r"volumetric flow rate into the tank", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"\geq 0"})
            self.vars_and_params.append({"Group" : 2, "Symbol": r"Q^\mathrm{out}", "Description": r"volumetric flow rate out of the tank", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"\geq 0"})
            if self.has_filter:
                self.vars_and_params.append({"Group" : 2, "Symbol": r"Q^\mathrm{filter}", "Description": r"volumetric flow rate out of the tank (solvent only)", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"\geq 0"})
        else:
            self.vars_and_params.append({"Group" : 0, "Symbol": r"z", "Description": r"axial cylinder coordinate", "Unit": r"m", "Dependence": r"\text{independent variable}", "Property": r"\in (0, L)"})
            self.vars_and_params.append({"Group" : -1, "Symbol": r"L", "Description": r"length of cylinder", "Unit": r"m", "Dependence": r"\text{constant}", "Property": r" > 0"})
            if not without_pores_:
                self.vars_and_params.append({"Group" : 5, "Symbol": r"u", "Description": r"interstitial velocity", "Unit": r"\frac{m}{s}", "Dependence": param_deps, "Property": r"> 0"})
            else:
                self.vars_and_params.append({"Group" : 5, "Symbol": r"\tilde{u}", "Description": r"apparent interstitial velocity", "Unit": r"\frac{m}{s}", "Dependence": param_deps, "Property": r"> 0"})

        if self.resolution == "2D":
            self.vars_and_params.append({"Group" : 0, "Symbol": r"\rho", "Description": r"radial cylinder coordinate", "Unit": r"m", "Dependence": r"\text{independent variable}", "Property": r"\in (0, R^{\mathrm{c}})"})
            self.vars_and_params.append({"Group" : -1, "Symbol": r"R^{\mathrm{c}}", "Description": r"cylinder radius", "Unit": r"m", "Dependence": r"\text{constant}", "Property": r" > 0"})

        if self.resolution == "3D":
            self.vars_and_params.append({"Group" : 0, "Symbol": r"\phi", "Description": r"angular cylinder coordinate", "Unit": r"m", "Dependence": r"\text{independent variable}", "Property": r"\in (0, 2\pi)"})

        if self.N_p > 0:
            if not without_pores_:
                self.vars_and_params.append({"Group" : 0.1, "Symbol": r"R^\mathrm{p}", "Description": r"particle radius", "Unit": r"-", "Dependence": r"-", "Property": r"> 0"})
            if not without_pores_:
                self.vars_and_params.append({"Group" : 4, "Symbol": r"\varepsilon^{\mathrm{c}}", "Description": r"column porosity", "Unit": r"-", "Dependence": re.sub("t, ", "", state_deps), "Property": r"\in (0, 1)"})
            else:
                self.vars_and_params.append({"Group" : 4, "Symbol": r"\varepsilon^{\mathrm{t}}", "Description": r"total porosity", "Unit": r"-", "Dependence": re.sub("t, ", "", state_deps), "Property": r"\in (0, 1)"})
            if not self.nonlimiting_filmDiff:
                symbol_name_ = r"k^\mathrm{f}_{i}" if self.N_p<=1 else r"k^\mathrm{f}_{j,i}"
                self.vars_and_params.append({"Group" : 7, "Symbol": symbol_name_, "Description": r"film diffusion coefficient", "Unit": r"\frac{m}{s}", "Dependence": r"i" if self.N_p<=1 else r"j,i", "Property": r"> 0"})

        for var_ in self.vars_and_params:
            var_["Symbol"] = rerender_variables(var_["Symbol"], var_format_)
            
        self.vars_and_params = sorted(self.vars_and_params, key=lambda x: x['Group'])

    def interstitial_volume_equation(self):

        without_pores_ = self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D"       

        if self.resolution == "0D":

            filter_str = r" - Q_{\mathrm{filter}}" if self.has_filter else ""

            equation = r"""
    \frac{\mathrm{d}V^{\b}}{\mathrm{d}t} &= Q_{\mathrm{in}} - Q_{\mathrm{out}}""" + filter_str + r""",
    \\
    \frac{\mathrm{d}}{\mathrm{d} t} \left( V^{\b} c^{\b}_i \right)"""

            if self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D":
                equation += r" + V^{\p} \varepsilon^{\mathrm{p}} \frac{\partial c^{b}_i}{\partial t}"
                if self.req_binding:
                    equation += r" + V^{\p} \left( 1 - \varepsilon^{\mathrm{p}} \right) \frac{\partial c^{\s}_i}{\partial t}"

            equation += r"&= Q_{\mathrm{in}} c^{\b}_{\mathrm{in},i} - Q_{\mathrm{out}} c^{\b}_i"

            if self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D" and not self.req_binding:
                equation += r" - V^{\p} \left( 1 - \varepsilon^{\mathrm{p}} \right) \frac{\partial c^{\s}_i}{\partial t}"

        else:
            equation = eq.bulk_time_derivative(r"\varepsilon^{\mathrm{c}}") if not without_pores_ else eq.bulk_time_derivative()
            if without_pores_:
                equation += r" + " + eq.solid_time_derivative(r"\varepsilon^{\mathrm{t}}")

            equation += " = " + eq.axial_convection() if without_pores_ else " = " + eq.axial_convection(r"\varepsilon^{\mathrm{c}}")

            if self.has_axial_dispersion:
                equation += " + " + eq.axial_dispersion(r"\varepsilon^{\mathrm{c}}")
            if self.has_radial_dispersion:
                equation += " + " + eq.radial_dispersion(r"\varepsilon^{\mathrm{c}}")
            if self.has_angular_dispersion:
                equation += " + " + eq.angular_dispersion(r"\varepsilon^{\mathrm{c}}")

            if self.N_p == 0:  # remove occurencies of porosity, which is just constant one in this case
                equation = re.sub(r"\\varepsilon^{\\mathrm{c}}", "", re.sub(
                    r"\\left\( \\varepsilon^{\\mathrm{c}} c^{\\l}_i \\right\)", r"c^{\\l}_i", equation))

        # if self.nonlimiting_filmDiff and 1Dparticle # entscheidende faktoren sind particle resolution und filmDiffMode. the following loop has thus to change
        par_added = 0
        for par_uniq in self.par_unique_intV_contribution_counts.keys():

            if dev_mode_:
                equation += eq.int_filmDiff_term(
                    Particle(
                        self.particle_models[0].geometry, self.particle_models[0].has_core, self.particle_models[0].resolution
                    ),
                    1 + par_added, par_added + self.par_unique_intV_contribution_counts[par_uniq], self.N_p == 1, self.nonlimiting_filmDiff, self.particle_models[0].has_surfDiff
                )
            else:
                equation += eq.int_filmDiff_term(
                    Particle(
                        self.particle_models[0].geometry, self.particle_models[0].has_core, self.particle_models[0].resolution
                    ),
                    1, r"N^{\mathrm{p}}", self.N_p == 1, self.nonlimiting_filmDiff, self.particle_models[0].has_surfDiff
                )

            par_added += self.par_unique_intV_contribution_counts[par_uniq]

        if self.resolution == "0D":
            equation = re.sub(
                r"\\left\(1 - \\varepsilon^{\\mathrm{c}} \\right\)", r"V^s", equation)

        equation = r"""\begin{align}
""" + equation + r""",
\end{align}"""

        return equation

    def interstitial_volume_bc(self):
        if not self.resolution == "0D":
            return eq.int_vol_BC(self.resolution, self.has_axial_dispersion)
        else:
            return None

    def particle_equations(self):

        eqs = {}
        boundary_conditions = {}

        for par_type in self.par_type_counts.keys():

            eqs[par_type] = eq.particle_transport(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=self.nonlimiting_filmDiff,
                                                  has_surfDiff=self.has_surfDiff, has_binding=self.has_binding, req_binding=self.req_binding, has_mult_bnd_states=self.has_mult_bnd_states, PTD=self.PTD)
            eqs[par_type] = eqs[par_type]

            boundary_conditions[par_type] = eq.particle_boundary(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=self.nonlimiting_filmDiff,
                                                                 has_surfDiff=self.has_surfDiff, has_binding=self.has_binding, req_binding=self.req_binding, has_mult_bnd_states=self.has_mult_bnd_states)
            boundary_conditions[par_type] = boundary_conditions[par_type]

            if self.N_p == 1:
                eqs[par_type] = re.sub(",j", "", eqs[par_type])
                eqs[par_type] = re.sub("j,", "", eqs[par_type])
                eqs[par_type] = re.sub("_{j}", "", eqs[par_type])
                boundary_conditions[par_type] = re.sub(",j", "", boundary_conditions[par_type])
                boundary_conditions[par_type] = re.sub("j,", "", boundary_conditions[par_type])
                boundary_conditions[par_type] = re.sub("_{j}", "", boundary_conditions[par_type])

        return eqs, boundary_conditions

    def domain_interstitial(self, with_time_domain=True):
        return r"$" + eq.int_vol_domain(self.resolution, with_time_domain=with_time_domain) + r"$"

    def domain_particle(self):
        if self.N_p > 0:
            return eq.full_particle_conc_domain(self.resolution, self.particle_models[0].resolution, self.particle_models[0].has_core, False, False)
        else:
            return ""
        
    def model_name(self):

        if self.resolution == "0D":
            
            if self.N_p > 0:
                model_name = "Finite Bath"
                
                if self.particle_models[0].resolution == "0D":
                    model_name += " without Pores" if self.nonlimiting_filmDiff else " with Pores"
                else:
                    model_name = "General " + model_name
                
                return model_name
            else:
                return "Continuously Stirred Tank"

        if self.has_angular_coordinate:
            model_name = "3D "
        elif self.has_radial_coordinate:
            model_name = "2D "
        # elif self.has_axial_coordinate: # default case, no name prefix !
        else:
            model_name = ""

        if self.N_p > 0:

            if self.particle_models[0].resolution == "1D":
                model_name += "General Rate Model"

            else:
                model_name += "Lumped Rate Model"

                model_name += " without Pores" if self.nonlimiting_filmDiff else " with Pores"

            if self.N_p > 1:
                if not any(par.geometry != "Sphere" for par in self.particle_models):
                    model_name += " with particle-size distribution"  # particle-size distribution
                else:
                    # particle-type distribution # TODO use when different kinds of geometry or binding
                    model_name += " with particle-type distribution"
        else:
            if self.has_axial_dispersion or self.has_radial_dispersion or self.has_angular_dispersion:
                model_name += "Dispersive "
            model_name += "Plug Flow"  # Reactor if we have reactions

        return model_name

    def model_assumptions(self):

        asmpts = {
            "General model assumptions": eq.HRM_asmpt(self.N_p, self.nonlimiting_filmDiff, self.has_binding, self.has_surfDiff, self.resolution),
            "Specific model assumptions": eq.int_vol_continuum_asmpt(self.resolution, self.N_p, self.nonlimiting_filmDiff) +
            (eq.particle_asmpt(
                self.particle_models[0].resolution, self.has_surfDiff) if self.N_p > 0 else [])
        }

        if self.nonlimiting_filmDiff:
            asmpts["Specific model assumptions"].append(
                r"the film around the particles does not limit mass transfer. That is, we assume $k^{\mathrm{f}}_i = \infty$)")
        if self.req_binding:
            asmpts["Specific model assumptions"].append(
                r"adsorption and desorption happen on a much faster time scale than the other mass transfer processes (e.g., convection, diffusion). Hence, we consider them to be equilibrated instantly, that is, to always be in (local) equilibrium")

        for idx in range(0, len(asmpts["Specific model assumptions"])):
            
            asmpts["Specific model assumptions"][idx] = rerender_variables(asmpts["Specific model assumptions"][idx], var_format_)
            
        for idx in range(0, len(asmpts["General model assumptions"])):
            
            asmpts["General model assumptions"][idx] = rerender_variables(asmpts["General model assumptions"][idx], var_format_)
            
        return asmpts

    def vars_params_description(self):

        description_ = ""

        idx_ = 1
        num_VP = len(self.vars_and_params)

        for thing in self.vars_and_params:

            if thing.get("Group", -1) < 0: # dont print symbols with negative group no.
                num_VP -= 1
                continue

            if not idx_ == 1:
                description_ += ", " if idx_ < num_VP else ", and "
            description_ += r"$" + thing["Symbol"]

            if not thing.get("Domain", "-") == "-":
                description_ += r"\colon " + re.sub(r"\$", "", thing["Domain"]) + r" \mapsto \mathbb{R}"

            description_ += thing.get("Property", "") + r"$"
            
            description_ += " is the " + thing["Description"]

            idx_ += 1
                
        return description_ + "."

# %% Streamlit UI

st.logo("images/logo_CADET.png", size="large", link=None, icon_image=None)

st.set_page_config(
    page_title=r"CADET-Equations",
    page_icon=r":material/biotech:", #":material/modeling:",
)

st.sidebar.title(
    "CADET-Equations: Packed-Bed Chromatography Model Equation Generator")
st.sidebar.write(
    "Configure a chromatography model to get the corresponding governing equations.")

# File uploader for JSON file
uploaded_file = st.sidebar.file_uploader(
    "Upload a configuration file (optional, see documentation)", type=["json", "h5"])

if uploaded_file is not None:

    uploaded_file_name = uploaded_file.name.lower()
    
    config = None

    if uploaded_file_name.endswith(".json"):

        config = json.load(uploaded_file)

    elif uploaded_file_name.endswith(".h5"):

        config = load_CADET_h5.get_config_from_CADET_h5(uploaded_file,
        str(st.sidebar.number_input("Unit index in CADET file", key=r"h5_input_unit_index", min_value=-1, max_value=999, step=1, value=-1)).zfill(3)
        )

    if config is not None:
        # Update Streamlit session state
        for key, value in config.items():

            st.session_state[key] = value

        st.sidebar.success("Uploaded configuration applied!")
    else:
        st.sidebar.error("Uploaded configuration could not be read!")

# User configuration of the model

var_format_ = st.sidebar.selectbox("Select format (e.g. $c^s$ or $q$ as the solid phase concentration)", [
                                   "CADET", "Legacy"], key=r"var_format")

advanced_mode_ = st.sidebar.selectbox("Advanced setup options (enables e.g. particle size distribution)", [
                                      "Off", "On"], key=r"advanced_mode") == "On"
if advanced_mode_:
    dev_mode_ = st.sidebar.selectbox("Developer setup options (not tested! Enables e.g. particle type distribution)", [
                                     "Off", "On"], key=r"dev_mode") == "On"
    if dev_mode_:
        advanced_mode_ = True
else:
    dev_mode_ = False

column_model = Column(dev_mode=dev_mode_, advanced_mode=advanced_mode_)

# %% Display equations

file_content = []  # used to export model to files

if st.toggle("Show Model Assumptions", key=r"model_assumptions"):

    asmpts = column_model.model_assumptions()

    for key in asmpts.keys():

        st.write(key + ":\n" + "\n".join(f"- {item}" for item in asmpts[key]))

        file_content.append(
            key + r""":
\begin{itemize}
""" + "\n".join(f"\\item {item}" for item in asmpts[key]) + r"""
\end{itemize}

"""
        )

if st.toggle("Show symbol table", key=r"sym_table"):

    df = pd.DataFrame(column_model.vars_and_params)
    if column_model.N_p > 0:
        df_par = pd.DataFrame(column_model.particle_models[0].vars_and_params, columns=["Group", 'Symbol', 'Description', "Dependence", 'Unit'])
        df = pd.concat([df, df_par], ignore_index=True)
    
    df[['Symbol', "Dependence", 'Unit']] = df[['Symbol', "Dependence", 'Unit']].map(lambda x: f"${x}$" if isinstance(x, str) else x)
    
    df = df.sort_values(by=r"Group").reset_index()
    
    st.table(df[['Symbol', 'Description', 'Dependence', 'Unit']])

interstitial_volume_eq = column_model.interstitial_volume_equation()

nComp_list = r"$i\in\{" + ", ".join(str(i) for i in range(1, column_model.N_c + 1)) + \
    r"\}$" if column_model.N_c > 0 else r"$i\in\{1, \dots, N^{\mathrm{c}} \}$"

show_eq_description = st.toggle("Show equation description", key=r"show_eq_description", value=True)

# The following function is used to both print the output and collect it to later generate and export output files
def write_and_save(output: str, as_latex: bool = False):

    output = rerender_variables(output, var_format_)

    if output is not None:

        file_content.append(output)

        if as_latex:
            st.latex(output)
        else:
            st.write(output)


st.write("### " + column_model.model_name())
file_content.append(r"\section*{" + column_model.model_name() + r"}")

if column_model.resolution == "0D":
    intro_str = r"Consider a continuously stirred tank "
else:
    intro_str = r"Consider a cylindrical column of length $L > 0$ "
    if column_model.resolution == "2D" or column_model.resolution == "3D":
        intro_str += r" and radius $R^{\mathrm{c}} > 0$ "

if column_model.N_p == 0:
    write_and_save(intro_str + r"filled with a liquid phase, and observed over a time interval $(0, T^{\mathrm{end}})$.")
elif column_model.N_p == 1:
    # TODO particle geometries?
    write_and_save(intro_str + r"packed with spherical particles, and observed over a time interval $(0, T^{\mathrm{end}})$.")
else:
    if column_model.resolution == "0D":
        d_j_def = r"$d_j \in [0, 1]$"
    else:
        d_j_def = r"$d_j \colon " + re.sub(r"\$", "", column_model.domain_interstitial(with_time_domain=False)) + r" \to [0, 1]$"
    write_and_save(intro_str + r"packed with $N^{\mathrm{p}}$ different particle-sizes indexed by $j \in \{1, \dots, N^{\mathrm{p}}\}$ and distributed according to the volume fractions " + d_j_def + r", which satisfy")

    if column_model.resolution == "0D":
        d_j_dep = r""
        d_j_dep2 = r""
    else:
        if column_model.resolution == "1D":
            d_j_dep = r"z"
        elif column_model.resolution == "2D":
            d_j_dep = r"z, \rho"
        elif column_model.resolution == "3D":
            d_j_dep = r"z, \rho, \phi"

        d_j_dep2 = r", \quad \forall """ + d_j_dep + r""" \in """ + re.sub(r"\$", "", column_model.domain_interstitial(with_time_domain=False))
        d_j_dep = "(" + d_j_dep + ")"

    write_and_save(r"""
    \begin{equation*}
	    \sum_{j=1}^{N_{\mathrm{p}}} d_j""" + d_j_dep + r" = 1 " + d_j_dep2 + r""".
    \end{equation*}

    """, as_latex=True)


if column_model.resolution == "0D":
    write_and_save(
        r"The evolution of the liquid volume $V^{\l}\colon (0, T^{\mathrm{end}}) \to \mathbb{R}$ and the concentrations $c_i\colon (0, T^{\mathrm{end}}) \to \mathbb{R}$ of the components in the tank is governed by")
    write_and_save(interstitial_volume_eq, as_latex=True)
else:
    eq_type = "convection" # first order derivative
    if not column_model.resolution == "1D" or column_model.has_axial_dispersion: # second order derivative
        eq_type += "-diffusion"
    if column_model.N_p > 0 and not column_model.nonlimiting_filmDiff: # i.e. film diffusion term (reaction type)
        eq_type += "-reaction"

    write_and_save(r"In the interstitial volume, mass transfer is governed by the following " + eq_type +
                   " equations in " + column_model.domain_interstitial() + r" and for all components " + nComp_list)
    write_and_save(interstitial_volume_eq, as_latex=True)
    write_and_save("with boundary conditions")
    write_and_save(column_model.interstitial_volume_bc(), as_latex=True)

if show_eq_description:
    write_and_save("Here, " + column_model.vars_params_description())

if column_model.N_p > 0:

    particle_eq, particle_bc = column_model.particle_equations()

    cur_par_count = 0
    for par_type in column_model.par_type_counts.keys():

        # in this case, we dont have a particle model. this configuration is still allowed for educational purpose.
        if not column_model.has_binding and column_model.nonlimiting_filmDiff and par_type.resolution == "0D":
            break

        if dev_mode_:
            nPar_list = ', '.join(str(j) for j in range(cur_par_count, column_model.par_type_counts[par_type] + 1))
        else:
            nPar_list = r"$j\in\{1, \dots, N^{\mathrm{p}}\}$"

        eq_type_ = "reaction" if column_model.particle_models[0].resolution == "0D" else "diffusion-reaction"
        
        tmp_str = r" and all particle sizes " + nPar_list if column_model.N_p > 1 else r""
        write_and_save(
            "In the particles, mass transfer is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True) + r" and for all components" + tmp_str)

        write_and_save(particle_eq[par_type], as_latex=True)
        cur_par_count += column_model.par_type_counts[par_type]

        if not particle_bc[par_type] == "":
            write_and_save("with boundary conditions")
            write_and_save(particle_bc[par_type], as_latex=True)

        if show_eq_description:
            write_and_save("Here, " + column_model.particle_models[0].vars_params_description())

write_and_save("Consistent initial values for all solution variables (concentrations) are defined at $t = 0$.")

st.session_state.latex_string = [
    r"""\documentclass{article}
""",
    r"""\usepackage{amssymb,amsmath,mleftright}
""",
    r"""\begin{document}
"""
]
st.session_state.latex_string.extend(file_content + [r"\end{document}"])
st.session_state.latex_string = "\n".join(st.session_state.latex_string)

st.session_state.latex_string = str(st.session_state.latex_string) # for testing purposes

st.download_button("Download .tex", st.session_state.latex_string, "model.tex", "text/plain")

if st.button("Generate PDF", key=r"generate_pdf"):
    with tempfile.TemporaryDirectory() as temp_dir:
        tex_path = f"{temp_dir}/model.tex"
        pdf_path = f"{temp_dir}/model.pdf"

        # Write LaTeX content to a temporary file
        with open(tex_path, "w") as f:
            f.write(st.session_state.latex_string)

        result = subprocess.run(
            ["pdflatex", "-output-directory", temp_dir, tex_path],
            capture_output=True,
            text=True,
            cwd=temp_dir
        )

        if result.returncode != 0:
            st.error(
                "PDF generation failed. Please check the LaTeX compilation output:\n" + result.stderr)
        else:
            with open(pdf_path, "rb") as pdf_file:
                st.download_button("Download PDF", pdf_file,
                                   "model.pdf", "application/pdf")

if st.button("Generate configuration file", key=r"generate_config"):

    def sort_session_states(key):
        # a proper sorting is required for the tests, where we cannot apply a
        # global session state but must loop over the keys
        if key == "advanced_mode":
            return 0
        elif key == "dev_mode":
            return 1
        elif key == "column_resolution":
            return 2
        elif key == "add_particles":
            return 3
        elif key == "has_binding":
            return 4
        return 10

    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        json_path = f"{temp_dir}/model.json"

        config = {key: st.session_state[key]
                  for key in st.session_state if key not in ui_config._NO_MODEL_CONFIG_STATES_}

        config = OrderedDict(sorted(config.items(), key=lambda x: (sort_session_states(x[0]), x[0])))

        config_json = json.dumps(config, indent=4)

        st.write("Current configuration:\n", config_json)

        with open(json_path, "w") as json_file:
            json_file.write(config_json)

        st.download_button(
            "Download Configuration (as json)",
            config_json,
            "config.json",
            "application/json"
        )
