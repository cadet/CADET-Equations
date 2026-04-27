# -*- coding: utf-8 -*-
"""
This script implements the `Column` dataclass and related logic.
"""

from dataclasses import dataclass, field
from collections import Counter
from typing import List, Literal, Optional
import re
import streamlit as st

from src import equations as eq
from src.utils import format_variables
from src.model_particle import Particle


@dataclass
class Column:
    """Minimal Description of a packed-bed column model used by the UI.

    The class collects configuration state, produces equation strings
    and exposes small helpers used by the Streamlit generator.
    """

    dev_mode: bool  # dev mode including untested, unstable and wip features
    advanced_mode: bool  # Ask for detailed parameter inputs
    var_format: Literal
    resolution: Literal = None
    N_c: int = -1
    N_p: int = -1

    column_type: Literal["Axial", "Radial", "Frustum"] = "Axial"

    has_axial_coordinate: bool = False
    has_radial_coordinate: bool = False
    has_angular_coordinate: bool = False
    has_axial_dispersion: bool = False
    has_radial_dispersion: bool = False
    has_angular_dispersion: bool = False

    nonlimiting_filmDiff: bool = False

    # Per-component configuration (used in dev_mode when N_c > 0)
    req_binding_per_comp: Optional[List[bool]] = None
    nonlimiting_filmDiff_per_comp: Optional[List[bool]] = None
    has_surfDiff_per_comp: Optional[List[bool]] = None
    has_mult_bnd_states_per_comp: Optional[List[bool]] = None

    # Per-particle-type configuration (used when N_p > 1 and N_c <= 0)
    nonlimiting_filmDiff_per_partype: Optional[List[bool]] = None
    has_surfDiff_per_partype: Optional[List[bool]] = None

    particle_models: Optional[List[Particle]] = None
    # counts per unique particle type (geometry, has_core, resolution)
    par_type_counts: Counter[Particle] = field(default_factory=Counter)
    # puts particle types together that have a similar contribution to the interstitial volume equation and counts them
    par_unique_intV_contribution_counts: Counter[Particle] = field(
        default_factory=Counter)
    has_surfDiff: bool = False
    has_binding: bool = True  # and thus solid phase
    req_binding: bool = False
    has_mult_bnd_states: bool = False
    # nCompReq: int
    # nCompKin: int
    # noPoresButReqBinding: bool

    has_filter: bool = False

    binding_model: str = "Arbitrary"

    has_reaction_bulk: bool = False
    has_reaction_particle_liquid: bool = False
    has_reaction_particle_solid: bool = False
    req_reaction_bulk: bool = False
    req_reaction_particle_liquid: bool = False
    req_reaction_particle_solid: bool = False

    vars_and_params: List[dict] = field(default_factory=list)

    def __post_init__(self):

        # Configure interstitial column transport
        st.sidebar.write("Configure interstitial volume model")
        col1, col2 = st.columns(2)

        with col1:
            column_type_label = st.sidebar.selectbox("Column geometry", [
                "Axial flow cylinder", "Radial flow cylinder", "Frustum"], key=r"column_type")
            self.column_type = {"Axial flow cylinder": "Axial", "Radial flow cylinder": "Radial", "Frustum": "Frustum"}[column_type_label]

            resolution_options = {
                "Axial": ["1D (axial coordinate)", "0D (Homogeneous Tank)", "2D (axial and radial coordinate)", "3D (axial, radial and angular coordinate)"],
                "Radial": ["1D (radial coordinate)"],
                "Frustum": ["1D (axial coordinate)"],
            }
            self.resolution = re.search(r'\dD', st.sidebar.selectbox("Column resolution",
                                        resolution_options[self.column_type], key=r"column_resolution")).group()

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

        if self.dev_mode:
            n_c_choice = st.sidebar.selectbox(
                "Number of components (enables per-component configuration)",
                ["Arbitrary"] + list(range(1, 11)),
                key=r"N_c_choice")
            self.N_c = -1 if n_c_choice == "Arbitrary" else int(n_c_choice)
        else:
            self.N_c = -1

        if self.has_axial_coordinate:
            with col2:
                dispersion_label = "Add radial Dispersion" if self.column_type == "Radial" else "Add axial Dispersion"
                self.has_axial_dispersion = st.sidebar.selectbox(
                    dispersion_label, ["No", "Yes"], key=r"has_axial_dispersion") == "Yes"

        if self.dev_mode:
            if self.has_radial_coordinate:
                self.has_radial_dispersion = st.sidebar.selectbox(
                    "Add radial Dispersion", ["Yes", "No"], key=r"has_radial_dispersion") == "Yes"
            if self.has_angular_coordinate:
                self.has_angular_dispersion = st.sidebar.selectbox(
                    "Add angular Dispersion", ["Yes", "No"], key=r"has_angular_dispersion") == "Yes"

        # Configure particle model
        with st.sidebar.expander("Configure particles", expanded=True):

            if self.dev_mode:
                self.N_p = st.number_input("Number of particle types", key=r"N^\mathrm{p}", min_value=0, step=1)
                
            elif self.advanced_mode:
                par_config = st.selectbox("Add particles", ["No", "Yes", "Particle size distribution"], key=r"PSD")
                self.N_p = 1 if par_config == "Yes" else 0 if par_config == "No" else 2
                
            else:
                self.N_p = int(st.selectbox("Add particles", ["No", "Yes"], key=r"add_particles") == "Yes")


            self.particle_models = []

            if self.N_p == 1 or not self.dev_mode:
                
                configure_particle_type(typeCounter=-1)
                
            else:
                for j in range(self.N_p):
                    
                    st.expander(f"Configure particle type {j + 1}")
                    
                    configure_particle_type(typeCounter=j)

            # We need to count and thus sort particle_models:
            #  Particle types need individual particle equations, which is why we count them
            #  Only specific differences lead to changes in the interstitial volume equations: geometry if kinetic film diffusion. else geometry + resolution
            #  Sorting by (geometry, resolution, nonlimiting_filmDiff, has_surfDiff) ensures we are using the same indices across the interstitial volume and particle equations.
            self.particle_models = sorted(self.particle_models, key=lambda particle: (
                particle.geometry, particle.resolution, particle.nonlimiting_filmDiff, particle.has_surfDiff))
            self.par_type_counts = Counter(self.particle_models)
            if self.nonlimiting_filmDiff:
                self.par_unique_intV_contribution_counts = Counter(
                    (particle.geometry, particle.resolution, particle.nonlimiting_filmDiff) for particle in self.particle_models)
            else:
                self.par_unique_intV_contribution_counts = Counter(
                    (particle.geometry, particle.nonlimiting_filmDiff) for particle in self.particle_models)
                

        with st.sidebar.expander("Configure reactions"):
        
            if not self.dev_mode:
                
                self.has_reaction_bulk = st.selectbox(
                    "Add bulk liquid reaction (kinetic)", ["No", "Yes"], key=r"has_reaction_bulk") == "Yes"
            else:
                
                self.has_reaction_bulk = st.selectbox(
                    "Add bulk liquid reaction", ["No", "Yes"], key=r"has_reaction_bulk") == "Yes"
            
            if self.has_reaction_bulk and self.dev_mode:
                
                self.req_reaction_bulk = st.selectbox(
                    "Bulk reaction kinetics mode", ["Kinetic", "Rapid-equilibrium"], key=r"req_reaction_bulk") == "Rapid-equilibrium"
            
            if self.N_p > 0 and self.dev_mode:
            
                self.has_reaction_particle_liquid = st.selectbox(
                    "Add particle liquid reaction", ["No", "Yes"], key=r"has_reaction_particle_liquid") == "Yes"
                
                if self.has_reaction_particle_liquid:
                
                    self.req_reaction_particle_liquid = st.selectbox(
                        "Particle liquid reaction kinetics mode", ["Kinetic", "Rapid-equilibrium"], key=r"req_reaction_particle_liquid") == "Rapid-equilibrium"
                
                if self.has_binding:
                    self.has_reaction_particle_solid = st.selectbox(
                        "Add particle solid reaction", ["No", "Yes"], key=r"has_reaction_particle_solid") == "Yes"
                    
                    if self.has_reaction_particle_solid:
                        self.req_reaction_particle_solid = st.selectbox(
                            "Particle solid reaction kinetics mode", ["Kinetic", "Rapid-equilibrium"], key=r"req_reaction_particle_solid") == "Rapid-equilibrium"

        self.fill_vars_and_params()

    def configure_particle_type(self, typeCounter):

        typeDiff = typeCounter >= 0
        j = typeCounter

        keyPrefix = "parType_{j+1}_" if typeDiff else "particle_"

        col1, col2 = st.columns(2)

        with col1:

            resolution = re.search(r'\dD', st.selectbox(f"Spatial resolution", [
                                   "1D (radial coordinate)", "0D (homogeneous)"], key=keyPrefix+"resolution")).group()
            
            has_core = st.selectbox(f"Add impenetrable core-shell to particles (i.e. " + r"$R^\mathrm{pc} > 0$)", [
                                           "No", "Yes"], key=keyPrefix+r"has_core") == "Yes" if (resolution == "1D" and self.advanced_mode) else False
            
            if self.dev_mode:
                geometry = st.selectbox(
                    f"Geometry", ["Sphere", "Cylinder", "Slab"], key=keyPrefix+r"geometry"
                    )
            else:
                geometry = "Sphere"

        with col2:
            if self.N_c <= 0:  # particle film diffusion (hidden when per-component is active)
                self.nonlimiting_filmDiff = st.selectbox(
                    "Infinite film diffusion rate", ["No", "Yes"], key=keyPrefix+r"nonlimiting_filmDiff") == "Yes"
            # else: Todo : here goes the per-component config of film diffusion

        # todo: continue from here

        # # Configure binding model
        # st.write("Configure binding model")

        # self.has_binding = st.selectbox(
        #     "Add binding", ["No", "Yes"], key=r"has_binding") == "Yes"

        # if self.has_binding:

        #     if self.N_c <= 0:
        #         # Global options (shown when per-component is not active)
        #         self.req_binding = st.selectbox("Binding kinetics mode", [
        #                                                 "Kinetic", "Rapid-equilibrium"], key=r"req_binding") == "Rapid-equilibrium"
        #     self.binding_model = st.selectbox("Binding model", eq.BINDING_MODELS, key=r"binding_model")
        #     if self.N_c <= 0:
        #         self.has_mult_bnd_states = st.selectbox("Add multiple bound states", [
        #                                                         "No", "Yes"], key=r"has_mult_bnd_states") == "Yes" if self.advanced_mode else False
        #         if self.N_p <= 1:
        #             self.has_surfDiff = st.selectbox("Add surface diffusion", [
        #                                                      "No", "Yes"], key=r"has_surfDiff") == "Yes" if resolution == "1D" else False

        # # Per-particle-type film diffusion and surface diffusion
        # if self.N_p > 1:
        #     self.nonlimiting_filmDiff_per_partype = []
        #     self.has_surfDiff_per_partype = []
        #     for jj in range(self.N_p):
        #         with st.expander(f"Particle type {jj + 1}"):
        #             self.nonlimiting_filmDiff_per_partype.append(
        #                 st.selectbox("Infinite film diffusion rate",
        #                              ["No", "Yes"],
        #                              key=f"nonlimiting_filmDiff_partype_{jj}") == "Yes"
        #             )
        #             if self.has_binding and resolution == "1D":
        #                 self.has_surfDiff_per_partype.append(
        #                     st.selectbox("Surface diffusion",
        #                                  ["No", "Yes"],
        #                                  key=f"has_surfDiff_partype_{jj}") == "Yes"
        #                 )
        #             else:
        #                 self.has_surfDiff_per_partype.append(False)
        #     self.nonlimiting_filmDiff = all(self.nonlimiting_filmDiff_per_partype)
        #     self.has_surfDiff = any(self.has_surfDiff_per_partype)

        # # Per-component configuration (independent of binding)
        # # When N_p > 1, film diffusion and surface diffusion are handled
        # # per particle type above; per-component only configures binding settings
        # if self.dev_mode and self.N_c > 0:
        #     st.write("Per-component configuration")
        #     self.req_binding_per_comp = []
        #     self.nonlimiting_filmDiff_per_comp = []
        #     self.has_surfDiff_per_comp = []
        #     self.has_mult_bnd_states_per_comp = []
        #     for comp_i in range(self.N_c):
        #         with st.expander(f"Component {comp_i + 1}"):
        #             if self.N_p <= 1:
        #                 self.nonlimiting_filmDiff_per_comp.append(
        #                     st.selectbox(f"Infinite film diffusion rate",
        #                                  ["No", "Yes"],
        #                                  key=f"nonlimiting_filmDiff_comp_{comp_i}") == "Yes"
        #                 )
        #             else:
        #                 self.nonlimiting_filmDiff_per_comp.append(self.nonlimiting_filmDiff)
        #             if self.has_binding:
        #                 self.req_binding_per_comp.append(
        #                     st.selectbox(f"Binding kinetics mode",
        #                                  ["Kinetic", "Rapid-equilibrium"],
        #                                  key=f"req_binding_comp_{comp_i}") == "Rapid-equilibrium"
        #                 )
        #                 if self.N_p <= 1 and resolution == "1D":
        #                     self.has_surfDiff_per_comp.append(
        #                         st.selectbox(f"Surface diffusion",
        #                                      ["No", "Yes"],
        #                                      key=f"has_surfDiff_comp_{comp_i}") == "Yes"
        #                     )
        #                 else:
        #                     self.has_surfDiff_per_comp.append(self.has_surfDiff if self.N_p > 1 else False)
        #                 self.has_mult_bnd_states_per_comp.append(
        #                     st.selectbox(f"Multiple bound states",
        #                                  ["No", "Yes"],
        #                                  key=f"has_mult_bnd_states_comp_{comp_i}") == "Yes"
        #                 )
        #             else:
        #                 self.req_binding_per_comp.append(False)
        #                 self.has_surfDiff_per_comp.append(False)
        #                 self.has_mult_bnd_states_per_comp.append(False)

        # nonlimiting_filmDiff_j = (
        #     self.nonlimiting_filmDiff_per_partype[j]
        #     if self.nonlimiting_filmDiff_per_partype is not None
        #     else self.nonlimiting_filmDiff
        # )
        # has_surfDiff_j = (
        #     self.has_surfDiff_per_partype[j]
        #     if self.has_surfDiff_per_partype is not None
        #     else self.has_surfDiff
        # )

        self.particle_models.append(
            Particle(
                geometry=geometry,
                var_format=self.var_format,
                resolution=resolution,
                has_core=has_core,
                has_binding=self.has_binding,
                req_binding=self.req_binding,
                has_mult_bnd_states=self.has_mult_bnd_states,
                has_surfDiff=has_surfDiff_j,
                nonlimiting_filmDiff=nonlimiting_filmDiff_j,
                interstitial_volume_resolution=self.resolution,
                column_type=self.column_type,
                single_partype=(self.N_p == 1),
                binding_model=self.binding_model,
                has_reaction_liquid=self.has_reaction_particle_liquid,
                has_reaction_solid=self.has_reaction_particle_solid,
                req_reaction_liquid=self.req_reaction_particle_liquid,
                req_reaction_solid=self.req_reaction_particle_solid
            )
        )

    def available_CADET_Core(self):
        """
        Return the availability status of the model in CADET-Core.
    
        Returns
        -------
        int
            -1 if model is not present,
             0 if model can be approximated,
             1 if model is present.
        """
        if self.has_angular_coordinate:
            return -1

        availability = 1

        par1D = None
        
        if self.N_p > 0:
            
            par1D = False

            for p in self.particle_models:

                availability = int(availability and p.available_CADET_Core())
                par1D = par1D or p.resolution == "1D"
                
            # no particles with pore diffusion but no film diffusion (yet)
            if par1D and self.nonlimiting_filmDiff:
                return 0
        
        # tank model
        if not self.has_axial_coordinate:
            # only 0D particles with non limiting film diffusion available
            if self.N_p > 0:
                if not par1D and self.nonlimiting_filmDiff:
                    return 1
                elif par1D and not self.nonlimiting_filmDiff: # can be approximated using 1 FV cell
                    return 0
                else:
                    return -1
            else:
                return 1

        return availability

    def available_CADET_Process(self):
        """
        Return the availability status of the model in CADET-Process.
    
        Returns
        -------
        int
            -1 if model is not present,
             0 if model can be approximated,
             1 if model is present.
        """

        # All CADET-Core models supported except multiple particle type and 2D models
        
        core_availability = self.available_CADET_Core()
        
        if core_availability == -1:
            return -1
        
        if self.N_p > 1:
            return -1

        if self.has_radial_coordinate:
            return -1

        if self.column_type in ["Radial", "Frustum"]:
            return -1

        else:
            return core_availability

    def fill_vars_and_params(self):

        without_pores_ = self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D"

        state_deps = r"t"
        param_deps = r"\text{constant}"
        if self.resolution == "1D":
            if self.column_type == "Radial":
                state_deps = r"t, \rho"
            elif self.column_type == "Frustum":
                state_deps = r"t, x"
            else:
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
            {"Group" : 1, "Symbol": r"c^{\b}_i", "Description": r"bulk liquid concentration", "Unit": r"\frac{mol}{m^3}", "Dependence" : state_deps, "Domain": eq.int_vol_domain(self.resolution, column_type=self.column_type)},
            {"Group" : -1, "Symbol": r"T^{\mathrm{end}}", "Description": r"process end time", "Unit": r"s", "Dependence": r"\text{constant}", "Property": r" > 0"},
            {"Group" : -0.2, "Symbol": r"N^\mathrm{c}", "Description": r"number of components", "Unit": r"-", "Dependence": r"-", "Property": r"\in \mathbb{N}"},
            {"Group" : -0.1, "Symbol": r"i", "Description": r"component index", "Unit": r"s", "Dependence": r"-", "Property": r"-"},
            ]

        if self.has_axial_dispersion:
            if self.column_type == "Radial":
                disp_symbol = r"D^\mathrm{rad}_i"
                disp_desc = r"radial dispersion coefficient"
                disp_symbol_apparent = r"\tilde{D}^\mathrm{rad}_i"
                disp_desc_apparent = r"apparent radial dispersion coefficient"
            else:
                disp_symbol = r"D^\mathrm{ax}_i"
                disp_desc = r"axial dispersion coefficient"
                disp_symbol_apparent = r"\tilde{D}^\mathrm{ax}_i"
                disp_desc_apparent = r"apparent axial dispersion coefficient"
            if not without_pores_:
                self.vars_and_params.append({"Group" : 6, "Symbol": disp_symbol, "Description": disp_desc, "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})
            else:
                self.vars_and_params.append({"Group" : 6, "Symbol": disp_symbol_apparent, "Description": disp_desc_apparent, "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})
        if self.has_radial_dispersion:
            self.vars_and_params.append({"Group" : 6, "Symbol": r"D^\mathrm{rad}_i", "Description": r"radial dispersion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})
        if self.has_angular_dispersion:
            self.vars_and_params.append({"Group" : 6, "Symbol": r"D^\mathrm{ang}_i", "Description": r"angular dispersion coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": param_deps_comp, "Property": r"\geq 0"})

        if self.resolution == "0D":
            self.vars_and_params.append({"Group" : 2, "Symbol": r"Q^\mathrm{in}", "Description": r"volumetric flow rate into the tank", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"\geq 0"})
            self.vars_and_params.append({"Group" : 2, "Symbol": r"Q^\mathrm{out}", "Description": r"volumetric flow rate out of the tank", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"\geq 0"})
            if self.has_filter:
                self.vars_and_params.append({"Group" : 2, "Symbol": r"Q^\mathrm{filter}", "Description": r"volumetric flow rate out of the tank (solvent only)", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"\geq 0"})
        elif self.column_type == "Radial":
            self.vars_and_params.append({"Group" : 0, "Symbol": r"\rho", "Description": r"radial coordinate", "Unit": r"m", "Dependence": r"\text{independent variable}", "Property": r"\in (R^{\mathrm{in}}, R^{\mathrm{out}})"})
            self.vars_and_params.append({"Group" : -1, "Symbol": r"R^{\mathrm{in}}", "Description": r"inner cylinder radius", "Unit": r"m", "Dependence": r"\text{constant}", "Property": r" > 0"})
            self.vars_and_params.append({"Group" : -1, "Symbol": r"R^{\mathrm{out}}", "Description": r"outer cylinder radius", "Unit": r"m", "Dependence": r"\text{constant}", "Property": r" > R^{\mathrm{in}}"})
            self.vars_and_params.append({"Group" : 2, "Symbol": r"Q", "Description": r"volumetric flow rate", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"> 0"})
            self.vars_and_params.append({"Group" : 5, "Symbol": r"v", "Description": r"velocity coefficient", "Unit": r"\frac{m^2}{s}", "Dependence": r"\text{constant}", "Property": r":= \frac{Q}{2 \pi L}"})
        elif self.column_type == "Frustum":
            self.vars_and_params.append({"Group" : 0, "Symbol": r"x", "Description": r"axial coordinate", "Unit": r"m", "Dependence": r"\text{independent variable}", "Property": r"\in (0, L)"})
            self.vars_and_params.append({"Group" : -1, "Symbol": r"L", "Description": r"length of column", "Unit": r"m", "Dependence": r"\text{constant}", "Property": r" > 0"})
            self.vars_and_params.append({"Group" : -1, "Symbol": r"R^0", "Description": r"column radius at inlet", "Unit": r"m", "Dependence": r"\text{constant}", "Property": r" > 0"})
            self.vars_and_params.append({"Group" : -1, "Symbol": r"R^L", "Description": r"column radius at outlet", "Unit": r"m", "Dependence": r"\text{constant}", "Property": r" > 0"})
            self.vars_and_params.append({"Group" : 3, "Symbol": r"R", "Description": r"column radius function", "Unit": r"m", "Dependence": r"x", "Property": r"(x) = R^0 + \frac{R^L - R^0}{L} x"})
            self.vars_and_params.append({"Group" : 2, "Symbol": r"Q", "Description": r"volumetric flow rate", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r"> 0"})
            self.vars_and_params.append({"Group" : 5, "Symbol": r"v", "Description": r"velocity coefficient", "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}", "Property": r":= \frac{Q}{\pi}"})
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
                self.vars_and_params.append({"Group" : 7, "Symbol": symbol_name_, "Description": r"film diffusion coefficient", "Unit": r"\frac{m}{s}", "Dependence": r"i" if self.N_p<=1 else r"j,i", "Property": r"\geq 0"})

        if self.has_reaction_bulk and not self.req_reaction_bulk:
            self.vars_and_params.append({"Group" : 8, "Symbol": r"f^{\mathrm{react},\b}_{i}", "Description": r"bulk liquid phase reaction function", "Unit": r"\frac{mol}{m^3 \cdot s}", "Dependence": r"\vec{c}^{\b}; i"})
            self.vars_and_params.append({"Group" : 8.1, "Symbol": r"\vec{c}^{\b}", "Description": r"bulk liquid components vector", "Unit": r"[\frac{mol}{m^3}]", "Dependence": state_deps})

        if self.has_reaction_bulk and self.req_reaction_bulk:
            self.vars_and_params.append({"Group" : 8, "Symbol": r"g^{\mathrm{react,eq},\b}_{k}", "Description": r"bulk liquid phase equilibrium constraint function", "Unit": r"\frac{mol}{m^3}", "Dependence": r"\vec{c}^{\b}; k"})
            self.vars_and_params.append({"Group" : 8.1, "Symbol": r"N^{\mathrm{react,eq},\b}", "Description": r"number of rapid-equilibrium bulk reactions", "Unit": r"-", "Dependence": r"-"})
            self.vars_and_params.append({"Group" : 8.1, "Symbol": r"M^{\b}", "Description": r"conserved moiety matrix for bulk reactions", "Unit": r"-", "Dependence": r"-"})
            self.vars_and_params.append({"Group" : 8.1, "Symbol": r"\vec{c}^{\b}", "Description": r"bulk liquid components vector", "Unit": r"[\frac{mol}{m^3}]", "Dependence": state_deps})

        for var_ in self.vars_and_params:
            var_["Symbol"] = format_variables(var_["Symbol"], self.var_format)
            
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

            if self.has_reaction_bulk and not self.req_reaction_bulk:
                equation += " + " + eq.bulk_reaction_term()

        else:
            equation = eq.bulk_time_derivative(r"\varepsilon^{\mathrm{c}}") if not without_pores_ else eq.bulk_time_derivative()
            if without_pores_:
                equation += r" + " + eq.solid_time_derivative(r"\varepsilon^{\mathrm{t}}")

            needs_linebreak = self.has_radial_dispersion or self.has_angular_dispersion
            eq_sign = " &= " if needs_linebreak else " = "

            # Select convection and dispersion terms based on column geometry
            if self.column_type == "Radial":
                convection_func = eq.radial_flow_convection
                dispersion_func = eq.radial_flow_dispersion
            elif self.column_type == "Frustum":
                convection_func = eq.frustum_convection
                dispersion_func = eq.frustum_dispersion
            else:
                convection_func = eq.axial_convection
                dispersion_func = eq.axial_dispersion

            equation += eq_sign + convection_func() if without_pores_ else eq_sign + convection_func(r"\varepsilon^{\mathrm{c}}")

            if self.has_axial_dispersion:
                equation += " + " + dispersion_func(r"\varepsilon^{\mathrm{c}}")
            if self.has_radial_dispersion:
                equation += r" \nonumber \\ & + " + eq.radial_dispersion(r"\varepsilon^{\mathrm{c}}")
            if self.has_angular_dispersion:
                equation += r" \nonumber \\ & + " + eq.angular_dispersion(r"\varepsilon^{\mathrm{c}}")

            if self.N_p == 0:  # remove occurencies of porosity, which is just constant one in this case
                equation = re.sub(r"\\varepsilon^{\\mathrm{c}}", "", re.sub(
                    r"\\left\( \\varepsilon^{\\mathrm{c}} c^{\\l}_i \\right\)", r"c^{\\l}_i", equation))

        par_added = 0
        for par_uniq in self.par_unique_intV_contribution_counts.keys():

            if self.dev_mode:
                equation += eq.int_filmDiff_term(
                    Particle(
                        self.particle_models[par_added].geometry, self.particle_models[par_added].has_core, self.var_format, self.particle_models[par_added].resolution
                    ),
                    1 + par_added, par_added + self.par_unique_intV_contribution_counts[par_uniq], self.N_p == 1, self.particle_models[par_added].nonlimiting_filmDiff, self.particle_models[par_added].has_surfDiff
                )
            else:
                equation += eq.int_filmDiff_term(
                    Particle(
                        self.particle_models[0].geometry, self.particle_models[0].has_core, self.var_format, self.particle_models[0].resolution
                    ),
                    1, r"N^{\mathrm{p}}", self.N_p == 1, self.particle_models[0].nonlimiting_filmDiff, self.particle_models[0].has_surfDiff
                )

            par_added += self.par_unique_intV_contribution_counts[par_uniq]

        if self.has_reaction_bulk and not self.req_reaction_bulk and self.resolution != "0D":
            equation += " + " + eq.bulk_reaction_term()

        if self.resolution == "0D":
            equation = re.sub(
                r"\\left\(1 - \\varepsilon^{\\mathrm{c}} \\right\)", r"V^s", equation)

        equation = r"""\begin{align}
""" + equation + r""",
\end{align}"""

        return equation

    def interstitial_volume_bc(self):
        if not self.resolution == "0D":
            return eq.int_vol_BC(self.resolution, self.has_axial_dispersion, self.column_type)
        else:
            return None

    def particle_equations(self):

        eqs = {}
        boundary_conditions = {}

        for par_type in self.par_type_counts.keys():

            eqs[par_type] = eq.particle_transport(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=par_type.nonlimiting_filmDiff,
                                                  has_surfDiff=par_type.has_surfDiff, has_binding=self.has_binding, req_binding=self.req_binding, has_mult_bnd_states=self.has_mult_bnd_states,
                                                  has_reaction_liquid=self.has_reaction_particle_liquid and not self.req_reaction_particle_liquid,
                                                  has_reaction_solid=self.has_reaction_particle_solid and not self.req_reaction_particle_solid,
                                                  binding_model=self.binding_model)

            boundary_conditions[par_type] = eq.particle_boundary(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=par_type.nonlimiting_filmDiff,
                                                                 has_surfDiff=par_type.has_surfDiff, has_binding=self.has_binding, req_binding=self.req_binding, has_mult_bnd_states=self.has_mult_bnd_states)

            if self.N_p == 1:
                eqs[par_type] = re.sub(",j", "", eqs[par_type])
                eqs[par_type] = re.sub("j,", "", eqs[par_type])
                eqs[par_type] = re.sub("_{j}", "", eqs[par_type])
                boundary_conditions[par_type] = re.sub(",j", "", boundary_conditions[par_type])
                boundary_conditions[par_type] = re.sub("j,", "", boundary_conditions[par_type])
                boundary_conditions[par_type] = re.sub("_{j}", "", boundary_conditions[par_type])

        return eqs, boundary_conditions

    def has_per_component_config(self):
        """Return True if per-component configuration is active and settings differ across components."""
        return self.req_binding_per_comp is not None and self.N_c > 0

    def component_groups(self):
        """Group components by their shared per-component settings.

        Returns a list of dicts, each containing:
          - 'components': list of 1-based component indices
          - 'req_binding', 'nonlimiting_filmDiff', 'has_surfDiff', 'has_mult_bnd_states': the shared settings
        """
        if not self.has_per_component_config():
            return None

        groups = {}
        for i in range(self.N_c):
            key = (
                self.req_binding_per_comp[i],
                self.nonlimiting_filmDiff_per_comp[i],
                self.has_surfDiff_per_comp[i],
                self.has_mult_bnd_states_per_comp[i],
            )
            groups.setdefault(key, []).append(i + 1)  # 1-based index

        result = []
        for (req_b, nlf, sd, mbs), comps in groups.items():
            result.append({
                'components': comps,
                'req_binding': req_b,
                'nonlimiting_filmDiff': nlf,
                'has_surfDiff': sd,
                'has_mult_bnd_states': mbs,
            })
        return result

    def particle_equations_for_group(self, group):
        """Generate particle equations using per-component group settings."""
        eqs = {}
        boundary_conditions = {}

        for par_type in self.par_type_counts.keys():

            eqs[par_type] = eq.particle_transport(
                par_type, singleParticle=self.N_p == 1,
                nonlimiting_filmDiff=group['nonlimiting_filmDiff'],
                has_surfDiff=group['has_surfDiff'],
                has_binding=self.has_binding,
                req_binding=group['req_binding'],
                has_mult_bnd_states=group['has_mult_bnd_states'],
                has_reaction_liquid=self.has_reaction_particle_liquid,
                has_reaction_solid=self.has_reaction_particle_solid,
                binding_model=self.binding_model)

            boundary_conditions[par_type] = eq.particle_boundary(
                par_type, singleParticle=self.N_p == 1,
                nonlimiting_filmDiff=group['nonlimiting_filmDiff'],
                has_surfDiff=group['has_surfDiff'],
                has_binding=self.has_binding,
                req_binding=group['req_binding'],
                has_mult_bnd_states=group['has_mult_bnd_states'])

            if self.N_p == 1:
                eqs[par_type] = re.sub(",j", "", eqs[par_type])
                eqs[par_type] = re.sub("j,", "", eqs[par_type])
                eqs[par_type] = re.sub("_{j}", "", eqs[par_type])
                boundary_conditions[par_type] = re.sub(",j", "", boundary_conditions[par_type])
                boundary_conditions[par_type] = re.sub("j,", "", boundary_conditions[par_type])
                boundary_conditions[par_type] = re.sub("_{j}", "", boundary_conditions[par_type])

        return eqs, boundary_conditions

    @staticmethod
    def format_component_set(components):
        """Format a list of component indices as a LaTeX set string."""
        if len(components) == 1:
            return r"$i = " + str(components[0]) + r"$"
        else:
            return r"$i \in \{" + ", ".join(str(c) for c in components) + r"\}$"

    def has_per_partype_config(self):
        """Return True if per-particle-type configuration is active."""
        return self.nonlimiting_filmDiff_per_partype is not None and self.N_p > 1

    def partype_indices(self, par_type):
        """Return 1-based indices of particle types matching the given Particle."""
        return [j + 1 for j, p in enumerate(self.particle_models) if p == par_type]

    @staticmethod
    def format_partype_set(indices):
        """Format a list of particle type indices as a LaTeX set string."""
        if len(indices) == 1:
            return r"$j = " + str(indices[0]) + r"$"
        else:
            return r"$j \in \{" + ", ".join(str(j) for j in indices) + r"\}$"

    def domain_interstitial(self, with_time_domain=True):
        return r"$" + eq.int_vol_domain(self.resolution, with_time_domain=with_time_domain, column_type=self.column_type) + r"$"

    def domain_particle(self):
        if self.N_p > 0:
            return eq.full_particle_conc_domain(self.resolution, self.particle_models[0].resolution, self.particle_models[0].has_core, False, False, column_type=self.column_type)
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
            else:
                return "Continuously Stirred Tank"

        else:

            if self.has_angular_coordinate:
                model_name = "3D "
            elif self.has_radial_coordinate:
                model_name = "2D "
            # elif self.has_axial_coordinate: # default case, no name prefix !
            else:
                model_name = ""

            if self.N_p > 0:

                if self.column_type == "Radial":
                    model_name += "Radial Flow "
                elif self.column_type == "Frustum":
                    model_name += "Frustum "

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
                if self.column_type == "Radial":
                    model_name += "Radial "
                elif self.column_type == "Frustum":
                    model_name += "Frustum "
                if self.has_axial_dispersion or self.has_radial_dispersion or self.has_angular_dispersion:
                    model_name += "Dispersive "
                model_name += "Plug Flow"  # Reactor if we have reactions

        if self.has_binding and self.binding_model != "Arbitrary":
            model_name += f" with {self.binding_model} binding"

        return model_name

    def model_assumptions(self):

        asmpts = {
            "General model assumptions": eq.HRM_asmpt(self.N_p, self.nonlimiting_filmDiff, self.has_binding, self.has_surfDiff, self.resolution),
            "Specific model assumptions": eq.int_vol_continuum_asmpt(self.resolution, self.N_p, self.nonlimiting_filmDiff, self.column_type) +
            (eq.particle_asmpt(
                self.particle_models[0].resolution, self.has_surfDiff) if self.N_p > 0 else [])
        }

        if self.nonlimiting_filmDiff:
            asmpts["Specific model assumptions"].append(
                r"the film around the particles does not limit mass transfer. That is, we assume $k^{\mathrm{f}}_i = \infty$)")
        if self.req_binding:
            asmpts["Specific model assumptions"].append(
                r"adsorption and desorption happen on a much faster time scale than the other mass transfer processes (e.g., convection, diffusion). Hence, we consider them to be equilibrated instantly, that is, to always be in (local) equilibrium")

        if self.req_reaction_bulk or self.req_reaction_particle_liquid or self.req_reaction_particle_solid:
            phases = []
            if self.req_reaction_bulk:
                phases.append("bulk liquid")
            if self.req_reaction_particle_liquid:
                phases.append("particle liquid")
            if self.req_reaction_particle_solid:
                phases.append("particle solid")
            asmpts["Specific model assumptions"].append(
                r"reactions in the " + ", ".join(phases) + r" phase happen on a much faster time scale than the other mass transfer processes (e.g., convection, diffusion). Hence, the reaction equilibria are attained instantaneously and the system is reduced through conserved moieties")

        for idx in range(0, len(asmpts["Specific model assumptions"])):
            
            asmpts["Specific model assumptions"][idx] = format_variables(asmpts["Specific model assumptions"][idx], self.var_format)
            
        for idx in range(0, len(asmpts["General model assumptions"])):
            
            asmpts["General model assumptions"][idx] = format_variables(asmpts["General model assumptions"][idx], self.var_format)
            
        if not self.binding_model == "Arbitrary":
            
            asmpts.update({"Binding model assumptions": eq.binding_model_assumptions(self.binding_model)})
            
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
