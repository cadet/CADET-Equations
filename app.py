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

from typing import List
from typing import Literal

import equations


# %%

def rerender_variables(input_str: str, var_format: int):
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

# %% src


# immutable & hashable dataclass -> unique particle type counter
@dataclass(frozen=True)
class Particle:
    geometry: Literal["Sphere", "Cylinder", "Slab"]
    hasCore: bool
    resolution: Literal["1D", "0D"]
    surface_volume_ratio: float = None
    # volume fraction ?
    # binding -> is_kinetic, nBound

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
    # counts per unique particle type (geometry, hasCore, resolution)
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

    def __post_init__(self):

        # Configure interstitial column transport
        st.sidebar.write("Configure interstitial volume model")
        col1, col2 = st.columns(2)

        with col1:
            self.resolution = re.search(r'\dD', st.sidebar.selectbox("Column resolution", [
                                        "1D (axial coordinate)", "0D (Homogeneous Tank)", "2D (axial and radial coordinate)", "3D (axial, radial and angular coordinate)"], key="column_resolution")).group()

        valid_resolutions = {"3D", "2D", "1D", "0D"}

        if self.resolution not in valid_resolutions:
            raise ValueError(
                f"Invalid resolution: {self.resolution}. Must be one of {valid_resolutions}.")
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
            "Number of components", key="N_c", min_value=1, step=1) if advanced_mode_ else -1

        if self.has_axial_coordinate:
            with col2:
                self.has_axial_dispersion = st.sidebar.selectbox(
                    "Add axial Dispersion", ["No", "Yes"], key="has_axial_dispersion") == "Yes"

        if self.dev_mode:
            if self.has_radial_coordinate:
                self.has_radial_dispersion = st.sidebar.selectbox(
                    "Add radial Dispersion", ["Yes", "No"], key="has_radial_dispersion") == "Yes"
            if self.has_angular_coordinate:
                self.has_angular_dispersion = st.sidebar.selectbox(
                    "Add angular Dispersion", ["Yes", "No"], key="has_angular_dispersion") == "Yes"

        # Configure particle model
        st.sidebar.write("Configure particle model")

        self.N_p = st.sidebar.number_input("Number of particle types", key="N_p", min_value=0, step=1) if dev_mode_ else int(
            st.sidebar.selectbox("Add particles", ["No", "Yes"], key="add_particles") == "Yes")

        if self.N_p > 0:

            self.particle_models = []

            for j in range(self.N_p):

                if self.N_p > 1:
                    st.sidebar.write(f"Configure particle type {j + 1}")

                col1, col2 = st.columns(2)

                with col1:

                    if dev_mode_:  # multiple particle types
                        resolution = re.search(r'\dD', st.sidebar.selectbox(f"Select spatial resolution of particle type {j + 1}", [
                                               "1D (radial coordinate)", "0D (homogeneous)"], key=f"parType_{j+1}_resolution")).group()
                        hasCore = st.sidebar.selectbox(f"Choose if particle type {j + 1} is a core-shell particle (i.e. " + r"$R_\mathrm{pc} > 0$)", [
                                                       "No core-shell", "Has core-shell"], key=f"parType_{j+1}_hasCore") == "Has core-shell"
                        geometry = st.sidebar.selectbox(f"Select geometry of particle type {j + 1}", [
                                                        "Sphere", "Cylinder", "Slab"], key=f"parType_{j+1}__geometry")
                    else:
                        resolution = re.search(r'\dD', st.sidebar.selectbox(f"Select spatial resolution of particles", [
                                               "1D (radial coordinate)", "0D (homogeneous)"], key="particle_resolution")).group()
                        hasCore = st.sidebar.selectbox(f"Add impenetrable core-shell (i.e. " + r"$R_\mathrm{pc} > 0$)", [
                                                       "No", "Yes"], key=f"parType_{j+1}_hasCore") == "Yes" if (resolution == "1D" and advanced_mode_) else False
                        geometry = "Sphere"

                with col2:
                    if j == 0:  # todo make this configurable for every particle type
                        self.nonlimiting_filmDiff = st.sidebar.selectbox(
                            "Non-limiting film diffusion", ["No", "Yes"], key="nonlimiting_filmDiff") == "Yes"

                if j == 0:  # todo make this configurable for every particle type
                    # Configure binding model
                    st.sidebar.write("Configure binding model")

                    self.has_binding = st.sidebar.selectbox(
                        "Add binding", ["No", "Yes"], key="has_binding") == "Yes"

                    if self.has_binding:

                        col1, col2 = st.columns(2)

                        with col1:
                            self.req_binding = st.sidebar.selectbox("Binding kinetics mode", [
                                                                    "Kinetic", "Rapid-equilibrium"], key="req_binding") == "Rapid-equilibrium"
                            self.has_mult_bnd_states = st.sidebar.selectbox("Add multiple bound states", [
                                                                            "No", "Yes"], key="has_mult_bnd_states") == "Yes" if advanced_mode_ else False
                        with col2:
                            self.has_surfDiff = st.sidebar.selectbox("Add surface diffusion", [
                                                                     "No", "Yes"], key="has_surfDiff") == "Yes" if resolution == "1D" else False

                self.particle_models.append(
                    Particle(
                        geometry=geometry,
                        resolution=resolution,
                        hasCore=hasCore
                    )
                )

            # We need to count and thus sort particle_models:
            #  Particle types need individual particle equations, which is why we count them
            #  Only specific differences lead to changes in the interstitial volume equations: geometry if kinetic film diffusion. else geometry + resolution
            #  Sorting by (geometry, resolution) ensures we are using the same indices across the interstitial volume and particle equations.
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

    def interstitial_volume_equation(self):

        if self.resolution == "0D":
            equation = r"""
    \frac{\mathrm{d}V^{\b}}{\mathrm{d}t} &= Q_{\mathrm{in}} - Q_{\mathrm{out}} - Q_{\mathrm{filter}},
    \\
    \frac{\mathrm{d}}{\mathrm{d} t} \left( V^{\b} c^{\b}_i \right)"""

            if self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D":
                equation += r" + V^{\p} \varepsilon_\mathrm{p} \frac{\partial c^{b}_i}{\partial t}"
                if self.req_binding:
                    equation += r" + V^{\p} \left( 1 - \varepsilon_\mathrm{p} \right) \frac{\partial c^{\s}_i}{\partial t}"

            equation += r"&= Q_{\mathrm{in}} c^{\b}_{\mathrm{in},i} - Q_{\mathrm{out}} c^{\b}_i"

            if self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D" and not self.req_binding:
                equation += r" - V^{\p} \left( 1 - \varepsilon_\mathrm{p} \right) \frac{\partial c^{\s}_i}{\partial t}"

        else:
            equation = equations.bulk_time_derivative_eps
            if self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D" and self.req_binding:
                equation += r" + " + equations.solid_time_derivative_eps

            equation += " = " + equations.axial_convection_eps

            if self.has_axial_dispersion:
                equation += " + " + equations.axial_dispersion_eps
            if self.has_radial_dispersion:
                equation += " + " + equations.radial_dispersion_eps
            if self.has_angular_dispersion:
                equation += " + " + equations.angular_dispersion_eps

            if self.N_p == 0:  # remove occurencies of porosity, which is just constant one in this case
                equation = re.sub(r"\\varepsilon_{\\mathrm{c}}", "", re.sub(
                    r"\\left\( \\varepsilon_{\\mathrm{c}} c^{\\l}_i \\right\)", r"c^{\\l}_i", equation))

        # if self.nonlimiting_filmDiff and 1Dparticle # entscheidende faktoren sind particle resolution und filmDiffMode. the following loop has thus to change
        par_added = 0
        for par_uniq in self.par_unique_intV_contribution_counts.keys():
            if self.nonlimiting_filmDiff:
                equation += equations.int_filmDiff_term(Particle(par_uniq[0], False, par_uniq[1]), 1 + par_added, par_added +
                                                        self.par_unique_intV_contribution_counts[par_uniq], self.N_p == 1, self.nonlimiting_filmDiff)
            else:
                equation += equations.int_filmDiff_term(Particle(par_uniq, False, "0D"), 1 + par_added, par_added +
                                                        self.par_unique_intV_contribution_counts[par_uniq], self.N_p == 1, self.nonlimiting_filmDiff)

            par_added += self.par_unique_intV_contribution_counts[par_uniq]

        if self.resolution == "0D":
            equation = re.sub(
                r"\\left\(1 - \\varepsilon_{\\mathrm{c}} \\right\)", r"V^s", equation)

        if self.nonlimiting_filmDiff and self.has_binding and self.particle_models[0].resolution == "0D" and self.req_binding:
            equation = re.sub(
                r"\\varepsilon_{\\mathrm{c}}", r"\\varepsilon_{\\mathrm{t}}", equation)

        equation = r"""\begin{align}
""" + equation + r""",
\end{align}"""

        return equation

    def interstitial_volume_bc(self):
        if not self.resolution == "0D":
            return equations.int_vol_BC(self.resolution, self.has_axial_dispersion)
        else:
            return None

    def particle_equations(self):

        eqs = {}
        boundary_conditions = {}

        for par_type in self.par_type_counts.keys():

            eqs[par_type] = equations.particle_transport(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=self.nonlimiting_filmDiff,
                                                         has_surfDiff=self.has_surfDiff, has_binding=self.has_binding, req_binding=self.req_binding, has_mult_bnd_states=self.has_mult_bnd_states)
            eqs[par_type] = eqs[par_type]

            boundary_conditions[par_type] = equations.particle_boundary(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=self.nonlimiting_filmDiff,
                                                                        has_surfDiff=self.has_surfDiff, has_binding=self.has_binding, req_binding=self.req_binding, has_mult_bnd_states=self.has_mult_bnd_states)
            boundary_conditions[par_type] = boundary_conditions[par_type]

        return eqs, boundary_conditions

    def model_name(self):

        if self.resolution == "0D":
            return "Finite Bath" if self.N_p > 0 else "Continuously Stirred Tank"

        if self.has_angular_coordinate:
            model_name = "3D"
        elif self.has_radial_coordinate:
            model_name = "2D"
        # elif self.has_axial_coordinate: # default case, no name prefix !
        else:
            model_name = ""

        if self.N_p > 0:

            if self.N_p > 1:
                if not any(par.geometry != "Sphere" for par in self.particle_models):
                    model_name += " PSD"  # particle-size distribution
                else:
                    # particle-type distribution # TODO use when different kinds of geometry or binding
                    model_name += " PTD"

            if self.particle_models[0].resolution == "1D":
                model_name += "General Rate Model"

                if self.has_surfDiff:
                    model_name += "  with surface diffusion"

            else:
                model_name += "Lumped Rate Model"

                if self.nonlimiting_filmDiff:
                    model_name += " without Pores"
        else:
            if self.has_axial_dispersion or self.has_radial_dispersion or self.has_angular_dispersion:
                model_name += "Dispersive "
            model_name += "Plug Flow"  # Reactor if we have reactions

        return model_name

    def model_assumptions(self):

        asmpts = {
            "General model assumptions": equations.HRM_asmpt(self.N_p, self.nonlimiting_filmDiff, self.has_binding, self.has_surfDiff, self.resolution),
            "Specific model assumptions": equations.int_vol_continuum_asmpt(self.resolution, self.N_p, self.nonlimiting_filmDiff) +
            (equations.particle_asmpt(
                self.particle_models[0].resolution, self.has_surfDiff) if self.N_p > 0 else [])
        }

        if self.nonlimiting_filmDiff:
            asmpts["Specific model assumptions"].append(
                r"the film around the particles does not limit mass transfer (i.e. $k_{\mathrm{f},i} \to \infty$)")
        if self.req_binding:
            asmpts["Specific model assumptions"].append(
                r"adsorption and desorption happen on a much faster time scale than the other mass transfer processes (e.g., convection, diffusion). Hence, we consider them to be equilibrated instantly, that is, to always be in (local) equilibrium")

        return asmpts


# %% Streamlit UI
st.sidebar.title(
    "CADET-Equations: Packed-Bed Chromatography Model Equation Generator")
st.sidebar.write(
    "Configure a chromatography model to get the corresponding governing equations.")

# File uploader for JSON file
uploaded_file = st.sidebar.file_uploader(
    "Define the model configuration below or upload configuration file", type="json")

if uploaded_file is not None:
    # Load JSON data
    config = json.load(uploaded_file)

    # Update Streamlit session state
    for key, value in config.items():

        st.session_state[key] = value

    st.sidebar.success("Configuration applied from JSON file!")

# User configuration of the model

var_format_ = st.sidebar.selectbox("Select format (e.g. $c^s$ or $q$ as the solid phase concentration)", [
                                   "CADET", "Legacy"], key="var_format")

advanced_mode_ = st.sidebar.selectbox("Advanced setup options (enables e.g. multiple bound states)", [
                                      "Off", "On"], key="advanced_mode") == "On"
if advanced_mode_:
    dev_mode_ = st.sidebar.selectbox("Developer setup options (not tested! Enables e.g. multiple particle types)", [
                                     "Off", "On"], key="dev_mode") == "On"
    advanced_mode_ = True if dev_mode_ else False
else:
    dev_mode_ = False

column_model = Column(dev_mode=dev_mode_, advanced_mode=advanced_mode_)

# %% Display equations

file_content = []  # used to export model to files

if st.toggle("Show Model Assumptions", key="model_assumptions"):

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

interstitial_volume_eq = column_model.interstitial_volume_equation()

nComp_list = r"$i\in\{" + ", ".join(str(i) for i in range(1, column_model.N_c + 1)) + \
    r"\}$" if column_model.N_c > 0 else r"$i\in\{1, \dots, N_c\}$"
nPar_list = ', '.join(str(j) for j in range(1, column_model.N_p + 1))

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
    intro_str = r"Consider a cylindrical column "

if column_model.N_p == 0:
    write_and_save(intro_str + r"filled with liquid phase.")
elif column_model.N_p == 1:
    # TODO particle geometries?
    write_and_save(intro_str + r"packed with spherical particles.")
else:
    write_and_save(intro_str + r"""packed with $N_{\mathrm{p}}\geq 0$ different particle types indexed by $j \in \{1, \dots, N_{\mathrm{p}}\}$ and distributed according to the volume fractions $d_j \colon (0, R_{\mathrm{c}}) \times (0,L) \to [0, 1]$.
	For all $(\rho,z,\varphi) \in (0, R_{\mathrm{c}}) \times (0,L) \times [0,2\pi)$ the volume fractions satisfy""")
    write_and_save(r"""
	\begin{equation*}
		\sum_{j=1}^{N_{\mathrm{p}}} d_j(\rho, z, \varphi) = 1.
	\end{equation*}

""", as_latex=True)


if column_model.resolution == "0D":
    write_and_save(
        r"The evolution of the liquid volume $V^{\l}\colon (0, T_{\mathrm{end}}) \to \Reals$ and the concentrations $c_i\colon (0, T_{\mathrm{end}}) \to \Reals$ of the components in the tank is governed by")
    write_and_save(interstitial_volume_eq, as_latex=True)
else:
    eq_type = "convection-diffusion-reaction"
    write_and_save(r"In the interstitial volume, mass transfer is governed by the following " + eq_type +
                   " equations in " + equations.int_vol_domain[column_model.resolution] + r" and for all components " + nComp_list)
    write_and_save(interstitial_volume_eq, as_latex=True)
    write_and_save("with boundary conditions")
    write_and_save(column_model.interstitial_volume_bc(), as_latex=True)


# add a variable and parameter description for the interstitial volume equation
diff_vars = r"D_{\mathrm{ax},i}" if column_model.has_axial_dispersion else ""
diff_vars += r"D_{\mathrm{rad},i}" if column_model.has_radial_dispersion else ""
diff_vars += r"D_{\mathrm{ang},i}" if column_model.has_angular_dispersion else ""

if column_model.resolution == "1D":
    diff_string = r", $D_{\mathrm{ax},i} \geq 0$ is the lumped axial diffusion coefficient" if column_model.has_axial_dispersion else ""
elif column_model.resolution == "2D":
    tmp_str = "are the lumped diffusion coefficients in axial and radial direction" if column_model.has_axial_dispersion else "is the lumped radial diffusion coefficient"
    diff_string = r", " + diff_vars + tmp_str
elif column_model.resolution == "3D":
    diff_string = r", " + diff_vars + \
        r"are the lumped diffusion coefficients in axial, radial and angular direction"
elif column_model.resolution == "0D":
    diff_string = ""

# if column_model.N_p == 0:
sol_vars_int_vol_eq = r"c^{\l}_i \colon " + re.sub(
    r"\$", "", equations.int_vol_domain[column_model.resolution]) + r" \to \mathbb{R}"
# else:# todo req bnd mit 0D par
if column_model.particle_models is None or column_model.particle_models[0].resolution == "0D":
    sol_vars_int_vol_eq = r"c^{\l}_i, c^{\p}_i \colon " + re.sub(
        r"\$", "", equations.int_vol_domain[column_model.resolution]) + r" \to \mathbb{R}"
else:
    sol_vars_int_vol_eq += r", c^{\p}_i \colon " + re.sub(r"\$", "", equations.particle_domain(
        column_model.resolution, column_model.particle_models[0].resolution, column_model.particle_models[0].hasCore, with_par_index=column_model.N_p, with_time_domain=True)) + r" \to \mathbb{R}"

sol_vars_int_vol_eq = sol_vars_int_vol_eq
sol_vars_int_vol_eq_names = "is the liquid concentration" if column_model.N_p == 0 else "are the bulk (interstitial volume) and particle liquid concentrations"
# if column_model.N_p > 0 and column_model.particle_models[0].resolution == "0D" and column_model.has_req_binding: # TODO add for req. binding
#     sol_vars_int_vol_eq_names = "are the liquid and solid concentrations"

filmDiff_str = r", $k_{\mathrm{f},i}\geq 0$ is the film diffusion coefficient" if column_model.N_p > 0 else ""

porosity_str = r", $\varepsilon_{\mathrm{c}} \colon " + re.sub(r"\(0, T_\\mathrm\{end\}\) \\times", "", re.sub(
    r"\$", "", equations.int_vol_domain[column_model.resolution])) + r" \mapsto (0, 1)$ is the interstitial column porosity" if column_model.N_p > 0 else ""
if column_model.nonlimiting_filmDiff and column_model.has_binding and column_model.particle_models[0].resolution == "0D" and column_model.req_binding:
    porosity_str = re.sub(
        r"\\varepsilon_{\\mathrm{c}}", r"\\varepsilon_{\\mathrm{t}}", porosity_str)
    porosity_str = re.sub("interstitial column porosity",
                          "total porosity", porosity_str)

u_string = r"$u>0$ is the interstitial velocity" if not column_model.resolution == "0D" else ""
if column_model.N_p > 0:
    write_and_save(
        r"Here, $" + sol_vars_int_vol_eq + r"$, " + sol_vars_int_vol_eq_names + r" of component $i \in \{1, \dots, N_{\mathrm{c}}\}$, $T_{\mathrm{end}} > 0$ is the simulation end time, " + u_string +
        diff_string + filmDiff_str + porosity_str +
        r", $c_{\mathrm{in},i}\colon " +
        equations.int_vol_inlet_domain[column_model.resolution] +
        r" \to \mathbb{R}$ is a given inlet concentration profile.",
    )

if column_model.N_p > 0:
    tmp = r"and particle types $j\in\{" + str(nPar_list) + r"\}$"
else:
    tmp = r""

if column_model.N_p > 0:

    particle_eq, particle_bc = column_model.particle_equations()

    # in this case, we dont have a particle model. this configuration is still allowed for educational purpose.
    if not (not column_model.has_binding and column_model.nonlimiting_filmDiff and column_model.particle_models[0].resolution == "0D"):
        write_and_save(
            r"In the particle models, mass transfer is governed by the following diffusion-reaction equations,")

    tmp = 0
    for par_type in column_model.par_type_counts.keys():

        # in this case, we dont have a particle model. this configuration is still allowed for educational purpose.
        if not column_model.has_binding and column_model.nonlimiting_filmDiff and par_type.resolution == "0D":
            break

        tmp_textblock = "in " + equations.particle_domain(column_model.resolution, par_type.resolution,
                                                          par_type.hasCore, with_par_index=True, with_time_domain=True) + " for "

        if column_model.N_p > 1:

            tmp_textblock += r"particle types $j\in\{" + \
                nPar_list + r"\}$ and components " + nComp_list
        else:
            tmp_textblock += r"components " + nComp_list

        write_and_save(tmp_textblock)
        write_and_save(particle_eq[par_type], as_latex=True)
        tmp += column_model.par_type_counts[par_type]

        if column_model.has_binding:

            if par_type.resolution == "0D":
                diffusion_description = "."
            else:
                diffusion_description = r", $D_i^{\p}$ is the particle diffusion coefficient for component $i$" if column_model.has_surfDiff else r", $D_i^{\p}, D_i^{\s}\geq 0$ are the pore and surface diffusion coefficients for component $i$."

            write_and_save(
                r"Here, $c^{\s}_i \colon " + re.sub(r"\$", "", equations.particle_domain(column_model.resolution, column_model.particle_models[0].resolution, column_model.particle_models[
                                                    0].hasCore, with_par_index=column_model.N_p, with_time_domain=True)) + r" \to \mathbb{R}$ is the solid phase concentration" + diffusion_description,
            )

        if not particle_bc[par_type] == "":
            write_and_save("The boundary conditions are given by")
            write_and_save(particle_bc[par_type], as_latex=True)

write_and_save(
    "Initial values for all solution variables (concentrations) are defined at $t = 0$.")

latex_string = [
    r"""\documentclass{article}
""",
    r"""\usepackage{amssymb,amsmath,mleftright}
""",
    r"""\begin{document}
"""
]
latex_string.extend(file_content + [r"\end{document}"])
latex_string = "\n".join(latex_string)

st.download_button("Download .tex", latex_string, "model.tex", "text/plain")

if st.button("Generate PDF", key="generate_pdf"):
    with tempfile.TemporaryDirectory() as temp_dir:
        tex_path = f"{temp_dir}/model.tex"
        pdf_path = f"{temp_dir}/model.pdf"

        # Write LaTeX content to a temporary file
        with open(tex_path, "w") as f:
            f.write(latex_string)

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

if st.button("Generate configuration file", key="generate_config"):

    no_config_state = ["generate_pdf", "generate_config"]

    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        json_path = f"{temp_dir}/model.json"

        config = {key: st.session_state[key]
                  for key in st.session_state if key not in no_config_state}

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
