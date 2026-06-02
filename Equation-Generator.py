# -*- coding: utf-8 -*-
"""
Main script for the Streamlit app that generates equations for packed-bed chromatography models based on user configuration.
The app allows users to specify model parameters, assumptions, and configurations, and then generates the corresponding governing equations in LaTeX format.
Users can also download the generated equations as a .tex file or a PDF, and export their configuration as a JSON file.
"""

import streamlit as st
import re
from collections import Counter
from dataclasses import dataclass, field
import json
import subprocess
import tempfile
from pathlib import Path
import pandas as pd

from typing import List
from typing import Literal
from typing import Optional
from collections import OrderedDict

from src import equations as eq
from src import load_CADET_h5
from src import ui_config
from src.utils import format_variables
from src.renderer import availability_badge_html, write_and_save as renderer_write_and_save
from src.model_column import Column
from src.model_crystallization import Crystallization
from src.generate_template import generate_unit_operation_script, generate_crystallization_script
from src.handler_cite import load_bibliography, cite_html, cite, render_references

# %% Streamlit UI

st.logo("images/logo_CADET.png", size="large", link=None, icon_image=None)

st.set_page_config(
    page_title=r"CADET-Equations",
    page_icon=r":material/biotech:", #":material/modeling:",
    layout="wide",
)

# Custom CSS to fix equation tag overlap and improve responsiveness
st.markdown("""
<style>
    /* Prevent equation numbers from overlapping with long equations */
    .katex-display {
        overflow-x: auto !important;
        overflow-y: hidden !important;
        text-align: center !important;
        padding: 0.5em 3em !important; /* Safety padding for the tag at the far right */
        position: relative !important;
        display: flex !important;
        justify-content: center !important;
        align-items: center !important;
    }
    .katex-display > .katex {
        display: inline-block !important;
        white-space: nowrap !important;
        position: relative !important;
    }
    /* Position the tag at the right-most point of the display container */
    .katex .tag {
        position: absolute !important;
        right: -2.5em !important; /* Pull it slightly inside the padding, or pin to edge of formula */
        top: 50% !important;
        transform: translateY(-50%) !important;
    }
</style>
""", unsafe_allow_html=True)

st.sidebar.title(
    "CADET-Equations: Model Equation Generator")
st.sidebar.write(
    "Configure a model to get the corresponding governing equations.")

# File uploader for JSON file
uploaded_file = st.sidebar.file_uploader(
    "Upload a configuration file (optional, see documentation)", type=["json", "h5"])

if uploaded_file is not None:

    uploaded_file_name = uploaded_file.name.lower()
    
    config = None

    if uploaded_file_name.endswith(".json"):

        json_data = json.load(uploaded_file)
        if isinstance(json_data, dict) and ('column_resolution' in json_data or 'model_type' in json_data):
            config = json_data
            st.sidebar.success("Uploaded configuration applied!")
        else:
            st.sidebar.error("Uploaded JSON is not a valid CADET-Equations configuration.")

    elif uploaded_file_name.endswith(".h5"):

        config = load_CADET_h5.get_config_from_CADET_h5(uploaded_file,
        str(st.sidebar.number_input("Unit index in CADET file", key=r"h5_input_unit_index", min_value=-1, max_value=999, step=1, value=-1)).zfill(3)
        )

    if config is not None:
        # Update Streamlit session state
        for key, value in config.items():

            st.session_state[key] = value

# ---------------------------------
# User configuration of the model
# ---------------------------------

import streamlit as st

# %% Sidebar selection of model family

if "model_type" not in st.session_state:
    st.session_state["model_type"] = "Chromatography"

# Custom button styling
st.markdown("""
<style>

/* Primary button */
div.stButton > button[kind="primary"] {
    background-color: #023d6b;
    color: white;
    border: none;
}

/* Primary button hover */
div.stButton > button[kind="primary"]:hover {
    background-color: #145a86;
    color: white;
}

</style>
""", unsafe_allow_html=True)


#%% Developer mode toggle enabling untested features

if "dev_mode" not in st.session_state:
    st.session_state["dev_mode"] = False

if st.sidebar.button(
    "Developer mode",
    key="dev_mode_button",
    use_container_width=True,
    type="primary" if st.session_state["dev_mode"] else "secondary",
):
    st.session_state["dev_mode"] = not st.session_state["dev_mode"]
    st.rerun()

st.sidebar.caption(
    (
        "Developer mode enables untested features such as Crystallization and further chromatography model variants."
    )
)

dev_mode_ = st.session_state["dev_mode"]

if dev_mode_:
    st.warning(
        "**Warning:** Models in developer mode are "
        "work in progress and have not been thoroughly verified. Please "
        "double-check all equations and parameters before relying on them.",
        icon="⚠️"
    )

if dev_mode_:
    
    st.sidebar.markdown("### Model Family")

    col1, col2 = st.sidebar.columns(2)

    if col1.button(
        "Chromatography",
        key="model_type_chromatography_button",
        use_container_width=True,
        type="primary" if st.session_state["model_type"] == "Chromatography" else "secondary",
    ):
        st.session_state["model_type"] = "Chromatography"
        st.rerun()

    if col2.button(
        "Crystallization",
        key="model_type_crystallization_button",
        use_container_width=True,
        type="primary" if st.session_state["model_type"] == "Crystallization" else "secondary",
    ):
        st.session_state["model_type"] = "Crystallization"
        st.rerun()

    st.write("Selected:", st.session_state["model_type"])

    model_type_ = st.session_state["model_type"]

    st.sidebar.caption(
        (
            "Chromatography: Convection, dispersion, fixed-bed particles, binding, and reactions."
            if model_type_ == "Chromatography"
            else "Crystallization: Population balance model, nucleation, growth, aggregation and breakage."
        )
    )

else:
    st.session_state["model_type"] = "Chromatography"
    model_type_ = st.session_state["model_type"]

# %% Variable format: CADET vs. Legacy

var_format_ = st.sidebar.selectbox("Select parameter format", [
                                   "CADET", "Legacy"], key=r"var_format")

# %% Display equations

file_content = []  # used to export model to files

bibliography_entries = load_bibliography("CITATION.bib")
used_citation_keys = []

# The following function is used to both print the output and collect it to later generate and export output files
def write_and_save(output: str, as_latex: bool = False):

    renderer_write_and_save(output, var_format_, file_content, as_latex)


def write_html_and_save(html_output: str, export_output: str):
    """Render HTML in Streamlit and save plain text for LaTeX export."""
    export_output = format_variables(export_output, var_format_)
    if export_output is not None:
        file_content.append(export_output)
    st.markdown(html_output, unsafe_allow_html=True)


if model_type_ == "Crystallization":

    cry_model = Crystallization(var_format=var_format_)

    show_eq_description = st.toggle("Show equation description", key=r"show_eq_description", value=True)

    if st.toggle("Show Model Assumptions", key=r"model_assumptions"):
        asmpts = cry_model.model_assumptions()
        for key in asmpts.keys():
            st.write(key + ":\n" + "\n".join(f"- {item}" for item in asmpts[key]))
            file_content.append(
                key + r""":
\begin{itemize}
""" + "\n".join(f"\\item {item}" for item in asmpts[key]) + r"""
\end{itemize}

""")

    if st.toggle("Show symbol table", key=r"sym_table"):
        df = pd.DataFrame(cry_model.vars_and_params)
        df[['Symbol', "Dependence", 'Unit']] = df[['Symbol', "Dependence", 'Unit']].map(
            lambda x: f"${x}$" if isinstance(x, str) else x)
        df = df.sort_values(by=r"Group").reset_index()
        st.table(df[['Symbol', 'Description', 'Dependence', 'Unit']])

    # Title
    st.write("### " + cry_model.model_name())
    file_content.append(r"\section*{" + cry_model.model_name() + r"}")

    html = availability_badge_html("CADET-Core", cry_model.available_CADET_Core())
    html += availability_badge_html("CADET-Process", cry_model.available_CADET_Process())
    st.markdown(html, unsafe_allow_html=True)

    has_primary = cry_model.has_primary_formation
    has_agg = cry_model.has_aggregation
    has_frag = cry_model.has_fragmentation

    write_html_and_save(
        "The crystallization models stated in the following were established in "
        + cite_html("ZHANG2024108612", bibliography_entries, used_citation_keys)
        + " and "
        + cite_html("ZHANG2025108860", bibliography_entries, used_citation_keys)
        + ".",
        "The crystallization models stated in the following were established in "
        + cite("ZHANG2024108612", bibliography_entries, used_citation_keys)
        + " and "
        + cite("ZHANG2025108860", bibliography_entries, used_citation_keys)
        + "."
    )

    if cry_model.column_type == "CSTR":
        write_and_save(r"Consider a continuously stirred tank reactor (CSTR) with volume $V(t)$, "
                       r"observed over a time interval $(0, T^{\mathrm{end}})$. "
                       r"The particle population is described by the number density $n(t, x)$ "
                       r"over the internal coordinate (particle size) $x\in(x_0, \infty)$.")
    else:
        disp_str = " with axial dispersion" if cry_model.has_axial_dispersion else ""
        write_and_save(r"Consider a dispersive plug flow reactor (DPFR) of length $L > 0$" + disp_str +
                       r", observed over a time interval $(0, T^{\mathrm{end}})$. "
                       r"The particle population is described by the number density $n(t, x, z)$ "
                       r"over the internal coordinate (particle size) $x\in(x_0, \infty)$ and axial position $z$.")

    # PBE
    write_and_save(r"The evolution of the number density is governed by the population balance equation")
    if cry_model.column_type == "CSTR":
        write_and_save(eq.cry_pbe_cstr(has_primary, cry_model.has_growth_dispersion,
                                       has_agg, has_frag), as_latex=True)
    else:
        write_and_save(eq.cry_pbe_dpfr(has_primary, cry_model.has_axial_dispersion,
                                       cry_model.has_growth_dispersion,
                                       has_agg, has_frag), as_latex=True)

    # Boundary conditions for PBE (internal coordinate)
    if has_primary:
        bc_internal = eq.cry_pbe_bc_internal(has_primary, cry_model.has_growth_dispersion)
        if bc_internal:
            write_and_save("with nucleation kinetics and regularity boundary conditions in the internal coordinate")
            write_and_save(bc_internal, as_latex=True)

    # Boundary conditions for PBE (external coordinate, DPFR only)
    if cry_model.column_type == "DPFR":
        write_and_save("Danckwerts boundary conditions are employed in the external (axial) coordinate")
        write_and_save(eq.cry_pbe_bc_external_dpfr(cry_model.has_axial_dispersion), as_latex=True)

    # Mass balance
    write_and_save(r"The solute mass balance is given by")
    if cry_model.column_type == "CSTR":
        write_and_save(eq.cry_mass_balance_cstr(has_primary), as_latex=True)
        write_and_save(r"Evolution of the reactor’s volume is governed by")
        write_and_save(eq.cry_volume_cstr(), as_latex=True)
    else:
        write_and_save(eq.cry_mass_balance_dpfr(has_primary, cry_model.has_axial_dispersion), as_latex=True)
        write_and_save("with Danckwerts boundary conditions")
        write_and_save(eq.cry_solute_bc_dpfr(cry_model.has_axial_dispersion), as_latex=True)

    # Constitutive equations
    if has_primary:
        write_and_save(r"The constitutive relations for growth and nucleation are defined as follows. "
                       r"The relative supersaturation is")
        write_and_save(r"\begin{align}" + eq.cry_supersaturation() + r". \end{align}", as_latex=True)

        write_and_save("The growth rate is")
        write_and_save(r"\begin{align}" + eq.cry_growth_rate(cry_model.size_dependent_growth) + r". \end{align}", as_latex=True)

        write_and_save("The primary nucleation rate is")
        write_and_save(r"\begin{align}" + eq.cry_primary_nucleation() + r". \end{align}", as_latex=True)

        if cry_model.has_secondary_nucleation:
            write_and_save("The secondary nucleation rate is")
            write_and_save(r"\begin{align}" + eq.cry_secondary_nucleation() + r". \end{align}", as_latex=True)
            write_and_save("where the suspension density is")
            write_and_save(r"\begin{align}" + eq.cry_suspension_density() + r". \end{align}", as_latex=True)

        write_and_save("The total nucleation rate is")
        if cry_model.has_secondary_nucleation:
            write_and_save(r"\begin{align}" + eq.cry_total_nucleation() + r". \end{align}", as_latex=True)
        else:
            write_and_save(r"\begin{align} B_0 = B_p. \end{align}", as_latex=True)

    # Aggregation details
    if has_agg:
        write_and_save(r"The aggregation birth and death terms are")
        write_and_save(eq.cry_aggregation_birth_death(), as_latex=True)
        write_and_save("The aggregation kernel is defined by a " + cry_model.aggregation_kernel_name() + " kernel")
        write_and_save(eq.cry_aggregation_kernel(cry_model.aggregation_kernel_index), as_latex=True)

    # Fragmentation details
    if has_frag:
        write_and_save(r"The fragmentation birth and death terms are")
        write_and_save(eq.cry_fragmentation_birth_death(), as_latex=True)
        write_and_save("The selection (fragmentation rate) function is")
        write_and_save(r"\begin{align}" + eq.cry_selection_function() + r". \end{align}", as_latex=True)
        write_and_save("The breakage probability density function is")
        write_and_save(r"\begin{align}" + eq.cry_breakage_function() + r". \end{align}", as_latex=True)

    if show_eq_description:
        write_and_save("Here, " + cry_model.vars_params_description())

    write_and_save("Consistent initial values for all solution variables are defined at $t = 0$.")


else: # Chromatography model family

    if dev_mode_:
        advanced_mode_ = True
    else:
        advanced_mode_ = st.sidebar.selectbox("Advanced options (enables e.g. particle size distribution)", [
                                            "Off", "On"], key=r"advanced_mode") == "On"

    column_model = Column(dev_mode=dev_mode_, advanced_mode=advanced_mode_, var_format=var_format_)

    show_eq_description = st.toggle("Show equation description", key=r"show_eq_description", value=True)

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

    # Title
    st.write("### " + column_model.model_name())
    file_content.append(r"\section*{" + column_model.model_name() + r"}")

    #%% CADET model availability badge

    html = availability_badge_html(
        "CADET-Core",
        column_model.available_CADET_Core()
    )

    html += availability_badge_html(
        "CADET-Process",
        column_model.available_CADET_Process()
    )

    st.markdown(html, unsafe_allow_html=True)

    #%% Continue with model

    if column_model.resolution == "0D":
        intro_str = r"Consider a continuously stirred tank "
    elif column_model.column_type == "Radial":
        intro_str = r"Consider a hollow cylindrical column with inner radius $R^{\mathrm{in}} > 0$ and outer radius $R^{\mathrm{out}} > R^{\mathrm{in}}$ "
    elif column_model.column_type == "Frustum":
        intro_str = r"Consider a conical frustum column of length $L > 0$ with inlet radius $R^0 > 0$ and outlet radius $R^L > 0$ "
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
        component_groups_intV = column_model.component_groups()
        if column_model.interstitial_groups_differ_in_film_diff(component_groups_intV):
            # Per-component mode: film diffusion settings differ, show separate bulk equations per group
            for group in component_groups_intV:
                nlf = group['nonlimiting_filmDiff_per_partype'][0]
                comp_set_str = column_model.format_component_set(group['components'])

                eq_type = "convection"
                if not column_model.resolution == "1D" or column_model.has_axial_dispersion:
                    eq_type += "-diffusion"
                if column_model.N_p > 0 and not nlf:
                    eq_type += "-reaction"

                write_and_save(
                    "For component(s) " + comp_set_str +
                    r", in the interstitial volume, mass transfer is governed by the following " + eq_type +
                    " equations in " + column_model.domain_interstitial())
                write_and_save(column_model.interstitial_volume_equation_for_group(group), as_latex=True)

            write_and_save("with boundary conditions")
            write_and_save(column_model.interstitial_volume_bc(), as_latex=True)
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

    if column_model.has_reaction_bulk and column_model.req_reaction_bulk:
        write_and_save(
            r"The bulk liquid phase reactions are in rapid equilibrium. "
            r"The resulting overdetermined system of $N^{\mathrm{c}}$ component transport equations is reduced through conserved moieties. "
            r"Specifically, $N^{\mathrm{react,eq},\b}$ of the component equations are replaced by algebraic equilibrium constraints")
        write_and_save(r"""
\begin{align}
""" + eq.req_reaction_bulk_constraint() + r""", \quad k = 1, \ldots, N^{\mathrm{react,eq},\b}.
\end{align}
""", as_latex=True)
        write_and_save(
            r"The remaining $N^{\mathrm{c}} - N^{\mathrm{react,eq},\b}$ equations are replaced by conserved moiety equations, "
            r"obtained by left-multiplying the component transport equations with the conserved moiety matrix "
            r"$M^{\b} \in \mathbb{R}^{(N^{\mathrm{c}} - N^{\mathrm{react,eq},\b}) \times N^{\mathrm{c}}}$, "
            r"whose rows span the left null space of the stoichiometric matrix")
        write_and_save(r"""
\begin{align}
""" + eq.conserved_moiety_equation_bulk() + r""", \quad l = 1, \ldots, N^{\mathrm{c}} - N^{\mathrm{react,eq},\b}.
\end{align}
""", as_latex=True)
        write_and_save(r"By construction, the reaction terms cancel in the conserved moiety equations.")

    if column_model.N_p > 0:

        component_groups = column_model.component_groups()

        if component_groups is not None and len(component_groups) > 1:
            # Per-component mode: generate separate equations for each group of components
            # that share the same settings
            for group in component_groups:
                comp_set_str = column_model.format_component_set(group['components'])

                for par_type in column_model.par_type_counts.keys():

                    orig_idx = column_model._original_partype_indices[par_type][0] - 1
                    if not column_model.has_binding and group['nonlimiting_filmDiff_per_partype'][orig_idx] and par_type.resolution == "0D":
                        break

                    if dev_mode_:
                        par_indices = column_model.partype_indices(par_type)
                        nPar_list = column_model.format_partype_set(par_indices)
                    else:
                        nPar_list = r"$j\in\{1, \dots, N^{\mathrm{p}}\}$"

                    eq_type_ = "reaction" if column_model.particle_models[0].resolution == "0D" else "diffusion-reaction"
                    tmp_str = r" and all particle types " + nPar_list if column_model.N_p > 1 else r""

                    whatComp = eq.primary_binding_eq_what_comps(column_model.binding_model)

                    write_and_save(
                        "For component(s) " + comp_set_str + ", mass transfer in the particles is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True, column_type=column_model.column_type) + tmp_str)

                    group_eq, group_bc = column_model.particle_equations_for_group(group)

                    write_and_save(group_eq[par_type], as_latex=True)

                    if not group_bc[par_type] == "":
                        write_and_save("with boundary conditions")
                        write_and_save(group_bc[par_type], as_latex=True)
                    
                    if column_model.binding_model not in ["Arbitrary", None] and eq.binding_model_references(column_model.binding_model, bibliography_entries, used_citation_keys) is not None:
                        write_html_and_save(
                            "Further details on the binding model can be found in " + eq.binding_model_references(column_model.binding_model, bibliography_entries, used_citation_keys) + ".",
                            "Further details on the binding model can be found in " + eq.binding_model_references(column_model.binding_model, bibliography_entries, used_citation_keys) + "."
                        )

                    if show_eq_description:
                        write_and_save("Here, " + column_model.particle_models[0].vars_params_description())

                    # Rapid-equilibrium reactions in the particle phase
                    if column_model.has_reaction_particle_liquid and column_model.req_reaction_particle_liquid:
                        write_and_save(
                            r"The particle liquid phase reactions are in rapid equilibrium. "
                            r"The system is reduced through conserved moieties: "
                            r"$N^{\mathrm{react,eq},\p}$ algebraic equilibrium constraints")
                        write_and_save(r"""
\begin{align}
""" + eq.req_reaction_particle_liquid_constraint(column_model.N_p == 1) + r""", \quad k = 1, \ldots, N^{\mathrm{react,eq},\p}
\end{align}
""", as_latex=True)
                        write_and_save(
                            r"replace $N^{\mathrm{react,eq},\p}$ of the component equations. "
                            r"The remaining $N^{\mathrm{c}} - N^{\mathrm{react,eq},\p}$ equations are conserved moiety equations")
                        write_and_save(r"""
\begin{align}
""" + eq.conserved_moiety_equation_particle_liquid(column_model.N_p == 1) + r""", \quad l = 1, \ldots, N^{\mathrm{c}} - N^{\mathrm{react,eq},\p}.
\end{align}
""", as_latex=True)

                    if column_model.has_reaction_particle_solid and column_model.req_reaction_particle_solid:
                        write_and_save(
                            r"The particle solid phase reactions are in rapid equilibrium. "
                            r"The system is reduced through conserved moieties: "
                            r"$N^{\mathrm{react,eq},\s}$ algebraic equilibrium constraints")
                        write_and_save(r"""
\begin{align}
""" + eq.req_reaction_particle_solid_constraint(column_model.N_p == 1) + r""", \quad k = 1, \ldots, N^{\mathrm{react,eq},\s}
\end{align}
""", as_latex=True)
                        write_and_save(
                            r"replace $N^{\mathrm{react,eq},\s}$ of the component equations. "
                            r"The remaining $N^{\mathrm{c}} - N^{\mathrm{react,eq},\s}$ equations are conserved moiety equations")
                        write_and_save(r"""
\begin{align}
""" + eq.conserved_moiety_equation_particle_solid(column_model.N_p == 1) + r""", \quad l = 1, \ldots, N^{\mathrm{c}} - N^{\mathrm{react,eq},\s}.
\end{align}
""", as_latex=True)

                    # SMA additional equations
                    if column_model.binding_model == "SMA" and column_model.has_binding:

                        PTD_ = column_model.PTD and column_model.N_p > 1
                        write_and_save(r"The number of available binding sites $\bar{q}_0$ is given by")
                        write_and_save(r"\begin{align}" + eq.sma_free_binding_sites(PTD=PTD_) + r".\end{align}", as_latex=True)

                        write_and_save(
                            "For the salt component, mass transfer is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True, column_type=column_model.column_type) + tmp_str)

                        # The salt component is in rapid equilibrium binding mode
                        salt_group = dict(group)
                        salt_group['req_binding_per_partype'] = [True] * len(group['req_binding_per_partype'])
                        tmpBndModel = column_model.binding_model
                        column_model.binding_model = "SMA_salt"
                        particle_eq_salt, particle_bc_salt = column_model.particle_equations_for_group(salt_group)
                        column_model.binding_model = tmpBndModel

                        if PTD_:
                            write_and_save(re.sub(r"_{i}", r"_{0}", particle_eq_salt[par_type]), as_latex=True)
                        else:
                            write_and_save(re.sub(r"_{j,i}", r"_{j,0}", particle_eq_salt[par_type]), as_latex=True)

                        if not particle_bc_salt[par_type] == "":
                            write_and_save(r"where the counter-ion concentration $c^{\s}_0$ satisfies the electroneutrality constraint. Boundary conditions are")
                            write_and_save(particle_bc_salt[par_type], as_latex=True)
                        else:
                            write_and_save(r"where the counter-ion concentration $c^{\s}_0$ satisfies the electroneutrality constraint.")


        else:
            # Standard mode: single set of equations for all components
            particle_eq, particle_bc = column_model.particle_equations()

            # Build per-particle-type index mapping for display
            has_mixed_partypes = column_model.has_per_partype_config() and len(column_model.par_type_counts) > 1

            for par_type in column_model.par_type_counts.keys():

                # in this case, we dont have a particle model. this configuration is still allowed for educational purpose.
                if not column_model.has_binding and par_type.nonlimiting_filmDiff and par_type.resolution == "0D":
                    break

                if dev_mode_:
                    par_indices = column_model.partype_indices(par_type)
                    nPar_list = column_model.format_partype_set(par_indices)
                else:
                    nPar_list = r"$j\in\{1, \dots, N^{\mathrm{p}}\}$"

                eq_type_ = "reaction" if column_model.particle_models[0].resolution == "0D" else "diffusion-reaction"

                # Add per-particle-type label when particle types have different settings
                partype_label = ""
                if has_mixed_partypes:
                    partype_label = " for particle type(s) " + nPar_list
                    tmp_str = ""
                elif column_model.N_p > 1:
                    tmp_str = r" and all particle types " + nPar_list
                else:
                    tmp_str = ""

                whatComp = eq.primary_binding_eq_what_comps(par_type.binding_model)

                write_and_save(
                    "In the particles, mass transfer is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True, column_type=column_model.column_type) + r" and for " + whatComp + " components" + tmp_str + partype_label)

                write_and_save(particle_eq[par_type], as_latex=True)

                if not particle_bc[par_type] == "":
                    write_and_save("with boundary conditions")
                    write_and_save(particle_bc[par_type], as_latex=True)

                if column_model.binding_model not in ["Arbitrary", None] and eq.binding_model_references(column_model.binding_model, bibliography_entries, used_citation_keys) is not None:
                    write_html_and_save(
                        "Further details on the binding model can be found in " + eq.binding_model_references(column_model.binding_model, bibliography_entries, used_citation_keys) + ".",
                        "Further details on the binding model can be found in " + eq.binding_model_references(column_model.binding_model, bibliography_entries, used_citation_keys) + "."
                    )

                if show_eq_description:
                    write_and_save("Here, " + column_model.particle_models[0].vars_params_description())

                # Rapid-equilibrium reactions in the particle phase
                if column_model.has_reaction_particle_liquid and column_model.req_reaction_particle_liquid:
                    write_and_save(
                        r"The particle liquid phase reactions are in rapid equilibrium. "
                        r"The system is reduced through conserved moieties: "
                        r"$N^{\mathrm{react,eq},\p}$ algebraic equilibrium constraints")
                    write_and_save(r"""
\begin{align}
""" + eq.req_reaction_particle_liquid_constraint(column_model.N_p == 1) + r""", \quad k = 1, \ldots, N^{\mathrm{react,eq},\p}
\end{align}
""", as_latex=True)
                    write_and_save(
                        r"replace $N^{\mathrm{react,eq},\p}$ of the component equations. "
                        r"The remaining $N^{\mathrm{c}} - N^{\mathrm{react,eq},\p}$ equations are conserved moiety equations")
                    write_and_save(r"""
\begin{align}
""" + eq.conserved_moiety_equation_particle_liquid(column_model.N_p == 1) + r""", \quad l = 1, \ldots, N^{\mathrm{c}} - N^{\mathrm{react,eq},\p}.
\end{align}
""", as_latex=True)

                if column_model.has_reaction_particle_solid and column_model.req_reaction_particle_solid:
                    write_and_save(
                        r"The particle solid phase reactions are in rapid equilibrium. "
                        r"The system is reduced through conserved moieties: "
                        r"$N^{\mathrm{react,eq},\s}$ algebraic equilibrium constraints")
                    write_and_save(r"""
\begin{align}
""" + eq.req_reaction_particle_solid_constraint(column_model.N_p == 1) + r""", \quad k = 1, \ldots, N^{\mathrm{react,eq},\s}
\end{align}
""", as_latex=True)
                    write_and_save(
                        r"replace $N^{\mathrm{react,eq},\s}$ of the component equations. "
                        r"The remaining $N^{\mathrm{c}} - N^{\mathrm{react,eq},\s}$ equations are conserved moiety equations")
                    write_and_save(r"""
\begin{align}
""" + eq.conserved_moiety_equation_particle_solid(column_model.N_p == 1) + r""", \quad l = 1, \ldots, N^{\mathrm{c}} - N^{\mathrm{react,eq},\s}.
\end{align}
""", as_latex=True)

                # Some more complicated binding models require additional equations
                if par_type.binding_model == "SMA" and par_type.has_binding:

                    PTD_ = column_model.PTD and column_model.N_p > 1
                    write_and_save(r"The number of available binding sites $\bar{q}_0$ is given by")
                    write_and_save(r"\begin{align}" + eq.sma_free_binding_sites(PTD=PTD_) + r".\end{align}", as_latex=True)

                    write_and_save(
                        "For the salt component, mass transfer is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True, column_type=column_model.column_type) + tmp_str)

                    # The salt component is in rapid equilibrium binding mode
                    salt_eq, salt_bc = column_model.particle_salt_equations(par_type)

                    if PTD_:
                        write_and_save(re.sub(r"_{i}", r"_{0}", salt_eq), as_latex=True)
                    else:
                        write_and_save(re.sub(r"_{j,i}", r"_{j,0}", salt_eq), as_latex=True)

                    if salt_bc != "":
                        write_and_save(r"where the counter-ion concentration $c^{\s}_0$ satisfies the electroneutrality constraint. Boundary conditions are")
                        write_and_save(salt_bc, as_latex=True)
                    else:
                       write_and_save(r"where the counter-ion concentration $c^{\s}_0$ satisfies the electroneutrality constraint.")


    # %% Reaction model definition section
    if column_model.reaction_model != "Arbitrary":
        has_kinetic_bulk = column_model.has_reaction_bulk and not column_model.req_reaction_bulk
        has_kinetic_par_liq = column_model.has_reaction_particle_liquid and not column_model.req_reaction_particle_liquid
        has_kinetic_par_sol = column_model.has_reaction_particle_solid and not column_model.req_reaction_particle_solid

        if has_kinetic_bulk or has_kinetic_par_liq or has_kinetic_par_sol:

            phase_map = [
                (has_kinetic_bulk, "bulk", "bulk liquid"),
                (has_kinetic_par_liq, "particle_liquid", "particle liquid"),
                (has_kinetic_par_sol, "particle_solid", "particle solid"),
            ]
            for active, phase_key, phase_label in phase_map:
                if not active:
                    continue
                definition = eq.reaction_model_definition(column_model.reaction_model, phase_key)
                if definition is None:
                    continue
                main_eq, flux_eq = definition
                write_html_and_save(
                r"The " + phase_label + " phase reaction term is defined by " + column_model.reaction_model + " type reactions ("
                + eq.reaction_model_references(column_model.reaction_model, bibliography_entries, used_citation_keys) + ").",
                r"The " + phase_label + " phase reaction term is defined by " + column_model.reaction_model + " type reactions ("
                + eq.reaction_model_references(column_model.reaction_model, bibliography_entries, used_citation_keys) + ")."
                )
                write_and_save(r"""
\begin{align}
""" + main_eq + r"""
\end{align}
""", as_latex=True)
                write_and_save("where")
                write_and_save(r"""
\begin{align}
""" + flux_eq + r""".
\end{align}
""", as_latex=True)

    write_and_save("Consistent initial values for all solution variables (concentrations) are defined at $t = 0$.")

render_references(bibliography_entries, used_citation_keys, file_content)

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

if model_type_ == "Chromatography":
    st.session_state.template_script = generate_unit_operation_script(column_model)
    st.download_button("Download CADET-Python template (.py file)", st.session_state.template_script, "unit_operation.py", "text/x-python")
elif model_type_ == "Crystallization":
    st.session_state.template_script = generate_crystallization_script(cry_model)
    st.download_button("Download CADET-Python template (.py file)", st.session_state.template_script, "unit_operation.py", "text/x-python")

if st.button("Generate PDF", key=r"generate_pdf"):
    with tempfile.TemporaryDirectory() as temp_dir:
        tex_path = f"{temp_dir}/model.tex"
        pdf_path = f"{temp_dir}/model.pdf"

        # Write LaTeX content to a temporary file
        with open(tex_path, "w", encoding="utf-8") as f:
            f.write(st.session_state.latex_string)

        try:
            result = subprocess.run(
                [
                    "pdflatex",
                    "-interaction=nonstopmode",
                    "-halt-on-error",
                    "-file-line-error",
                    "-output-directory",
                    temp_dir,
                    tex_path,
                ],
                capture_output=True,
                text=True,
                cwd=temp_dir,
                timeout=90,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            partial = (exc.stderr or exc.stdout or "").strip()
            if partial:
                partial = partial[-4000:]
            st.error(
                "PDF generation timed out after 90 seconds. "
                "The LaTeX process likely waited on an error prompt.\n"
                + partial
            )
        else:
            if result.returncode != 0:
                log = (result.stderr or result.stdout or "").strip()
                if log:
                    log = log[-4000:]
                st.error(
                    "PDF generation failed. LaTeX output:\n" + log
                )
            elif not Path(pdf_path).exists():
                st.error(
                    "PDF generation failed: pdflatex finished but no PDF file was produced."
                )
            else:
                with open(pdf_path, "rb") as pdf_file:
                    st.download_button(
                        "Download PDF",
                        pdf_file,
                        "model.pdf",
                        "application/pdf",
                    )

if st.button("Generate configuration file", key=r"generate_config"):

    def sort_session_states(key):
        # a proper sorting is required for the tests, where we cannot apply a
        # global session state but must loop over the keys
        if key == "model_type":
            return -1
        elif key == "advanced_mode":
            return 0
        elif key == "dev_mode":
            return 1
        elif key == "column_type":
            return 1.5
        elif key == "column_resolution":
            return 2
        elif key in ("add_particles", "PSD"):
            return 3
        elif key == "has_binding":
            return 3.3
        elif key == "particle_resolution":
            return 3.5
        elif key in ("has_reaction_bulk", "has_reaction_particle_liquid", "has_reaction_particle_solid"):
            return 3.6
        elif key == "reaction_model":
            return 3.9
        elif key.startswith("cry_"):
            return 0.5
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
