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

column_model = Column(dev_mode=dev_mode_, advanced_mode=advanced_mode_, var_format=var_format_)

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
    
    renderer_write_and_save(output, var_format_, file_content, as_latex)

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

            cur_par_count = 0
            for par_type in column_model.par_type_counts.keys():

                if not column_model.has_binding and group['nonlimiting_filmDiff'] and par_type.resolution == "0D":
                    break

                if dev_mode_:
                    nPar_list = ', '.join(str(j) for j in range(cur_par_count, column_model.par_type_counts[par_type] + 1))
                else:
                    nPar_list = r"$j\in\{1, \dots, N^{\mathrm{p}}\}$"

                eq_type_ = "reaction" if column_model.particle_models[0].resolution == "0D" else "diffusion-reaction"
                tmp_str = r" and all particle sizes " + nPar_list if column_model.N_p > 1 else r""

                whatComp = eq.primary_binding_eq_what_comps(column_model.binding_model)

                write_and_save(
                    "For component(s) " + comp_set_str + ", mass transfer in the particles is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True) + tmp_str)

                group_eq, group_bc = column_model.particle_equations_for_group(group)

                write_and_save(group_eq[par_type], as_latex=True)

                if not group_bc[par_type] == "":
                    write_and_save("with boundary conditions")
                    write_and_save(group_bc[par_type], as_latex=True)

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
                        "For the salt component, mass transfer is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True) + tmp_str)

                    # The salt component is in rapid equilibrium binding mode
                    salt_group = dict(group)
                    salt_group['req_binding'] = True
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

                cur_par_count += column_model.par_type_counts[par_type]

    else:
        # Standard mode: single set of equations for all components
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

            whatComp = eq.primary_binding_eq_what_comps(column_model.binding_model)

            write_and_save(
                "In the particles, mass transfer is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True) + r" and for " + whatComp + " components" + tmp_str)

            write_and_save(particle_eq[par_type], as_latex=True)

            if not particle_bc[par_type] == "":
                write_and_save("with boundary conditions")
                write_and_save(particle_bc[par_type], as_latex=True)

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
            if column_model.binding_model == "SMA" and column_model.has_binding:

                PTD_ = column_model.PTD and column_model.N_p > 1
                write_and_save(r"The number of available binding sites $\bar{q}_0$ is given by")
                write_and_save(r"\begin{align}" + eq.sma_free_binding_sites(PTD=PTD_) + r".\end{align}", as_latex=True)

                write_and_save(
                    "For the salt component, mass transfer is governed by " + eq_type_ + " equations in " + eq.full_particle_conc_domain(column_model.resolution, par_type.resolution, par_type.has_core, with_par_index=False, with_time_domain=True) + tmp_str)

                # The salt component is in rapid equilibrium binding mode
                tmpReqBnd = column_model.req_binding
                column_model.req_binding = True
                tmpBndModel = column_model.binding_model
                column_model.binding_model = "SMA_salt"
                particle_eq_salt, particle_bc_salt = column_model.particle_equations()
                column_model.req_binding = tmpReqBnd
                column_model.binding_model = tmpBndModel

                if PTD_:
                    write_and_save(re.sub(r"_{i}", r"_{0}", particle_eq_salt[par_type]), as_latex=True)
                    re.sub(r"_{i}", r"_{0}", particle_bc_salt[par_type])
                else:
                    write_and_save(re.sub(r"_{j,i}", r"_{j,0}", particle_eq_salt[par_type]), as_latex=True)
                    re.sub(r"_{j,i}", r"_{j,0}", particle_bc_salt[par_type])

                if not particle_bc_salt[par_type] == "":
                    write_and_save(r"where the counter-ion concentration $c^{\s}_0$ satisfies the electroneutrality constraint. Boundary conditions are")
                    write_and_save(particle_bc_salt[par_type], as_latex=True)
                else:
                   write_and_save(r"where the counter-ion concentration $c^{\s}_0$ satisfies the electroneutrality constraint.")

            cur_par_count += column_model.par_type_counts[par_type]

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
