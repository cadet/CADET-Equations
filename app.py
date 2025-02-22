# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import streamlit as st
import re
from collections import Counter
from dataclasses import dataclass, field

import subprocess
import tempfile

from typing import List
from typing import Literal


#%% Get the equations
exec(open("equations.py").read())

#%% helper functionality

@dataclass(frozen=True)  # immutable & hashable dataclass -> unique particle type counter
class Particle:
    geometry: Literal["Sphere", "Cylinder", "Slab"]
    hasCore: bool
    resolution: Literal["1D", "0D"]
    hasSurfDiff: bool = False
    surface_volume_ratio: float = None
    # volume fraction ?
    # binding -> is_kinetic, nBound

    def __post_init__(self):
        valid_geometries = {"Sphere", "Cylinder", "Slab"}
        valid_resolutions = {"1D", "0D"}

        if self.geometry not in valid_geometries:
            raise ValueError(f"Invalid geometry: {self.geometry}. Must be one of {valid_geometries}.")
        if self.resolution not in valid_resolutions:
            raise ValueError(f"Invalid resolution: {self.resolution}. Must be one of {valid_resolutions}.")
        
        if self.geometry == "Sphere":
            object.__setattr__(self, 'surface_volume_ratio', 3)
        elif self.geometry == "Cylinder":
            object.__setattr__(self, 'surface_volume_ratio', 2)
        elif self.geometry == "Slab":
            object.__setattr__(self, 'surface_volume_ratio', 1)

@dataclass
class Column:

    dev_mode: bool # dev mode including untested, unstable and wip features
    advanced_mode: bool # Ask for detailed parameter inputs
    resolution: Literal["1D", "2D", "3D"]
    N_c: int
    N_p: int

    has_axial_coordinate: bool = False
    has_radial_coordinate: bool = False
    has_angular_coordinate: bool = False
    has_axial_dispersion: bool = False
    has_radial_dispersion: bool = False
    has_angular_dispersion: bool = False

    nonlimiting_filmDiff: bool = False # List[bool] = None # TODO make this possible per component and particle type in advanced mode

    particle_models: List[Particle] = None
    par_type_counts: Counter[Particle] = field(default_factory=Counter) # counts per unique particle type (geometry, hasCore, resolution)
    par_unique_intV_contribution_counts: Counter[Particle] = field(default_factory=Counter) # puts particle types together that have a similar contribution to the interstitial volume equation and counts them
    has_surfDiff: bool = False
    has_binding: bool = True # and thus solid phase
    has_mult_bnd_states: bool = False
    # nCompReq: int
    # nCompKin: int
    # noPoresButReqBinding: bool

    def __post_init__(self):
        valid_resolutions = {"3D", "2D", "1D"}

        if self.resolution not in valid_resolutions:
            raise ValueError(f"Invalid resolution: {self.resolution}. Must be one of {valid_resolutions}.")
        if int(re.search("\\d", self.resolution).group()) > 0:
             self.has_axial_coordinate = True
             self.has_axial_dispersion = True
        if int(re.search("\\d", self.resolution).group()) > 1:
             self.has_radial_coordinate = True
             self.has_radial_dispersion = True # TODO make corresponding BC optional
        if int(re.search("\\d", self.resolution).group()) > 2:
             self.has_angular_coordinate = True
             self.has_angular_dispersion = True # TODO make corresponding BC optional
        
        if self.has_axial_coordinate:
            self.has_axial_dispersion = st.selectbox("Add axial Dispersion", ["No", "Yes"], key="has_axial_dispersion") == "Yes"

        if self.dev_mode:
             if self.has_radial_coordinate:
                self.has_radial_dispersion = st.selectbox("Add radial Dispersion", ["Yes", "No"], key="has_radial_dispersion") == "Yes"
             if self.has_angular_coordinate:
                self.has_angular_dispersion = st.selectbox("Add angular Dispersion", ["Yes", "No"], key="has_angular_dispersion") == "Yes"

        if self.N_p > 0:

            self.nonlimiting_filmDiff = st.selectbox("Non-limiting film diffusion", ["No", "Yes"], key="nonlimiting_filmDiff") == "Yes"

            self.has_binding = st.selectbox("Add binding", ["Yes", "No"]) == "Yes"
            self.has_mult_bnd_states = st.selectbox("Add multiple bound states", ["No", "Yes"], key="has_mult_bnd_states") == "Yes" if advanced_mode_ else False

            self.particle_models = []

            for j in range(self.N_p):

                self.particle_models.append(
                    Particle(
                        geometry=st.selectbox(f"Select geometry of particle type {j + 1}", ["Sphere", "Cylinder", "Slab"], key=f"parType_{j+1}__geometry") if self.dev_mode else "Sphere",
                        resolution=re.search(r'\dD', st.selectbox(f"Select spatial resolution of particle type {j + 1}", ["1D (radial coordinate)", "0D (homogeneous)"], key=f"parType_{j+1}_resolution")).group() if self.dev_mode else re.search(r'\dD', st.selectbox(f"Select spatial resolution of particles", ["1D (radial coordinate)", "0D (homogeneous)"], key="particle_resolution")).group(),
                        hasCore=st.selectbox(f"Choose if particle type {j + 1} is a core-shell particle (i.e." + r"$R_\mathrm{pc} > 0$)", ["No core-shell", "Has core-shell"], key=f"parType_{j+1}_hasCore") == "Has core-shell" if self.dev_mode else False
                        )
                )

            if any(particle.resolution == "1D" for particle in self.particle_models):
                self.has_surfDiff = st.selectbox("Add surface diffusion", ["No", "Yes"], key="has_surfDiff") == "Yes"
            else:
                self.has_surfDiff = False
            
            # We need to count and thus sort particle_models:
            #  Particle types need individual particle equations, which is why we count them
            #  Only specific differences lead to changes in the interstitial volume equations: geometry if kinetic film diffusion. else geometry + resolution
            #  Sorting by (geometry, resolution) ensures we are using the same indices across the interstitial volume and particle equations.
            self.particle_models = sorted(self.particle_models, key=lambda particle: (particle.geometry, particle.resolution))
            self.par_type_counts = Counter(self.particle_models)
            if self.nonlimiting_filmDiff:
                self.par_unique_intV_contribution_counts = Counter((particle.geometry, particle.resolution) for particle in self.particle_models)
            else:
                self.par_unique_intV_contribution_counts = Counter(particle.geometry for particle in self.particle_models)


    def interstitial_volume_equation(self):

        equation = bulk_time_derivative_eps + " = " + axial_convection_eps
        if self.has_axial_dispersion:
            equation += " + " + axial_dispersion_eps
        if self.has_radial_dispersion:
            equation += " + " + radial_dispersion_eps
        if self.has_angular_dispersion:
            equation += " + " + angular_dispersion_eps

        # if self.nonlimiting_filmDiff and 1Dparticle # entscheidende faktoren sind particle resolution und filmDiffMode. the following loop has thus to change
        par_added = 0
        for par_uniq in self.par_unique_intV_contribution_counts.keys():
                if self.nonlimiting_filmDiff:
                    equation += int_filmDiff_term(Particle(par_uniq[0], False, par_uniq[1]), 1 + par_added, par_added + self.par_unique_intV_contribution_counts[par_uniq], self.N_p == 1, self.nonlimiting_filmDiff)
                else:
                    equation += int_filmDiff_term(Particle(par_uniq, False, "0D"), 1 + par_added, par_added + self.par_unique_intV_contribution_counts[par_uniq], self.N_p == 1, self.nonlimiting_filmDiff)

                par_added += self.par_unique_intV_contribution_counts[par_uniq]

        equation = r"\begin{align}" + equation + r",\end{align}"

        return rerender_variables(equation)
    
    def interstitial_volume_bc(self):
        return rerender_variables(int_vol_BC(self.resolution, self.has_axial_dispersion))
    
    def particle_equations(self):

        equations = {}
        boundary_conditions = {}

        for par_type in self.par_type_counts.keys():

            equations[par_type] = particle_transport(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=self.nonlimiting_filmDiff, has_surfDiff=self.has_surfDiff, has_binding=self.has_binding, has_mult_bnd_states=self.has_mult_bnd_states)
            equations[par_type] = rerender_variables(equations[par_type])

            boundary_conditions[par_type] = particle_boundary(par_type, singleParticle=self.N_p == 1, nonlimiting_filmDiff=self.nonlimiting_filmDiff, has_surfDiff=self.has_surfDiff, has_binding=self.has_binding, has_mult_bnd_states=self.has_mult_bnd_states)
            boundary_conditions[par_type] = rerender_variables(boundary_conditions[par_type])

        return equations, boundary_conditions
    

#%% Streamlit UI
st.title("CADET-Equations: Packed-Bed Chromatography Model Equation Generator")
st.write("Configure a chromatography model to get the corresponding governing equations.")

# User configuration of the model

# Do not change the label of the dev mode, otherwise it will be included in the CI
dev_mode_=st.selectbox("Developer setup options (not tested! Enables e.g. multiple particle types)", ["Off", "On"], key="dev_mode") == "On"

advanced_mode_=st.selectbox("Advanced setup options (enables e.g. multiple bound states)", ["Off", "On"], key="advanced_mode") == "On"

column_model = Column(
    dev_mode=dev_mode_, advanced_mode=advanced_mode_,
    resolution=re.search(r'\dD', st.selectbox("Column resolution", ["1D (axial coordinate)", "2D (axial and radial coordinate)", "3D (axial, radial and angular coordinate)"], key="column_resolution")).group(),
    N_c = st.number_input("Number of components", key="N_c", min_value=1, step=1),
    N_p = st.number_input("Number of particle types", key="N_p", min_value=0, step=1) if dev_mode_ else int(st.selectbox("Add particles", ["No", "Yes"], key="add_particles") == "Yes")
    )


#%% Display equations

interstitial_volume_eq = column_model.interstitial_volume_equation()

nComp_list = ', '.join(str(i) for i in range(1, column_model.N_c + 1))
nPar_list = ', '.join(str(j) for j in range(1, column_model.N_p + 1))

st.write("MODEL NAME BASED ON INPUT")

file_content = []

def write_and_save(output:str, as_latex:bool=False):

    file_content.append(output)

    if as_latex:
        st.latex(output)
    else:
        st.write(output)


write_and_save(r"In the interstitial volume, mass transfer is governed by the following convection-diffusion-reaction equations in " + int_vol_domain[column_model.resolution] + r" and for all components $i\in\{" + str(nComp_list) + r"\}$")
write_and_save(interstitial_volume_eq, as_latex=True)
write_and_save("with boundary conditions")
write_and_save(column_model.interstitial_volume_bc(), as_latex=True)
if column_model.N_p >0:
     tmp = r"and particle types $j\in\{" + str(nPar_list) + r"\}$"
else:
     tmp = r""

if column_model.N_p > 0:

    particle_eq, particle_bc = column_model.particle_equations()

    write_and_save(r"In the particle models, mass transfer is governed by the following diffusion-reaction equations,")

    tmp = 0
    for par_type in column_model.par_type_counts.keys():

        tmp_textblock = "in " + particle_domain(column_model.resolution, par_type.resolution, par_type.hasCore, with_par_index=True, with_time_domain=True) + " for "

        if column_model.N_p > 1:
            
            tmp_textblock += r"particle types $j\in\{" + nPar_list + r"\}$ and components $i\in\{" + str(nComp_list) + r"\}$"
        else:
            tmp_textblock += r"components $i\in\{" + str(nComp_list) + r"\}$"

        write_and_save(tmp_textblock)
        write_and_save(particle_eq[par_type], as_latex=True)
        tmp += column_model.par_type_counts[par_type]

        if not particle_bc[par_type] == "":
            write_and_save("The boundary conditions are given by")
            write_and_save(particle_bc[par_type], as_latex=True)

write_and_save("Initial values for all solution variables (concentrations) are defined at $t = 0$.")

latex_string = [
    r"\documentclass{article}",
    r"\usepackage{amsmath}",
    r"\begin{document}",
    r"\section*{Chromatography Model}"
    ]
latex_string.extend(file_content + [r"\end{document}"])
latex_string = "\n".join(latex_string)

st.download_button("Download .tex", latex_string, "output.tex", "text/plain")

with tempfile.TemporaryDirectory() as temp_dir:
    tex_path = f"{temp_dir}/output.tex"
    pdf_path = f"{temp_dir}/output.pdf"

    # Write LaTeX content to a temp file
    with open(tex_path, "w") as f:
        f.write(latex_string)

    result = subprocess.run(
        ["pdflatex", "-output-directory", temp_dir, tex_path],
        capture_output=True,
        text=True,
        cwd=temp_dir
    )

    # Print LaTeX output for debugging
    print(result.stdout)
    print(result.stderr)

    # Check if pdflatex was successful
    if result.returncode != 0:
        print("LaTeX compilation failed.")
        print("Error output:", result.stderr)

    # Read and provide the PDF for download
    with open(pdf_path, "rb") as pdf_file:
        st.download_button("Download PDF", pdf_file, "output.pdf", "application/pdf")
