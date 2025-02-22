# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import re

def rerender_variables(input_str:str):
     input_str = re.sub(r"\\l(?![a-zA-Z])", r"\\mathrm{\\ell}", input_str)
     input_str = re.sub(r"\\p(?![a-zA-Z])", r"\\mathrm{p}", input_str)
     input_str = re.sub(r"\\s(?![a-zA-Z])", r"\\mathrm{s}", input_str)
     return input_str

#%% Equation definitions


# Interstitial volume transport terms including spatially variable porosity (epsilon)
bulk_time_derivative_eps = r"\varepsilon_{\mathrm{c}} \frac{\partial c^{\l}_i}{\partial t}"
axial_convection_eps = r"- u \frac{\partial \left( \varepsilon_{\mathrm{c}} c^{\l}_i \right)}{\partial z}"
axial_dispersion_eps = r"\frac{\partial}{\partial z} \left( \varepsilon_{\mathrm{c}} D_{\mathrm{ax},i} \frac{\partial c^{\l}_i}{\partial z} \right)"
radial_dispersion_eps = r"\frac{1}{\rho} \frac{\partial}{\partial \rho} \left( \rho \varepsilon_{\mathrm{c}} D_{\mathrm{rad},i}  \frac{\partial c^{\l}_i}{\partial \rho} \right)"
angular_dispersion_eps = r"\frac{1}{\rho} \frac{\partial}{\partial \varphi} \left( \varepsilon_{\mathrm{c}} D_{\mathrm{ang},i}  \frac{\partial c^{\l}_i}{\partial \varphi} \right)"


# Film diffusion in the interstitial volume
def int_filmDiff_term(particle, numIdxBegin, numIdxEnd, singleParticle=False, nonLimitingFilmDiff=False):
    
    if singleParticle:
        term = r"- \left(1 - \varepsilon_{\mathrm{c}} \right) \frac{" + str(particle.surface_volume_ratio) + r"}{R_{\mathrm{p}}} k_{\mathrm{f},i} \left(c^{\l}_i - \left. c^{\p}_{i} \right|_{r = R_{\mathrm{p}}} \right)"
    else:
        term = r"- \left(1 - \varepsilon_{\mathrm{c}} \right) \sum_{j=" + str(numIdxBegin) + r"}^{" + str(numIdxEnd) + r"} \frac{" + str(particle.surface_volume_ratio) + r"d_j}{R_{\mathrm{p},j}} k_{\mathrm{f},i,j} \left(c^{\l}_i - \left. c^{\p}_{i,j} \right|_{r = R_{\mathrm{p},j}} \right)"

    if nonLimitingFilmDiff:
        
        if particle.resolution == "0D":
             return ""
        else:
             # substitute boundary condition into equation
             substitute = r"\\varepsilon_{\\mathrm{p},j} \\left(D^{\\p}_{i,j} \\left. \\frac{\\partial c^{\\p}_{i,j}}{\\partial r} \\right)\\right|_{r = R_{\\mathrm{p},j}}"
             term = re.sub("k_.*?right.$", substitute, term)
             if singleParticle:
                  term = re.sub(",j", "", term)

    return term

# Boundary conditions of the interstitial volume equations

def int_vol_BC(resolution:str, hasAxialDispersion:bool):

     # handle domain the BC is defined on
     ax_bc_domain = r"(0, T_{\mathrm{end}})"
     rad_bc_domain = ""
     ang_bc_domain = ""

     if resolution == "2D":
          ax_bc_domain += r"\times (0, R_{\mathrm{c}})"
          rad_bc_domain = r"(0, T_{\mathrm{end}}) \times (0, L)"

     elif resolution == "3D":
          ax_bc_domain += r"\times (0, R_{\mathrm{c}}) \times [0,2\pi)"
          rad_bc_domain = r"(0, T_{\mathrm{end}}) \times (0, L) \times [0,2\pi)"
          ang_bc_domain = r"(0, T_{\mathrm{end}}) \times (0, L) \times (0, R_{\mathrm{c}})"

     # define single equations
     inflow_bc = r"u c_{\mathrm{in},i} &= \left.\left( u c^{\l}_i - D_{\mathrm{ax},i} \frac{\partial c^{\l}_i}{\partial z} \right)\right|_{z=0} & &\qquad\text{on }" + ax_bc_domain
     outflow_bc = r"0 &= - D_{\mathrm{ax},i} \left. \frac{\partial c^{\l}_i}{\partial z} \right|_{z=L} & &\qquad\text{on }" +  ax_bc_domain

     rad_wall_bc = r"0 &= - \left(D_{\mathrm{rad},i} \left. \frac{\partial c^{\l}_i}{\partial \rho} \right) \right|_{\rho=R_{\mathrm{c}}} & &\qquad\text{on }" + rad_bc_domain
     rad_inner_bc = r"0 &= - \left(D_{\mathrm{rad},i} \left. \frac{\partial c^{\l}_i}{\partial \rho} \right) \right|_{\rho=0} & &\qquad\text{on }" + rad_bc_domain

     ang_periodic_bc = r"0 &= D_{\mathrm{ang},i} \, c^{\l}_i \Big|_{\varphi=0} - D_{\mathrm{ang},i} \, c^{\l}_i \Big|_{\varphi=2\pi} & &\quad \text{on }" + ang_bc_domain

     # collect the required equations for full BC
     boundary_conditions = inflow_bc
     
     if hasAxialDispersion:
          boundary_conditions += r""",\\
               """ + outflow_bc
     
     if resolution in ["2D", "3D"]:
          boundary_conditions += r""",\\
""" + rad_wall_bc + r""",\\
""" + rad_inner_bc
     
     if resolution == "3D":
          boundary_conditions += r""",\\
""" + ang_periodic_bc
          
     return r"""
\begin{alignat}{2}
""" + boundary_conditions + r""".
\end{alignat}"""

def int_vol_initial(resolution:str, includeParLiquid:bool):
     
     if resolution == "1D":

          bulk_liquid_eq = r"\left. c^{\l}_i \right|_{t = 0} &= c^{\l}_{\mathrm{init},i} & & \qquad\text{in } (0, L)"
          par_liquid_eq = r"\left. c^{\p}_{i,j} \right|_{t = 0, r = R_{\mathrm{p},j}} &= c^{\p}_{\mathrm{init},i,j} & & \qquad\text{in } (0, L)"

     if resolution == "2D":
          bulk_liquid_eq = r"\left. c^{\l}_i \right|_{t = 0} &= c^{\l}_{\mathrm{init},i} & & \qquad\text{in } (0, R_{\mathrm{c}}) \times (0, L)"
          par_liquid_eq = r"\left. c^{\p}_{i,j} \right|_{t = 0, r = R_{\mathrm{p},j}} &= c^{\p}_{\mathrm{init},i,j} & & \qquad\text{in } (0, R_{\mathrm{c}}) \times (0, L)"

     if resolution == "3D":
          bulk_liquid_eq = r"\left. c^{\l}_i \right|_{t = 0} &= c^{\l}_{\mathrm{init},i} & & \qquad\text{in } (0, R_{\mathrm{c}}) \times (0, L)\times [0,2\pi)"
          par_liquid_eq = r"\left. c^{\p}_{i,j} \right|_{t = 0, r = (0, R_{\mathrm{p},j})} &= c^{\p}_{\mathrm{init},i,j} & & \qquad\text{in } (0, R_{\mathrm{c}}) \times (0, L)\times [0,2\pi)"

     if includeParLiquid:
          equation = r"""
          \begin{align}
          """ + bulk_liquid_eq + r""",
          \\""" + par_liquid_eq + r""".
          \end{align}
          """
     else:
          equation = r"""
          \begin{align}
          """ + bulk_liquid_eq + r""".
          \end{align}
          """

     return equation


int_vol_domain = {
     "1D": r"$(0, T_\mathrm{end}) \times (0, L)$",
     "2D": r"$(0, T_\mathrm{end}) \times (0, L) \times (0, R_\mathrm{c})$",
     "3D": r"$(0, T_\mathrm{end}) \times (0, L) \times (0, R_\mathrm{c}) \times (0, 2\pi)$"
}

# Particle transport terms
def particle_transport(particle, singleParticle:bool, nonlimiting_filmDiff:bool, has_surfDiff:bool, has_binding:bool, has_mult_bnd_states:bool):
     
     ret_term = ""          

     if particle.resolution == "0D":
          if nonlimiting_filmDiff and has_binding:
               ret_term = r"""\begin{align}""" + re.sub(r"{\\p}", r"{\\l}", particle_transport_homogeneous_solid(has_mult_bnd_states)) + r""". \end{align}"""
          elif not has_binding:
               ret_term = r"""\begin{align}""" + re.sub(r"{\\p}", r"{\\l}", particle_transport_homogeneous_liquid(has_mult_bnd_states)) + r""". \end{align}"""
          else:
               ret_term = particle_transport_homogeneous(has_mult_bnd_states)
     else:
          ret_term = particle_transport_radial(particle.geometry, has_surfDiff, has_binding, has_mult_bnd_states)

     if singleParticle:
          ret_term = re.sub(",j", "", ret_term)

     return ret_term

def particle_transport_homogeneous_liquid(has_mult_bnd_states:bool):

     bnd_term = r"f_{\mathrm{bind},i,j,k} \left( \vec{c}^{\p}, \vec{c}^{\s} \right)"
     if has_mult_bnd_states: # add sum
          bnd_term = r"\sum_{k=1}^{N_{\mathrm{b},i}}"

     return r"""
		\varepsilon_{\mathrm{p},j} \frac{\partial c^{\p}_{i,j}}{\partial t}
        &= \frac{3}{R_{\mathrm{p},j}} k_{\mathrm{f},i,j} \left( c^{\l}_{i} - c^{\p}_{i,j} \right) -  \left( 1 - \varepsilon_{\mathrm{p},j} \right)""" + bnd_term

def particle_transport_homogeneous_solid(has_mult_bnd_states:bool):
     
     bnd_term = r"f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right)"
     if has_mult_bnd_states: # add bound state index
          bnd_term = re.sub("i,j", "i,j,k", bnd_term)

     return r"""
          \left( 1- \varepsilon_{\mathrm{p},j} \right) \frac{\partial c^{\s}_{i,j}}{\partial t}
          &=  \left( 1 - \varepsilon_{\mathrm{p},j} \right) """ + bnd_term

def particle_transport_homogeneous(has_mult_bnd_states:bool):
     return r"""
	\begin{align}
     """ + particle_transport_homogeneous_liquid(has_mult_bnd_states) + r""",
          \\""" + particle_transport_homogeneous_solid(has_mult_bnd_states) + r""".
	\end{align}
     """

def particle_transport_radial(geometry:str, has_surfDiff:bool, has_binding:bool, has_mult_bnd_states:bool):
     
     surfDiffTerm = ""

     if geometry == "Sphere": 

          if has_surfDiff:
               surfDiffTerm = r"\left( 1 - \varepsilon_{\mathrm{p},j} \right) \frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i,j}^{\s} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) +"
          
          particle_solid = r"\left( 1 - \varepsilon_{\mathrm{p},j} \right) \frac{\partial c^{\s}_{i,j}}{\partial t} &= " + surfDiffTerm + r" \left( 1 - \varepsilon_{\mathrm{p},j} \right) f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right)"

          binding_term = r"f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right)"

          if has_mult_bnd_states:
               particle_solid = re.sub("i,j", "i,j,k", particle_solid)
               binding_term = re.sub("i,j", "i,j,k", binding_term)
               binding_term = r"\sum_{k=1}^{N_{\mathrm{b},i}} " + binding_term

          particle_liquid = r"\varepsilon_{\mathrm{p},j} \frac{\partial c^{\p}_{i,j}}{\partial t} &= \varepsilon_{\mathrm{p},j} \frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i,j}^{\p} \frac{\partial c^{\p}_{i,j}}{\partial r} \right) - \left( 1 - \varepsilon_{\mathrm{p},j} \right)" + binding_term

          if has_binding:
               return r"""
               \begin{align}
               """ + particle_liquid + r""", \\
               """ + particle_solid + r""".
               \end{align}
               """
          else:
               return r"""
               \begin{align}
               """ + particle_liquid + r""".
               \end{align}
               """
          
     if geometry == "Cylinder": 

          if has_surfDiff:
               surfDiffTerm = r"\left( 1 - \varepsilon_{\mathrm{p},j} \right) \frac{1}{r} \frac{\partial }{\partial r} \left( r D_{i,j}^{\s} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) + "
     
          return r"""
               \begin{align}    
               \varepsilon_{\mathrm{p},j} \frac{\partial c^{\p}_{i,j}}{\partial t}
               &=
               \varepsilon_{\mathrm{p},j} \frac{1}{r} \frac{\partial }{\partial r} \left( r D_{i,j}^{\p} \frac{\partial c^{\p}_{i,j}}{\partial r} \right) - \left( 1 - \varepsilon_{\mathrm{p},j} \right) f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right), \\
                    \left( 1 - \varepsilon_{\mathrm{p},j} \right) \frac{\partial c^{\s}_{i,j}}{\partial t}
               &=
               """ + surfDiffTerm + r"""\left( 1 - \varepsilon_{\mathrm{p},j} \right) f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right).
               \end{align}
               """
     if geometry == "Slab":

          if has_surfDiff:
               surfDiffTerm = r"\left( 1 - \varepsilon_{\mathrm{p},j} \right) \frac{\partial }{\partial r} \left( D_{i,j}^{\s} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) + "

          return r"""
               \begin{align}    
               \varepsilon_{\mathrm{p},j} \frac{\partial c^{\p}_{i,j}}{\partial t}
               &=
               \varepsilon_{\mathrm{p},j} \frac{\partial }{\partial r} \left( D_{i,j}^{\p} \frac{\partial c^{\p}_{i,j}}{\partial r} \right) - \left( 1 - \varepsilon_{\mathrm{p},j} \right) f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right), \\
                    \left( 1 - \varepsilon_{\mathrm{p},j} \right) \frac{\partial c^{\s}_{i,j}}{\partial t}
               &=
               """ + surfDiffTerm + r"""\left( 1 - \varepsilon_{\mathrm{p},j} \right) f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right).
               \end{align}
               """

def particle_boundary(particle, singleParticle:bool, nonlimiting_filmDiff:bool, has_surfDiff:bool, has_binding:bool, has_mult_bnd_states:bool):

     if particle.resolution == "0D":
          return ""
     
     if nonlimiting_filmDiff:
          outerLiquidBC = r"\left. c^{\p}_{i,j} \right|_{r = R_{\mathrm{p},j}} &= c^{\l}_i"
     else:
          outerLiquidBC = r"""\varepsilon_{\mathrm{p},j} \left(D^{\p}_{i,j} \left. \frac{\partial c^{\p}_{i,j}}{\partial r} \right)\right|_{r = R_{\mathrm{p},j}}
          &= k_{\mathrm{f},i,j} \left( c^{\l}_i - \left. c^{\p}_{i,j} \right|_{r = R_{\mathrm{p},j}} \right)"""

     inner_boundary = r"R_{\mathrm{pc},j}" if particle.hasCore else r"0"

     particleLiquidBC =  r"""
          -\varepsilon_{\mathrm{p},j} \left( \left. D^{\p}_{i,j} \frac{\partial c^{\p}_{i,j}}{\partial r} \right) \right|_{r=""" + inner_boundary + r"""}
          &= 0, \\
          """ + outerLiquidBC
     
     if has_surfDiff and has_binding:

          particleSolidBC = r""",\\
          -\left( 1 - \varepsilon_{\mathrm{p},j} \right) \left( \left. D^{\s}_{i,j} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) \right|_{r=""" + inner_boundary + r"""}
          &= 0, \\
          \left(1 - \varepsilon_{\mathrm{p},j} \right) \left( \left. D^{\s}_{i,j} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) \right|_{r = R_{\mathrm{p},j}}
          &= 0.
          """
          if has_mult_bnd_states:
               particleSolidBC = re.sub("i,j", "i,j,k", particleSolidBC)

          boundary_condition = particleLiquidBC + particleSolidBC
     else:
          boundary_condition = particleLiquidBC + "."
     
     boundary_condition = r"\begin{align}" + boundary_condition + r"\end{align}"
     
     if singleParticle:
          boundary_condition = re.sub(",j", "", boundary_condition)

     return boundary_condition

def particle_initial(domain:str, singleParticle:bool, includeParLiquid:bool):
     
     if includeParLiquid:
          initial_condition = r"""
          \begin{alignat}{2}
          \left. c^{\p}_{i,j} \right|_{t = 0} &= c^{\p}_{\mathrm{init},i,j} & & \qquad\text{in }""" + re.sub("\\$", "", domain) + r""",\\
          \left. c^{\s}_{i,j} \right|_{t = 0} &= c^{\s}_{\mathrm{init},i,j} & & \qquad\text{in }""" + re.sub("\\$", "", domain) + r""".
          \end{alignat}
          """
     else:
          initial_condition = r"""
          \begin{alignat}{2}
          \left. c^{\s}_{i,j} \right|_{t = 0} &= c^{\s}_{\mathrm{init},i,j} & & \qquad\text{in }""" + re.sub("\\$", "", domain) + r""".
          \end{alignat}
          """

     if singleParticle:
          initial_condition = re.sub(",j", "", initial_condition)

     return initial_condition

def particle_domain(column_resolution:str, particle_resolution:str, hasCore:bool, with_par_index=False, with_time_domain=True):
    
    domain = r"$ (0, T_\mathrm{end}) \times (0, L)" if with_time_domain else r"$ \times (0, L)"

    if column_resolution in ["2D", "3D"]:
        domain += r"\times (0, R_\mathrm{c})"
    if column_resolution == "3D":
        domain += r"\times (0, 2\pi)"
         
    par1D_domain = r"\times (R_{\mathrm{pc},j}, R_{\mathrm{p},j})$" if hasCore else r"\times (0, R_{\mathrm{p},j})$"
    if not with_par_index:
         par1D_domain = re.sub(r",j", "", par1D_domain)

    return domain + r"$" if particle_resolution == "0D" else domain + par1D_domain
