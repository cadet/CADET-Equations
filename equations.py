# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import re

#%% Model assumptions

def HRM_asmpt(N_p:int, nonlimiting_filmDiff:bool, has_binding:bool, has_surfDiff:bool, column_resolution:str):
    
    if int(re.search("\\d", column_resolution).group()) == 0:
        device_str = "tank"
    else:
        device_str = "column"

    asmpts= [
        r"the fluid is incompressible;",
        r"the fluid flow is laminar;",
        r"the fluid density and viscosity are constant;",
        r"the " + device_str + " is operated under constant conditions (e.g., temperature, flow rate);",
        r"the process is isothermal (i.e., there are no thermal effects);",
        r"diffusion does not depend on the concentration of the components, viscosity of the fluid, or pressure;"if N_p > 0 or int(re.search(r'\d?', column_resolution).group()) > 0 else "",
        r"the fluid in the particles is stagnant (i.e., there is no convective flow);" if N_p > 0 else "there is no solid phase in the " + device_str + " (and thus no adsorption);",
        r"the porous particles are spherical, homogeneous, rigid, of uniform porosity, and do not move;" if N_p > 0 else "",
        r"the partial molar volumes are the same in mobile and solid phase;" if has_binding else "",
        r"the binding model parameters do not depend on pressure and are constant along the column;"if has_binding else "",
    ]

    if N_p > 0:
        asmpts.append(r"the solvent is not adsorbed;" if has_binding else "both the solvent and the solutes are not adsorbed (i.e. there is no adsorption);")

    return [asmpt for asmpt in asmpts if not asmpt == ""]

def int_vol_3DContinuum_asmpt(N_p:int, nonlimiting_filmDiff:bool):

    volume_string = "interstitial volume" if N_p > 0 else "column"

    asmpts = [
	r"the particles form a continuum inside the column (i.e., there is interstitial and particle volume at every point in the column);" if N_p > 0 else "",
	r"there is a small number $N_{\mathrm{p}}\geq 1$ of (representative) particle radii $R_{\mathrm{p},j}$, $j \in \{ 1, \dots, N_{\mathrm{p}} \}$;"if N_p > 0 else "",
    ]

    return [asmpt for asmpt in asmpts if not asmpt == ""]

def int_vol_2DContinuum_asmpt(N_p:int, nonlimiting_filmDiff:bool):

    return int_vol_3DContinuum_asmpt(N_p, nonlimiting_filmDiff) + [
	"the column is radially symmetric (i.e., concentration profiles and parameters only depend on the axial and radial position);"
    ]

def int_vol_1DContinuum_asmpt(N_p:int, nonlimiting_filmDiff:bool):

    return int_vol_2DContinuum_asmpt(N_p, nonlimiting_filmDiff) + [
	r"the column is radially homogeneous (i.e., concentration profiles and parameters only depend on the axial position);"
    ]

def int_vol_0DContinuum_asmpt(N_p:int, nonlimiting_filmDiff:bool):
    
    volume_string = "interstitial volume" if N_p > 0 else "tank"

    return [r"the tank is spatially homogeneous (i.e., the spatial position inside the " + volume_string + " does not matter);"]

def int_vol_continuum_asmpt(resolution:str, N_p:int, nonlimiting_filmDiff:bool):
    
    asmpts = []
    if int(re.search(r'\d?', resolution).group()) > 0:
        asmpts.append(r"the fluid only flows in the axial direction of the column (i.e., there is no flow in the radial and angular direction);")

    if resolution == "3D":
        return int_vol_3DContinuum_asmpt(N_p, nonlimiting_filmDiff)
    elif resolution == "2D":
        return int_vol_2DContinuum_asmpt(N_p, nonlimiting_filmDiff)
    elif resolution == "1D":
        return int_vol_1DContinuum_asmpt(N_p, nonlimiting_filmDiff)
    elif resolution == "0D":
        return int_vol_0DContinuum_asmpt(N_p, nonlimiting_filmDiff)

def particle_1D_asmpt(has_surfDiff:bool):

     asmpts = [
          r"the interstitial liquid phase concentration is spatially constant on the particle surface;"
     ]

     if has_surfDiff:
          asmpts.append(r"there is no surface diffusion (i.e. adsorbed molecules are spatially constant)")

     return asmpts

def particle_0D_asmpt():

     return particle_1D_asmpt(False) + [
          r"the pore and surface diffusion are infinitely fast. That is, we assume $D_{i}^{\p} = D_{i}^{\s} = \infty$;"
     ]

def particle_asmpt(resolution:str, has_surfDiff:bool):

     if resolution == "0D":
          return particle_0D_asmpt()
     elif resolution == "1D":
          return particle_1D_asmpt(has_surfDiff)

#%% Equation definitions


# Interstitial volume transport terms including spatially variable porosity (epsilon)
bulk_time_derivative_eps = r"\varepsilon_{\mathrm{c}} \frac{\partial c^{\b}_i}{\partial t}"
solid_time_derivative_eps = r"\left( 1 - \varepsilon_{\mathrm{c}} \right) \frac{\partial c^{\s}_i}{\partial t}"
axial_convection_eps = r"- u \frac{\partial \left( \varepsilon_{\mathrm{c}} c^{\b}_i \right)}{\partial z}"
axial_dispersion_eps = r"\frac{\partial}{\partial z} \left( \varepsilon_{\mathrm{c}} D_{\mathrm{ax},i} \frac{\partial c^{\b}_i}{\partial z} \right)"
radial_dispersion_eps = r"\frac{1}{\rho} \frac{\partial}{\partial \rho} \left( \rho \varepsilon_{\mathrm{c}} D_{\mathrm{rad},i}  \frac{\partial c^{\b}_i}{\partial \rho} \right)"
angular_dispersion_eps = r"\frac{1}{\rho} \frac{\partial}{\partial \varphi} \left( \varepsilon_{\mathrm{c}} D_{\mathrm{ang},i}  \frac{\partial c^{\b}_i}{\partial \varphi} \right)"


# Film diffusion in the interstitial volume
def int_filmDiff_term(particle, numIdxBegin, numIdxEnd, singleParticle=False, nonLimitingFilmDiff=False):
    
    if singleParticle:
        term = r"- \left(1 - \varepsilon_{\mathrm{c}} \right) \frac{" + str(particle.surface_volume_ratio) + r"}{R_{\mathrm{p}}} k_{\mathrm{f},i} \left(c^{\b}_i - \left. c^{\p}_{i} \right|_{r = R_{\mathrm{p}}} \right)"
    else:
        term = r"- \left(1 - \varepsilon_{\mathrm{c}} \right) \sum_{j=" + str(numIdxBegin) + r"}^{" + str(numIdxEnd) + r"} \frac{" + str(particle.surface_volume_ratio) + r"d_j}{R_{\mathrm{p},j}} k_{\mathrm{f},i,j} \left(c^{\b}_i - \left. c^{\p}_{i,j} \right|_{r = R_{\mathrm{p},j}} \right)"

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
     ax_diff_term = r"- D_{\mathrm{ax},i} \frac{\partial c^{\b}_i}{\partial z}" if hasAxialDispersion else ""
     inflow_bc = r"u c_{\mathrm{in},i} &= \left.\left( u c^{\b}_i " + ax_diff_term + r"\right)\right|_{z=0} & &\qquad\text{on }" + ax_bc_domain
     outflow_bc = r"0 &= - D_{\mathrm{ax},i} \left. \frac{\partial c^{\b}_i}{\partial z} \right|_{z=L} & &\qquad\text{on }" +  ax_bc_domain

     rad_wall_bc = r"0 &= - \left(D_{\mathrm{rad},i} \left. \frac{\partial c^{\b}_i}{\partial \rho} \right) \right|_{\rho=R_{\mathrm{c}}} & &\qquad\text{on }" + rad_bc_domain
     rad_inner_bc = r"0 &= - \left(D_{\mathrm{rad},i} \left. \frac{\partial c^{\b}_i}{\partial \rho} \right) \right|_{\rho=0} & &\qquad\text{on }" + rad_bc_domain

     ang_periodic_bc = r"0 &= D_{\mathrm{ang},i} \, c^{\b}_i \Big|_{\varphi=0} - D_{\mathrm{ang},i} \, c^{\b}_i \Big|_{\varphi=2\pi} & &\quad \text{on }" + ang_bc_domain

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

          bulk_liquid_eq = r"\left. c^{\b}_i \right|_{t = 0} &= c^{\b}_{\mathrm{init},i} & & \qquad\text{in } (0, L)"
          par_liquid_eq = r"\left. c^{\p}_{i,j} \right|_{t = 0, r = R_{\mathrm{p},j}} &= c^{\p}_{\mathrm{init},i,j} & & \qquad\text{in } (0, L)"

     if resolution == "2D":
          bulk_liquid_eq = r"\left. c^{\b}_i \right|_{t = 0} &= c^{\b}_{\mathrm{init},i} & & \qquad\text{in } (0, R_{\mathrm{c}}) \times (0, L)"
          par_liquid_eq = r"\left. c^{\p}_{i,j} \right|_{t = 0, r = R_{\mathrm{p},j}} &= c^{\p}_{\mathrm{init},i,j} & & \qquad\text{in } (0, R_{\mathrm{c}}) \times (0, L)"

     if resolution == "3D":
          bulk_liquid_eq = r"\left. c^{\b}_i \right|_{t = 0} &= c^{\b}_{\mathrm{init},i} & & \qquad\text{in } (0, R_{\mathrm{c}}) \times (0, L)\times [0,2\pi)"
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
     "0D": r"$(0, T_\mathrm{end})$",
     "1D": r"$(0, T_\mathrm{end}) \times (0, L)$",
     "2D": r"$(0, T_\mathrm{end}) \times (0, L) \times (0, R_\mathrm{c})$",
     "3D": r"$(0, T_\mathrm{end}) \times (0, L) \times (0, R_\mathrm{c}) \times (0, 2\pi)$"
}
int_vol_inlet_domain = {
     "0D": r"(0, T_{\mathrm{end}})",
     "1D": r"(0, T_{\mathrm{end}})",
     "2D": r"(0, T_{\mathrm{end}}) \times (0, R_\mathrm{c})",
     "3D": r"(0, T_{\mathrm{end}}) \times (0, R_\mathrm{c}) \times (0, 2\pi)"
}
int_vol_vars = {
     "1D": r"z",
     "2D": r"z, \rho",
     "3D": r"z, \rho, \varphi"
}

# Particle transport terms
def particle_transport(particle, singleParticle:bool, nonlimiting_filmDiff:bool, has_surfDiff:bool, has_binding:bool, req_binding:bool, has_mult_bnd_states:bool):
     
     ret_term = ""          

     if particle.resolution == "0D":
          if nonlimiting_filmDiff and has_binding:
               ret_term = r"""\begin{align}""" + re.sub(r"{\\p}", r"{\\l}", particle_transport_homogeneous_solid(req_binding, has_mult_bnd_states)) + r""". \end{align}"""
          elif not has_binding:
               ret_term = r"""\begin{align}""" + re.sub(r"{\\p}", r"{\\l}", particle_transport_homogeneous_liquid(req_binding, has_mult_bnd_states)) + r""". \end{align}"""
          else:
               ret_term = particle_transport_homogeneous(req_binding, has_mult_bnd_states)
     else:
          ret_term = particle_transport_radial(particle.geometry, has_surfDiff, has_binding, req_binding, has_mult_bnd_states)

     if singleParticle:
          ret_term = re.sub(",j", "", ret_term)

     return ret_term

def particle_transport_homogeneous_liquid(req_binding:bool, has_mult_bnd_states:bool):

     bnd_term = r"f_{\mathrm{bind},i,j} \left( \vec{c}^{\p}, \vec{c}^{\s} \right) "

     lhs_term = r"\varepsilon_{\mathrm{p},j} \frac{\partial c^{\p}_{i,j}}{\partial t} "
     rhs_term = r"\frac{3}{R_{\mathrm{p},j}} k_{\mathrm{f},i,j} \left( c^{\b}_{i} - c^{\p}_{i,j} \right) "

     if has_mult_bnd_states: # add sum over bound states
          if req_binding:
               lhs_term += r"\sum_{k=1}^{N_{\mathrm{b},i}} "
          else:
               bnd_term = r"\sum_{k=1}^{N_{\mathrm{b},i}} " + bnd_term

     if req_binding:
          lhs_term += r" + \left( 1 - \varepsilon_{\mathrm{p},j}\right) \frac{\partial c^{\s}_{i,j}}{\partial t} "
     else:
          rhs_term += r"- \left( 1 - \varepsilon_{\mathrm{p},j} \right) " + bnd_term

     if has_mult_bnd_states: # add sum and bound state indices
          bnd_term = r"\sum_{k=1}^{N_{\mathrm{b},i}}"
          lhs_term = re.sub("i,j", "i,j,k", lhs_term)
          rhs_term = re.sub("i,j", "i,j,k", rhs_term)

     return r"""
          """ + lhs_term + r"""
          &=""" + rhs_term

def particle_transport_homogeneous_solid(req_binding:bool, has_mult_bnd_states:bool):
     
     bnd_term = r"f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right)"

     if req_binding:
          lhs_term = r"0"
     else:
          lhs_term = r"\frac{\partial c^{\s}_{i,j}}{\partial t}"

     if has_mult_bnd_states: # add bound state index
          lhs_term = re.sub("i,j", "i,j,k", lhs_term)
          bnd_term = re.sub("i,j", "i,j,k", bnd_term)

     return r"""
          """ + lhs_term + r"""
          &= """ + bnd_term

def particle_transport_homogeneous(req_binding:bool, has_mult_bnd_states:bool):
     return r"""
\begin{align}
""" + particle_transport_homogeneous_liquid(req_binding, has_mult_bnd_states) + r""",\\
""" + particle_transport_homogeneous_solid(req_binding, has_mult_bnd_states) + r""".
\end{align}
"""

def particle_transport_radial(geometry:str, has_surfDiff:bool, has_binding:bool, req_binding:bool, has_mult_bnd_states:bool):
     
     surfDiffTerm = ""
     binding_term = r"f_{\mathrm{bind},i,j}\left( \vec{c}^{\p}, \vec{c}^{\s} \right) "
     if has_mult_bnd_states:
          binding_term = re.sub("i,j", "i,j,k", binding_term)

     if geometry == "Sphere": 

          if has_surfDiff:
               surfDiffTerm = r" \frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i,j}^{\s} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) "
          
          liquid_lhs = r"\frac{\partial c^{\p}_{i,j}}{\partial t} "
          liquid_rhs = r"\frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i,j}^{\p} \frac{\partial c^{\p}_{i,j}}{\partial r} \right)"

          if req_binding:
               solid_lhs = r"0 "
               solid_rhs = binding_term
               liquid_lhs += r"+ "

               if has_mult_bnd_states:
                    liquid_lhs += r"\sum_{k=1}^{N_{\mathrm{b},i}} "
               
               liquid_lhs += r"\frac{1 - \varepsilon_{\mathrm{p},j}}{\varepsilon_{\mathrm{p},j}} \frac{\partial c^{\s}_{i,j}}{\partial t}"
               
               if has_surfDiff:
                    liquid_rhs += r" - \frac{1 - \varepsilon_{\mathrm{p},j}}{\varepsilon_{\mathrm{p},j}}"

                    if has_mult_bnd_states:
                         liquid_rhs += r"\sum_{k=1}^{N_{\mathrm{b},i}} "

                    liquid_rhs += surfDiffTerm

          else:
               solid_lhs = r"\frac{\partial c^{\s}_{i,j}}{\partial t} "
               solid_rhs = surfDiffTerm + " + " if has_surfDiff else ""
               solid_rhs += binding_term
               liquid_rhs += r" - \frac{1 - \varepsilon_{\mathrm{p},j}}{\varepsilon_{\mathrm{p},j}}" + binding_term

          solid_eq = solid_lhs + r"&= " + solid_rhs
          liquid_eq = liquid_lhs + r"&= " + liquid_rhs

          if has_mult_bnd_states:
               solid_eq = re.sub("i,j", "i,j,k", solid_eq)
               liquid_eq = re.sub("i,j", "i,j,k", liquid_eq)

          if has_binding:
               return r"""
\begin{align}
""" + liquid_eq + r""", \\
""" + solid_eq + r""".
\end{align}
"""
          else:
               return r"""
\begin{align}
""" + liquid_eq + r""".
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

def particle_boundary(particle, singleParticle:bool, nonlimiting_filmDiff:bool, has_surfDiff:bool, has_binding:bool, req_binding:bool, has_mult_bnd_states:bool):

     if particle.resolution == "0D":
          return ""
     
     if nonlimiting_filmDiff:
          outerLiquidBC = r"\left. c^{\p}_{i,j} \right|_{r = R_{\mathrm{p},j}} &= c^{\b}_i"
     else:
          if not req_binding:
               outerLiquidBC = r""" \varepsilon_p \left. \left( D^{\p}_{i,j} \frac{\partial c^{\p}_{i,j}}{\partial r} \right)\right|_{r = R_{\mathrm{p},j}}
               &= k_{\mathrm{f},i,j} \left. \left( c^{\b}_i - c^{\p}_{i,j} \right|_{r = R_{\mathrm{p},j}} \right)"""
          else:
               outerLiquidBC = r""" \left. \left( \varepsilon_p  D^{\p}_{i,j} \frac{\partial c^{\p}_{i,j}}{\partial r} + (1 - \varepsilon_p ) D^{\s}_{i,j} \frac{\partial c^{\s}_{i,j}}{\partial r} \right)\right|_{r = R_{\mathrm{p},j}}
               &= k_{\mathrm{f},i,j} \left. \left( c^{\b}_i - c^{\p}_{i,j} \right|_{r = R_{\mathrm{p},j}} \right)"""

     inner_boundary = r"R_{\mathrm{pc},j}" if particle.hasCore else r"0"

     if req_binding and has_surfDiff:
          innerLiquidBC = r"- \left. \left( \varepsilon_p D^{\p}_{i,j} \frac{\partial c^{\p}_{i,j}}{\partial r} + (1 - \varepsilon_p) D^{\s}_{i,j} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) \right|_{r=" + inner_boundary + r"}"
     else:
          innerLiquidBC = r"- \left. \left( D^{\p}_{i,j} \frac{\partial c^{\p}_{i,j}}{\partial r} \right) \right|_{r=" + inner_boundary + r"}"
     particleLiquidBC =  r"""
""" + innerLiquidBC + r"""
&= 0, \\
""" + outerLiquidBC
     
     if has_surfDiff and has_binding and not req_binding:

          particleSolidBC = r""",\\
-\left( \left. D^{\s}_{i,j} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) \right|_{r=""" + inner_boundary + r"""}
&= 0, \\
\left( \left. D^{\s}_{i,j} \frac{\partial c^{\s}_{i,j}}{\partial r} \right) \right|_{r = R_{\mathrm{p},j}}
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
    
    if not column_resolution == "0D":
          domain = r"$ (0, T_\mathrm{end}) \times (0, L)" if with_time_domain else r"$ \times (0, L)"
    else:
          domain = r"$ (0, T_\mathrm{end})" if with_time_domain else r"$ "

    if column_resolution in ["2D", "3D"]:
        domain += r"\times (0, R_\mathrm{c})"
    if column_resolution == "3D":
        domain += r"\times (0, 2\pi)"
         
    par1D_domain = r"\times (R_{\mathrm{pc},j}, R_{\mathrm{p},j})$" if hasCore else r"\times (0, R_{\mathrm{p},j})$"
    if not with_par_index:
         par1D_domain = re.sub(r",j", "", par1D_domain)

    return domain + r"$" if particle_resolution == "0D" else domain + par1D_domain
