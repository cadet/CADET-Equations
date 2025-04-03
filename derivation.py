import streamlit as st
import sympy as sp


dim_int_vol_eq = 1
dim_par_eq = 1
# Independent variables
t = sp.symbols('t')
indep_vars = (t)

if dim_int_vol_eq > 0:

    z = sp.symbols('z')
    indep_vars += (z,)

    if dim_int_vol_eq > 1:

        rho = sp.symbols('rho')
        indep_vars += (rho,)

    if dim_int_vol_eq > 2:
        
        phi = sp.symbols('phi')
        indep_vars += (phi,)

if dim_par_eq > 0:

    r = sp.symbols('r')
    indep_vars += (r,)

t, z, rho, phi, r = sp.symbols('t z rho phi r')

# Concentration functions
cb_i = sp.Function(r'c^{\mathrm{b}}_i')(t, z, rho, phi)  # Bulk concentration
cp_ji = sp.Function(r'c^{\mathrm{p}}_{j,i}')(t, z, rho, phi, r)  # Particle concentration at surface
cs_ji = sp.Function(r'c^{\mathrm{s}}_{j,i}')(t, z, rho, phi, r)  # Particle concentration at surface

# PDE derivatives
cb_t = sp.diff(cb_i, t)  # ∂c^b_i/∂t
cp_t = sp.diff(cp_ji, t)  # ∂c^b_i/∂t
cs_t = sp.diff(cs_ji, t)  # ∂c^b_i/∂t
cb_z = sp.diff(cb_i, z)  # ∂c^b_i/∂z
cb_zz = sp.diff(cb_i, z, z)  # ∂²c^b_i/∂z²
cb_rr = sp.diff(cb_i, rho, rho)  # ∂²c^b_i/∂rho²
cb_aa = sp.diff(cb_i, phi, phi)  # ∂²c^b_i/∂phi²
cp_rr = sp.diff(cb_i, r, r)  # ∂²c^p_i/∂r²
cs_rr = sp.diff(cb_i, r, r)  # ∂²c^s_i/∂r²

# Parameters
eps_c = sp.Function(r'\varepsilon^{\mathrm{c}}')(rho, phi)  # Column porosity, cannot depend on z due to incompressibility
u = sp.Function('u')(rho, phi)  # Interstitial velocity, cannot depend on z due to incompressibility
D_ax_i = sp.Function(r'D^{\mathrm{ax}}_{i}')(z, rho, phi)  # Axial dispersion coefficient
D_rad_i = sp.Function(r'D^{\mathrm{rad}}_{i}')(z, rho, phi)  # Radial dispersion coefficient
D_ang_i = sp.Function(r'D^{\mathrm{ang}}_{i}')(z, rho, phi)  # Angular dispersion coefficient
d_j = sp.Function(r'd_{j}')(z, rho, phi) # Particle type volume fraction
kf_ji = sp.Function(r'k^{\mathrm{f}}_{j,i}')(z, rho, phi)  # Film difusion coefficient
Rp_j = sp.Symbol(r'R^{\mathrm{p}}_{j}')  # Particle radius
Np = sp.Symbol(r'N^{\mathrm{p}}', integer=True)  # Number of particle types

# Binding
fbind_j = sp.Function(r'f^{\mathrm{bind}_j}')(cp_ji, cs_ji) # todo: vectors?!

# Indices
i = sp.Symbol('i', integer=True) # Component index
j = sp.Symbol('j', integer=True) # Particle type index
k = sp.Symbol('k', integer=True) # Bound state index



# Film diffusion sum
fd_sum = sp.Sum((3 * d_j / Rp_j) * kf_ji * (cb_i - cp_ji.subs(r, Rp_j)), (j, 1, Np))

# Define the PDE equation
pde = sp.Eq(
    cb_t,
    -u * cb_z + D_ax_i * cb_zz + D_rad_i * cb_rr + D_ang_i * cb_aa - (1 - eps_c) / eps_c * fd_sum
)

# Model derivation and display of it

# Set the page layout to wide
st.set_page_config(layout="wide")
st.write("Starting point: 3D GRM")
st.latex(sp.latex(pde))

pde = pde.subs(rho, 0)

st.latex(sp.latex(pde))


