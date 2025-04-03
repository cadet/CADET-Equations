import streamlit as st
import sympy as sp


st.set_page_config(
    page_title="CADET-Equations",
    page_icon=":material/calculate:",
)



# Indices
i = sp.Idx('i') # Component index

# Independent variables
t = sp.symbols('t')
z = sp.symbols('z')
r = sp.symbols('r')
indep_vars = (t, z, r)

t, z, r = sp.symbols('t z r')

# Concentration functions
cb_i = sp.Function(r'c^{\mathrm{b}}_i')(t, z)  # Bulk concentration
cp_i = sp.Function(r'c^{\mathrm{p}}_{i}')(t, z, r)  # Particle concentration at surface
cs_i = sp.Function(r'c^{\mathrm{s}}_{i}')(t, z, r)  # Particle concentration at surface

# PDE derivatives
cb_dt = sp.diff(cb_i, t)  # ∂c^b_i/∂t
cp_i_dt = sp.diff(cp_i, t)  # ∂c^b_i/∂t
cs_i_dt = sp.diff(cs_i, t)  # ∂c^b_i/∂t
cb_dz = sp.diff(cb_i, z)  # ∂c^b_i/∂z
cb_dzz = sp.diff(cb_i, z, z)  # ∂²c^b_i/∂z²
cp_i_drr = sp.diff(cp_i, r, r)  # ∂²c^p_i/∂r²
cs_i_drr = sp.diff(cs_i, r, r)  # ∂²c^s_i/∂r²

# Parameters
eps_c = sp.symbols(r'\varepsilon^{\mathrm{c}}')  # Column porosity
eps_p = sp.symbols(r'\varepsilon^{\mathrm{p}}')  # Particle porosity
u = sp.symbols('u') # Interstitial velocity, cannot depend on z due to incompressibility
D_ax_i = sp.symbols(r'D^{\mathrm{ax}}_{i}')  # Axial dispersion coefficient
kf = sp.IndexedBase(r'k^{\mathrm{f}}')  # Film diffusion coefficient
Rp = sp.Symbol(r'R^{\mathrm{p}}')  # Particle radius
Np = sp.Symbol(r'N^{\mathrm{p}}', integer=True)  # Number of particle types

# Binding
fbind_j = sp.Function(r'f^{\mathrm{bind}}_{j}')(cp_i, cs_i)

# Film diffusion term with Np = 1
fd_term = 3 / Rp * kf[i] * (cb_i - cp_i.subs(r, Rp))

# Define the PDE equation
lrmp_bulk = sp.Eq(
    cb_dt,
    -u * cb_dz + D_ax_i * cb_dzz - (1 - eps_c) / eps_c * fd_term
)
lrmp_particle_liquid = sp.Eq(
    cp_i_dt,
    fd_term / eps_p - (1 - eps_p) / eps_p * fbind_j
)
lrmp_particle_solid = sp.Eq(
    cs_i_dt,
    fbind_j
)


# Model derivation

st.write("Starting point: LRMP with $N^p=1$")
st.latex(sp.latex(lrmp_bulk))
st.latex(sp.latex(lrmp_particle_liquid))
st.latex(sp.latex(lrmp_particle_solid)) # TODO remove component index from binding function

st.write("Next, we remove $k^{\mathrm{f}}$ from the interstitial volume mass balance by substituting the corresponding expression from the particle liquid equation")
to_be_subs = 3 / Rp * kf[i] * (cb_i - cp_i.subs(r, Rp))
substitute = sp.solve(lrmp_particle_liquid, to_be_subs)[0]
lrm_bulk = lrmp_bulk.subs(to_be_subs, substitute)
st.latex(sp.latex(lrm_bulk))

st.write("Substituting the solid phase derivate into this equation, we get")
lrm_bulk = lrm_bulk.subs(fbind_j, cs_i_dt)
st.latex(sp.latex(lrm_bulk))

st.write("We rearrangethe equation and get")
lrm_bulk = sp.Eq(lrm_bulk - rhs, lrm_bulk - )
