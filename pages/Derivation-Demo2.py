import streamlit as st
import sympy as sp

sp.init_printing()

# Variables
t, z, r = sp.symbols("t z r")
cb, cp, cs = sp.symbols(r"c^\mathrm{b} c^\mathrm{p} c^\mathrm{s}", cls=sp.Function)
cb = cb(t, z)
cp = cp(t, z)
cs = cs(t, z)

# Parameters
eps_c, eps_p = sp.symbols(r"\varepsilon^{c} \varepsilon^{p}")
u, D_ax, Rp, kf = sp.symbols("u D^\mathrm{ax} R^\mathrm{p} k^\mathrm{f}")
fbind = sp.Function(r"f^{\mathrm{bind}}")(cb, cs)

# Derivatives
cb_t = sp.Derivative(cb, t)
cb_z = sp.Derivative(cb, z)
cb_zz = sp.Derivative(cb, z, z)
cp_t = sp.Derivative(cp, t)
cs_t = sp.Derivative(cs, t)

# Original system
eq_bulk = sp.Eq(eps_c * cb_t, -u * sp.Derivative(eps_c * cb, z) + sp.Derivative(eps_c * D_ax * cb_z, z) - (1 - eps_c) * (3 / Rp) * kf * (cb - cp))
eq_cp = sp.Eq(cp_t, (3 / (eps_p * Rp)) * kf * (cb - cp) - (1 - eps_p) / eps_p * fbind)
eq_cs = sp.Eq(cs_t, fbind)

# Step container
steps = []

# Starting point: LRMP
steps.append(("Lumped Rate Model with Pores (LRMP)", (eq_bulk, eq_cp, eq_cs)))

# Step 1: Substitute the solid concentration into the particle liquid equation
eq_cp = eq_cp.subs(fbind, cs_t)
steps.append(("Substitute the solid concentration into the particle liquid equation", (eq_bulk, eq_cp, eq_cs)))

# Step 2: Substitute cp_t from eq_cp into eq_bulk
# Rearranged cp_t for substitution
replaced_cp = sp.solve(eq_cp, cb - cp)[0]
fd_term = (3 / Rp) * kf * (cb - cp)
fd_term_sub = sp.simplify(fd_term.subs(cb - cp, replaced_cp))

# Substituted into eq_bulk
rhs_bulk2 = eq_bulk.rhs.subs(fd_term, fd_term_sub)
eq_bulk2 = sp.Eq(eq_bulk.lhs, rhs_bulk2)

steps.append(("Remove $k^\mathrm{f}$ from bulk equation by substitution with the particle liquid equation", (eq_bulk2, eq_cp, eq_cs)))

# Step 3: Take k_f -> ∞ (i.e., cb = cp)
eq_cp2 = sp.Eq(eq_cp.lhs / kf, eq_cp.rhs / kf).simplify()
lhs_lim = sp.limit(eq_cp2.lhs, kf, sp.oo)
rhs_lim = sp.limit(eq_cp2.rhs, kf, sp.oo)
eq_cp3 = sp.Eq(lhs_lim, rhs_lim).simplify()
eq_cp4 = sp.Eq(eq_cp3.lhs * (Rp * eps_p / 3), eq_cp3.rhs * (Rp * eps_p / 3)).simplify()

steps.append((r"Divide the particle liquid eq. by $k^\mathrm{f}$, then take the limit $k^\mathrm{f} \to \infty$, then divide by $\frac{R^\mathrm{p} \varepsilon^\mathrm{p}}{3}$", (eq_bulk2, eq_cp4, eq_cs)))

eq_bulk3 = eq_bulk2.subs(cp, cb).simplify()
steps.append((r"Substitute $c^\mathrm{p}$ with $c^\mathrm{b}$ in the bulk equation", (eq_bulk3, eq_cp4, eq_cs)))

# Step 4: Move time derivatives to the left side
rhs_time_der = (1 - eps_c) * (eps_p * (cb_t - cs_t) + cs_t)
eq_bulk4 = sp.Eq(eq_bulk3.lhs + rhs_time_der, eq_bulk3.rhs + rhs_time_der).simplify()
steps.append((r"Move time derivatives to RHS", (eq_bulk4, eq_cp4, eq_cs)))


# Step 5: Introduce total porosity
eps_t = sp.Symbol(r"\varepsilon^{t}")


# Render the derivation steps
for title, (eq1, eq2, eq3) in steps:
    st.markdown(f"### {title}")
    st.latex(sp.latex(eq1))
    if eq2 is not None:
        st.latex(sp.latex(eq2))
    if eq3 is not None:
        st.latex(sp.latex(eq3))
    st.markdown("---")
