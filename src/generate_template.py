# -*- coding: utf-8 -*-
"""
Generate a downloadable Python template script from a Column or Crystallization model.

The generated script defines a ``get_unit_operation()`` function that returns
a dictionary with the CADET unit operation configuration corresponding to the
currently specified equations, excluding numerics/discretization.
"""

_DEFAULTS = {
    'col_length': 0.014,
    'col_porosity': 0.37,
    'velocity': 5.75e-4,
    'col_dispersion': 5.75e-8,
    'col_radius_inner': 0.01,
    'col_radius_outer': 0.024,
    'velocity_coeff': 5.75e-4,
    'init_volume': 1e-3,
    'par_radius': 4.5e-5,
    'par_coreradius': 1e-5,
    'par_porosity': 0.75,
    'film_diffusion': 6.9e-6,
    'pore_diffusion': 6.07e-11,
    'surface_diffusion': 1e-11,
}

_BINDING_DEFAULTS = {
    'Linear': {
        'cadet_name': 'LINEAR',
        'per_comp': {'lin_ka': 3.55, 'lin_kd': 0.1},
        'scalar': {},
    },
    'Langmuir': {
        'cadet_name': 'MULTI_COMPONENT_LANGMUIR',
        'per_comp': {'mcl_ka': 1.14, 'mcl_kd': 0.002, 'mcl_qmax': 4.88},
        'scalar': {},
    },
    'SMA': {
        'cadet_name': 'STERIC_MASS_ACTION',
        'per_comp': {'sma_ka': 35.5, 'sma_kd': 1000.0, 'sma_nu': 4.7, 'sma_sigma': 11.83},
        'salt': {'sma_ka': 0.0, 'sma_kd': 0.0, 'sma_nu': 0.0, 'sma_sigma': 0.0},
        'scalar': {'sma_lambda': 1200.0},
    },
}

_PAR_GEOM_MAP = {
    'Sphere': 'SPHERE',
    'Cylinder': 'CYLINDER',
    'Slab': 'SLAB',
}


def _get_unit_type(column_model):
    if column_model.resolution == "0D":
        return "CSTR"
    prefix = ""
    if column_model.column_type == "Radial":
        prefix = "RADIAL_"
    elif column_model.column_type == "Frustum":
        prefix = "FRUSTUM_"
    return f"{prefix}COLUMN_MODEL_{column_model.resolution}"


def _per_comp(ncomp_fixed, default_value):
    if ncomp_fixed is not None:
        return repr([default_value] * ncomp_fixed)
    return f"[{default_value}] * ncomp"


def _ncomp_var(ncomp_fixed):
    if ncomp_fixed is not None:
        return str(ncomp_fixed)
    return "ncomp"


def _indent(lines, level=1):
    prefix = "    " * level
    return [prefix + line for line in lines]


def _binding_params_lines(binding_model, ncomp_fixed, req_binding):
    lines = []
    is_kinetic = 0 if req_binding else 1

    if binding_model not in _BINDING_DEFAULTS:
        lines.append(f"par['adsorption_model'] = 'NONE'  # TODO: specify binding model")
        lines.append(f"par['adsorption'] = {{")
        lines.append(f"    'is_kinetic': {is_kinetic},")
        lines.append(f"}}")
        return lines

    info = _BINDING_DEFAULTS[binding_model]
    lines.append(f"par['adsorption_model'] = '{info['cadet_name']}'")
    lines.append(f"par['adsorption'] = {{")
    lines.append(f"    'is_kinetic': {is_kinetic},")

    if binding_model == "SMA":
        salt = info['salt']
        comp = info['per_comp']
        for key in info['per_comp']:
            if ncomp_fixed is not None:
                vals = [salt[key]] + [comp[key]] * (ncomp_fixed - 1)
                lines.append(f"    '{key}': {vals},")
            else:
                lines.append(f"    '{key}': [{salt[key]}] + [{comp[key]}] * (ncomp - 1),")
    else:
        for key, val in info['per_comp'].items():
            lines.append(f"    '{key}': {_per_comp(ncomp_fixed, val)},")

    for key, val in info['scalar'].items():
        lines.append(f"    '{key}': {val},")

    lines.append(f"}}")
    return lines


def _particle_lines(particle, par_idx, ncomp_fixed, binding_model, req_binding,
                     has_binding, has_mult_bnd_states):
    lines = []
    var = f"par_{par_idx}" if par_idx > 0 else "par"
    idx_str = str(par_idx).zfill(3)

    lines.append(f"")
    lines.append(f"{var} = {{}}")
    lines.append(f"{var}['par_geom'] = '{_PAR_GEOM_MAP.get(particle.geometry, 'SPHERE')}'")
    lines.append(f"{var}['par_radius'] = {_DEFAULTS['par_radius']}")

    if particle.has_core:
        lines.append(f"{var}['par_coreradius'] = {_DEFAULTS['par_coreradius']}")

    lines.append(f"{var}['par_porosity'] = {_DEFAULTS['par_porosity']}")

    if not particle.nonlimiting_filmDiff:
        lines.append(f"{var}['has_film_diffusion'] = 1")
        lines.append(f"{var}['film_diffusion'] = {_per_comp(ncomp_fixed, _DEFAULTS['film_diffusion'])}")
    else:
        lines.append(f"{var}['has_film_diffusion'] = 0")

    if particle.resolution == "1D":
        lines.append(f"{var}['has_pore_diffusion'] = 1")
        lines.append(f"{var}['pore_diffusion'] = {_per_comp(ncomp_fixed, _DEFAULTS['pore_diffusion'])}")
        if particle.has_surfDiff:
            lines.append(f"{var}['has_surface_diffusion'] = 1")
            lines.append(f"{var}['surface_diffusion'] = {_per_comp(ncomp_fixed, _DEFAULTS['surface_diffusion'])}")
        else:
            lines.append(f"{var}['has_surface_diffusion'] = 0")
    else:
        lines.append(f"{var}['has_pore_diffusion'] = 0")
        lines.append(f"{var}['has_surface_diffusion'] = 0")

    nbound_val = 2 if has_mult_bnd_states else 1
    lines.append(f"{var}['nbound'] = {_per_comp(ncomp_fixed, nbound_val)}")
    lines.append(f"{var}['init_cp'] = {_per_comp(ncomp_fixed, 0.0)}")
    lines.append(f"{var}['init_cs'] = {_per_comp(ncomp_fixed, 0.0)}")

    if has_binding:
        lines.append(f"")
        lines.extend(_binding_params_lines(binding_model, ncomp_fixed, req_binding))

    lines.append(f"")
    lines.append(f"unit['particle_type_{idx_str}'] = {var}")
    return lines


def generate_unit_operation_script(column_model):
    """Return a Python script string for the current Column configuration."""

    model_name = column_model.model_name()
    ncomp_fixed = column_model.N_c if column_model.N_c > 0 else None
    default_ncomp = 4 if (column_model.binding_model == "SMA" and column_model.has_binding) else 1

    lines = []

    # Module docstring
    lines.append(f'"""')
    lines.append(f"CADET unit operation template: {model_name}")
    lines.append(f"Generated by CADET-Equations")
    lines.append(f'"""')
    lines.append("")
    lines.append("")

    # Function signature
    if ncomp_fixed is not None:
        lines.append("def get_unit_operation():")
    else:
        lines.append(f"def get_unit_operation(ncomp={default_ncomp}):")

    # Docstring
    body = []
    body.append(f'"""Return a CADET unit operation dictionary for a {model_name}.')
    if ncomp_fixed is None:
        body.append("")
        body.append("Parameters")
        body.append("----------")
        body.append("ncomp : int")
        body.append(f"    Number of components (default: {default_ncomp}).")
    body.append("")
    body.append("Returns")
    body.append("-------")
    body.append("dict")
    body.append("    Unit operation configuration dictionary.")
    body.append('"""')
    body.append("")

    # Unit dict
    unit_type = _get_unit_type(column_model)
    body.append("unit = {}")
    body.append("")
    body.append(f"unit['UNIT_TYPE'] = '{unit_type}'")
    body.append(f"unit['ncomp'] = {_ncomp_var(ncomp_fixed)}")

    # Column geometry parameters
    if column_model.resolution == "0D":
        body.append(f"unit['init_volume'] = {_DEFAULTS['init_volume']}")
        if column_model.has_filter:
            body.append(f"unit['flowrate_filter'] = 0.0")
    elif column_model.column_type == "Radial":
        body.append(f"unit['col_radius_inner'] = {_DEFAULTS['col_radius_inner']}")
        body.append(f"unit['col_radius_outer'] = {_DEFAULTS['col_radius_outer']}")
        body.append(f"unit['velocity_coeff'] = {_DEFAULTS['velocity_coeff']}")
    elif column_model.column_type == "Frustum":
        body.append(f"unit['col_length'] = {_DEFAULTS['col_length']}")
        body.append(f"unit['col_radius_inlet'] = {_DEFAULTS['col_radius_inner']}")
        body.append(f"unit['col_radius_outlet'] = {_DEFAULTS['col_radius_outer']}")
        body.append(f"unit['velocity'] = {_DEFAULTS['velocity']}")
    else:
        body.append(f"unit['col_length'] = {_DEFAULTS['col_length']}")
        body.append(f"unit['velocity'] = {_DEFAULTS['velocity']}")

    # Porosity
    if column_model.N_p > 0:
        if column_model.nonlimiting_filmDiff and column_model.particle_models[0].resolution == "0D":
            total_por = round(
                _DEFAULTS['col_porosity'] + (1.0 - _DEFAULTS['col_porosity']) * _DEFAULTS['par_porosity'], 4
            )
            body.append(f"unit['total_porosity'] = {total_por}")
        else:
            body.append(f"unit['col_porosity'] = {_DEFAULTS['col_porosity']}")
    elif column_model.resolution != "0D":
        body.append(f"unit['col_porosity'] = 1.0")
        body.append(f"unit['total_porosity'] = 1.0")

    # Dispersion
    if column_model.has_axial_dispersion:
        body.append(f"unit['col_dispersion'] = {_DEFAULTS['col_dispersion']}")
    elif column_model.resolution != "0D":
        body.append(f"unit['col_dispersion'] = 0.0")

    if column_model.has_radial_dispersion:
        body.append(f"unit['col_dispersion_radial'] = {_DEFAULTS['col_dispersion']}")

    # Initial conditions
    body.append(f"unit['init_c'] = {_per_comp(ncomp_fixed, 0.0)}")

    # Particle configuration
    if column_model.N_p > 0:
        body.append("")
        body.append(f"unit['npartype'] = {column_model.N_p}")
        if column_model.N_p == 1:
            body.append(f"unit['par_type_volfrac'] = [1.0]")
        else:
            volfrac = [round(1.0 / column_model.N_p, 4)] * column_model.N_p
            body.append(f"unit['par_type_volfrac'] = {volfrac}")

        for j, particle in enumerate(column_model.particle_models):
            bnd_model = particle.binding_model if column_model.N_p > 1 else column_model.binding_model
            req_bnd = particle.req_binding if column_model.N_p > 1 else column_model.req_binding
            mult_bnd = particle.has_mult_bnd_states if column_model.N_p > 1 else column_model.has_mult_bnd_states

            body.extend(_particle_lines(
                particle, j, ncomp_fixed, bnd_model, req_bnd,
                column_model.has_binding, mult_bnd,
            ))

    body.append("")
    body.append("return unit")

    lines.extend(_indent(body))
    lines.append("")

    return "\n".join(lines)


_CRY_DEFAULTS = {
    'nbins': 50,
    'x_min': 1e-6,
    'x_max': 1e-3,
    'init_volume': 1e-3,
    'col_length': 0.47,
    'col_dispersion': 4.2e-5,
    'cross_section_area': 3.066e-5,
    'growth_rate_constant': 5e-6,
    'growth_dispersion_rate': 2.5e-15,
    'primary_nucleation_rate': 5.0,
    'secondary_nucleation_rate': 4e8,
    'nuclei_mass_density': 1200.0,
    'vol_shape_factor': 0.524,
    'aggregation_rate_constant': 5e-13,
    'fragmentation_rate_constant': 6e10,
    'fragmentation_kernel_gamma': 2,
    'fragmentation_selection_alpha': 1,
}

def generate_crystallization_script(cry_model):
    """Return a Python script string for the current Crystallization configuration."""

    model_name = cry_model.model_name()

    lines = []

    lines.append(f'"""')
    lines.append(f"CADET unit operation template: {model_name}")
    lines.append(f"Generated by CADET-Equations")
    lines.append(f'"""')
    lines.append("")
    lines.append("import numpy as np")
    lines.append("")
    lines.append("")

    nbins = _CRY_DEFAULTS['nbins']
    ncomp = nbins + 1

    lines.append(f"def get_unit_operation(nbins={nbins}):")

    body = []
    body.append(f'"""Return a CADET unit operation dictionary for a {model_name}.')
    body.append("")
    body.append("Parameters")
    body.append("----------")
    body.append("nbins : int")
    body.append(f"    Number of crystal size bins (default: {nbins}).")
    body.append("")
    body.append("Returns")
    body.append("-------")
    body.append("dict")
    body.append("    Unit operation configuration dictionary.")
    body.append('"""')
    body.append("")

    body.append("ncomp = nbins + 1")
    body.append("")
    body.append("unit = {}")
    body.append("")

    if cry_model.column_type == "CSTR":
        body.append("unit['UNIT_TYPE'] = 'CSTR'")
        body.append("unit['ncomp'] = ncomp")
        body.append(f"unit['init_volume'] = {_CRY_DEFAULTS['init_volume']}")
    else:
        body.append("unit['UNIT_TYPE'] = 'LUMPED_RATE_MODEL_WITHOUT_PORES'")
        body.append("unit['ncomp'] = ncomp")
        body.append(f"unit['col_length'] = {_CRY_DEFAULTS['col_length']}")
        body.append(f"unit['cross_section_area'] = {_CRY_DEFAULTS['cross_section_area']}")
        body.append("unit['total_porosity'] = 1.0")
        if cry_model.has_axial_dispersion:
            body.append(f"unit['col_dispersion'] = {_CRY_DEFAULTS['col_dispersion']}")
        else:
            body.append("unit['col_dispersion'] = 0.0")

    body.append("unit['adsorption_model'] = 'NONE'")
    body.append("unit['reaction_model'] = 'CRYSTALLIZATION'")
    body.append("unit['init_c'] = [0.0] * ncomp")
    body.append("")

    if cry_model.column_type == "CSTR":
        reaction_var = "reaction_bulk"
    else:
        reaction_var = "reaction"

    body.append(f"cry = {{}}")
    body.append(f"cry['cry_bins'] = np.linspace({_CRY_DEFAULTS['x_min']}, {_CRY_DEFAULTS['x_max']}, nbins).tolist()")
    body.append("")

    has_primary = cry_model.has_primary_formation
    has_agg = cry_model.has_aggregation
    has_frag = cry_model.has_fragmentation

    cry_mode = (1 if has_primary else 0) | (2 if has_agg else 0) | (4 if has_frag else 0)
    body.append(f"cry['cry_mode'] = {cry_mode}")

    if has_primary:
        body.append(f"cry['cry_growth_rate_constant'] = {_CRY_DEFAULTS['growth_rate_constant']}")
        body.append("cry['cry_a'] = 1.0")
        body.append("cry['cry_g'] = 1.0")
        body.append("cry['cry_k'] = 1.0")
        body.append(f"cry['cry_primary_nucleation_rate'] = {_CRY_DEFAULTS['primary_nucleation_rate']}")
        body.append("cry['cry_u'] = 10.0")
        body.append("cry['cry_b'] = 2.0")
        body.append(f"cry['cry_nuclei_mass_density'] = {_CRY_DEFAULTS['nuclei_mass_density']}")
        body.append(f"cry['cry_vol_shape_factor'] = {_CRY_DEFAULTS['vol_shape_factor']}")

        if cry_model.size_dependent_growth:
            body.append("cry['cry_p'] = 1.0")
            body.append("cry['cry_growth_constant'] = 1.0")
        else:
            body.append("cry['cry_p'] = 0.0")
            body.append("cry['cry_growth_constant'] = 0.0")

        if cry_model.has_growth_dispersion:
            body.append(f"cry['cry_growth_dispersion_rate'] = {_CRY_DEFAULTS['growth_dispersion_rate']}")
        else:
            body.append("cry['cry_growth_dispersion_rate'] = 0.0")

        if cry_model.has_secondary_nucleation:
            body.append(f"cry['cry_secondary_nucleation_rate'] = {_CRY_DEFAULTS['secondary_nucleation_rate']}")
        else:
            body.append("cry['cry_secondary_nucleation_rate'] = 0.0")

    if has_agg:
        body.append(f"cry['cry_aggregation_rate_constant'] = {_CRY_DEFAULTS['aggregation_rate_constant']}")
        body.append(f"cry['cry_aggregation_index'] = {cry_model.aggregation_kernel_index}")

    if has_frag:
        body.append(f"cry['cry_fragmentation_rate_constant'] = {_CRY_DEFAULTS['fragmentation_rate_constant']}")
        body.append(f"cry['cry_fragmentation_kernel_gamma'] = {_CRY_DEFAULTS['fragmentation_kernel_gamma']}")
        body.append(f"cry['cry_fragmentation_selection_function_alpha'] = {_CRY_DEFAULTS['fragmentation_selection_alpha']}")

    body.append("")
    body.append(f"unit['{reaction_var}'] = cry")
    body.append("")
    body.append("return unit")

    lines.extend(_indent(body))
    lines.append("")

    return "\n".join(lines)
