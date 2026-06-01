# -*- coding: utf-8 -*-
"""
Helpers to extract a generator configuration from CADET HDF5 files.

This module provides minimal mapping functions used when a user
uploads a CADET HDF5 file to pre-populate the generator UI.
"""

import streamlit as st
import h5py

import re

CADET_binding_model_map = {
    'LINEAR': 'Linear',
    'MULTI_COMPONENT_LANGMUIR': 'Langmuir',
    'STERIC_MASS_ACTION': 'SMA',
}


CADET_column_unit_types = [
                    'CSTR',
                    'COLUMN_MODEL_1D', 'COLUMN_MODEL_2D',
                    'GENERAL_RATE_MODEL', 'LUMPED_RATE_MODEL_WITHOUT_PORES', 'LUMPED_RATE_MODEL_WITH_PORES',
                    'GENERAL_RATE_MODEL_2D',
                    'FRUSTUM_COLUMN_MODEL_1D',
                    'FRUSTUM_GENERAL_RATE_MODEL', 'FRUSTUM_LUMPED_RATE_MODEL_WITHOUT_PORES', 'FRUSTUM_LUMPED_RATE_MODEL_WITH_PORES',
                    'RADIAL_COLUMN_MODEL_1D',
                    'RADIAL_GENERAL_RATE_MODEL', 'RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES', 'RADIAL_LUMPED_RATE_MODEL_WITH_PORES',
                    # old interface
                    'GENERAL_RATE_MODEL_DG', 'LUMPED_RATE_MODEL_WITHOUT_PORES_DG', 'LUMPED_RATE_MODEL_WITH_PORES_DG',
                ]


def is_v6_interface(unit_type, h5_unit_group):
    """Detect a v6-style CADET HDF5 unit group.

    Returns True when the file structure matches the newer v6 layout.
    """
    if re.search("COLUMN_MODEL", unit_type):
        return True
    nParType = get_h5_value(h5_unit_group, 'NPARTYPE')
    if nParType is not None and nParType >= 1 and 'particle_type_000' in h5_unit_group:
        return True
    return False


def get_h5_value(unit_group, key:str, firstEntryIfList=True):
    """Read a value from an HDF5 group and return a native Python type.

    The function decodes bytes and returns the first element of array-like
    values when appropriate.
    """

    value = unit_group.get(key)

    if value is None:
        return None

    value = value[()]

    if isinstance(value, bytes):

        value = value.decode('utf-8')

    # If the value is a NumPy array or similar, get the first element
    if hasattr(value, '__len__') and not isinstance(value, (str, bytes)):
        if len(value) > 0 and firstEntryIfList:
            value = value[0]

    return value


def map_unit_type_to_column_geometry(cadet_unit_type):
    """Map a CADET unit type string to the generator's column geometry label."""

    if re.search("RADIAL", cadet_unit_type):
        return "Radial flow cylinder"
    elif re.search("FRUSTUM", cadet_unit_type):
        return "Frustum"
    else:
        return "Axial flow cylinder"


def map_unit_type_to_column_model(cadet_unit_type):
    """Map a CADET unit type string to the generator's column resolution label."""

    if re.search("3D", cadet_unit_type):
        return "3D (axial, radial and angular coordinate)"
    elif re.search("2D", cadet_unit_type):
        return "2D (axial and radial coordinate)"
    elif re.search("CSTR", cadet_unit_type):
        return "0D (Homogeneous Tank)"
    elif re.search("RADIAL", cadet_unit_type):
        return "1D (radial coordinate)"
    elif cadet_unit_type in CADET_column_unit_types:
        return "1D (axial coordinate)"
    else:
        raise ValueError(
            f"Invalid unit type: {cadet_unit_type}. Must be one of {CADET_column_unit_types}.")


def map_unit_to_particle_model(cadet_unit_type, h5_unit_group):
    """Return a textual particle model mapping for a CADET unit type/group."""

    if is_v6_interface(cadet_unit_type, h5_unit_group):
        return _map_v6_particle_model(h5_unit_group)

    if re.search("WITHOUT_PORES", cadet_unit_type) or re.search("CSTR", cadet_unit_type):

        if get_h5_value(h5_unit_group, 'TOTAL_POROSITY') == 1.0:
            return None
        if get_h5_value(h5_unit_group, 'CONST_SOLID_VOLUME') == 0.0:
            return None

    elif get_h5_value(h5_unit_group, 'COL_POROSITY') == 1.0:
        return None

    if re.search("GENERAL_RATE", cadet_unit_type):
        return "1D (radial coordinate)"
    elif re.search("LUMPED_RATE", cadet_unit_type) or re.search("CSTR", cadet_unit_type):
        return "0D (homogeneous)"
    elif cadet_unit_type in CADET_column_unit_types:
        return None
    else:
        raise ValueError(
            f"Invalid unit type: {cadet_unit_type}. Must be one of {CADET_column_unit_types}.")


def _map_v6_particle_model(h5_unit_group):
    """Determine particle model for v6 interface using HAS_* flags in particle_type_000."""

    nParType = get_h5_value(h5_unit_group, 'NPARTYPE')
    if nParType is None or nParType < 1:
        return None

    pt_group = h5_unit_group.get('particle_type_000')
    if pt_group is None:
        return None

    has_pore_diff = get_h5_value(pt_group, 'HAS_PORE_DIFFUSION')
    has_surf_diff = get_h5_value(pt_group, 'HAS_SURFACE_DIFFUSION')

    if has_pore_diff or has_surf_diff:
        return "1D (radial coordinate)"

    has_film_diff = get_h5_value(pt_group, 'HAS_FILM_DIFFUSION')
    if has_film_diff:
        return "0D (homogeneous)"

    # Equilibrium particle (no film/pore/surface diffusion) — LRM-like
    return "0D (homogeneous)"


def extract_config_data_from_unit(unit_type, h5_unit_group):

    config = {}

    config['advanced_mode'] = "Off"
    config['var_format'] = "CADET"
    config['show_eq_description'] = True
    config['model_assumptions'] = True

    config['column_type'] = map_unit_type_to_column_geometry(unit_type)
    config['column_resolution'] = map_unit_type_to_column_model(unit_type)

    if re.search("0D", config['column_resolution']):
        flow_filter = get_h5_value(h5_unit_group, 'FLOWRATE_FILTER')
        config['has_filter'] = "No"
        if flow_filter is not None:
            config['has_filter'] = "Yes" if flow_filter > 0.0 else "No"

    if re.search("2D", config['column_resolution']):

        Dax = get_h5_value(h5_unit_group, 'COL_DISPERSION')
        if Dax is not None:
            config['has_axial_dispersion'] = "No" if Dax < 1E-20 else "Yes"
        
        Drad = get_h5_value(h5_unit_group, 'COL_DISPERSION_RADIAL')
        if Drad is not None:
            config['has_radial_dispersion'] = "No" if Drad < 1E-20 else "Yes"

    elif re.search("1D", config['column_resolution']):

        Dax = get_h5_value(h5_unit_group, 'COL_DISPERSION')
        if Dax is not None:
            config['has_axial_dispersion'] = "No" if Dax < 1E-20 else "Yes"
    
    par_model = map_unit_to_particle_model(unit_type, h5_unit_group)

    if par_model is not None:

        config['add_particles'] = "Yes"

        config['particle_resolution'] = par_model

        if is_v6_interface(unit_type, h5_unit_group):
            _extract_v6_particle_config(config, h5_unit_group, par_model)
        else:
            _extract_v5_particle_config(config, unit_type, h5_unit_group, par_model)

    else:

        config['add_particles'] = "No"

    if config['advanced_mode'] == "On":
        # Advanced mode uses single PSD selectbox with 3 options
        add_par = config.pop('add_particles', "No")
        if add_par == "Yes":
            if config.get('PSD') == "Yes":
                config['PSD'] = "Particle size distribution"
            else:
                config['PSD'] = "Yes"
        else:
            config['PSD'] = "No"
    else:
        config.pop('dev_mode', None)
        config.pop('PSD', None)
        config.pop('particle_has_core', None)
        config.pop('has_radial_dispersion', None)
        config.pop('has_mult_bnd_states', None)

    return config


def _extract_v5_particle_config(config, unit_type, h5_unit_group, par_model):
    """Extract particle configuration from v5 interface (particle info at unit level)."""

    nonlimiting = bool(re.search("WITHOUT_PORES", unit_type))

    nParType = get_h5_value(h5_unit_group, 'NPARTYPE')
    nParType = 1 if nParType is None else nParType

    if nParType > 1:
        config['advanced_mode'] = "On"
        config['PSD'] = "Yes"

    config['particle_nonlimiting_filmDiff'] = "Yes" if nonlimiting else "No"

    binding_model = get_h5_value(h5_unit_group, 'ADSORPTION_MODEL', firstEntryIfList=False)

    config['has_binding'] = "No"

    if binding_model is not None:

        if not isinstance(binding_model, str):
            if len(binding_model) > 1:
                binding_model = binding_model[0]

        if binding_model != "NONE":

            config['has_binding'] = "Yes"

            config['binding_model'] = CADET_binding_model_map.get(binding_model, "Arbitrary")
            if binding_model not in CADET_binding_model_map:
                st.sidebar.warning(f"Binding model {binding_model} not implemented in CADET-Equations, default to arbitrary binding")

            config['has_mult_bnd_states'] = "No"

            if nParType > 1:
                ads_group = h5_unit_group['adsorption_000']
            else:
                ads_group = h5_unit_group['adsorption']

            config['req_binding'] = "Kinetic" if get_h5_value(ads_group, 'IS_KINETIC') else "Rapid-equilibrium"
            if get_h5_value(ads_group, 'NBOUND') is not None:
                config['has_mult_bnd_states'] = "Yes" if get_h5_value(ads_group, 'NBOUND') > 1 else "No"

            if par_model == "1D (radial coordinate)":

                surfDiff = get_h5_value(h5_unit_group, 'PAR_SURFDIFFUSION')
                config['particle_has_surfDiff'] = "Yes" if surfDiff is not None and surfDiff > 0.0 else "No"

    if par_model == "1D (radial coordinate)":
        _extract_particle_core_config(config, h5_unit_group)


def _extract_v6_particle_config(config, h5_unit_group, par_model):
    """Extract particle configuration from v6 interface (particle info in particle_type_xxx subgroups)."""

    nParType = get_h5_value(h5_unit_group, 'NPARTYPE')
    nParType = 1 if nParType is None else nParType

    if nParType > 1:
        config['advanced_mode'] = "On"
        config['PSD'] = "Yes"

    pt_group = h5_unit_group['particle_type_000']

    has_film_diff = get_h5_value(pt_group, 'HAS_FILM_DIFFUSION')
    config['particle_nonlimiting_filmDiff'] = "No" if has_film_diff else "Yes"

    binding_model = get_h5_value(pt_group, 'ADSORPTION_MODEL', firstEntryIfList=False)

    config['has_binding'] = "No"

    if binding_model is not None:

        if not isinstance(binding_model, str):
            if len(binding_model) > 1:
                binding_model = binding_model[0]

        if binding_model != "NONE":

            config['has_binding'] = "Yes"

            config['binding_model'] = CADET_binding_model_map.get(binding_model, "Arbitrary")
            if binding_model not in CADET_binding_model_map:
                st.sidebar.warning(f"Binding model {binding_model} not implemented in CADET-Equations, default to arbitrary binding")

            config['has_mult_bnd_states'] = "No"

            ads_group = pt_group['adsorption']
            config['req_binding'] = "Kinetic" if get_h5_value(ads_group, 'IS_KINETIC') else "Rapid-equilibrium"
            if get_h5_value(pt_group, 'NBOUND') is not None:
                config['has_mult_bnd_states'] = "Yes" if get_h5_value(pt_group, 'NBOUND') > 1 else "No"

            if par_model == "1D (radial coordinate)":
                has_surf_diff = get_h5_value(pt_group, 'HAS_SURFACE_DIFFUSION')
                config['particle_has_surfDiff'] = "Yes" if has_surf_diff else "No"

    if par_model == "1D (radial coordinate)":
        _extract_particle_core_config(config, pt_group)


def _extract_particle_core_config(config, group):
    """Extract particle core radius config. Shared between v5 (unit group) and v6 (particle_type group)."""
    config['particle_has_core'] = "No"
    parCore = get_h5_value(group, 'PAR_CORERADIUS')
    if parCore is not None:
        if parCore > 0.0:
            config['particle_has_core'] = "Yes"
            config['advanced_mode'] = "On"


CADET_crystallization_unit_types = ['CSTR', 'LUMPED_RATE_MODEL_WITHOUT_PORES']

AGGREGATION_KERNEL_MAP = {
    0: "Constant",
    1: "Brownian",
    2: "Smoluchowski",
    3: "Golovin",
    4: "Differential force",
}


def _is_crystallization_unit(h5_unit_group):
    reaction_model = get_h5_value(h5_unit_group, 'REACTION_MODEL')
    return reaction_model == 'CRYSTALLIZATION'


def _get_crystallization_reaction_group(h5_unit_group):
    for name in ('reaction_bulk', 'reaction'):
        grp = h5_unit_group.get(name)
        if grp is not None:
            return grp
    return None


def extract_crystallization_config(unit_type, h5_unit_group):

    config = {}
    config['model_type'] = 'Crystallization'
    config['dev_mode'] = True
    config['var_format'] = 'CADET'
    config['show_eq_description'] = True
    config['model_assumptions'] = True

    if unit_type == 'CSTR':
        config['cry_column_type'] = 'CSTR'
    else:
        config['cry_column_type'] = 'DPFR'

    reaction_group = _get_crystallization_reaction_group(h5_unit_group)

    if reaction_group is None:
        config['cry_has_primary_formation'] = 'Yes'
        config['cry_has_aggregation'] = 'No'
        config['cry_has_fragmentation'] = 'No'
        return config

    cry_mode = get_h5_value(reaction_group, 'CRY_MODE')
    if cry_mode is None:
        cry_mode = 1

    has_primary = bool(cry_mode & 1)
    has_agg = bool(cry_mode & 2)
    has_frag = bool(cry_mode & 4)

    config['cry_has_primary_formation'] = 'Yes' if has_primary else 'No'

    if has_primary:
        gd_rate = get_h5_value(reaction_group, 'CRY_GROWTH_DISPERSION_RATE')
        config['cry_has_growth_dispersion'] = 'Yes' if (gd_rate is not None and gd_rate > 0) else 'No'

        sec_rate = get_h5_value(reaction_group, 'CRY_SECONDARY_NUCLEATION_RATE')
        config['cry_has_secondary_nucleation'] = 'Yes' if (sec_rate is not None and sec_rate > 0) else 'No'

        cry_p = get_h5_value(reaction_group, 'CRY_P')
        config['cry_size_dependent_growth'] = 'Yes' if (cry_p is not None and cry_p != 0) else 'No'

    config['cry_has_aggregation'] = 'Yes' if has_agg else 'No'

    if has_agg:
        agg_idx = get_h5_value(reaction_group, 'CRY_AGGREGATION_INDEX')
        if agg_idx is None:
            agg_idx = 0
        config['cry_aggregation_kernel'] = AGGREGATION_KERNEL_MAP.get(agg_idx, 'Constant')

    config['cry_has_fragmentation'] = 'Yes' if has_frag else 'No'

    if config['cry_column_type'] == 'DPFR':
        col_disp = get_h5_value(h5_unit_group, 'COL_DISPERSION')
        config['cry_has_axial_dispersion'] = 'Yes' if (col_disp is not None and col_disp > 1E-20) else 'No'

    return config


def get_config_from_CADET_h5(h5_filename, unit_idx):

    with h5py.File(h5_filename, 'r') as f:

        model_group = f['input/model']

        if unit_idx == "-01":
        
            unit_keys = [k for k in model_group.keys() if re.match(r'^unit_\d{3}$', k)]
            unit_keys.sort(key=lambda x: int(x.split('_')[1]))

            for unit_key in unit_keys:

                unit_group = model_group[unit_key]

                unit_type = get_h5_value(unit_group, 'UNIT_TYPE')

                if unit_type is not None:

                    if unit_type in CADET_column_unit_types:

                        st.sidebar.success(unit_type + " was found in " + re.sub(r"input/model/", "", unit_key) + " and is applied!")

                        if _is_crystallization_unit(unit_group):
                            return extract_crystallization_config(unit_type, unit_group)

                        return extract_config_data_from_unit(unit_type, unit_group)

        else:

            if 'unit_' + unit_idx in model_group.keys(): 
                unit_group = model_group['unit_' + unit_idx]
                unit_type = get_h5_value(unit_group, 'UNIT_TYPE')
            else:
                st.error(f"unit_{unit_idx} does not exist!")
                return None

            if unit_type in CADET_column_unit_types:

                st.sidebar.success(f"Unit type {unit_type} is applied.")

                if _is_crystallization_unit(unit_group):
                    return extract_crystallization_config(unit_type, unit_group)

                return extract_config_data_from_unit(unit_type, unit_group)

            else:
                st.error(f"Equations for {unit_type} are not available.")
            
        return None
