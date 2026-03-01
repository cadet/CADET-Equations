# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import streamlit as st
import h5py

import re

CADET_column_unit_types_v5 = [
                    'GENERAL_RATE_MODEL', 'LUMPED_RATE_MODEL_WITHOUT_PORES', 'LUMPED_RATE_MODEL_WITH_PORES',
                    'GENERAL_RATE_MODEL_DG', 'LUMPED_RATE_MODEL_WITHOUT_PORES_DG', 'LUMPED_RATE_MODEL_WITH_PORES_DG',
                    'GENERAL_RATE_MODEL_2D', 'CSTR'
                ]

CADET_column_unit_types_v6 = [
                    'COLUMN_MODEL_1D', 'COLUMN_MODEL_2D', 'RADIAL_COLUMN_MODEL_1D', 
                    'FRUSTUM_COLUMN_MODEL_1D', 'COLUMN_MODEL_3D'
                ]

CADET_column_unit_types = CADET_column_unit_types_v5 + CADET_column_unit_types_v6



def get_h5_value(unit_group, key:str, firstEntryIfList=True):

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

def map_unit_type_to_column_model(cadet_unit_type):

    if re.search("3D", cadet_unit_type):
        return "3D (axial, radial and angular coordinate)"
    elif re.search("2D", cadet_unit_type):
        return "2D (axial and radial coordinate)"
    elif re.search("CSTR", cadet_unit_type):
        return "0D (Homogeneous Tank)"
    elif cadet_unit_type in CADET_column_unit_types:
        return "1D (axial coordinate)"
    else:
        raise ValueError(
            f"Invalid unit type: {cadet_unit_type}. Must be one of {CADET_column_unit_types}.")


def map_unit_to_particle_model(cadet_unit_type, h5_unit_group, is_v6=False):

    if re.search("WITHOUT_PORES", cadet_unit_type) or re.search("CSTR", cadet_unit_type):
        
        if get_cadet_unit_value(h5_unit_group, 'TOTAL_POROSITY', is_v6) == 1.0:
            return None
        if get_cadet_unit_value(h5_unit_group, 'CONST_SOLID_VOLUME', is_v6) == 0.0:
            return None
        
    elif get_cadet_unit_value(h5_unit_group, 'COL_POROSITY', is_v6) == 1.0:
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


def get_cadet_unit_value(unit_group, key, is_v6=False, firstEntryIfList=True):
    if is_v6:
        for sub in ['model', 'discretization']:
            if sub in unit_group:
                val = get_h5_value(unit_group[sub], key, firstEntryIfList)
                if val is not None:
                    return val
    return get_h5_value(unit_group, key, firstEntryIfList)


def get_cadet_subgroup(unit_group, subgroup_path, is_v6=False):
    if is_v6:
        # Most subgroups like 'adsorption' should be under 'model' in v6
        v6_path = 'model/' + subgroup_path
        if v6_path in unit_group:
            return unit_group[v6_path]
    return unit_group[subgroup_path] if subgroup_path in unit_group else None


def extract_config_data_from_unit(unit_type, h5_unit_group, is_v6=False):

    config = {}

    config['advanced_mode'] = "Off"
    config['dev_mode'] = "Off"
    config['var_format'] = "CADET"
    config['show_eq_description'] = True
    config['model_assumptions'] = True

    config['column_resolution'] = map_unit_type_to_column_model(unit_type)

    if re.search("0D", config['column_resolution']):
        flow_filter = get_cadet_unit_value(h5_unit_group, 'FLOWRATE_FILTER', is_v6)
        config['has_filter'] = "No"
        if flow_filter is not None:
            config['has_filter'] = "Yes" if flow_filter > 0.0 else "No"

    if re.search("2D", config['column_resolution']):

        Dax = get_cadet_unit_value(h5_unit_group, 'COL_DISPERSION', is_v6)
        if Dax is not None:
            config['has_axial_dispersion'] = "No" if Dax < 1E-20 else "Yes"
        
        Drad = get_cadet_unit_value(h5_unit_group, 'COL_DISPERSION_RADIAL', is_v6)
        if Drad is not None:
            config['has_radial_dispersion'] = "No" if Drad < 1E-20 else "Yes"

    elif re.search("1D", config['column_resolution']):

        Dax = get_cadet_unit_value(h5_unit_group, 'COL_DISPERSION', is_v6)
        if Dax is not None:
            config['has_axial_dispersion'] = "No" if Dax < 1E-20 else "Yes"
    
    par_model = map_unit_to_particle_model(unit_type, h5_unit_group, is_v6)

    if par_model is not None:

        config['add_particles'] = "Yes"

        config['particle_resolution'] = par_model

        config['nonlimiting_filmDiff'] = "Yes" if re.search("WITHOUT_PORES", unit_type) else "No"

        nParType = get_cadet_unit_value(h5_unit_group, 'NPARTYPE', is_v6)
        nParType = 1 if nParType is None else nParType

        if nParType > 1:
            
            config['advanced_mode'] = "On"

            config['PSD'] = "Yes"

        binding_model = get_cadet_unit_value(h5_unit_group, 'ADSORPTION_MODEL', is_v6, firstEntryIfList=False)
        
        config['has_binding'] = "No"

        if binding_model is not None:
            
            config['PTD'] = "No"
            if not isinstance(binding_model, str):
                if len(binding_model) > 1:
                    return binding_model
                    if len(set(binding_model)) > 1:
                        config['PTD'] = "Yes"
                    binding_model = binding_model[0]
            
            if binding_model != "NONE":

                config['has_binding'] = "Yes"
                
                config['has_mult_bnd_states'] = "No"

                if nParType > 1:
                    ads_group = get_cadet_subgroup(h5_unit_group, 'adsorption_000', is_v6)
                    config['req_binding'] = "Kinetic" if get_h5_value(ads_group, 'IS_KINETIC') else "Rapid-equilibrium"
                    if get_h5_value(ads_group, 'NBOUND') is not None:
                        config['has_mult_bnd_states'] = "Yes" if get_h5_value(ads_group, 'NBOUND') > 1 else "No"
                else:
                    ads_group = get_cadet_subgroup(h5_unit_group, 'adsorption', is_v6)
                    config['req_binding'] = "Kinetic" if get_h5_value(ads_group, 'IS_KINETIC') else "Rapid-equilibrium"
                    if get_h5_value(ads_group, 'NBOUND') is not None:
                        config['has_mult_bnd_states'] = "Yes" if get_h5_value(ads_group, 'NBOUND') > 1 else "No"

                if par_model == "1D (radial coordinate)":

                    config['has_surfDiff'] = "No"
                    surfDiff = get_cadet_unit_value(h5_unit_group, 'PAR_SURFDIFFUSION', is_v6)
                    if surfDiff is not None:
                        config['has_surfDiff'] = "Yes" if surfDiff > 0.0 else "No"

                    config['particle_has_core'] = "No"
                    parCore = get_cadet_unit_value(h5_unit_group, 'PAR_CORERADIUS', is_v6)
                    if parCore is not None:
                        if parCore > 0.0:
                            config['particle_has_core'] = "Yes"
                            config['advanced_mode'] = "On"

    else:

        config['add_particles'] = "No"

    if not config['advanced_mode'] == "On":
        config.pop('dev_mode', None)
        config.pop('particle_has_core', None)
        config.pop('has_radial_dispersion', None)
        config.pop('has_mult_bnd_states', None)
        config.pop('PTD', None)

    return config


def get_config_from_CADET_h5(h5_filename, unit_idx):

    with h5py.File(h5_filename, 'r') as f:

        # Check for CADET version
        is_v6 = False
        if 'meta' in f:
            sim_version = get_h5_value(f['meta'], 'SIMULATOR_VERSION')
            if sim_version and str(sim_version).startswith('6'):
                is_v6 = True
        
        model_group = f['input/model']

        if unit_idx == "-01":
        
            unit_keys = [k for k in model_group.keys() if re.match(r'^unit_\d{3}$', k)]
            unit_keys.sort(key=lambda x: int(x.split('_')[1]))

            for unit_key in unit_keys:

                unit_group = model_group[unit_key]

                unit_type = get_h5_value(unit_group, 'UNIT_TYPE')

                if unit_type is not None:

                    # Optional: infer v6 from unit_type if not explicitly set in meta
                    if not is_v6 and unit_type in CADET_column_unit_types_v6:
                        is_v6 = True

                    if unit_type in CADET_column_unit_types:
                        
                        st.sidebar.success(unit_type + " was found in " + re.sub(r"input/model/", "", unit_key) + " and is applied!" + (" (CADET v6 detected)" if is_v6 else ""))

                        return extract_config_data_from_unit(unit_type, unit_group, is_v6)

        else:

            if 'unit_' + unit_idx in model_group.keys(): 
                unit_group = model_group['unit_' + unit_idx]
                unit_type = get_h5_value(unit_group, 'UNIT_TYPE')
            else:
                st.error(f"unit_{unit_idx} does not exist!")
                return None
            
            if unit_type is not None:
                if not is_v6 and unit_type in CADET_column_unit_types_v6:
                    is_v6 = True

                if unit_type in CADET_column_unit_types:

                    st.sidebar.success(f"Unit type {unit_type} is applied." + (" (CADET v6 detected)" if is_v6 else ""))

                    return extract_config_data_from_unit(unit_type, unit_group, is_v6)

                else:
                    st.error(f"Equations for {unit_type} are not available.")
            
        return None
