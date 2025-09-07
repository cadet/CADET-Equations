# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import streamlit as st
import h5py

import re

CADET_column_unit_types = [
                    'GENERAL_RATE_MODEL', 'LUMPED_RATE_MODEL_WITHOUT_PORES', 'LUMPED_RATE_MODEL_WITH_PORES',
                    'GENERAL_RATE_MODEL_DG', 'LUMPED_RATE_MODEL_WITHOUT_PORES_DG', 'LUMPED_RATE_MODEL_WITH_PORES_DG',
                    'GENERAL_RATE_MODEL_2D', 'CSTR'
                ]


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


def map_unit_to_particle_model(cadet_unit_type, h5_unit_group):

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


def extract_config_data_from_unit(unit_type, h5_unit_group):

    config = {}

    config['advanced_mode'] = "Off"
    config['dev_mode'] = "Off"
    config['var_format'] = "CADET"
    config['show_eq_description'] = True
    config['model_assumptions'] = True

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

        config['nonlimiting_filmDiff'] = "Yes" if re.search("WITHOUT_PORES", unit_type) else "No"

        nParType = get_h5_value(h5_unit_group, 'NPARTYPE')
        nParType = 1 if nParType is None else nParType

        if nParType > 1:
            
            config['advanced_mode'] = "On"

            config['PSD'] = "Yes"

        binding_model = get_h5_value(h5_unit_group, 'ADSORPTION_MODEL', firstEntryIfList=False)
        
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
                    config['req_binding'] = "Kinetic" if get_h5_value(h5_unit_group['adsorption_000'], 'IS_KINETIC') else "Rapid-equilibrium"
                    if get_h5_value(h5_unit_group['adsorption_000'], 'NBOUND') is not None:
                        config['has_mult_bnd_states'] = "Yes" if get_h5_value(h5_unit_group['adsorption_000'], 'NBOUND') > 1 else "No"
                else:
                    config['req_binding'] = "Kinetic" if get_h5_value(h5_unit_group['adsorption'], 'IS_KINETIC') else "Rapid-equilibrium"
                    if get_h5_value(h5_unit_group['adsorption'], 'NBOUND') is not None:
                        config['has_mult_bnd_states'] = "Yes" if get_h5_value(h5_unit_group['adsorption'], 'NBOUND') > 1 else "No"

                if par_model == "1D (radial coordinate)":

                    config['has_surfDiff'] = "No"
                    surfDiff = get_h5_value(h5_unit_group, 'PAR_SURFDIFFUSION')
                    if surfDiff is not None:
                        config['has_surfDiff'] = "Yes" if surfDiff > 0.0 else "No"

                    config['particle_has_core'] = "No"
                    parCore = get_h5_value(h5_unit_group, 'PAR_CORERADIUS')
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

                        return extract_config_data_from_unit(unit_type, unit_group)

        else:

            if 'unit_' + unit_idx in model_group.keys(): 
                unit_group = model_group['unit_' + unit_idx]
                unit_type = get_h5_value(unit_group, 'UNIT_TYPE')
            else:
                st.error(f"unit_{unit_idx} does not exist!")
                return None

            if unit_type in CADET_column_unit_types:

                st.sidebar.success(f"{unit_type} unit_{unit_idx} is applied.")

                return extract_config_data_from_unit(unit_type, unit_group)

            else:
                st.error(f"Equations for {unit_type} are not available.")
            
        return None
