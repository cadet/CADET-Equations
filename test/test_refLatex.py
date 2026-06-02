# -*- coding: utf-8 -*-
"""
This script implements tests comparing generated latex output to reference data
"""

from streamlit.testing.v1 import AppTest
import json
import pytest

from src import load_CADET_h5

# The following test compares the latex output with a reference output file
# At the time of writing, streamlit.testing did not provide functionality for upload_file,
# which is circumvened by applying the configuration variables explicitly
# Also, testing the download_button is not supported yet, which is why this is circumvened
# using the session_state variable latex_string

def read_tex_file(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()
        return content
    except FileNotFoundError:
        print(f"Error: The file {file_path} does not exist.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def _config_key_order(key):
    """Sort config keys so that dependent widgets are set after their parents."""
    order = {
        "dev_mode": -2, "model_type": -1,
        "advanced_mode": 0,
        "column_type": 1.5, "column_resolution": 2,
        "add_particles": 3, "PSD": 3.2, "has_binding": 3.3, "particle_resolution": 3.5,
        "has_reaction_bulk": 3.6, "has_reaction_particle_liquid": 3.7,
        "has_reaction_particle_solid": 3.8, "reaction_model": 3.9,
        "cry_column_type": 4,
        "cry_has_axial_dispersion": 4.5,
        "cry_has_primary_formation": 5,
        "cry_size_dependent_growth": 6,
        "cry_has_growth_dispersion": 6.1,
        "cry_has_secondary_nucleation": 6.2,
        "cry_has_aggregation": 7,
        "cry_aggregation_kernel": 7.1,
        "cry_has_fragmentation": 8,
    }
    return order.get(key, 10)


def apply_model_from_config(at, model_config):

    for config in sorted(model_config.keys(), key=_config_key_order):

        if config == "dev_mode":
            if model_config[config]:
                at.button(key="dev_mode_button").click().run()
        elif config == "model_type":
            if model_config[config] == "Crystallization":
                at.button(key="model_type_crystallization_button").click().run()
        elif config in [box.key for box in at.toggle]:
            at.toggle(key=config).set_value(model_config[config]).run()
        elif config in [box.key for box in at.selectbox]:
            at.selectbox(key=config).set_value(model_config[config]).run()
        elif config in [box.key for box in at.button]:
            if model_config[config]:
                at.button(key=config).click().run()
        else:
            raise ValueError(f"Error: {config} is neither a toggle nor a selectbox nor a button in the current session state")

    return model_config

@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("model_name", ["CSTR", "Plug_Flow", "LRM_dynLin", "LRMP_dynLin", "LRMP_reqLin", "GRMsd_dynLin", "GRM", "GRM_dynLin", "GRMsd_PSD_dynLin", "GRMsd2D_dynLin", "GRMsd_nonLimFD_parCore", "GRMsd_nonLimFD_reqBnd", "GeneralFiniteBath", "Radial_Plug_Flow", "Radial_Dispersive_Plug_Flow", "Radial_GRM", "Radial_LRM", "Radial_LRMP", "Frustum_Plug_Flow", "Frustum_GRM", "cry_CSTR_primaryGrowth", "cry_CSTR_Agg_Frag", "cry_DPFR_aggregation", "cry_DPFR_primarySecondaryNucleation", "GRMsd_MassActionLaw", "GRMsd_MichaelisMenten"])
def test_json_config_output_against_latex_reference(model_name, test_dir):

    at = AppTest.from_file("../Equation-Generator.py")
    
    at.run()
    assert not at.exception

    at.toggle(key="model_assumptions").set_value(True).run()

    json_config_dir = test_dir + '/data/EquationGenerator_configs/'
    ref_latex_dir = test_dir + '/data/ref_latex/'

    with open(json_config_dir + model_name + ".json", 'r') as file:
        
        model_config = json.load(file)

        apply_model_from_config(at, model_config)

    latex_string = at.session_state.latex_string

    ref_string = read_tex_file(ref_latex_dir + model_name + ".tex")

    assert latex_string == ref_string


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("model_name", ["CSTR", "Plug_Flow", "LRM_dynLin", "LRMP_dynLin", "LRMP_reqLin", "GRM", "GRM_dynLin", "GRMsd_dynLin", "GRMsd_PSD_dynLin", "GRMsd2D_dynLin", "GRMsd_nonLimFD_parCore", "GRMsd_nonLimFD_reqBnd"])
def test_CADET_config_output_against_latex_reference(model_name, test_dir):

    CADET_file_ref = {
        "CSTR": None,
        "Plug_Flow": r"PlugFlow_1comp_benchmark1_DG_P3Z1",
        "LRM_dynLin": r"LRM_dynLin_1comp_benchmark1_DG_P3Z1",
        "LRMP_dynLin": r"LRMP_dynLin_1comp_benchmark1_DG_P3Z1",
        "LRMP_reqLin": r"LRMP_reqLin_1comp_benchmark1_DG_P3Z1",
        "GRM": r"GRM_arbBnd_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1",
        "GRM_dynLin": r"GRM_dynLin_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1",
        "GRMsd_dynLin": r"GRMsd_dynLin_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1",
        "GRMsd_PSD_dynLin": r"GRMsdParType2_dynLin_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1",
        "GRMsd2D_dynLin": r"2DGRMsd3Zone_dynLin_1Comp_FV_axZ4radZ3parZ3",
        "GRMsd_nonLimFD_parCore": None, # does not exist in cadet-core
        "GRMsd_nonLimFD_reqBnd": None # does not exist in cadet-core
    }

    ref_name = CADET_file_ref[model_name]

    if ref_name is not None:

        at = AppTest.from_file("../Equation-Generator.py")
    
        at.run()
        assert not at.exception

        at.toggle(key="model_assumptions").set_value(True).run()

        config_file = test_dir + "/data/CADET_configs/" + ref_name + ".h5"

        model_config = load_CADET_h5.get_config_from_CADET_h5(config_file, "-01")

        ###########  test against config json  ##############

        json_config_dir = test_dir + '/data/EquationGenerator_configs/'

        with open(json_config_dir + model_name + ".json", 'r') as file:
        
            model_config_json = json.load(file)

            assert model_config == model_config_json

        ###########  test against latex (redundant when
        # test_json_config_output_against_latex_reference is executed) ##############

        apply_model_from_config(at, model_config)

        latex_string = at.session_state.latex_string

        ref_latex_dir = test_dir + '/data/ref_latex/'

        ref_string = read_tex_file(ref_latex_dir + model_name + ".tex")

        assert latex_string == ref_string


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("model_name", ["CSTR", "Plug_Flow", "LRM_dynLin", "LRMP_dynLin", "LRMP_reqLin", "GRM", "GRM_dynLin", "GRMsd_dynLin", "GRMsd_PSD_dynLin", "GRMsd2D_dynLin"])
def test_CADET_v6_config_output_against_latex_reference(model_name, test_dir):

    CADET_v6_file_ref = {
        "CSTR": "v6_CSTR",
        "Plug_Flow": "v6_PlugFlow_1comp",
        "LRM_dynLin": "v6_LRM_dynLin_1comp",
        "LRMP_dynLin": "v6_LRMP_dynLin_1comp",
        "LRMP_reqLin": "v6_LRMP_reqLin_1comp",
        "GRM": "v6_GRM_arbBnd_1comp",
        "GRM_dynLin": "v6_GRM_dynLin_1comp",
        "GRMsd_dynLin": "v6_GRMsd_dynLin_1comp",
        "GRMsd_PSD_dynLin": "v6_GRMsd_PSD_dynLin_1comp",
        "GRMsd2D_dynLin": "v6_GRMsd2D_dynLin_1comp",
    }

    ref_name = CADET_v6_file_ref[model_name]

    at = AppTest.from_file("../Equation-Generator.py")

    at.run()
    assert not at.exception

    at.toggle(key="model_assumptions").set_value(True).run()

    config_file = test_dir + "/data/CADET_configs/" + ref_name + ".h5"

    model_config = load_CADET_h5.get_config_from_CADET_h5(config_file, "-01")

    ###########  test against config json  ##############

    json_config_dir = test_dir + '/data/EquationGenerator_configs/'

    with open(json_config_dir + model_name + ".json", 'r') as file:

        model_config_json = json.load(file)

        assert model_config == model_config_json

    ###########  test against latex  ##############

    apply_model_from_config(at, model_config)

    latex_string = at.session_state.latex_string

    ref_latex_dir = test_dir + '/data/ref_latex/'

    ref_string = read_tex_file(ref_latex_dir + model_name + ".tex")

    assert latex_string == ref_string


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.debug
def test_CADET_config_manual_unitIdx(test_dir):
    
    ref_latex_dir = test_dir + '/data/ref_latex/'

    at = AppTest.from_file("../Equation-Generator.py")
    
    at.run()
    assert not at.exception
    
    config_file = test_dir + "/data/CADET_configs/ref_acyclicSystem1_LRMP_linBnd_1comp.h5"
    model_config = load_CADET_h5.get_config_from_CADET_h5(config_file, "-01")
    apply_model_from_config(at, model_config)
    latex_string = at.session_state.latex_string

    ref_string = read_tex_file(ref_latex_dir + "LRMP_dynLin.tex")

    assert latex_string == ref_string

    model_config = load_CADET_h5.get_config_from_CADET_h5(config_file, "003")
    apply_model_from_config(at, model_config)
    latex_string = at.session_state.latex_string

    ref_string = read_tex_file(ref_latex_dir + "LRMP_reqLin.tex")

    assert latex_string == ref_string


@pytest.mark.ci
@pytest.mark.reference
def test_crystallization_h5_manual_unitIdx(test_dir):

    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception
    at.toggle(key="model_assumptions").set_value(True).run()

    config_file = test_dir + "/data/CADET_configs/ref_cry_CSTR_PBM_primaryNucleationGrowthGrowthRateDispersion_benchmark1.h5"
    model_config = load_CADET_h5.get_config_from_CADET_h5(config_file, "001")

    json_config_dir = test_dir + '/data/EquationGenerator_configs/'
    with open(json_config_dir + "cry_CSTR_primaryGrowth.json", 'r') as file:
        model_config_json = json.load(file)
        assert model_config == model_config_json


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("model_name", ["cry_CSTR_primaryGrowth", "cry_CSTR_Agg_Frag", "cry_DPFR_aggregation", "cry_DPFR_primarySecondaryNucleation"])
def test_crystallization_h5_config_output(model_name, test_dir):

    CADET_file_ref = {
        "cry_CSTR_primaryGrowth": "ref_cry_CSTR_PBM_primaryNucleationGrowthGrowthRateDispersion_benchmark1",
        "cry_CSTR_Agg_Frag": "ref_cry_CSTR_PBM_Agg_Frag_benchmark1",
        "cry_DPFR_aggregation": "ref_cry_DPFR_PBM_aggregation_benchmark1",
        "cry_DPFR_primarySecondaryNucleation": "ref_cry_DPFR_PBM_primarySecondaryNucleationGrowth_benchmark1",
    }

    ref_name = CADET_file_ref[model_name]

    at = AppTest.from_file("../Equation-Generator.py")

    at.run()
    assert not at.exception

    at.toggle(key="model_assumptions").set_value(True).run()

    config_file = test_dir + "/data/CADET_configs/" + ref_name + ".h5"

    model_config = load_CADET_h5.get_config_from_CADET_h5(config_file, "-01")

    json_config_dir = test_dir + '/data/EquationGenerator_configs/'

    with open(json_config_dir + model_name + ".json", 'r') as file:

        model_config_json = json.load(file)

        assert model_config == model_config_json

    apply_model_from_config(at, model_config)

    latex_string = at.session_state.latex_string

    ref_latex_dir = test_dir + '/data/ref_latex/'

    ref_string = read_tex_file(ref_latex_dir + model_name + ".tex")

    assert latex_string == ref_string
