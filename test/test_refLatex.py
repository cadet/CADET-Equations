# -*- coding: utf-8 -*-
"""
@author: jmbr
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

def apply_model_from_config(at, model_config):

    for config in model_config.keys():

        if config in [box.key for box in at.toggle]:
            at.toggle(key=config).set_value(model_config[config]).run()
        elif config in [box.key for box in at.selectbox]:
            at.selectbox(key=config).set_value(model_config[config]).run()
        else:
            raise ValueError(f"Error: {config} is neither a toggle nor a selectbox in the current session state")

    return model_config

@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("model_name", ["CSTR", "Plug_Flow", "LRM", "LRMP", "GRMsd", "GRM", "GRMsd_PSD", "GRMsd2D", "GRMsd_nonLimFD_parCore", "GRMsd_nonLimFD_reqBnd"])
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
@pytest.mark.parametrize("model_name", ["CSTR", "Plug_Flow", "LRM", "LRMP", "GRM", "GRMsd", "GRMsd_PSD", "GRMsd2D", "GRMsd_nonLimFD_parCore", "GRMsd_nonLimFD_reqBnd"])
def test_CADET_config_output_against_latex_reference(model_name, test_dir):

    CADET_file_ref = {
        "CSTR": "CSTR",
        "Plug_Flow": r"PlugFlow_1comp_benchmark1_DG_P3Z1",
        "LRM": r"LRM_dynLin_1comp_benchmark1_DG_P3Z1",
        "LRMP": r"LRMP_dynLin_1comp_benchmark1_DG_P3Z1",
        "GRM": r"GRM_dynLin_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1",
        "GRMsd": r"GRMsd_dynLin_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1",
        "GRMsd_PSD": r"GRMsdParType2_dynLin_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1",
        "GRMsd2D": r"2DGRMsd3Zone_dynLin_1Comp_FV_axZ4radZ3parZ3",
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

    ref_string = read_tex_file(ref_latex_dir + "LRMP.tex")

    assert latex_string == ref_string

    model_config = load_CADET_h5.get_config_from_CADET_h5(config_file, "003")
    apply_model_from_config(at, model_config)
    latex_string = at.session_state.latex_string

    ref_string = read_tex_file(ref_latex_dir + "LRMP_reqBnd.tex")

    assert latex_string == ref_string
