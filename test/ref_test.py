# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

from streamlit.testing.v1 import AppTest
import json
import pytest

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

def apply_model_from_config_json(at, config_path:str):

    with open(config_path, 'r') as file:
        
        model_config = json.load(file)

        for config in model_config.keys():

            if config in [box.key for box in at.toggle]:
                at.toggle(key=config).set_value(model_config[config]).run()
            elif config in [box.key for box in at.selectbox]:
                at.selectbox(key=config).set_value(model_config[config]).run()
            else:
                raise ValueError(f"Error: {config} is neither a toggle nor a selectbox")

    return model_config

@pytest.mark.ci
def test_latex_model_output_with_reference():

    at = AppTest.from_file("../app.py")
    
    at.run()
    assert not at.exception

    at.toggle(key="model_assumptions").set_value(True).run()

    model_list = ["Plug_Flow", "GRM", "LRMP"]

    # Test models
    for mode_name in model_list:
        apply_model_from_config_json(at, "test/data/" + mode_name + ".json")

        latex_string = at.session_state.latex_string

        ref_string = read_tex_file("test/data/" + mode_name + ".tex")

        assert latex_string == ref_string


















