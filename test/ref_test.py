# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

from streamlit.testing.v1 import AppTest

# The following test compares the latex output with a reference output file
# At the time of writing, streamlit.testing did not provide functionality for upload_file
# which is why the tested configuration is explicitly set in the test via the input widgets
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



def test_latex_model_output_with_reference():

    at = AppTest.from_file("../app.py")
    
    at.run()
    assert not at.exception

    latex_string = at.session_state.latex_string

    ref_string = read_tex_file("test/data/Plug_Flow.tex")

    assert latex_string == ref_string
