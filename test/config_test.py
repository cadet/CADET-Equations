# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

from streamlit.testing.v1 import AppTest


def test_streamlit_app():
    
    at = AppTest.from_file("../app.py")
    
    at.run()
    assert not at.exception

    assert at.selectbox(key="column_resolution").value == "1D (axial coordinate)"

    # We test every configuration by recursively iteration through all combinations
    def config_recursion(boxies, counter):

        boxy = at.selectbox(key=boxies.pop(0)) # TODO allow other boxes. if number input define min max options

        for opt in boxy.options:

            boxy.set_value(opt).run()
            counter += 1

            assert not at.exception

            if len(boxies) > 0:
                counter = config_recursion(boxies, counter)

        return counter
    
    # Some boxes are conditional, e.g. film_diffusion can only be configured when particles are present.
    # We thus first identify the boxes that change the number of boxes and thus combinations.
    # Their keys are: dev_mode, advanced_mode, add_particles
    # dev_mode is not tested, TODO: advanced_mode
    critical_variables = ["dev_mode", "advanced_mode",
                          "add_particles"]


    ########## 1st iteration: add_particles = False
    # check for the critical variables
    assert at.selectbox(key="dev_mode").value == "Off"
    assert at.selectbox(key="advanced_mode").value == "Off"
    assert at.selectbox(key="add_particles").value == "No"

    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]

    # check if variables were added or removed
    assert len(config_keys) == 2
    assert len(at.get("number_input")) == 1
    assert len(at.get("selectbox") + at.get("number_input")) == 6

    # test all configs with add_particles = "No"
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 5

    # ########## 2nd iteration: add_particles = True # TODO : noBinding aber fragt trotzdem nach surfDIff
    # at.selectbox(key="add_particles").set_value("Yes").run()
    # # check for the critical variables
    # assert at.selectbox(key="dev_mode").value == "Off"
    # assert at.selectbox(key="advanced_mode").value == "Off"
    # assert at.selectbox(key="add_particles").value == "Yes"

    # config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    # assert len(config_keys) == 6
    # assert len(at.get("number_input")) == 1
    # assert len(at.get("selectbox") + at.get("number_input")) == 10

    # # test all configs with add_particles = "Yes"
    # num_configs = config_recursion(config_keys, 0)
    # assert num_configs == 0










    # selectboxes = at.selectbox

    # for box in selectboxes:

    #     # for option in at.selectbox(key="column_resolution").options:
    #     for option in box.options:

    #         # at.selectbox(key="column_resolution").set_value(option).run()
    #         # boxy.set_value(option).run() # fails, proving that its doing the loops
    #         box.set_value(option).run()

