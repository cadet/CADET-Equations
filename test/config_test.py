# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

from streamlit.testing.v1 import AppTest


untested_variables = ["dev_mode", "advanced_mode"] # dev_mode is not tested, TODO: test the advanced_mode

# Some boxes are conditional, e.g. film_diffusion can only be configured when particles are present.
# We thus first identify the boxes that change the number of boxes and thus combinations (critical_variables).
# We then start the recursion for every combination of the aforementioned critical keys.
critical_variables = untested_variables + ["add_particles", "has_binding", "particle_resolution"]

def test_streamlit_app():
    
    at = AppTest.from_file("../app.py")
    
    at.run()
    assert not at.exception

    assert at.selectbox(key="column_resolution").value == "1D (axial coordinate)"

    # We test every configuration by recursively iteratiing through all combinations
    def config_recursion(boxies, counter, firstCall:bool=True):

        if boxies[0] in at.session_state:

            boxy = at.selectbox(key=boxies.pop(0)) # TODO allow other boxes. if number input define min max options

            for opt in [opti for opti in boxy.options]: # Note: some values will be run repeatedly since the current value of that box is applied again. Excluding this value here would make the recursion incomplete though

                boxy.set_value(opt).run()
                counter += 1
                assert not at.exception

                for box in at.selectbox:
                    if box.key not in untested_variables and box.key not in critical_variables:
                        print("box: ", box.key, ", value: ", box.value)

                if len(boxies) > 0:
                    counter = config_recursion(boxies, counter, firstCall=False)

            boxies.append(boxy.key)

        return counter

    # check if the untested variables are set correctly
    assert at.selectbox(key="dev_mode").value == "Off"
    assert at.selectbox(key="advanced_mode").value == "Off"

    # 1) test all configs with add_particles = "No" -> no has_binding and no particle_resolution
    assert at.selectbox(key="add_particles").value == "No"
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 6 + 3 # Note: three repetitive configs due to the setup of the recursion

    # 2) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "Yes"
    at.selectbox(key="add_particles").set_value("Yes").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 93 # Note: some repetitive configs due to the setup of the recursion

    # 3) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "Yes"
    at.selectbox(key="particle_resolution").set_value("0D (homogeneous)").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 45 # Note: nine repetitive configs due to the setup of the recursion

    # 4) + 5) test all configs with add_particles = "Yes", "has_binding" = "No" -> no particle_resolution
    at.selectbox(key="has_binding").set_value("No").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 12 + 9 # Note: nine repetitive configs due to the setup of the recursion

