# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

from streamlit.testing.v1 import AppTest


untested_variables = ["dev_mode", "advanced_mode"] # dev_mode is not tested, TODO: test the advanced_mode

# Some boxes are conditional, e.g. film_diffusion can only be configured when particles are present.
# We thus first identify the boxes that change the number of boxes and thus combinations (critical_variables).
# We then start the recursion for every combination of the aforementioned critical keys.
# Sometimes, only specific options can be critical (e.g. 0D tank). Such options can be added here and handled independently later
critical_variables = untested_variables + ["add_particles", "has_binding", "particle_resolution", "0D (Homogeneous Tank)"]

def test_streamlit_app():
    
    at = AppTest.from_file("../app.py")
    
    at.run()
    assert not at.exception

    assert at.selectbox(key="column_resolution").value == "1D (axial coordinate)"

    # We test every configuration by recursively iteratiing through all combinations
    def config_recursion(boxies, counter):

        if boxies[0] in at.session_state:

            boxy = at.selectbox(key=boxies.pop(0)) # TODO allow other boxes. if number input define min max options

            for opt in [opti for opti in boxy.options if opti not in critical_variables]: # Note: some values will be run repeatedly since the current value of that box is applied again. Excluding this value here would make the recursion incomplete though

                boxy.set_value(opt).run()
                counter += 1
                assert not at.exception

                for box in at.selectbox:
                    if box.key not in untested_variables and box.key not in critical_variables:
                        print("box: ", box.key, ", value: ", box.value)

                if len(boxies) > 0:
                    counter = config_recursion(boxies, counter)

            boxies.append(boxy.key)

        return counter

    # check if we are in the mode we support
    assert at.selectbox(key="advanced_mode").value == "Off"
    at.selectbox(key="advanced_mode").set_value("On").run()
    assert at.selectbox(key="dev_mode").value == "Off"

    # 1) test all configs with add_particles = "No" -> no has_binding and no particle_resolution
    assert at.selectbox(key="add_particles").value == "No"
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 6 + 3 # Note: three repetitive configs due to the setup of the recursion

    # 2) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "No"
    at.selectbox(key="add_particles").set_value("Yes").run()
    assert at.selectbox(key="particle_resolution").value == "1D (radial coordinate)"
    assert at.selectbox(key="has_binding").value == "No"
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 45 # Note: some repetitive configs due to the setup of the recursion

    # 3) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "No"
    at.selectbox(key="particle_resolution").set_value("0D (homogeneous)").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 12 + 9 # Note: nine repetitive configs due to the setup of the recursion

    # 4) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "Yes"
    at.selectbox(key="has_binding").set_value("Yes").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 45 # Note: repetitive configs due to the setup of the recursion

    # 4) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "Yes"
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 93 # Note: repetitive configs due to the setup of the recursion

    # 6) test all configs for the tank
    # enable the tank and remove column resolution from the test set
    critical_variables.remove("0D (Homogeneous Tank)")
    critical_variables.append("column_resolution")
    at.selectbox(key="column_resolution").set_value("0D (Homogeneous Tank)").run()

    # reset variables
    at.selectbox(key="add_particles").set_value("No").run()

    # 6a) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "No"
    at.selectbox(key="add_particles").set_value("Yes").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    print(config_keys)
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 6 # Note: repetitive configs due to the setup of the recursion

    # 6b) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "No"
    at.selectbox(key="particle_resolution").set_value("0D (homogeneous)").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 2 # Note: repetitive configs due to the setup of the recursion

    # 6c) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "No"
    at.selectbox(key="has_binding").set_value("Yes").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 6 # Note: repetitive configs due to the setup of the recursion

    # 6d) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "Yes"
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    config_keys = [box.key for box in at.selectbox if box.key not in critical_variables]
    num_configs = config_recursion(config_keys, 0)
    assert num_configs == 14 # Note: repetitive configs due to the setup of the recursion
    