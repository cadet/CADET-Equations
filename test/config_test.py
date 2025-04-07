# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

from streamlit.testing.v1 import AppTest
import pytest

untested_variables = ["dev_mode", "advanced_mode", "var_format", "param_table", "show_eq_description"] # dev_mode is not tested, TODO: test the advanced_mode

# Some boxes are conditional, e.g. film_diffusion can only be configured when particles are present.
# We thus first identify the boxes that change the number of boxes and thus combinations (critical_variables).
# We then start the recursion for every combination of the aforementioned critical keys.
# Sometimes, only specific options can be critical (e.g. 0D tank). Such options can be added here and handled independently later
critical_variables = untested_variables + ["add_particles", "has_binding", "particle_resolution", "0D (Homogeneous Tank)"]

# We test every configuration by recursively iteratiing through all combinations
def config_recursion(at, widgies, counter, test_file_generator_buttons:bool):

    if widgies[0] in at.session_state:

        widgyKey = widgies.pop(0)
        if any(key == widgyKey for key in [box.key for box in at.selectbox]):
            widgyKind = "selectbox"
            widgy = at.selectbox(key=widgyKey)
        elif any(key == widgyKey for key in [box.key for box in at.toggle]):
            widgyKind = "toggle"
            widgy = at.toggle(key=widgyKey)
        else:
            raise ValueError(f"Error: {widgyKey} is neither a toggle nor a dropdown")

        options = [opti for opti in widgy.options if opti not in critical_variables] if widgyKind == "selectbox" else [True, False]

        for opt in options: # Note: some values will be run repeatedly since the current value of that widget is applied again. Excluding this value here would make the recursion incomplete though

            widgy.set_value(opt).run()

            if test_file_generator_buttons:

                at.button(key="generate_pdf").click().run()
                at.button(key="generate_config").click().run()

            counter += 1
            assert not at.exception

            for box in at.selectbox:
                if box.key not in untested_variables and box.key not in critical_variables:
                    print("box: ", box.key, ", value: ", box.value)

            if len(widgies) > 0:
                counter = config_recursion(at, widgies, counter, test_file_generator_buttons)

        widgies.append(widgy.key)

    return counter


def run_configs(at, toBeFixed:dict, configKeys:list, numConfigsToBeChecked:int):

    # Check if we have the correct initial configuration
    for fixKey in toBeFixed.keys():

        if any(key == fixKey for key in [box.key for box in at.selectbox]):
            widgyKind = "selectbox"
            widgy = at.selectbox(key=fixKey)
        elif any(key == fixKey for key in [box.key for box in at.toggle]):
            widgyKind = "toggle"
            widgy = at.toggle(key=fixKey)
        else:
            raise ValueError(f"Error: {fixKey} is neither a toggle nor a dropdown")

        if widgy.value != toBeFixed[fixKey]:
            widgy.set_value(toBeFixed[fixKey]).run()

    numConfigs = config_recursion(at, configKeys, 0, True)
    assert numConfigs == numConfigsToBeChecked # Note: repetitive configs due to the setup of the recursion


def getInputToBeTested(at):

    return [box.key for box in at.selectbox if box.key not in critical_variables] + [box.key for box in at.toggle if box.key not in critical_variables]


@pytest.mark.ci
@pytest.mark.all_config_runs
def test_streamlit_app():
    
    at = AppTest.from_file("../Equation-Generator.py")
    
    at.run()
    assert not at.exception

    assert at.selectbox(key="column_resolution").value == "1D (axial coordinate)"
    assert at.toggle(key="show_eq_description")

    # check if we are in the mode we support
    assert at.selectbox(key="advanced_mode").value == "Off"
    at.selectbox(key="advanced_mode").set_value("On").run()
    assert at.selectbox(key="dev_mode").value == "Off"

    # 1) test all configs with add_particles = "No" -> no has_binding and no particle_resolution
    assert at.selectbox(key="add_particles").value == "No"
    run_configs(at, {"add_particles" : "No"}, getInputToBeTested(at), 6 * 2 + 9) # Note: 3 * 2 + 3 repetitive configs due to the setup of the recursion

    # 2) test all configs with "add_particles": "Yes", "particle_resolution": "1D (radial coordinate)", "has_binding": "No"
    preSet = {"add_particles": "Yes", "particle_resolution": "1D (radial coordinate)", "has_binding": "No"}
    at.selectbox(key="add_particles").set_value("Yes").run()
    assert at.selectbox(key="particle_resolution").value == "1D (radial coordinate)"
    assert at.selectbox(key="has_binding").value == "No"
    run_configs(at, preSet, getInputToBeTested(at), 45 * 2 + 3) # Note: some repetitive configs due to the setup of the recursion

    # 3) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "No"
    preSet = {"add_particles": "Yes", "particle_resolution": "0D (homogeneous)", "has_binding": "No"}
    at.selectbox(key="particle_resolution").set_value("0D (homogeneous)").run()
    run_configs(at, preSet, getInputToBeTested(at), 21 * 2 + 3) # Note: nine repetitive configs due to the setup of the recursion

    # 4) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "Yes"
    preSet = {"add_particles": "Yes", "particle_resolution": "0D (homogeneous)", "has_binding": "Yes"}
    at.selectbox(key="has_binding").set_value("Yes").run()
    run_configs(at, preSet, getInputToBeTested(at), 90 * 2 + 9) # Note: some repetitive configs due to the setup of the recursion

    # 5) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "Yes"
    preSet = {"add_particles": "Yes", "particle_resolution": "1D (radial coordinate)", "has_binding": "Yes"}
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    run_configs(at, preSet, getInputToBeTested(at), 765) # Note: some repetitive configs due to the setup of the recursion

    # 6) test all configs for the tank
    # enable the tank and remove column resolution from the test set
    critical_variables.remove("0D (Homogeneous Tank)")
    critical_variables.append("column_resolution")
    at.selectbox(key="column_resolution").set_value("0D (Homogeneous Tank)").run()

    # reset variables
    at.selectbox(key="add_particles").set_value("No").run()

    # 6a) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "No"
    preSet = {"add_particles": "Yes", "particle_resolution": "1D (radial coordinate)", "has_binding": "No"}
    at.selectbox(key="add_particles").set_value("Yes").run()
    run_configs(at, preSet, getInputToBeTested(at), 14) # Note: some repetitive configs due to the setup of the recursion

    # 6b) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "No"
    preSet = {"add_particles": "Yes", "particle_resolution": "0D (homogeneous)", "has_binding": "No"}
    at.selectbox(key="particle_resolution").set_value("0D (homogeneous)").run()
    run_configs(at, preSet, getInputToBeTested(at), 4 + 2) # Note: some repetitive configs due to the setup of the recursion

    # 6c) test all configs with add_particles = "Yes", "particle_resolution" = "0D (homogeneous)", "has_binding" = "Yes"
    preSet = {"add_particles": "Yes", "particle_resolution": "0D (homogeneous)", "has_binding": "Yes"}
    at.selectbox(key="has_binding").set_value("Yes").run()
    run_configs(at, preSet, getInputToBeTested(at), 30) # Note: some repetitive configs due to the setup of the recursion

    # 6d) test all configs with add_particles = "Yes", "particle_resolution" = "1D (radial coordinate)", "has_binding" = "Yes"
    preSet = {"add_particles": "Yes", "particle_resolution": "1D (radial coordinate)", "has_binding": "Yes"}
    at.selectbox(key="particle_resolution").set_value("1D (radial coordinate)").run()
    run_configs(at, preSet, getInputToBeTested(at), 126) # Note: some repetitive configs due to the setup of the recursion
    