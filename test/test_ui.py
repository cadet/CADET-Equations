"""
This script implements tests iterating through all configurable states of the App and checking if an error is thrown.

The LaTeX output of every visited configuration is collected and compiled with pdflatex in a single combined document
per configuration stage, which verifies that the generated LaTeX is valid without paying the cost of one pdflatex
run per configuration.
"""

import subprocess
import tempfile
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

untested_variables = [
    "dev_mode",
    "advanced_mode",
    "var_format",
    "sym_table",
    "show_eq_description",
    "PSD",
    "PTD",
    "has_filter",
    "binding_model",
    "reaction_model",
    "has_reaction_bulk",
    "has_reaction_particle_liquid",
    "has_reaction_particle_solid",
    "N_c_choice",
    "model_type",
]
# dev_mode is not tested, TODO: test the advanced_mode (which includes PSD, PTD). test filter

# Some boxes are conditional, e.g. film_diffusion can only be configured when particles are present.
# We thus first identify the boxes that change the number of boxes and thus combinations (critical_variables).
# We then start the recursion for every combination of the aforementioned critical keys.
# Sometimes, only specific options can be critical (e.g. 0D tank). Such options can be added here and handled independently later
critical_variables = untested_variables + ["PSD", "has_binding", "particle_resolution", "Mixed tank", "column_type"]


# We test every configuration by recursively iterating through all combinations
def config_recursion(at, widgies, counter, criticals, latex_documents):

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

        options = (
            [opti for opti in widgy.options if opti not in criticals] if widgyKind == "selectbox" else [True, False]
        )

        for opt in options:  # Note: some values will be run repeatedly since the current value of that widget is applied again. Excluding this value here would make the recursion incomplete though
            widgy.set_value(opt).run()

            counter += 1
            assert not at.exception

            latex_documents.add(at.session_state.latex_string)

            if len(widgies) > 0:
                counter = config_recursion(at, widgies, counter, criticals, latex_documents)

        widgies.append(widgy.key)

    return counter


def run_configs(at, toBeFixed: dict, criticals: list, numConfigsToBeChecked: int, latex_documents: set):

    # Check if we have the correct initial configuration
    for fixKey in toBeFixed:
        if any(key == fixKey for key in [box.key for box in at.selectbox]):
            widgy = at.selectbox(key=fixKey)
        elif any(key == fixKey for key in [box.key for box in at.toggle]):
            widgy = at.toggle(key=fixKey)
        else:
            raise ValueError(f"Error: {fixKey} is neither a toggle nor a dropdown")

        if widgy.value != toBeFixed[fixKey]:
            widgy.set_value(toBeFixed[fixKey]).run()

    numConfigs = config_recursion(at, getInputToBeTested(at, criticals), 0, criticals, latex_documents)
    assert numConfigs == numConfigsToBeChecked  # Note: repetitive configs due to the setup of the recursion

    # Exercise the file generator buttons once per stage. Per-configuration LaTeX validity
    # is verified much faster by assert_latex_compiles on the collected documents.
    at.button(key="generate_pdf").click().run()
    assert not at.exception
    assert len(at.error) == 0, "PDF generation failed:\n" + "\n".join(str(err.value) for err in at.error)

    at.button(key="generate_config").click().run()
    assert not at.exception


def getInputToBeTested(at, criticals):

    return [box.key for box in at.selectbox if box.key not in criticals] + [
        box.key for box in at.toggle if box.key not in criticals
    ]


def assert_latex_compiles(latex_documents: set):
    """Compile all collected LaTeX documents in a single pdflatex run."""

    docs = sorted(latex_documents)
    preamble = docs[0].split(r"\begin{document}", 1)[0]

    lines = preamble.splitlines() + [r"\begin{document}"]
    doc_start_lines = []  # line number (1-based) where each collected document starts
    for doc in docs:
        doc_start_lines.append(len(lines) + 1)
        body = doc.split(r"\begin{document}", 1)[1].rsplit(r"\end{document}", 1)[0]
        lines.extend(body.splitlines())
        lines.append(r"\clearpage")
    lines.append(r"\end{document}")
    combined = "\n".join(lines) + "\n"

    with tempfile.TemporaryDirectory() as temp_dir:
        tex_path = Path(temp_dir) / "combined.tex"
        tex_path.write_text(combined, encoding="utf-8")

        result = subprocess.run(
            [
                "pdflatex",
                "-interaction=nonstopmode",
                "-halt-on-error",
                "-file-line-error",
                "-output-directory",
                temp_dir,
                str(tex_path),
            ],
            capture_output=True,
            text=True,
            cwd=temp_dir,
            timeout=900,
            check=False,
        )

    if result.returncode != 0:
        # Map the error line reported by -file-line-error back to the failing document
        message = f"pdflatex failed on the combined document of {len(docs)} configurations."
        for log_line in result.stdout.splitlines():
            if ".tex:" in log_line:
                try:
                    err_line = int(log_line.split(".tex:", 1)[1].split(":", 1)[0])
                except ValueError:
                    continue
                doc_idx = max(i for i, start in enumerate(doc_start_lines) if start <= err_line)
                snippet_start = doc_start_lines[doc_idx] - 1
                snippet_end = doc_start_lines[doc_idx + 1] - 1 if doc_idx + 1 < len(docs) else len(lines)
                message += f"\nFirst error: {log_line.strip()}\nFailing document (index {doc_idx}):\n" + "\n".join(
                    lines[snippet_start:snippet_end]
                )
                break
        message += "\npdflatex log tail:\n" + result.stdout[-2000:]
        raise AssertionError(message)


@pytest.mark.ci
@pytest.mark.all_config_runs
def test_streamlit_app_defaults():

    at = AppTest.from_file("../Equation-Generator.py")

    at.run()
    assert not at.exception

    assert at.selectbox(key="column_resolution").value == "1D (axial coordinate)"
    assert at.toggle(key="show_eq_description")

    # check if we are in the mode we support
    assert at.selectbox(key="advanced_mode").value == "Off"
    at.selectbox(key="advanced_mode").set_value("On").run()
    assert not at.button(key="dev_mode_button").value
    assert at.selectbox(key="PSD").value == "No"


# Every stage fixes one combination of the critical variables and recursively tests all
# combinations of the remaining widgets. The stages are independent tests so that they can run in parallel.
config_stages = [
    pytest.param("Axial", {"PSD": "No"}, 62, id="axial_noPSD"),
    pytest.param(
        "Axial",
        {"PSD": "Yes", "has_binding": "No", "particle_resolution": "1D (radial coordinate)"},
        269,
        id="axial_par1D_noBinding",
    ),
    pytest.param(
        "Axial",
        {"PSD": "Yes", "has_binding": "No", "particle_resolution": "0D (homogeneous)"},
        133,
        id="axial_par0D_noBinding",
    ),
    pytest.param(
        "Axial",
        {"PSD": "Yes", "has_binding": "Yes", "particle_resolution": "0D (homogeneous)"},
        549,
        id="axial_par0D_binding",
    ),
    pytest.param(
        "Axial",
        {"PSD": "Yes", "has_binding": "Yes", "particle_resolution": "1D (radial coordinate)"},
        2221,
        id="axial_par1D_binding",
    ),
    pytest.param(
        "Mixed tank",
        {"PSD": "Yes", "has_binding": "No", "particle_resolution": "1D (radial coordinate)"},
        45,
        id="tank_par1D_noBinding",
    ),
    pytest.param(
        "Mixed tank",
        {"PSD": "Yes", "has_binding": "No", "particle_resolution": "0D (homogeneous)"},
        21,
        id="tank_par0D_noBinding",
    ),
    pytest.param(
        "Mixed tank",
        {"PSD": "Yes", "has_binding": "Yes", "particle_resolution": "0D (homogeneous)"},
        93,
        id="tank_par0D_binding",
    ),
    pytest.param(
        "Mixed tank",
        {"PSD": "Yes", "has_binding": "Yes", "particle_resolution": "1D (radial coordinate)"},
        381,
        id="tank_par1D_binding",
    ),
]


@pytest.mark.ci
@pytest.mark.all_config_runs
@pytest.mark.parametrize("column_type,preSet,numConfigs", config_stages)
def test_streamlit_app_configs(column_type, preSet, numConfigs):

    at = AppTest.from_file("../Equation-Generator.py")

    at.run()
    assert not at.exception

    at.selectbox(key="advanced_mode").set_value("On").run()

    criticals = list(critical_variables)
    if column_type == "Mixed tank":
        # the tank has no spatial resolution to be configured
        at.selectbox(key="column_type").set_value("Mixed tank").run()
        criticals.remove("Mixed tank")
        criticals.append("column_resolution")

    latex_documents = set()
    run_configs(at, preSet, criticals, numConfigs, latex_documents)

    assert_latex_compiles(latex_documents)
