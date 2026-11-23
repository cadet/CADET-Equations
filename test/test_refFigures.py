"""
This script compares the generated model figures to reference data.

It mirrors test_refLatex.py: every configuration under
``data/EquationGenerator_configs`` is applied to the app and the SVG of each
figure it renders is compared byte for byte with
``data/ref_svg/<model>/<figure>.svg``, one file per figure so that each
reference is a valid document a viewer or an editor can open on its own.

Because the figures are display only and never enter the exported document,
this is the only test that pins their content. Regenerate the references with
``test/regenerate_ref_svg.py`` after an intended change, then check ``git
status``: the files that change are exactly the models the edit affected.
"""

import json
import os
import re

import pytest
from streamlit.testing.v1 import AppTest
from test_refLatex import apply_model_from_config

# Crystallization models have no figures yet, so they are not listed here.
FIGURE_MODELS = [
    # 0D (Tank) models
    "CSTR",
    "CSTR_filter",
    "FiniteBath_wo_pores",
    "FiniteBath_w_pores",
    "GeneralFiniteBath",
    # 1D Axial models
    "Plug_Flow",
    "Dispersive_Plug_Flow",
    "LRM",
    "LRM_reqBnd",
    "LRMP",
    "LRMP_noBnd",
    "GRM",
    "GRM_reqBnd",
    # Binding models with a figure of their own
    "GRM_dynLin",
    "GRM_Langmuir",
    "GRM_SMA",
    "GRMsd",
    "GRMsd_PSD",
    "GRMsd_parCore",
    "GRMsd_nonLimFD_parCore",
    "GRMsd_nonLimFD_reqBnd",
    "GRMsd_legacyFormat",
    # Particle geometries other than a sphere
    "GRMsd_cylinderPar",
    "GRMsd_slabPar",
    # Side-cavity models
    "GRM_sideCavities",
    "GRM_sideCavities_nonLimFD",
    # 2D and 3D models
    "2D_Plug_Flow",
    "2D_LRM",
    "GRMsd2D",
    "3D_GRM",
    # Radial flow models
    "Radial_Plug_Flow",
    "Radial_GRM",
    "Radial_LRMP",
    # Frustum models
    "Frustum_Plug_Flow",
    "Frustum_GRM",
    # Generalized (smoothly varying cross section) models
    "Generalized_Plug_Flow",
    "Generalized_GRM",
]


def figures_of(at) -> list[str]:
    """Return the SVG of every figure the app rendered, in document order."""
    return [str(block.value) for block in at.markdown if str(block.value).startswith("<svg")]


def name_of(figure: str) -> str:
    """Return the name a figure carries in its ``data-figure`` attribute."""
    return re.search(r'data-figure="([^"]+)"', figure).group(1)


def render_figures(model_name: str, test_dir: str) -> dict[str, str]:
    """Apply a stored configuration and return its figures by name."""
    at = AppTest.from_file("../Equation-Generator.py")
    # The first run of a process imports the whole app; three seconds, the
    # default, is not always enough for it on a cold cache.
    at.run(timeout=30)
    assert not at.exception

    with open(test_dir + "/data/EquationGenerator_configs/" + model_name + ".json") as file:
        apply_model_from_config(at, json.load(file))

    assert not at.exception
    return {name_of(figure): figure + "\n" for figure in figures_of(at)}


def reference_dir(model_name: str, test_dir: str) -> str:
    """Return the directory holding one file per figure of ``model_name``."""
    return os.path.join(test_dir, "data", "ref_svg", model_name)


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("model_name", FIGURE_MODELS)
def test_figures_against_reference(model_name, test_dir):
    directory = reference_dir(model_name, test_dir)
    assert os.path.isdir(directory), (
        f"No reference figures for {model_name}. Run test/regenerate_ref_svg.py to create them."
    )

    figures = render_figures(model_name, test_dir)
    stored = sorted(file[: -len(".svg")] for file in os.listdir(directory) if file.endswith(".svg"))
    assert stored == sorted(figures), (
        f"{model_name} now draws {sorted(figures)} but {stored} are on file. Run test/regenerate_ref_svg.py."
    )

    for name, figure in figures.items():
        with open(os.path.join(directory, name + ".svg"), encoding="utf-8") as file:
            assert figure == file.read(), f"{model_name}/{name}.svg"


@pytest.mark.ci
@pytest.mark.reference
def test_a_model_without_particles_has_only_the_column_figure(test_dir):
    assert len(figures_of(_app_for("Plug_Flow", test_dir))) == 1


@pytest.mark.ci
@pytest.mark.reference
def test_a_model_with_particles_also_has_a_particle_figure(test_dir):
    assert len(figures_of(_app_for("GRM", test_dir))) == 2


def _app_for(model_name: str, test_dir: str):
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    with open(test_dir + "/data/EquationGenerator_configs/" + model_name + ".json") as file:
        apply_model_from_config(at, json.load(file))
    return at
