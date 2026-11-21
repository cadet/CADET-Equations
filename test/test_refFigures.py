"""
This script compares the generated model figures to reference data.

It mirrors test_refLatex.py: every configuration under
``data/EquationGenerator_configs`` is applied to the app and the SVG of each
figure it renders is compared byte for byte with ``data/ref_svg/<name>.svg``.

Because the figures are display only and never enter the exported document,
this is the only test that pins their content. Regenerate the references with
``test/regenerate_ref_svg.py`` after an intended change, then check ``git
status``: the files that change are exactly the models the edit affected.
"""

import json
import os

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


def render_figures(model_name: str, test_dir: str) -> str:
    """Apply a stored configuration and return its figures as one document."""
    at = AppTest.from_file("../Equation-Generator.py")
    at.run()
    assert not at.exception

    with open(test_dir + "/data/EquationGenerator_configs/" + model_name + ".json") as file:
        apply_model_from_config(at, json.load(file))

    assert not at.exception
    return "\n".join(figures_of(at)) + "\n"


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("model_name", FIGURE_MODELS)
def test_figures_against_reference(model_name, test_dir):
    reference_path = test_dir + "/data/ref_svg/" + model_name + ".svg"
    assert os.path.exists(reference_path), (
        f"No reference figure for {model_name}. Run test/regenerate_ref_svg.py to create it."
    )

    with open(reference_path, encoding="utf-8") as file:
        reference = file.read()

    assert render_figures(model_name, test_dir) == reference


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
