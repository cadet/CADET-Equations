"""
Unit tests for src/figures.

The tests assert on the *scene* wherever possible rather than on the rendered
SVG text, so that moving a drawing around does not break a test about which
features a model's figure shows.  Only the tests that are explicitly about the
renderer look at the SVG itself.
"""

import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pytest

from src.figures import column_figure
from src.figures.column import build_column_scene
from src.figures.labels import NORMAL, SUBSCRIPT, SUPERSCRIPT, Run, label_runs, label_text
from src.figures.scene import Arrow, Ellipse, EllipseArc, Label, Scene, Style
from src.figures.svg import render_svg

# %% Helpers


def _make_column(**overrides):
    """Create a minimal Column-like mock object.

    Mirrors the approach of test_equations.py: the figure builders take only a
    few attributes, so the Streamlit-dependent dataclass is not needed.
    """
    defaults = dict(has_axial_coordinate=True, model_name=lambda: "General Rate Model")
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def _labels(scene):
    """Return the raw label texts of a scene."""
    return [item.text for item in scene.items if isinstance(item, Label)]


# %% Label parsing


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize(
    "text, var_format, expected",
    [
        (r"c^{\b}_i", "CADET", "cbi"),
        (r"c^{\s}_i", "CADET", "csi"),
        # Legacy renames the bound phase, exactly as it does in the equations.
        (r"c^{\s}_i", "Legacy", "qi"),
        (r"c^{\b}", "Legacy", "c"),
        (r"\Gamma^{\mathrm{wall}}", "CADET", "Γwall"),
        (r"\rho", "CADET", "ρ"),
        (r"L^{\b}", "CADET", "Lb"),
    ],
)
def test_label_text_expands_shorthands(text, var_format, expected):
    assert label_text(text, var_format) == expected


@pytest.mark.ci
@pytest.mark.unit_test
def test_label_runs_carry_shift_and_style():
    runs = label_runs(r"D^{\mathrm{ax}}_i")
    assert runs == [
        Run("D", NORMAL, True),
        Run("ax", SUPERSCRIPT, False),  # \mathrm is upright
        Run("i", SUBSCRIPT, True),
    ]


@pytest.mark.ci
@pytest.mark.unit_test
def test_capital_greek_is_upright_and_lowercase_is_italic():
    (capital,) = label_runs(r"\Omega")
    (lowercase,) = label_runs(r"\rho")
    assert not capital.italic
    assert lowercase.italic


@pytest.mark.ci
@pytest.mark.unit_test
def test_digits_are_upright():
    assert label_runs("0") == [Run("0", NORMAL, False)]


# %% Column scene


@pytest.mark.ci
@pytest.mark.unit_test
def test_column_scene_draws_a_cylinder_with_a_hidden_far_cap():
    scene = build_column_scene(_make_column())

    # One cap is drawn as two arcs so that its far half can be dashed.
    hidden = [item for item in scene.items if isinstance(item, EllipseArc) and item.style is Style.HIDDEN]
    assert len(hidden) == 1
    # The near cap is a full ellipse.
    assert any(isinstance(item, Ellipse) and item.style is Style.WALL for item in scene.items)


@pytest.mark.ci
@pytest.mark.unit_test
def test_column_scene_labels_the_boundaries():
    scene = build_column_scene(_make_column())
    labels = _labels(scene)
    assert r"\Gamma^{\mathrm{in}}" in labels
    assert r"\Gamma^{\mathrm{out}}" in labels
    assert r"\Gamma^{\mathrm{wall}}" in labels
    assert r"\Omega^{\b}" in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_column_scene_has_two_flow_arrows():
    scene = build_column_scene(_make_column())
    flow = [item for item in scene.items if isinstance(item, Arrow) and item.style is Style.FLOW]
    assert len(flow) == 2


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_axial_coordinate, expected", [(True, 1), (False, 0)])
def test_axial_coordinate_is_drawn_only_when_resolved(has_axial_coordinate, expected):
    """The resolution is shown by which coordinate axes the figure draws."""
    scene = build_column_scene(_make_column(has_axial_coordinate=has_axial_coordinate))
    axes = [item for item in scene.items if isinstance(item, Arrow) and item.style is Style.AXIS]
    assert len(axes) == expected
    assert ("z" in _labels(scene)) is has_axial_coordinate


# %% SVG rendering


@pytest.mark.ci
@pytest.mark.unit_test
def test_rendered_figure_is_well_formed_xml():
    ET.fromstring(column_figure(_make_column()))


@pytest.mark.ci
@pytest.mark.unit_test
def test_rendered_figure_follows_the_parameter_format():
    """Legacy format renames the bound phase in figures as well."""
    scene = Scene(name="test", xmin=0, ymin=0, xmax=1, ymax=1)
    scene.add(Label((0.5, 0.5), r"c^{\s}_i"))
    assert ">q<" in render_svg(scene, "Legacy")
    assert ">c<" in render_svg(scene, "CADET")


@pytest.mark.ci
@pytest.mark.unit_test
def test_colours_are_css_variables_so_the_figure_follows_the_theme():
    svg = column_figure(_make_column())
    assert "var(--cadet-fig-" in svg
    # A literal colour would not follow the light/dark appearance. The only
    # "#" a figure may contain is the one of a marker reference.
    assert re.search(r"#[0-9a-fA-F]{3,8}\b", svg) is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_coordinates_are_rounded_to_two_decimals():
    """Keeps cosmetic changes to the geometry out of the reference figures."""
    scene = Scene(name="test", xmin=0, ymin=0, xmax=10, ymax=10)
    scene.add(Ellipse((1 / 3, 2 / 3), 1 / 7, 1 / 7, Style.BULK))
    assert 'cx="0.33" cy="9.33"' in render_svg(scene)


@pytest.mark.ci
@pytest.mark.unit_test
def test_marker_ids_are_derived_from_the_scene_name():
    """Unique ids let several figures coexist on one page."""
    scene = Scene(name="widget", xmin=0, ymin=0, xmax=1, ymax=1)
    scene.add(Arrow((0, 0), (1, 1), Style.FLOW))
    svg = render_svg(scene)
    assert 'id="widget-head-flow"' in svg
    assert "url(#widget-head-flow)" in svg


@pytest.mark.ci
@pytest.mark.unit_test
def test_double_arrows_get_a_head_at_both_ends():
    scene = Scene(name="test", xmin=0, ymin=0, xmax=1, ymax=1)
    scene.add(Arrow((0, 0), (1, 1), Style.AXIS, double=True))
    svg = render_svg(scene)
    assert "marker-start=" in svg
    assert "marker-end=" in svg


@pytest.mark.ci
@pytest.mark.unit_test
def test_scene_y_axis_points_up():
    """Scene coordinates are maths-like; the renderer flips them for SVG."""
    scene = Scene(name="test", xmin=0, ymin=0, xmax=4, ymax=4)
    scene.add(Ellipse((1, 3), 1, 1, Style.BULK))  # near the top of the scene
    assert 'cy="1"' in render_svg(scene)  # which is a small y in SVG
