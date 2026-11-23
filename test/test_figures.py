"""
Unit tests for src/figures.

The tests assert on the *scene* wherever possible rather than on the rendered
SVG, so that moving a drawing around does not break a test about which features
a model's figure shows. Only the tests that are explicitly about the renderer
look at the SVG itself.
"""

import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from types import SimpleNamespace

import pytest

from src.equations import BINDING_MODELS
from src.figures import binding_figures, column_figure, particle_figures
from src.figures.binding import DRAWABLE, build_binding_scene
from src.figures.column import _AXIS_Y, _BODIES, build_column_scene
from src.figures.labels import NORMAL, SUBSCRIPT, SUPERSCRIPT, Run, label_runs, label_text
from src.figures.palette import css_rules
from src.figures.particle import build_particle_scene
from src.figures.scene import Arrow, Ellipse, EllipseArc, Label, Polyline, Rect, Scene, Style
from src.figures.svg import render_svg

# %% Helpers


@dataclass(frozen=True)
class FakeParticle:
    """Hashable stand-in for src.model_particle.Particle.

    Mirrors the approach of test_equations.py: the figure builders read only a
    few attributes, so the Streamlit-dependent dataclass is not needed. It must
    be hashable because Column counts particle types with a Counter.
    """

    geometry: str = "Sphere"
    resolution: str = "1D"
    has_core: bool = False
    has_binding: bool = True
    req_binding: bool = False
    has_surfDiff: bool = False
    nonlimiting_filmDiff: bool = False
    binding_model: str = "Arbitrary"


def _make_column(**overrides):
    """Create a minimal Column-like mock object."""
    particles = overrides.pop("particles", [FakeParticle()])
    defaults = dict(
        column_type="Axial",
        resolution="1D",
        has_axial_coordinate=True,
        has_radial_coordinate=False,
        has_angular_coordinate=False,
        has_axial_dispersion=True,
        has_filter=False,
        N_p=len(particles),
        particle_models=list(particles),
        par_type_counts={particle: 1 for particle in particles},
        model_name=lambda: "General Rate Model",
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def _labels(scene):
    return [item.text for item in scene.items if isinstance(item, Label)]


def _arrows(scene, style):
    return [item for item in scene.items if isinstance(item, Arrow) and item.style is style]


def _of_style(scene, style, kind=Ellipse):
    return [item for item in scene.items if isinstance(item, kind) and item.style is style]


def _flow_paths(scene):
    """Return the flow paths of a column figure, in drawing order."""
    return [item for item in scene.items if isinstance(item, Polyline) and item.style is Style.FLOW and item.arrow]


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
        (r"P^{L^{\b}}", "CADET", "PLb"),
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
def test_nested_scripts_ride_one_level_further_out():
    """P^{L^{b}} must not flatten both scripts onto one height."""
    runs = label_runs(r"P^{L^{\b}}")
    assert [run.shift for run in runs] == [NORMAL, SUPERSCRIPT, 2 * SUPERSCRIPT]


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


# %% Column bodies


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize(
    "column_type, resolution",
    [
        ("Axial", "1D"),
        ("Axial", "2D"),
        ("Axial", "3D"),
        ("Axial", "0D"),  # the mixed tank is an axial unit of zero dimensions
        ("Radial", "1D"),
        ("Frustum", "1D"),
        ("Generalized", "1D"),
    ],
)
def test_every_geometry_builds_a_scene(column_type, resolution):
    scene = build_column_scene(
        _make_column(
            column_type=column_type,
            resolution=resolution,
            has_axial_coordinate=resolution != "0D",
            has_radial_coordinate=resolution in ("2D", "3D"),
            has_angular_coordinate=resolution == "3D",
        )
    )
    assert scene.items
    assert scene.height > 0


@pytest.mark.ci
@pytest.mark.unit_test
def test_column_scene_draws_a_cylinder_with_a_hidden_far_cap():
    scene = build_column_scene(_make_column())
    assert len(_of_style(scene, Style.HIDDEN, EllipseArc)) == 1
    assert _of_style(scene, Style.WALL)


@pytest.mark.ci
@pytest.mark.unit_test
def test_column_scene_labels_the_boundaries():
    labels = _labels(build_column_scene(_make_column()))
    assert r"\Gamma^{\mathrm{in}}" in labels
    assert r"\Gamma^{\mathrm{out}}" in labels
    assert r"\Gamma^{\mathrm{wall}}" in labels
    assert r"\Omega^{\b}" in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_tank_has_no_wall_boundary_and_no_coordinate():
    """A homogeneous tank resolves nothing, so it carries no axis or wall."""
    scene = build_column_scene(_make_column(resolution="0D", has_axial_coordinate=False))
    assert r"\Gamma^{\mathrm{wall}}" not in _labels(scene)
    assert _arrows(scene, Style.AXIS) == []


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_tank_draws_a_filter_only_when_one_is_configured():
    without = build_column_scene(_make_column(resolution="0D", has_axial_coordinate=False))
    with_filter = build_column_scene(_make_column(resolution="0D", has_axial_coordinate=False, has_filter=True))
    assert "filter" not in _labels(without)
    assert "filter" in _labels(with_filter)


# %% Resolution: the coordinate arrows


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("resolution", ["1D", "2D", "3D"])
def test_the_resolution_decides_which_coordinates_are_drawn(resolution):
    """The resolution is shown by which coordinate axes the figure draws."""
    labels = _labels(
        build_column_scene(
            _make_column(
                resolution=resolution,
                has_radial_coordinate=resolution in ("2D", "3D"),
                has_angular_coordinate=resolution == "3D",
            )
        )
    )
    assert "z" in labels
    assert (r"\rho" in labels) is (resolution in ("2D", "3D"))
    assert (r"\varphi" in labels) is (resolution == "3D")


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize(
    "column_type, symbol",
    [("Axial", "z"), ("Frustum", "z"), ("Generalized", "x"), ("Radial", r"\rho")],
)
def test_the_flow_coordinate_matches_the_geometry(column_type, symbol):
    assert symbol in _labels(build_column_scene(_make_column(column_type=column_type)))


# %% Convection and dispersion


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize(
    "column_type, symbol",
    [
        ("Axial", r"D^{\mathrm{ax}}"),
        ("Frustum", r"D^{\mathrm{ax}}"),
        ("Radial", r"D^{\mathrm{rad}}"),
        ("Generalized", r"D^{\mathrm{x}}"),
    ],
)
def test_dispersion_is_named_after_the_geometry(column_type, symbol):
    scene = build_column_scene(_make_column(column_type=column_type, has_axial_dispersion=True))
    assert symbol in _labels(scene)


@pytest.mark.ci
@pytest.mark.unit_test
def test_pure_convection_draws_a_straight_flow_path():
    """Without dispersion the flow runs straight along the axis of the body."""
    scene = build_column_scene(_make_column(has_axial_dispersion=False))
    assert r"D^{\mathrm{ax}}" not in _labels(scene)

    path = _flow_paths(scene)[0]
    assert len(path.points) == 2
    assert {point[1] for point in path.points} == {_AXIS_Y}


@pytest.mark.ci
@pytest.mark.unit_test
def test_dispersion_bends_the_flow_path_around_the_packing():
    """Dispersion is drawn as a meander, not as a plot beside the column."""
    scene = build_column_scene(_make_column(has_axial_dispersion=True))

    path = _flow_paths(scene)[0]
    offsets = [point[1] - _AXIS_Y for point in path.points]
    assert max(offsets) > 0.2 and min(offsets) < -0.2, "the path should leave the axis both ways"
    # It starts and ends on the axis, so it joins the straight inlet and outlet.
    assert offsets[0] == offsets[-1] == 0.0


@pytest.mark.ci
@pytest.mark.unit_test
def test_the_flow_path_never_runs_through_a_bead():
    """The meander shares the body with the packing, so it needs its own lane."""
    for column_type in ("Axial", "Frustum", "Generalized", "Radial"):
        scene = build_column_scene(_make_column(column_type=column_type))
        beads = _of_style(scene, Style.PARTICLE)
        for path in _flow_paths(scene):
            for x, y in path.points:
                for bead in beads:
                    distance = math.hypot(x - bead.center[0], y - bead.center[1])
                    assert distance > bead.rx, (column_type, bead.center)


# %% Particle types in the bed


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_column_without_particles_draws_no_beads():
    scene = build_column_scene(_make_column(particles=[], N_p=0))
    assert not _of_style(scene, Style.PARTICLE)


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("nonlimiting, films", [(False, 8), (True, 0)])
def test_a_film_ring_marks_film_diffusion_limited_particles(nonlimiting, films):
    """Every bead of a film-limited type carries the ring, and no other does."""
    scene = build_column_scene(_make_column(particles=[FakeParticle(nonlimiting_filmDiff=nonlimiting)]))
    assert len(_of_style(scene, Style.FILM)) == films


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("nonlimiting, arrows", [(False, 8), (True, 0)])
def test_film_diffusion_is_drawn_as_a_flux_into_every_limited_bead(nonlimiting, arrows):
    """It is shown on the beads themselves, not explained in a key."""
    scene = build_column_scene(_make_column(particles=[FakeParticle(nonlimiting_filmDiff=nonlimiting)]))
    # The flow through the column is a path, so every flow arrow is a flux.
    assert len(_arrows(scene, Style.FLOW)) == arrows


@pytest.mark.ci
@pytest.mark.unit_test
def test_one_bead_of_each_type_is_named_after_its_film_regime():
    scene = build_column_scene(_make_column(particles=[FakeParticle(), FakeParticle(nonlimiting_filmDiff=True)]))
    labels = _labels(scene)
    assert r"j = 1: k^{\mathrm{f}}" in labels
    assert r"j = 2: c^{\b} = c^{\p}" in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_single_particle_type_needs_no_index_on_its_bead():
    assert r"k^{\mathrm{f}}" in _labels(build_column_scene(_make_column()))


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("column_type", ["Axial", "Frustum", "Generalized", "Radial"])
def test_beads_stay_inside_the_bed_of_every_geometry(column_type):
    """A tapered or annular bed is cramped, so bead centres must be clamped."""
    column = _make_column(column_type=column_type)
    probe = Scene(name="probe", xmin=0, ymin=0, xmax=1, ymax=1)
    body = _BODIES[column_type](probe, column)

    scene = build_column_scene(column)
    for bead in _of_style(scene, Style.PARTICLE):
        x, y = bead.center
        assert any(start <= x <= stop for start, stop in body.bed_spans), bead
        depth = abs(y - body.bed_axis)
        assert depth + bead.rx <= body.bed_profile(x) + 1e-9, bead
        assert depth - bead.rx >= body.corridor - 1e-9, bead


@pytest.mark.ci
@pytest.mark.unit_test
def test_the_frustum_narrows_along_the_flow():
    """Its inlet is the wide end, as in the TikZ source the figure follows."""
    probe = Scene(name="probe", xmin=0, ymin=0, xmax=1, ymax=1)
    body = _BODIES["Frustum"](probe, _make_column(column_type="Frustum"))
    start, stop = body.bed_spans[0]
    assert body.bed_profile(start) > body.bed_profile(stop)


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_radial_column_stands_upright_and_feeds_over_its_shell():
    """The inflow boundary of the radial model sits at the outer radius."""
    probe = Scene(name="probe", xmin=0, ymin=0, xmax=1, ymax=1)
    body = _BODIES["Radial"](probe, _make_column(column_type="Radial"))
    assert body.top - body.bottom > body.outlet_x - body.inlet_x, "the body is taller than it is wide"

    scene = build_column_scene(_make_column(column_type="Radial"))
    # Two mirrored paths run inwards, as int_vol_BC feeds the shell.
    centre_x = 0.5 * (body.inlet_x + body.outlet_x)
    paths = _flow_paths(scene)
    assert len(paths) == 2
    for path in paths:
        assert abs(path.points[0][0] - centre_x) > abs(path.points[-1][0] - centre_x)


@pytest.mark.ci
@pytest.mark.unit_test
def test_the_radial_channel_is_closed_at_the_base():
    """Both ends of the channel are drawn, the hidden one dashed."""
    scene = build_column_scene(_make_column(column_type="Radial"))
    hidden = _of_style(scene, Style.HIDDEN)
    assert len(hidden) == 1, "the base of the channel is the one hidden ellipse"
    # It is drawn as a whole ellipse: an arc of a full turn would be empty.
    assert "<ellipse" in render_svg(scene)


# %% Particle figure


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("geometry", ["Sphere", "Cylinder", "Slab"])
def test_every_particle_geometry_builds_a_scene(geometry):
    assert build_particle_scene(FakeParticle(geometry=geometry)).items


@pytest.mark.ci
@pytest.mark.unit_test
def test_the_particle_figure_shows_the_phases_of_the_model():
    labels = _labels(build_particle_scene(FakeParticle()))
    assert r"c^{\b}" in labels
    assert r"c^{\p}" in labels
    assert r"c^{\s}" in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_particle_without_binding_has_no_bound_phase():
    labels = _labels(build_particle_scene(FakeParticle(has_binding=False)))
    assert r"c^{\s}" not in labels
    assert r"f^{\mathrm{bind}}" not in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_rapid_equilibrium_binding_is_named_as_such():
    labels = _labels(build_particle_scene(FakeParticle(req_binding=True)))
    assert "rapid equilibrium" in labels
    assert r"f^{\mathrm{bind}} = 0" in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_non_limiting_film_diffusion_draws_an_equality_instead_of_a_flux():
    """Without a film resistance there is no film, and no transport step."""
    scene = build_particle_scene(FakeParticle(nonlimiting_filmDiff=True))
    labels = _labels(scene)
    assert "no film resistance" in labels
    assert "film diffusion" not in labels
    assert r"k^{\mathrm{f}}" not in labels
    assert not [item for item in scene.items if getattr(item, "style", None) is Style.FILM]


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_surfDiff", [True, False])
def test_surface_diffusion_is_drawn_only_when_configured(has_surfDiff):
    labels = _labels(build_particle_scene(FakeParticle(has_surfDiff=has_surfDiff)))
    assert (r"D^{\s}" in labels) is has_surfDiff


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_homogeneous_particle_resolves_no_interior_coordinate():
    labels = _labels(build_particle_scene(FakeParticle(resolution="0D")))
    assert "r" not in labels
    assert r"D^{\p}" not in labels
    assert "homogeneous particle" in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_particle_core_adds_an_inner_boundary_and_a_tick():
    scene = build_particle_scene(FakeParticle(has_core=True))
    assert [item for item in scene.items if getattr(item, "style", None) is Style.CORE]
    assert r"R^{\mathrm{c}}" in _labels(scene)


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("geometry", ["Sphere", "Cylinder", "Slab"])
def test_the_particle_radius_is_named_as_the_equations_name_it(geometry):
    """The symbol table writes r in (0, R^p), so the figure must not say P."""
    labels = _labels(build_particle_scene(FakeParticle(geometry=geometry)))
    assert r"R^{\mathrm{p}}" in labels
    assert not [label for label in labels if label.startswith("P^")]


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("geometry", ["Sphere", "Cylinder", "Slab"])
def test_the_particle_is_drawn_as_a_porous_matrix(geometry):
    """Pore diffusion and surface diffusion need a pore to happen in."""
    scene = build_particle_scene(FakeParticle(geometry=geometry))
    assert any(getattr(item, "style", None) is Style.MATRIX for item in scene.items), "the solid skeleton"
    assert len(_of_style(scene, Style.PORE)) >= 6, "the pores through it"


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("has_binding, bound", [(True, True), (False, False)])
def test_the_bound_phase_sits_on_the_pore_wall(has_binding, bound):
    scene = build_particle_scene(FakeParticle(has_binding=has_binding, has_surfDiff=has_binding))
    assert bool(_of_style(scene, Style.SOLID)) is bound


@pytest.mark.ci
@pytest.mark.unit_test
def test_one_particle_figure_is_produced_per_distinct_type():
    column = _make_column(particles=[FakeParticle(), FakeParticle(has_core=True)])
    figures = particle_figures(column)
    assert len(figures) == 2
    # Distinct scene names keep the arrow-head ids of the two figures apart.
    assert 'id="particle-1-head-flow"' in figures[0]
    assert 'id="particle-2-head-flow"' in figures[1]


# %% Binding models


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("binding_model", ["Linear", "Langmuir", "SMA"])
def test_every_drawable_binding_model_builds_a_scene(binding_model):
    assert build_binding_scene(binding_model).items


@pytest.mark.ci
@pytest.mark.unit_test
def test_an_arbitrary_isotherm_has_no_figure():
    """There is nothing characteristic to draw for an unspecified isotherm."""
    assert build_binding_scene("Arbitrary") is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_the_drawable_models_are_a_subset_of_the_models_on_offer():
    assert set(DRAWABLE) <= set(BINDING_MODELS)
    assert all(build_binding_scene(model) is not None for model in DRAWABLE)


@pytest.mark.ci
@pytest.mark.unit_test
def test_the_isotherms_are_plotted_against_the_two_concentrations():
    for binding_model in ("Linear", "Langmuir"):
        labels = _labels(build_binding_scene(binding_model))
        assert r"c^{\p}" in labels, binding_model
        assert r"c^{\s}" in labels, binding_model


@pytest.mark.ci
@pytest.mark.unit_test
def test_only_langmuir_shows_a_capacity():
    assert r"q^{\mathrm{max}}" in _labels(build_binding_scene("Langmuir"))
    assert r"q^{\mathrm{max}}" not in _labels(build_binding_scene("Linear"))


@pytest.mark.ci
@pytest.mark.unit_test
def test_sma_names_the_quantities_of_its_binding_term():
    """The cartoon must use the same symbols as equations.binding_term_sma."""
    labels = _labels(build_binding_scene("SMA"))
    for symbol in (r"\Lambda", r"\nu_{i}", r"\sigma_{i}", r"c^{\s}_{i}", r"c^{\p}_{i}"):
        assert symbol in labels, symbol
    # The displaced counter-ions go into the pore liquid, not onto the resin.
    assert r"c^{\p}_{0}" in labels


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_binding_figure_is_produced_only_for_a_particle_that_binds():
    binding = binding_figures(_make_column(particles=[FakeParticle(has_binding=False)]))
    assert binding == [None]


@pytest.mark.ci
@pytest.mark.unit_test
def test_binding_figures_run_parallel_to_the_particle_figures():
    column = _make_column(
        particles=[
            FakeParticle(binding_model="Langmuir"),
            FakeParticle(binding_model="Arbitrary", has_core=True),
        ]
    )
    figures = binding_figures(column)
    assert len(figures) == len(particle_figures(column)) == 2
    assert figures[0] is not None and 'id="binding-1-' in figures[0]
    assert figures[1] is None  # Arbitrary has no figure


# %% SVG rendering


@pytest.mark.ci
@pytest.mark.unit_test
def test_rendered_figure_is_well_formed_xml():
    ET.fromstring(column_figure(_make_column()))
    for figure in particle_figures(_make_column()):
        ET.fromstring(figure)


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
def test_a_figure_paints_itself_so_it_keeps_its_colours_outside_the_app():
    """Inkscape and most other drawing programs implement neither a CSS custom
    property nor its fallback, and a paint they cannot read leaves ``fill`` at
    its initial value: black."""
    for svg in _every_figure():
        assert "var(" not in svg
        for paint in re.findall(r'(?:fill|stroke)="([^"]+)"', svg):
            assert paint == "none" or re.fullmatch(r"#[0-9a-f]{6}", paint), paint


@pytest.mark.ci
@pytest.mark.unit_test
def test_every_painted_element_names_its_style_so_the_page_can_repaint_it():
    """The class is what lets the app override the light paint for dark mode."""
    rules = css_rules()
    for svg in _every_figure():
        painted = re.findall(r'(fill|stroke)="#[0-9a-f]{6}"', svg)
        named = [name for value in re.findall(r'class="([^"]*)"', svg) for name in value.split()]
        assert len([name for name in named if name != "cadet-figure"]) == len(painted)
        for name in named:
            assert name == "cadet-figure" or f".{name} {{" in rules


def _every_figure() -> list[str]:
    """Return one figure of every kind the app can draw."""
    column = _make_column()
    return [column_figure(column), *particle_figures(column), *(f for f in binding_figures(column) if f)]


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
def test_an_arc_can_carry_an_arrow_head():
    scene = Scene(name="test", xmin=0, ymin=0, xmax=2, ymax=2)
    scene.add(EllipseArc((1, 1), 0.5, 0.5, 0, 180, Style.AXIS, arrow=True))
    assert "marker-end=" in render_svg(scene)


@pytest.mark.ci
@pytest.mark.unit_test
def test_an_arc_of_a_full_turn_is_drawn_rather_than_dropped():
    """Its start and end coincide, where a single SVG arc draws nothing."""
    scene = Scene(name="test", xmin=0, ymin=0, xmax=2, ymax=2)
    scene.add(EllipseArc((1, 1), 0.5, 0.25, 0, 360, Style.HIDDEN))
    assert render_svg(scene).count("<path") == 2


@pytest.mark.ci
@pytest.mark.unit_test
def test_scene_y_axis_points_up():
    """Scene coordinates are maths-like; the renderer flips them for SVG."""
    scene = Scene(name="test", xmin=0, ymin=0, xmax=4, ymax=4)
    scene.add(Ellipse((1, 3), 1, 1, Style.BULK))  # near the top of the scene
    assert 'cy="1"' in render_svg(scene)  # which is a small y in SVG


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_rectangle_is_placed_by_its_lower_left_corner():
    scene = Scene(name="test", xmin=0, ymin=0, xmax=4, ymax=4)
    scene.add(Rect((1, 1), 2, 1, Style.BULK))
    # Lower-left (1, 1) with height 1 means an SVG top edge at y = 4 - 2 = 2.
    assert 'x="1" y="2" width="2" height="1"' in render_svg(scene)
