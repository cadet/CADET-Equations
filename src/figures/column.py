"""Build the scene of the interstitial volume (column) figure.

A column figure is assembled from four independent layers, each driven by a
different part of the model:

=====================  ====================================================
layer                  driven by
=====================  ====================================================
body                   ``column_type`` and, for the tank, ``resolution``
flow path              ``has_axial_dispersion``
coordinate arrows      which coordinates the resolution resolves
particle beads         ``par_type_counts`` and each type's film diffusion
=====================  ====================================================

Each body returns a :class:`_Body` describing where the other layers may
attach, so adding a geometry means writing one body function and registering
it in ``_BODIES`` -- the layers above it need no change.

The drawings are adapted from the TikZ sources of the CADET modelling papers,
so a figure in the app and a figure in a paper show the same object.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field

from src.figures.annotate import TICK, boundary, dimension
from src.figures.scene import (
    Anchor,
    Arrow,
    Ellipse,
    EllipseArc,
    Label,
    Point,
    Polyline,
    Rect,
    Scene,
    Style,
)

# %% Shared layout, in scene units

_LEFT, _RIGHT = 2.2, 10.2  # the straight part of a horizontal body
_AXIS_Y = 4.5  # height of the body's axis of symmetry
_BODY_RY = 1.95  # half-height of a cylindrical body
_PERSPECTIVE = 0.355  # how far an end cap bulges out, relative to its height
_SAMPLES = 48  # points used to approximate a smooth wall or a flow path

_WAVE = 0.46  # amplitude of a dispersed (meandering) flow path
_WAVELENGTH = 2.3  # distance it takes to complete one meander
_LEAD = 1.5  # how far a flow path reaches beyond the body
_DIAGONAL = 0.707  # cosine of 45 degrees, where a cap is named


@dataclass
class _Body:
    """Where the annotation layers may attach to a drawn body.

    ``inlet_x`` and ``outlet_x`` are where the flow path meets the silhouette
    and ``wave`` is the stretch of it that runs through the packing, so that a
    dispersed path meanders there and nowhere else.  ``bed_spans`` are the
    x-ranges the packing occupies and ``bed_profile`` gives its half-height
    about ``bed_axis``, which together keep the beads of the scatter inside a
    body of any shape.
    ``flow_symbol`` and ``flow_ticks`` name the coordinate that runs along the
    flow, and ``dispersion`` is the symbol of the flow-directional dispersion
    coefficient.
    """

    top: float
    bottom: float
    bed_spans: tuple[tuple[float, float], ...] = ()
    bed_profile: Callable[[float], float] = field(default=lambda x: 0.0)
    bed_axis: float = _AXIS_Y
    inlet_x: float = _LEFT
    outlet_x: float = _RIGHT
    wave: tuple[float, float] = (_LEFT, _RIGHT)
    flow_symbol: str = "z"
    inlet_ry: float = _BODY_RY  # half-height of the cap the flow enters over
    outlet_ry: float = _BODY_RY  # and of the one it leaves over
    flow_ticks: tuple[str, str] = ("0", r"L^{\b}")
    dispersion: str = r"D^{\mathrm{ax}}"
    draws_own_flow: bool = False
    has_wall: bool = True
    wall_at: Point = (_LEFT + 0.78 * (_RIGHT - _LEFT), _AXIS_Y + _BODY_RY)
    wall_side: float = 1.0  # -1 points the wall label at the body from below
    headroom: float = 1.65  # space above the body, for what is drawn there
    footroom: float = 1.35  # space below it, before the coordinate axis
    #: Where a bead's annotation sits, above and below the packing.
    label_y: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        if self.label_y is None:
            self.label_y = (self.top + 0.45, self.bottom - 0.4)

    bead_scale: float = 1.0  # shrinks the beads where the bed is cramped
    xmin: float = 0.0  # left edge of the frame
    xmax: float = 13.0  # right edge of the frame
    domain_at: Point = (_LEFT + 1.25, _AXIS_Y - 1.2)


def _cap_rx(ry: float) -> float:
    """Horizontal bulge of an end cap of half-height ``ry``."""
    return _PERSPECTIVE * ry


def _constant(value: float) -> Callable[[float], float]:
    """Return a bed profile of constant half-height."""
    return lambda x: value


# %% Bodies


def _axial_cylinder(scene: Scene, model) -> _Body:
    """A cylinder with the flow running along its axis, left to right."""
    del model
    rx = _cap_rx(_BODY_RY)
    top, bottom = _AXIS_Y + _BODY_RY, _AXIS_Y - _BODY_RY
    left_centre, right_centre = (_LEFT, _AXIS_Y), (_RIGHT, _AXIS_Y)

    scene.add(
        Ellipse(left_centre, rx, _BODY_RY, Style.BULK, stroked=False),
        Rect((_LEFT, bottom), _RIGHT - _LEFT, 2 * _BODY_RY, Style.BULK, stroked=False),
        Ellipse(right_centre, rx, _BODY_RY, Style.BULK_FACE, stroked=False),
        EllipseArc(left_centre, rx, _BODY_RY, 90, 270, Style.WALL),
        EllipseArc(left_centre, rx, _BODY_RY, -90, 90, Style.HIDDEN),
        Polyline(((_LEFT, top), (_RIGHT, top)), Style.WALL),
        Polyline(((_LEFT, bottom), (_RIGHT, bottom)), Style.WALL),
        Ellipse(right_centre, rx, _BODY_RY, Style.WALL, filled=False),
    )
    return _Body(
        top=top,
        bottom=bottom,
        bed_spans=((_LEFT + 1.8, _RIGHT - 0.85),),
        bed_profile=_constant(_BODY_RY - 0.22),
        inlet_x=_LEFT - rx,
        outlet_x=_RIGHT + rx,
        wave=(_LEFT + 0.3, _RIGHT - rx - 0.15),
        wall_at=(_LEFT + 0.78 * (_RIGHT - _LEFT), top),
        domain_at=(_LEFT + 1.2, _AXIS_Y - 1.3),
    )


def _frustum(scene: Scene, model) -> _Body:
    """A conical frustum: the cross section narrows along the flow."""
    del model
    inlet_ry, outlet_ry = 2.3, 1.45
    inlet_rx, outlet_rx = _cap_rx(inlet_ry), _cap_rx(outlet_ry)
    top, bottom = _AXIS_Y + inlet_ry, _AXIS_Y - inlet_ry
    left_centre, right_centre = (_LEFT, _AXIS_Y), (_RIGHT, _AXIS_Y)
    length = _RIGHT - _LEFT

    def radius_at(x: float) -> float:
        return inlet_ry + (outlet_ry - inlet_ry) * (x - _LEFT) / length

    silhouette = (
        (_LEFT, _AXIS_Y + inlet_ry),
        (_RIGHT, _AXIS_Y + outlet_ry),
        (_RIGHT, _AXIS_Y - outlet_ry),
        (_LEFT, _AXIS_Y - inlet_ry),
    )
    scene.add(
        Ellipse(left_centre, inlet_rx, inlet_ry, Style.BULK, stroked=False),
        Polyline(silhouette, Style.BULK, closed=True, filled=True, stroked=False),
        Ellipse(right_centre, outlet_rx, outlet_ry, Style.BULK_FACE, stroked=False),
        EllipseArc(left_centre, inlet_rx, inlet_ry, 90, 270, Style.WALL),
        EllipseArc(left_centre, inlet_rx, inlet_ry, -90, 90, Style.HIDDEN),
        Polyline(silhouette[:2], Style.WALL),
        Polyline(silhouette[2:][::-1], Style.WALL),
        Ellipse(right_centre, outlet_rx, outlet_ry, Style.WALL, filled=False),
    )

    # The two radii are what distinguishes a frustum from a cylinder. Each is
    # labelled above its end cap, where nothing else competes for the space.
    # A radius is an extent, so it is drawn without a head; the tick at the
    # axis is what says where it is measured from.
    scene.add(
        Polyline((left_centre, (_LEFT, _AXIS_Y + inlet_ry)), Style.ANNOTATION),
        Polyline(((_LEFT - TICK, _AXIS_Y), (_LEFT + TICK, _AXIS_Y)), Style.ANNOTATION),
        # Spelled as the symbol table spells them, upright script and all.
        Label((_LEFT, _AXIS_Y + inlet_ry + 0.12), r"P^{\mathrm{0}}", Anchor.ABOVE, size=0.9),
        Polyline((right_centre, (_RIGHT, _AXIS_Y + outlet_ry)), Style.ANNOTATION),
        Polyline(((_RIGHT - TICK, _AXIS_Y), (_RIGHT + TICK, _AXIS_Y)), Style.ANNOTATION),
        Label((_RIGHT, _AXIS_Y + outlet_ry + 0.12), r"P^{\mathrm{L^{\b}}}", Anchor.ABOVE, size=0.9),
    )
    return _Body(
        top=top,
        bottom=bottom,
        bed_spans=((_LEFT + 1.75, _RIGHT - 0.55),),
        bed_profile=lambda x: radius_at(x) - 0.22,
        bead_scale=0.86,
        inlet_x=_LEFT - inlet_rx,
        outlet_x=_RIGHT + outlet_rx,
        inlet_ry=inlet_ry,
        outlet_ry=outlet_ry,
        wave=(_LEFT + 0.35, _RIGHT - outlet_rx - 0.15),
        wall_at=(_LEFT + 0.55 * length, _AXIS_Y + radius_at(_LEFT + 0.55 * length)),
        domain_at=(_LEFT + 1.15, _AXIS_Y - 1.55),
    )


def _generalized_radius(fraction: float) -> float:
    """Half-height of the generalized body at ``fraction`` along the flow."""
    # A venturi-like waist makes "the cross section varies smoothly" obvious.
    return 2.05 - 0.75 * math.exp(-(((fraction - 0.45) / 0.23) ** 2))


def _generalized(scene: Scene, model) -> _Body:
    """A column whose cross-section area varies smoothly along the flow."""
    del model
    length = _RIGHT - _LEFT
    fractions = [index / _SAMPLES for index in range(_SAMPLES + 1)]
    wall = [(_LEFT + fraction * length, _generalized_radius(fraction)) for fraction in fractions]

    upper = tuple((x, _AXIS_Y + radius) for x, radius in wall)
    lower = tuple((x, _AXIS_Y - radius) for x, radius in wall)
    inlet_ry, outlet_ry = _generalized_radius(0.0), _generalized_radius(1.0)
    inlet_rx, outlet_rx = _cap_rx(inlet_ry), _cap_rx(outlet_ry)
    left_centre, right_centre = (_LEFT, _AXIS_Y), (_RIGHT, _AXIS_Y)

    scene.add(
        Ellipse(left_centre, inlet_rx, inlet_ry, Style.BULK, stroked=False),
        Polyline(upper + lower[::-1], Style.BULK, closed=True, filled=True, stroked=False),
        Ellipse(right_centre, outlet_rx, outlet_ry, Style.BULK_FACE, stroked=False),
        EllipseArc(left_centre, inlet_rx, inlet_ry, 90, 270, Style.WALL),
        EllipseArc(left_centre, inlet_rx, inlet_ry, -90, 90, Style.HIDDEN),
        Polyline(upper, Style.WALL),
        Polyline(lower, Style.WALL),
        Ellipse(right_centre, outlet_rx, outlet_ry, Style.WALL, filled=False),
    )

    # Mark the cross section at the waist, where the variation is clearest.
    waist_x = _LEFT + 0.45 * length
    waist_ry = _generalized_radius(0.45)
    scene.add(
        EllipseArc((waist_x, _AXIS_Y), _cap_rx(waist_ry), waist_ry, -90, 90, Style.HIDDEN),
        EllipseArc((waist_x, _AXIS_Y), _cap_rx(waist_ry), waist_ry, 90, 270, Style.HIDDEN),
        Label((waist_x, _AXIS_Y - waist_ry - 0.12), "A(x)", Anchor.BELOW, size=0.9),
    )
    top = _AXIS_Y + max(radius for _, radius in wall)
    return _Body(
        top=top,
        bottom=2 * _AXIS_Y - top,
        bed_spans=((_LEFT + 1.75, _RIGHT - 0.6),),
        bed_profile=lambda x: _generalized_radius((x - _LEFT) / length) - 0.24,
        bead_scale=0.84,
        inlet_x=_LEFT - inlet_rx,
        outlet_x=_RIGHT + outlet_rx,
        inlet_ry=inlet_ry,
        outlet_ry=outlet_ry,
        wave=(_LEFT + 0.35, _RIGHT - outlet_rx - 0.15),
        flow_symbol="x",
        flow_ticks=(r"x^{\mathrm{in}}", r"x^{\mathrm{out}}"),
        dispersion=r"D^{\mathrm{x}}",
        wall_at=(_LEFT + 0.82 * length, _AXIS_Y + _generalized_radius(0.82)),
        domain_at=(_LEFT + 1.1, _AXIS_Y - 1.5),
    )


# %% Radial flow column


#: Layout of the upright annular body, in scene units.
_RADIAL_OUTER, _RADIAL_INNER, _RADIAL_HALF_HEIGHT = 2.45, 0.5, 2.5


def _radial_cylinder(scene: Scene, model) -> _Body:
    """An annular column standing upright, as in ``cylinder-annular.tex``.

    The feed enters the central channel -- the model's inflow boundary sits at
    the inner radius -- and spreads radially outwards through the packing, so
    the two halves of the drawing are mirror images of one flow path.
    """
    del model
    centre_x = 0.5 * (_LEFT + _RIGHT)
    outer, inner, half_height = _RADIAL_OUTER, _RADIAL_INNER, _RADIAL_HALF_HEIGHT
    top, bottom = _AXIS_Y + half_height, _AXIS_Y - half_height
    cap_ry = _PERSPECTIVE * outer
    inner_cap_ry = cap_ry * inner / outer

    scene.add(
        Rect((centre_x - outer, bottom), 2 * outer, 2 * half_height, Style.BULK, stroked=False),
        Ellipse((centre_x, bottom), outer, cap_ry, Style.BULK, stroked=False),
        Ellipse((centre_x, top), outer, cap_ry, Style.BULK_FACE, stroked=False),
        # The channel is empty and opens at the top face.
        Ellipse((centre_x, top), inner, inner_cap_ry, Style.VOID, stroked=False),
        # Base: the near half of the rim is drawn solid, the far half dashed.
        EllipseArc((centre_x, bottom), outer, cap_ry, 180, 360, Style.WALL),
        EllipseArc((centre_x, bottom), outer, cap_ry, 0, 180, Style.HIDDEN),
        Ellipse((centre_x, bottom), inner, inner_cap_ry, Style.HIDDEN, filled=False),
        Polyline(((centre_x - outer, bottom), (centre_x - outer, top)), Style.WALL),
        Polyline(((centre_x + outer, bottom), (centre_x + outer, top)), Style.WALL),
        # The channel wall runs the full height, hidden by the packing.
        Polyline(((centre_x - inner, bottom), (centre_x - inner, top)), Style.HIDDEN),
        Polyline(((centre_x + inner, bottom), (centre_x + inner, top)), Style.HIDDEN),
        Ellipse((centre_x, top), outer, cap_ry, Style.WALL, filled=False),
        Ellipse((centre_x, top), inner, inner_cap_ry, Style.WALL, filled=False),
    )
    return _Body(
        top=top + cap_ry,
        bottom=bottom - cap_ry,
        bed_spans=(
            (centre_x - outer + 0.18, centre_x - inner - 0.18),
            (centre_x + inner + 0.18, centre_x + outer - 0.18),
        ),
        bed_profile=_constant(half_height - 0.75),
        bead_scale=0.74,
        inlet_x=centre_x - outer,
        outlet_x=centre_x + outer,
        flow_symbol=r"\rho",
        flow_ticks=(r"P^{\mathrm{inner}}", r"P^{\mathrm{outer}}"),
        dispersion=r"D^{\mathrm{rad}}",
        draws_own_flow=True,
        # Both end faces are impermeable; the wall is pointed out on the base
        # so that it is not confused with the domain on the open top face.
        wall_at=(centre_x - 0.55 * outer, bottom - 0.835 * cap_ry),
        wall_side=-1.0,
        headroom=1.35,  # the feed arrives from above the open channel
        footroom=1.45,  # and the wall is named below the base
        # The body is tall, so a bead's annotation leaves it sideways rather
        # than reaching all the way past the end faces.
        label_y=(_AXIS_Y + 2.1, _AXIS_Y - 2.1),
        xmin=1.8,
        xmax=10.7,
        domain_at=(centre_x - 0.6 * outer, top + 0.02),
    )


def _radial_flow(scene: Scene, body: _Body, dispersed: bool) -> None:
    """Draw the feed over the shell and the stream collected in the channel.

    The inflow boundary of the radial model sits at the outer radius, so the
    two mirrored paths run inwards and leave through the top of the channel.
    """
    centre_x = 0.5 * (body.inlet_x + body.outlet_x)
    outer, inner = _RADIAL_OUTER, _RADIAL_INNER
    cap_ry = _PERSPECTIVE * outer
    top = body.top - cap_ry
    mouth = top + cap_ry * inner / outer

    for side in (-1.0, 1.0):
        scene.add(
            Polyline(
                _flow_points(
                    (centre_x + side * (outer + _LEAD), _AXIS_Y),
                    (centre_x + side * inner, _AXIS_Y),
                    (centre_x + side * (outer - 0.1), centre_x + side * (inner + 0.2)),
                    dispersed,
                ),
                Style.FLOW,
                arrow=True,
            )
        )
    scene.add(
        Arrow((centre_x, mouth + 0.08), (centre_x, mouth + 1.5), Style.FLOW),
        Label((centre_x - 0.18, mouth + 1.0), r"\Gamma^{\mathrm{out}}", Anchor.RIGHT),
        Label((centre_x + outer + 0.78, _AXIS_Y + 0.24), r"\Gamma^{\mathrm{in}}", Anchor.ABOVE),
        Label((centre_x + outer + 0.78, _AXIS_Y - 0.3), "Flow", Anchor.BELOW, size=0.85, italic=False),
    )
    if dispersed:
        scene.add(Label((centre_x - outer - 0.78, _AXIS_Y - 0.3), body.dispersion, Anchor.BELOW, size=0.85))


# %% Tank


def _tank(scene: Scene, model) -> _Body:
    """A continuous stirred tank: no spatial coordinate, so no axes."""
    half_width, half_height = 2.35, 2.5
    centre_x = 0.5 * (_LEFT + _RIGHT)
    rim_ry = _cap_rx(half_width)
    top, bottom = _AXIS_Y + half_height, _AXIS_Y - half_height

    scene.add(
        Rect((centre_x - half_width, bottom), 2 * half_width, 2 * half_height, Style.BULK, stroked=False),
        Ellipse((centre_x, bottom), half_width, rim_ry, Style.BULK, stroked=False),
        Ellipse((centre_x, top), half_width, rim_ry, Style.BULK_FACE, stroked=False),
        EllipseArc((centre_x, bottom), half_width, rim_ry, 180, 360, Style.WALL),
        Polyline(((centre_x - half_width, bottom), (centre_x - half_width, top)), Style.WALL),
        Polyline(((centre_x + half_width, bottom), (centre_x + half_width, top)), Style.WALL),
        Ellipse((centre_x, top), half_width, rim_ry, Style.WALL, filled=False),
    )

    # Stirrer: the well-mixed assumption is what the tank model is about. It
    # sits low in the vessel so that the beads above it stay legible.
    blade_y, blade_half = bottom + 0.6, 0.7
    scene.add(
        Polyline(((centre_x, top + 0.95), (centre_x, blade_y)), Style.WALL),
        Polyline(((centre_x - blade_half, blade_y + 0.2), (centre_x + blade_half, blade_y - 0.2)), Style.WALL),
        Polyline(((centre_x - blade_half, blade_y - 0.2), (centre_x + blade_half, blade_y + 0.2)), Style.WALL),
        EllipseArc((centre_x, blade_y + 0.75), 1.3, 0.4, 15, 300, Style.ANNOTATION, arrow=True),
    )

    # Both nozzles sit on the front wall, below the rim and above the base.
    inlet_y, outlet_y = top - rim_ry - 0.4, bottom + 0.35
    scene.add(
        Arrow((centre_x - half_width - 1.5, inlet_y), (centre_x - half_width - 0.15, inlet_y), Style.FLOW),
        Arrow((centre_x + half_width + 0.15, outlet_y), (centre_x + half_width + 1.5, outlet_y), Style.FLOW),
        Label((centre_x - half_width - 0.85, inlet_y + 0.22), r"\Gamma^{\mathrm{in}}", Anchor.ABOVE),
        Label((centre_x + half_width + 1.4, outlet_y + 0.22), r"\Gamma^{\mathrm{out}}", Anchor.ABOVE),
    )

    if getattr(model, "has_filter", False):
        # A filter holds the particles back, so only liquid leaves the tank.
        screen_x = centre_x + half_width + 0.6
        for offset in (-0.16, 0.0, 0.16):
            scene.add(Polyline(((screen_x + offset, outlet_y - 0.4), (screen_x + offset, outlet_y + 0.4)), Style.WALL))
        scene.add(Label((screen_x, outlet_y - 0.5), "filter", Anchor.BELOW, size=0.8, italic=False))

    bed_low, bed_high = blade_y + 1.25, top - rim_ry - 0.4
    return _Body(
        top=top + 0.95,
        bottom=bottom - rim_ry,
        bed_spans=((centre_x - half_width + 0.5, centre_x + half_width - 0.5),),
        bed_profile=_constant(0.5 * (bed_high - bed_low)),
        bed_axis=0.5 * (bed_low + bed_high),
        # A suspension has no flow path to keep clear, so the beads may sit
        # anywhere in the bed, across its middle included.
        draws_own_flow=True,
        has_wall=False,  # a homogeneous tank has no wall boundary condition
        headroom=0.4,
        footroom=0.45,
        bead_scale=0.72,
        # The suspension sits in the middle of a tall vessel, so a bead's
        # annotation leaves past its side rather than past the rim.
        label_y=(bed_high + 1.35, bed_low - 0.45),
        xmin=2.0,
        xmax=11.0,
    )


#: Body builder per column geometry. The tank is selected by resolution.
_BODIES = {
    "Axial": _axial_cylinder,
    "Radial": _radial_cylinder,
    "Frustum": _frustum,
    "Generalized": _generalized,
}


# %% Annotation layers


def _flow_points(
    start: Point,
    end: Point,
    wave: tuple[float, float],
    dispersed: bool,
) -> tuple[Point, ...]:
    """Return the points of a flow path from ``start`` to ``end``.

    Pure convection carries a sharp front straight through the column; with
    dispersion the path meanders around the packing instead, which is the one
    difference between the two transport equations a drawing can show.
    """
    y = start[1]
    if not dispersed:
        return (start, end)

    wave_from, wave_to = wave
    span = wave_to - wave_from
    # A whole number of meanders, each about one wavelength long, so that the
    # path rejoins the straight stretch at the far end of the packing.
    meanders = max(1.0, round(abs(span) / _WAVELENGTH))
    amplitude = min(_WAVE, 0.22 * abs(span) / meanders)

    points = [start, (wave_from, y)]
    for index in range(1, _SAMPLES + 1):
        fraction = index / _SAMPLES
        # Tapered at both ends so the meander grows out of the straight path.
        envelope = min(1.0, 4.0 * min(fraction, 1.0 - fraction))
        points.append(
            (
                wave_from + fraction * span,
                y + amplitude * envelope * math.sin(2.0 * math.pi * meanders * fraction),
            )
        )
    points.append(end)
    return tuple(points)


def _flow_path(scene: Scene, body: _Body, dispersed: bool) -> None:
    """Draw the flow through a body with through-flow, inlet to outlet."""
    start = (body.inlet_x - _LEAD, _AXIS_Y)
    end = (body.outlet_x + _LEAD, _AXIS_Y)
    scene.add(
        Polyline(_flow_points(start, end, body.wave, dispersed), Style.FLOW, arrow=True),
        Label((body.inlet_x - 0.75, _AXIS_Y - 0.3), "Flow", Anchor.BELOW, size=0.85, italic=False),
    )
    # The two faces are named on the cap itself, rather than beside the flow
    # path: a label floating at the end of the path could as easily be naming
    # the stream as the face it crosses. The arrow lands where the cap is
    # widest and nothing else is written -- the far side of its upper quarter.
    for sign, x, ry, name in (
        (-1.0, _LEFT, body.inlet_ry, r"\Gamma^{\mathrm{in}}"),
        (1.0, _RIGHT, body.outlet_ry, r"\Gamma^{\mathrm{out}}"),
    ):
        corner = (x + sign * _cap_rx(ry) * _DIAGONAL, _AXIS_Y + ry * _DIAGONAL)
        boundary(scene, corner, (sign * 0.5, 0.87), name, reach=0.75)
    if dispersed:
        scene.add(Label((body.outlet_x + 0.75, _AXIS_Y - 0.3), body.dispersion, Anchor.BELOW, size=0.85))


def _wall_label(scene: Scene, body: _Body) -> None:
    """Point out the column wall from outside the body."""
    x, y = body.wall_at
    boundary(scene, (x, y), (0.0, body.wall_side), r"\Gamma^{\mathrm{wall}}", reach=0.8)


def _flow_coordinate(scene: Scene, body: _Body, baseline: float) -> None:
    """Draw the coordinate that runs along the flow direction."""
    start, end = body.flow_ticks
    dimension(scene, (_LEFT, baseline), (body.outlet_x, baseline), body.flow_symbol, start, end)


def _radial_flow_coordinate(scene: Scene, body: _Body, baseline: float) -> None:
    """Draw the radial coordinate of a radial flow column, below the body.

    It is measured from the axis of the column, so it starts under the centre
    of the body and its two ticks line up with the channel and the outer wall.
    """
    centre_x = 0.5 * (body.inlet_x + body.outlet_x)
    inner_tick, outer_tick = centre_x + _RADIAL_INNER, centre_x + _RADIAL_OUTER
    inner, outer = body.flow_ticks
    # Dashed guides tie the two ticks to the channel and to the outer shell.
    for guide in (inner_tick, outer_tick):
        scene.add(Polyline(((guide, body.bottom - 0.12), (guide, baseline + 0.52)), Style.HIDDEN))
    dimension(scene, (inner_tick, baseline), (outer_tick, baseline), body.flow_symbol, inner, outer)


def _radial_coordinate(scene: Scene, body: _Body) -> None:
    """Draw the radial coordinate a 2D or 3D model resolves.

    Both cross-sectional coordinates are drawn on the near end cap, where the
    circular cross section they parametrise is actually visible -- and where
    they do not run through the bed.
    """
    ry = body.top - _AXIS_Y
    dimension(
        scene,
        (_RIGHT, _AXIS_Y),
        (_RIGHT, _AXIS_Y + ry - 0.04),
        r"\rho",
        "0",
        r"P^{\mathrm{c}}",
        flip=True,
    )


def _angular_coordinate(scene: Scene, body: _Body) -> None:
    """Draw the angular coordinate a 3D model resolves, on the near cap."""
    ry = body.top - _AXIS_Y
    # Swept through the lower half of the cap, clear of the radial arrow.
    scene.add(
        EllipseArc((_RIGHT, _AXIS_Y), _cap_rx(ry) * 0.62, ry * 0.62, -170, -10, Style.AXIS, arrow=True),
        Label((_RIGHT + _cap_rx(ry) * 0.7, _AXIS_Y - ry * 0.66), r"\varphi", Anchor.BELOW_RIGHT),
    )


# %% Particle types


#: How many beads a bed holds at most, and the clearance between two of them.
#: A bed is a bed: the beads are scattered over the whole cross section rather
#: than lined up in rows, which is what a packing looks like and what keeps the
#: eye from reading the rows as a structure the model has.
_BEAD_COUNT = 11
_BEAD_GAP = 0.12

#: How many candidate positions the scatter may try before it gives up, so a
#: bed too small for ``_BEAD_COUNT`` beads ends the search rather than looping.
_BEAD_TRIES = 240


def _particle_types(model) -> list:
    """Return the distinct particle types of ``model``, in display order."""
    # Defined here rather than imported from the package root, which imports
    # this module.
    counts = getattr(model, "par_type_counts", None)
    if counts:
        return list(counts)
    return list(getattr(model, "particle_models", []) or [])


def _bead_radius(index: int, count: int) -> float:
    """Give each particle type its own size, smallest type first."""
    if count < 2:
        return 0.28
    return 0.21 + 0.13 * index / (count - 1)


def _beads(scene: Scene, body: _Body, model) -> None:
    """Fill the bed with beads, one visual class per particle type."""
    types = _particle_types(model)
    if not types or not body.bed_spans:
        return

    count = len(types)
    # Place every bead before drawing any, because the annotation goes on the
    # bead of each type that sits closest to the edge of the bed, where its
    # leader reaches the label without crossing the packing.
    beads = _scatter(body, count)

    for span_index, particle_index, centre, radius, reach in beads:
        del span_index
        particle = types[particle_index]
        side = 1.0 if centre[1] >= body.bed_axis else -1.0
        limiting = not getattr(particle, "nonlimiting_filmDiff", False)
        if limiting:
            scene.add(Ellipse(centre, reach, reach, Style.FILM))
            # Short arrows on small beads, so a cramped bed does not fill up.
            _film_flux(scene, centre, reach, side, _FLUX_LENGTH * radius / 0.28)
        scene.add(Ellipse(centre, radius, radius, Style.PARTICLE))

    # One bead of each type carries the annotation for that type: what the
    # arrows mean and, with more than one type, which bead is which. The types
    # take turns between the two sides of the bed so the labels do not meet.
    for particle_index in range(count):
        wanted = 1.0 if particle_index % 2 == 0 else -1.0
        mine = [bead for bead in beads if bead[1] == particle_index]
        if not mine:
            continue
        # On the type's own side of the bed and, where there is a choice, in
        # its upstream half: the far end of the bed is where the wall and the
        # outlet are named, and a leader across those would be unreadable.
        start, stop = body.bed_spans[0]
        wanted_side = [bead for bead in mine if (bead[2][1] - body.bed_axis) * wanted > 0]
        upstream = [bead for bead in wanted_side if bead[2][0] < start + 0.45 * (stop - start)]
        bead = max(upstream or wanted_side or mine, key=lambda bead: (bead[2][1] - body.bed_axis) * wanted)
        _, _, centre, radius, reach = bead
        side = 1.0 if centre[1] >= body.bed_axis else -1.0
        limiting = not getattr(types[particle_index], "nonlimiting_filmDiff", False)
        start = _along_flux(centre, side, reach + _FLUX_GAP + (_FLUX_LENGTH * radius / 0.28 if limiting else 0.0))
        _bead_label(scene, body, start, side, particle_index, count, limiting)


def _halton(index: int, base: int) -> float:
    """Return the ``index``-th point of the Halton sequence over ``base``.

    A scatter has to look unplanned and stay byte-identical from one run to
    the next, which rules out a random generator; this fills the bed evenly
    without ever repeating a pattern the eye can pick up.
    """
    fraction, result = 1.0, 0.0
    while index:
        fraction /= base
        result += fraction * (index % base)
        index //= base
    return result


def _scatter(body: _Body, count: int) -> list[tuple[int, int, Point, float, float]]:
    """Return where the beads of a bed go, as (span, type, centre, r, reach).

    Candidate positions are tried in turn and kept when the whole bead, film
    and all, fits inside the bed profile and clears every bead already placed.
    A bed too small for ``_BEAD_COUNT`` of them simply holds fewer.
    """
    placed: list[tuple[int, int, Point, float, float]] = []
    spans = body.bed_spans

    for try_index in range(1, _BEAD_TRIES + 1):
        if len(placed) == _BEAD_COUNT:
            break
        # The types take turns, so every one of them is in the bed; a bed with
        # several spans gets a bead in each in turn, an annulus on both sides.
        particle_index = len(placed) % count
        start, stop = spans[try_index % len(spans)]
        radius = _bead_radius(particle_index, count) * body.bead_scale
        # A film resistance is drawn as the liquid film it stands for, so the
        # bead needs room for it; with non-limiting film diffusion it does not.
        reach = radius * 1.34

        x = start + _halton(try_index, 2) * (stop - start)
        depth = body.bed_profile(x) - reach
        if depth <= 0.0:  # the bed is too shallow here to hold a whole bead
            continue
        centre = (x, body.bed_axis + (2.0 * _halton(try_index, 3) - 1.0) * depth)
        if any(
            math.hypot(centre[0] - other[0], centre[1] - other[1]) < reach + other_reach + _BEAD_GAP
            for _, _, other, _, other_reach in placed
        ):
            continue
        placed.append((try_index % len(spans), particle_index, centre, radius, reach))
    return placed


#: Where a film flux arrow comes from, how long it is, and how far its tip
#: stays clear of the film, in scene units.
_FLUX_DIRECTION, _FLUX_LENGTH, _FLUX_GAP = (-0.72, 0.69), 0.34, 0.05


def _along_flux(point: Point, side: float, distance: float) -> Point:
    """Return the point ``distance`` out from ``point`` along the flux arrow."""
    return (point[0] + _FLUX_DIRECTION[0] * distance, point[1] + _FLUX_DIRECTION[1] * side * distance)


def _film_flux(scene: Scene, centre: Point, reach: float, side: float, length: float) -> Point:
    """Draw the flux through the film of one bead and return where it starts.

    Film diffusion is what the ring around a bead stands for, so it is drawn
    on the bead itself rather than explained in a key beside the figure.
    """
    tip = _along_flux(centre, side, reach + _FLUX_GAP)
    tail = _along_flux(tip, side, length)
    # Double-headed: the film is crossed in whichever direction the gradient
    # points, and a single head would claim the bead only ever takes up. Drawn
    # as a flux rather than a flow, so the heads suit an arrow this short.
    scene.add(Arrow(tail, tip, Style.FLUX, double=True))
    return tail


def _bead_label(scene: Scene, body: _Body, start: Point, side: float, index: int, count: int, limiting: bool) -> None:
    """Name the film resistance, and the particle type when there are several.

    The label sits clear of the packing and is joined to the bead by a leader
    that continues the flux arrow, so the two read as one annotation.
    """
    edge = body.label_y[0] if side > 0 else body.label_y[1]
    run = (edge - start[1]) / (_FLUX_DIRECTION[1] * side)
    anchor_point = (start[0] + _FLUX_DIRECTION[0] * run, edge)
    scene.add(Polyline((start, anchor_point), Style.ANNOTATION))

    # Named rather than lettered: the arrows are the mechanism, and "film
    # diffusion" says what they are without the reader looking up k^f first.
    text = "film diffusion" if limiting else r"c^{\b} = c^{\p}"
    if count > 1:
        text = f"j = {index + 1}: " + text
    # Written away from the inlet, which is named on the cap just left of it.
    scene.add(Label(anchor_point, text, Anchor.ABOVE_RIGHT if side > 0 else Anchor.BELOW_RIGHT, size=0.85))


# %% Assembly


def build_column_scene(column_model) -> Scene:
    """Return the figure of ``column_model``'s interstitial volume.

    ``column_model`` is a :class:`~src.model_column.Column`; it is typed loosely
    to keep this module free of a dependency on the Streamlit-backed dataclass.
    """
    scene = Scene(name="column", xmin=0.0, ymin=0.0, xmax=13.0, ymax=0.0)

    is_tank = column_model.resolution == "0D"
    is_radial = column_model.column_type == "Radial" and not is_tank
    build_body = _tank if is_tank else _BODIES[column_model.column_type]
    body = build_body(scene, column_model)
    dispersed = bool(getattr(column_model, "has_axial_dispersion", False))

    if is_radial:
        _radial_flow(scene, body, dispersed)
    elif not body.draws_own_flow:
        _flow_path(scene, body, dispersed)
    # After the flow path: where the two cross, the bead is in front, so the
    # stream reads as passing between the beads rather than through one.
    _beads(scene, body, column_model)
    if body.has_wall:
        _wall_label(scene, body)
    if not is_tank:
        scene.add(Label(body.domain_at, r"\Omega^{\b}", Anchor.CENTER, size=1.1))

    top = body.top + body.headroom

    if getattr(column_model, "has_radial_coordinate", False):
        _radial_coordinate(scene, body)
    if getattr(column_model, "has_angular_coordinate", False):
        _angular_coordinate(scene, body)

    bottom = body.bottom - body.footroom
    if not is_tank:
        if is_radial:
            _radial_flow_coordinate(scene, body, bottom)
        else:
            _flow_coordinate(scene, body, bottom)
        bottom -= 0.75

    scene.xmin, scene.xmax = body.xmin, body.xmax
    scene.ymin, scene.ymax = bottom - 0.3, top
    scene.title = "Interstitial volume of the " + column_model.model_name()
    return scene


__all__ = ["build_column_scene"]
