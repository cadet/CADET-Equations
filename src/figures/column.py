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

_CORRIDOR = 0.62  # half-height of the bead-free lane the flow path runs in
_WAVE = 0.46  # amplitude of a dispersed (meandering) flow path
_WAVELENGTH = 2.3  # distance it takes to complete one meander
_LEAD = 1.5  # how far a flow path reaches beyond the body


@dataclass
class _Body:
    """Where the annotation layers may attach to a drawn body.

    ``inlet_x`` and ``outlet_x`` are where the flow path meets the silhouette
    and ``wave`` is the stretch of it that runs through the packing, so that a
    dispersed path meanders there and nowhere else.  ``bed_spans`` are the
    x-ranges the packing occupies and ``bed_profile`` gives its half-height
    about ``bed_axis``, which together keep beads inside a body of any shape.
    ``flow_symbol`` and ``flow_ticks`` name the coordinate that runs along the
    flow, and ``dispersion`` is the symbol of the flow-directional dispersion
    coefficient.
    """

    top: float
    bottom: float
    bed_spans: tuple[tuple[float, float], ...] = ()
    bed_profile: Callable[[float], float] = field(default=lambda x: 0.0)
    bed_axis: float = _AXIS_Y
    corridor: float = _CORRIDOR
    inlet_x: float = _LEFT
    outlet_x: float = _RIGHT
    wave: tuple[float, float] = (_LEFT, _RIGHT)
    flow_symbol: str = "z"
    flow_ticks: tuple[str, str] = ("0", r"L^{\b}")
    dispersion: str = r"D^{\mathrm{ax}}"
    draws_own_flow: bool = False
    has_wall: bool = True
    wall_at: Point = (_LEFT + 0.78 * (_RIGHT - _LEFT), _AXIS_Y + _BODY_RY)
    wall_side: float = 1.0  # -1 points the wall label at the body from below
    headroom: float = 1.65  # space above the body, for what is drawn there
    footroom: float = 0.95  # space below it, before the coordinate axis
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
    scene.add(
        Arrow(left_centre, (_LEFT, _AXIS_Y + inlet_ry), Style.ANNOTATION),
        # Spelled as the symbol table spells them, upright script and all.
        Label((_LEFT, _AXIS_Y + inlet_ry + 0.12), r"P^{\mathrm{0}}", Anchor.ABOVE, size=0.9),
        Arrow(right_centre, (_RIGHT, _AXIS_Y + outlet_ry), Style.ANNOTATION),
        Label((_RIGHT, _AXIS_Y + outlet_ry + 0.12), r"P^{\mathrm{L^{\b}}}", Anchor.ABOVE, size=0.9),
    )
    return _Body(
        top=top,
        bottom=bottom,
        bed_spans=((_LEFT + 1.75, _RIGHT - 0.55),),
        bed_profile=lambda x: radius_at(x) - 0.22,
        corridor=0.58,
        bead_scale=0.86,
        inlet_x=_LEFT - inlet_rx,
        outlet_x=_RIGHT + outlet_rx,
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
        corridor=0.58,
        bead_scale=0.84,
        inlet_x=_LEFT - inlet_rx,
        outlet_x=_RIGHT + outlet_rx,
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
        corridor=0.6,
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
        corridor=-0.45,
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
        Label((body.inlet_x - 0.75, _AXIS_Y + 0.24), r"\Gamma^{\mathrm{in}}", Anchor.ABOVE),
        Label((body.outlet_x + 0.75, _AXIS_Y + 0.24), r"\Gamma^{\mathrm{out}}", Anchor.ABOVE),
        Label((body.inlet_x - 0.75, _AXIS_Y - 0.3), "Flow", Anchor.BELOW, size=0.85, italic=False),
    )
    if dispersed:
        scene.add(Label((body.outlet_x + 0.75, _AXIS_Y - 0.3), body.dispersion, Anchor.BELOW, size=0.85))


def _wall_label(scene: Scene, body: _Body) -> None:
    """Point out the column wall from outside the body."""
    x, y = body.wall_at
    side = body.wall_side
    anchor = Anchor.ABOVE if side > 0 else Anchor.BELOW
    scene.add(
        Arrow((x, y + side * 0.8), (x, y + side * 0.1), Style.ANNOTATION),
        Label((x, y + side * 0.87), r"\Gamma^{\mathrm{wall}}", anchor),
    )


def _flow_coordinate(scene: Scene, body: _Body, baseline: float) -> None:
    """Draw the coordinate that runs along the flow direction."""
    tick = 0.13
    start, end = body.flow_ticks
    scene.add(
        Arrow((_LEFT, baseline), (body.outlet_x + 0.8, baseline), Style.AXIS),
        Polyline(((_LEFT, baseline - tick), (_LEFT, baseline + tick)), Style.AXIS),
        Polyline(((body.outlet_x, baseline - tick), (body.outlet_x, baseline + tick)), Style.AXIS),
        Label((_LEFT, baseline - tick - 0.05), start, Anchor.BELOW, size=0.9),
        Label((body.outlet_x, baseline - tick - 0.05), end, Anchor.BELOW, size=0.9),
        Label((body.outlet_x + 0.8, baseline + 0.1), body.flow_symbol, Anchor.ABOVE_RIGHT),
    )


def _radial_flow_coordinate(scene: Scene, body: _Body, baseline: float) -> None:
    """Draw the radial coordinate of a radial flow column, below the body.

    It is measured from the axis of the column, so it starts under the centre
    of the body and its two ticks line up with the channel and the outer wall.
    """
    centre_x = 0.5 * (body.inlet_x + body.outlet_x)
    inner_tick, outer_tick = centre_x + _RADIAL_INNER, centre_x + _RADIAL_OUTER
    inner, outer = body.flow_ticks
    tick = 0.13
    # Dashed guides tie the two ticks to the channel and to the outer shell.
    for guide in (inner_tick, outer_tick):
        scene.add(Polyline(((guide, body.bottom - 0.12), (guide, baseline + 0.28)), Style.HIDDEN))
    scene.add(
        Arrow((centre_x, baseline), (outer_tick + 1.0, baseline), Style.AXIS),
        Polyline(((inner_tick, baseline - tick), (inner_tick, baseline + tick)), Style.AXIS),
        Polyline(((outer_tick, baseline - tick), (outer_tick, baseline + tick)), Style.AXIS),
        Label((inner_tick, baseline - tick - 0.05), inner, Anchor.BELOW, size=0.9),
        Label((outer_tick, baseline - tick - 0.05), outer, Anchor.BELOW, size=0.9),
        Label((outer_tick + 1.0, baseline + 0.1), body.flow_symbol, Anchor.ABOVE_RIGHT),
    )


def _radial_coordinate(scene: Scene, body: _Body) -> None:
    """Draw the radial coordinate a 2D or 3D model resolves.

    Both cross-sectional coordinates are drawn on the near end cap, where the
    circular cross section they parametrise is actually visible -- and where
    they do not run through the bed.
    """
    ry = body.top - _AXIS_Y
    scene.add(
        Arrow((_RIGHT, _AXIS_Y), (_RIGHT, _AXIS_Y + ry - 0.04), Style.AXIS),
        Label((_RIGHT - 0.1, _AXIS_Y + 0.55 * ry), r"\rho", Anchor.RIGHT),
        Label((_RIGHT, body.top + 0.12), r"P^{\mathrm{c}}", Anchor.ABOVE, size=0.9),
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


#: How many beads a bed holds, and how deep in it they sit -- the depths run
#: through a cycle rather than alternating so that the packing does not read as
#: a grid. Both are fixed rather than sampled so that the rendered figure stays
#: byte-stable for the reference tests.
_BEAD_COUNT = 8
_ACROSS = (0.82, 0.2, 0.58, 0.34, 0.92, 0.46)


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
    spans = len(body.bed_spans)

    # Plan the beads before placing any, because where a bead goes along its
    # stretch of packing depends on how many others share that stretch. Beads
    # alternate between the two sides of the flow path and only then move on
    # to the next stretch, so a column with two separate beds -- an annulus --
    # shows every particle type on both sides of both of them.
    plan = [
        (spot % count, (spot // (2 * count)) % spans, 1.0 if (spot // count) % 2 == 0 else -1.0)
        for spot in range(_BEAD_COUNT)
    ]
    sizes: dict[tuple[int, float], int] = {}
    for _, span_index, side in plan:
        sizes[(span_index, side)] = sizes.get((span_index, side), 0) + 1

    placed: dict[tuple[int, float], int] = {}
    labelled: set[int] = set()
    for particle_index, span_index, side in plan:
        order = placed.get((span_index, side), 0)
        placed[(span_index, side)] = order + 1
        total = sizes[(span_index, side)]

        particle = types[particle_index]
        radius = _bead_radius(particle_index, count) * body.bead_scale
        # A film resistance is drawn as the liquid film it stands for, so the
        # bead needs room for it; with non-limiting film diffusion it does not.
        reach = radius * 1.34

        start, stop = body.bed_spans[span_index]
        usable = max(stop - start - 2.0 * reach, 0.0)
        if total < 3:
            # Too few to stagger without crowding: spread them to the edges.
            along = order / (total - 1) if total > 1 else 0.5
        else:
            # The two sides are offset by half a step so the packing does not
            # read as a grid of columns.
            along = order / (total - 1) if side > 0 else (order + 0.5) / total
        x = start + reach + along * usable

        near, far = body.corridor + reach, body.bed_profile(x) - reach
        if far <= near:  # the bed is too shallow here to hold a whole bead
            continue
        across = _ACROSS[(order + (0 if side > 0 else 3)) % len(_ACROSS)]
        centre = (x, body.bed_axis + side * (near + across * (far - near)))

        limiting = not getattr(particle, "nonlimiting_filmDiff", False)
        leader_from = _along_flux(centre, side, reach + _FLUX_GAP)
        if limiting:
            scene.add(Ellipse(centre, reach, reach, Style.FILM))
            # Short arrows on small beads, so a cramped bed does not fill up.
            leader_from = _film_flux(scene, centre, reach, side, _FLUX_LENGTH * radius / 0.28)
        scene.add(Ellipse(centre, radius, radius, Style.PARTICLE))

        # One bead of each type carries the annotation for that type: what the
        # arrows mean and, with more than one type, which bead is which. The
        # types take turns between the two sides so the labels do not meet.
        wanted_side = 1.0 if particle_index % 2 == 0 else -1.0
        if particle_index not in labelled and span_index == 0 and side == wanted_side:
            labelled.add(particle_index)
            _bead_label(scene, body, leader_from, side, particle_index, count, limiting)


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
    scene.add(Arrow(tail, tip, Style.FLOW))
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

    text = r"k^{\mathrm{f}}" if limiting else r"c^{\b} = c^{\p}"
    if count > 1:
        text = f"j = {index + 1}: " + text
    scene.add(Label(anchor_point, text, Anchor.ABOVE if side > 0 else Anchor.BELOW, size=0.85))


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

    _beads(scene, body, column_model)
    if is_radial:
        _radial_flow(scene, body, dispersed)
    elif not body.draws_own_flow:
        _flow_path(scene, body, dispersed)
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
