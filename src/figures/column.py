"""Build the scene of the interstitial volume (column) figure.

A column figure is assembled from four independent layers, each driven by a
different part of the model:

=====================  ====================================================
layer                  driven by
=====================  ====================================================
body                   ``column_type`` and, for the tank, ``resolution``
coordinate arrows      which coordinates the resolution resolves
dispersion inset       ``has_axial_dispersion`` (1D only)
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
from dataclasses import dataclass

from src.figures.scene import (
    Anchor,
    Arrow,
    Ellipse,
    EllipseArc,
    Label,
    Polyline,
    Rect,
    Scene,
    Style,
)

# %% Shared layout, in scene units

_LEFT, _RIGHT = 2.2, 10.2  # the straight part of a horizontal body
_AXIS_Y = 4.5  # height of the body's axis of symmetry
_BODY_RY = 1.75  # half-height of a cylindrical body
_PERSPECTIVE = 0.355  # how far an end cap bulges out, relative to its height
_SAMPLES = 48  # points used to approximate a smooth wall


@dataclass
class _Body:
    """Where the annotation layers may attach to a drawn body.

    ``inlet_x`` and ``outlet_x`` are where flow arrows meet the silhouette,
    ``beds`` are the rectangles (x0, y0, x1, y1) that beads may occupy -- a
    radial column's bed is an annulus and so contributes two of them --
    ``flow_symbol`` and ``flow_ticks`` name the coordinate that runs along the
    flow, and ``dispersion`` is the symbol of the flow-directional dispersion
    coefficient.
    """

    top: float
    bottom: float
    beds: tuple[tuple[float, float, float, float], ...]
    inlet_x: float = _LEFT
    outlet_x: float = _RIGHT
    flow_symbol: str = "z"
    flow_ticks: tuple[str, str] = ("0", r"L^{\b}")
    dispersion: str = r"D^{\mathrm{ax}}"
    draws_own_flow_arrows: bool = False
    has_wall: bool = True
    wall_x: float = _LEFT + 0.78 * (_RIGHT - _LEFT)
    bead_scale: float = 1.0  # shrinks the beads where the bed is cramped
    xmin: float = 0.0  # left edge of the frame, for annotations beside the body
    domain_at: tuple[float, float] | None = None  # where to write the bulk domain


def _cap_rx(ry: float) -> float:
    """Horizontal bulge of an end cap of half-height ``ry``."""
    return _PERSPECTIVE * ry


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
        beds=((_LEFT + 1.5, bottom + 0.3, _RIGHT - 0.3, top - 0.3),),
        inlet_x=_LEFT - rx,
        outlet_x=_RIGHT + rx,
    )


def _frustum(scene: Scene, model) -> _Body:
    """A conical frustum: the cross section grows along the flow."""
    del model
    inlet_ry, outlet_ry = 1.15, 1.95
    inlet_rx, outlet_rx = _cap_rx(inlet_ry), _cap_rx(outlet_ry)
    top, bottom = _AXIS_Y + outlet_ry, _AXIS_Y - outlet_ry
    left_centre, right_centre = (_LEFT, _AXIS_Y), (_RIGHT, _AXIS_Y)

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

    # The two radii are what distinguishes a frustum from a cylinder.
    scene.add(
        # The inlet radius points down, clear of the inlet boundary label.
        Arrow(left_centre, (_LEFT, _AXIS_Y - inlet_ry), Style.ANNOTATION),
        Label((_LEFT - 0.1, _AXIS_Y - 0.55 * inlet_ry), r"P^{0}", Anchor.RIGHT, size=0.9),
        Arrow(right_centre, (_RIGHT, _AXIS_Y + outlet_ry), Style.ANNOTATION),
        Label((_RIGHT + 0.12, _AXIS_Y + 0.55 * outlet_ry), r"P^{L^{\b}}", Anchor.LEFT, size=0.9),
    )
    return _Body(
        top=top,
        bottom=bottom,
        beds=((_LEFT + 1.4, _AXIS_Y - inlet_ry + 0.25, _RIGHT - 0.45, _AXIS_Y + inlet_ry - 0.25),),
        inlet_x=_LEFT - inlet_rx,
        outlet_x=_RIGHT + outlet_rx,
    )


def _generalized_radius(fraction: float) -> float:
    """Half-height of the generalized body at ``fraction`` along the flow."""
    # A venturi-like waist makes "the cross section varies smoothly" obvious.
    return 1.85 - 0.8 * math.exp(-(((fraction - 0.45) / 0.23) ** 2))


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
        Label((waist_x, _AXIS_Y - waist_ry - 0.1), "A(x)", Anchor.BELOW, size=0.9),
    )
    top = _AXIS_Y + max(radius for _, radius in wall)
    return _Body(
        top=top,
        bottom=2 * _AXIS_Y - top,
        beds=((_LEFT + 1.4, _AXIS_Y - 0.6, _RIGHT - 0.45, _AXIS_Y + 0.6),),
        bead_scale=0.85,
        inlet_x=_LEFT - inlet_rx,
        outlet_x=_RIGHT + outlet_rx,
        flow_symbol="x",
        flow_ticks=(r"x^{\mathrm{in}}", r"x^{\mathrm{out}}"),
        dispersion=r"D^{\mathrm{x}}",
    )


def _radial_cylinder(scene: Scene, model) -> _Body:
    """A hollow cylinder: the flow runs radially inwards into the channel."""
    del model
    rx, channel_ry = _cap_rx(_BODY_RY), 0.62
    channel_rx = _cap_rx(channel_ry)
    top, bottom = _AXIS_Y + _BODY_RY, _AXIS_Y - _BODY_RY
    left_centre, right_centre = (_LEFT, _AXIS_Y), (_RIGHT, _AXIS_Y)

    scene.add(
        Ellipse(left_centre, rx, _BODY_RY, Style.BULK, stroked=False),
        Rect((_LEFT, bottom), _RIGHT - _LEFT, 2 * _BODY_RY, Style.BULK, stroked=False),
        Ellipse(right_centre, rx, _BODY_RY, Style.BULK_FACE, stroked=False),
        # The collecting channel is empty and open at the near face.
        Ellipse(right_centre, channel_rx, channel_ry, Style.VOID, stroked=False),
        EllipseArc(left_centre, rx, _BODY_RY, 90, 270, Style.WALL),
        EllipseArc(left_centre, rx, _BODY_RY, -90, 90, Style.HIDDEN),
        Polyline(((_LEFT, top), (_RIGHT, top)), Style.WALL),
        Polyline(((_LEFT, bottom), (_RIGHT, bottom)), Style.WALL),
        # The channel wall runs the length of the column, hidden by the bed.
        Polyline(((_LEFT, _AXIS_Y + channel_ry), (_RIGHT, _AXIS_Y + channel_ry)), Style.HIDDEN),
        Polyline(((_LEFT, _AXIS_Y - channel_ry), (_RIGHT, _AXIS_Y - channel_ry)), Style.HIDDEN),
        Ellipse(right_centre, rx, _BODY_RY, Style.WALL, filled=False),
        Ellipse(right_centre, channel_rx, channel_ry, Style.WALL, filled=False),
    )

    # Feed enters over the whole outer surface and leaves through the channel.
    for fraction in (0.3, 0.52, 0.74):
        x = _LEFT + fraction * (_RIGHT - _LEFT)
        scene.add(
            Arrow((x, top + 1.0), (x, top + 0.12), Style.FLOW),
            Arrow((x, bottom - 1.0), (x, bottom - 0.12), Style.FLOW),
        )
    scene.add(
        Arrow((_RIGHT + rx + 0.15, _AXIS_Y), (_RIGHT + rx + 1.35, _AXIS_Y), Style.FLOW),
        Label((_LEFT + 0.3 * (_RIGHT - _LEFT) - 0.25, top + 0.55), r"\Gamma^{\mathrm{in}}", Anchor.RIGHT),
        Label((_RIGHT + rx + 1.45, _AXIS_Y + 0.22), r"\Gamma^{\mathrm{out}}", Anchor.ABOVE),
    )

    # The bed is the annulus, so beads go in the two bands around the channel.
    margin = 0.24
    return _Body(
        top=top + 1.0,
        bottom=bottom - 1.0,
        beds=(
            (_LEFT + 1.5, _AXIS_Y + channel_ry + margin, _RIGHT - 0.4, top - margin),
            (_LEFT + 1.5, bottom + margin, _RIGHT - 0.4, _AXIS_Y - channel_ry - margin),
        ),
        inlet_x=_LEFT - rx,
        outlet_x=_RIGHT + rx,
        flow_symbol=r"\rho",
        flow_ticks=(r"P^{\mathrm{inner}}", r"P^{\mathrm{outer}}"),
        dispersion=r"D^{\mathrm{rad}}",
        draws_own_flow_arrows=True,
        bead_scale=0.74,
        xmin=-1.6,
        # the channel is not part of the bed, so the label goes in the annulus
        domain_at=(_LEFT + 1.15, _AXIS_Y - channel_ry - 0.56),
    )


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
    blade_y, blade_half = bottom + 0.7, 0.7
    scene.add(
        Polyline(((centre_x, top + 0.95), (centre_x, blade_y)), Style.WALL),
        Polyline(((centre_x - blade_half, blade_y + 0.2), (centre_x + blade_half, blade_y - 0.2)), Style.WALL),
        Polyline(((centre_x - blade_half, blade_y - 0.2), (centre_x + blade_half, blade_y + 0.2)), Style.WALL),
        EllipseArc((centre_x, blade_y + 0.8), 1.3, 0.4, 15, 300, Style.ANNOTATION, arrow=True),
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

    return _Body(
        top=top + 0.95,
        bottom=bottom - rim_ry,
        beds=((centre_x - half_width + 0.4, blade_y + 1.35, centre_x + half_width - 0.4, top - rim_ry - 0.5),),
        draws_own_flow_arrows=True,
        has_wall=False,  # a homogeneous tank has no wall boundary condition
        bead_scale=0.78,
    )


#: Body builder per column geometry. The tank is selected by resolution.
_BODIES = {
    "Axial": _axial_cylinder,
    "Radial": _radial_cylinder,
    "Frustum": _frustum,
    "Generalized": _generalized,
}


# %% Annotation layers


def _flow_arrows(scene: Scene, body: _Body) -> None:
    """Draw the inlet and outlet arrows of a body with through-flow."""
    scene.add(
        Arrow((body.inlet_x - 1.45, _AXIS_Y), (body.inlet_x - 0.2, _AXIS_Y), Style.FLOW),
        Arrow((body.outlet_x + 0.2, _AXIS_Y), (body.outlet_x + 1.45, _AXIS_Y), Style.FLOW),
        Label((body.inlet_x - 0.8, _AXIS_Y + 0.22), r"\Gamma^{\mathrm{in}}", Anchor.ABOVE),
        Label((body.outlet_x + 0.8, _AXIS_Y + 0.22), r"\Gamma^{\mathrm{out}}", Anchor.ABOVE),
        Label((body.inlet_x - 0.8, _AXIS_Y - 0.28), "Flow", Anchor.BELOW, size=0.85, italic=False),
    )


def _wall_label(scene: Scene, body: _Body) -> None:
    """Point out the column wall from above the body."""
    scene.add(
        Arrow((body.wall_x, body.top + 0.8), (body.wall_x, body.top + 0.1), Style.ANNOTATION),
        Label((body.wall_x, body.top + 0.87), r"\Gamma^{\mathrm{wall}}", Anchor.ABOVE),
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


def _radial_flow_coordinate(scene: Scene, body: _Body) -> None:
    """Draw the radial coordinate of a radial flow column.

    The flow runs inwards, so unlike the axial geometries the coordinate is
    perpendicular to the body and is drawn on the body itself.
    """
    # Drawn beside the body rather than on it: the bed is full of beads, and
    # the two ticks line up with the channel wall and the column wall anyway.
    x = _LEFT - 1.15
    inner, outer = body.flow_ticks
    channel_y, wall_y = body.beds[0][1] - 0.24, body.beds[0][3] + 0.24
    tick = 0.13
    scene.add(
        Arrow((x, channel_y), (x, wall_y + 0.6), Style.AXIS),
        Polyline(((x - tick, channel_y), (x + tick, channel_y)), Style.AXIS),
        Polyline(((x - tick, wall_y), (x + tick, wall_y)), Style.AXIS),
        Label((x - tick - 0.08, channel_y), inner, Anchor.RIGHT, size=0.9),
        Label((x - tick - 0.08, wall_y), outer, Anchor.RIGHT, size=0.9),
        Label((x + 0.12, wall_y + 0.6), r"\rho", Anchor.ABOVE_RIGHT),
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


def _dispersion_inset(scene: Scene, body: _Body, dispersed: bool) -> None:
    """Sketch the concentration profile a 1D model transports.

    Pure convection moves a sharp front; dispersion smears it, which is the
    one feature of the transport that a drawing can show directly.
    """
    left, width = _LEFT + 0.6, 3.4
    base = body.top + 1.45
    height, front = 0.95, 0.45

    # The front travels to the right: high concentration behind it, zero ahead.
    if dispersed:
        profile = tuple(
            (
                left + fraction * width,
                base + height / (1.0 + math.exp((fraction - front) * 13.0)),
            )
            for fraction in (index / _SAMPLES for index in range(_SAMPLES + 1))
        )
    else:
        step = left + front * width
        profile = ((left, base + height), (step, base + height), (step, base), (left + width, base))

    scene.add(
        Arrow((left - 0.3, base), (left + width + 0.5, base), Style.AXIS),
        Arrow((left - 0.3, base - 0.12), (left - 0.3, base + height + 0.45), Style.AXIS),
        Label((left - 0.42, base + height * 0.6), r"c^{\b}", Anchor.RIGHT, size=0.85),
        Polyline(profile, Style.FLOW),
    )
    if dispersed:
        # The double arrow straddles the front, sitting clear above the curve.
        spread, middle = 0.6, left + front * width
        arrow_y = base + height + 0.3
        scene.add(
            Arrow((middle - spread, arrow_y), (middle + spread, arrow_y), Style.ANNOTATION, double=True),
            Label((middle, arrow_y + 0.14), body.dispersion, Anchor.ABOVE, size=0.85),
        )


# %% Particle types


#: Deterministic bead positions in the unit square, so that the rendered
#: figure stays byte-stable for the reference tests.
#: The heights run in a three-step cycle rather than alternating, so that with
#: two particle types the types do not end up sorted into rows.
_BEAD_SPOTS = tuple((0.06 + 0.127 * spot, (0.76, 0.26, 0.51)[spot % 3]) for spot in range(8))


def _particle_types(model) -> list:
    """Return the distinct particle types of ``model``, in display order."""
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
    if not types:
        return

    for spot, (fraction_x, fraction_y) in enumerate(_BEAD_SPOTS):
        index = spot % len(types)
        particle = types[index]
        radius = _bead_radius(index, len(types)) * body.bead_scale
        # Dividing by the type count keeps the band from tracking the type,
        # so an annular bed shows every particle type in both of its bands.
        x0, y0, x1, y1 = body.beds[(spot // len(types)) % len(body.beds)]
        # Clamped so that a thin band -- a radial column's annulus -- still
        # holds whole beads instead of letting them cross the wall.
        reach = radius * 1.34
        centre = (
            min(max(x0 + fraction_x * (x1 - x0), x0 + reach), x1 - reach),
            min(max(y0 + fraction_y * (y1 - y0), y0 + reach), y1 - reach),
        )
        # A film resistance is drawn as the liquid film it stands for; with
        # non-limiting film diffusion the bead sits directly in the bulk.
        if not getattr(particle, "nonlimiting_filmDiff", False):
            scene.add(Ellipse(centre, radius * 1.34, radius * 1.34, Style.FILM))
        scene.add(Ellipse(centre, radius, radius, Style.PARTICLE))


#: Legend entries laid out per row, and the horizontal pitch between them.
_LEGEND_COLUMNS, _LEGEND_PITCH, _LEGEND_ROW = 2, 5.6, 0.85


def _particle_legend(scene: Scene, model, baseline: float) -> float:
    """Explain the bead classes below the body; return the height used.

    The legend is what answers "which particle types are film-diffusion
    limited and which are in rapid equilibrium with the bulk".
    """
    types = _particle_types(model)
    if not types:
        return 0.0

    multiple = len(types) > 1
    for index, particle in enumerate(types):
        row, column = divmod(index, _LEGEND_COLUMNS)
        centre = (_LEFT + column * _LEGEND_PITCH, baseline - row * _LEGEND_ROW)
        radius = _bead_radius(index, len(types))
        limiting = not getattr(particle, "nonlimiting_filmDiff", False)
        if limiting:
            scene.add(Ellipse(centre, radius * 1.34, radius * 1.34, Style.FILM))
        scene.add(Ellipse(centre, radius, radius, Style.PARTICLE))

        prefix = f"j = {index + 1}: " if multiple else ""
        text = r"film diffusion k^{\mathrm{f}}" if limiting else r"c^{\b} = c^{\p}"
        scene.add(Label((centre[0] + 0.5, centre[1]), prefix + text, Anchor.LEFT, size=0.82))

    rows = -(-len(types) // _LEGEND_COLUMNS)
    return 0.35 + rows * _LEGEND_ROW


# %% Assembly


def build_column_scene(column_model) -> Scene:
    """Return the figure of ``column_model``'s interstitial volume.

    ``column_model`` is a :class:`~src.model_column.Column`; it is typed loosely
    to keep this module free of a dependency on the Streamlit-backed dataclass.
    """
    scene = Scene(name="column", xmin=0.0, ymin=0.0, xmax=13.0, ymax=0.0)

    is_tank = column_model.resolution == "0D"
    build_body = _tank if is_tank else _BODIES[column_model.column_type]
    body = build_body(scene, column_model)

    if not body.draws_own_flow_arrows:
        _flow_arrows(scene, body)
    if body.has_wall:
        _wall_label(scene, body)
    _beads(scene, body, column_model)
    if not is_tank:
        domain_at = body.domain_at or (body.beds[0][0] - 0.35, _AXIS_Y)
        scene.add(Label(domain_at, r"\Omega^{\b}", Anchor.CENTER, size=1.1))

    # Headroom for the wall label, whose superscript rises above its baseline.
    top = body.top + (1.65 if body.has_wall else 0.4)
    # The smearing of a front is a one-dimensional picture, so the inset is
    # drawn for the 1D models, where it is the whole of the transport story.
    if column_model.resolution == "1D":
        _dispersion_inset(scene, body, column_model.has_axial_dispersion)
        top = body.top + 4.3

    if getattr(column_model, "has_radial_coordinate", False):
        _radial_coordinate(scene, body)
    if getattr(column_model, "has_angular_coordinate", False):
        _angular_coordinate(scene, body)

    bottom = body.bottom - 0.5
    if column_model.column_type == "Radial":
        _radial_flow_coordinate(scene, body)
    elif not is_tank:
        bottom = body.bottom - 1.1
        _flow_coordinate(scene, body, bottom)
        bottom -= 0.8

    bottom -= _particle_legend(scene, column_model, bottom - 0.8)

    scene.xmin = body.xmin
    scene.ymin, scene.ymax = bottom - 0.3, top
    scene.title = "Interstitial volume of the " + column_model.model_name()
    return scene


__all__ = ["build_column_scene"]
