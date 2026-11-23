"""Build the scene of the particle figure.

The figure has two panels, which answer two different questions:

* the **particle** on the left shows *where* things happen -- the geometry, the
  optional impermeable core, the liquid film, and the transport that runs along
  the particle's radial coordinate;
* the **phase chain** on the right shows *what* happens -- which concentrations
  exist and how they are coupled, from the bulk through the film to the pore
  liquid and, if the model binds, on to the bound phase.

Every element is switched by one attribute of the particle, so the figure
cannot drift away from the equations: a model without surface diffusion simply
does not get the surface diffusion arrow.
"""

from __future__ import annotations

import math

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

# %% Layout, in scene units

_CENTRE = (3.3, 4.5)  # centre of the particle
_RADIUS = 2.05  # how far the particle reaches from its centre
_SLAB = 0.68  # half-thickness of a slab, relative to _RADIUS
_CAP = 0.3  # how far a cylinder's end face bulges, relative to _RADIUS
_FILM = 0.42  # thickness of the liquid film around it
_CORE = 0.42  # core radius, as a fraction of the particle radius

_PORE_STEP = 0.58  # pitch of the pores drawn through the solid matrix
_PORE_SIZE = 0.15  # their radius, before the small variation below
_CHANNEL_HALF = 0.68  # half-height of the one pore the transport is drawn in
_CHANNEL_RISE = 0.3  # how far above the centre that pore runs
_BOUND = 0.13  # radius of a bound molecule on a pore wall
#: Fixed offsets and size steps, so the matrix does not read as a lattice and
#: the rendered figure still stays byte-stable for the reference tests.
_JITTER = (0.0, 0.05, -0.04, 0.03, -0.054, 0.042, -0.02)

_CHAIN_X = 9.9  # centre of the phase chain
_CHAIN_HALF = 1.95  # half-width of a phase box
_BOX_HALF = 0.46  # half-height of a phase box
_ROWS = (7.0, 4.5, 2.0)  # bulk, pore liquid, bound phase


def _shape(scene: Scene, geometry: str, grow: float, style: Style) -> None:
    """Add the particle's outline, grown by ``grow``, in ``style``.

    One function for all three geometries keeps the film, the body and the core
    in step: each is the same shape at a different size.  All three reach the
    same distance up and down, so that the figure keeps its proportions
    whichever geometry the model uses, and all three are drawn with the
    particle's radial coordinate running to the right.
    """
    x, y = _CENTRE
    radius = _RADIUS + grow

    if geometry == "Slab":
        # Seen edge-on: resolved across its thickness, which is what the radial
        # coordinate runs along, and extended in the other direction.
        half_width = _SLAB * _RADIUS + grow
        scene.add(Rect((x - half_width, y - radius), 2 * half_width, 2 * radius, style))
    elif geometry == "Cylinder":
        # A squat cylinder standing on end, so that its radius -- the direction
        # the model resolves -- is the horizontal one, as for the other two.
        cap_ry = _CAP * radius
        half_length = radius - cap_ry
        scene.add(
            Rect((x - radius, y - half_length), 2 * radius, 2 * half_length, style, stroked=False),
            Ellipse((x, y - half_length), radius, cap_ry, style, stroked=False),
            Ellipse((x, y + half_length), radius, cap_ry, style, stroked=False),
            # Only the silhouette is stroked: the hidden half of the base would
            # otherwise show through whatever is drawn on top of this shape.
            EllipseArc((x, y - half_length), radius, cap_ry, 180, 360, style),
            Polyline(((x - radius, y - half_length), (x - radius, y + half_length)), style),
            Polyline(((x + radius, y - half_length), (x + radius, y + half_length)), style),
            Ellipse((x, y + half_length), radius, cap_ry, style, filled=False),
        )
    else:  # Sphere
        scene.add(Ellipse((x, y), radius, radius, style))


def _surface_point(geometry: str) -> tuple[float, float]:
    """Return a point on the particle surface for the boundary label's leader."""
    x, y = _CENTRE
    if geometry == "Slab":
        return x + _SLAB * _RADIUS, y + 0.72 * _RADIUS
    if geometry == "Cylinder":
        return x + _RADIUS, y + 0.3 * _RADIUS
    return x + 0.73 * _RADIUS, y + 0.68 * _RADIUS


def _surface_x(geometry: str, grow: float = 0.0) -> float:
    """Return the x where the radial coordinate leaves the particle."""
    x = _CENTRE[0]
    if geometry == "Slab":
        return x + _SLAB * _RADIUS + grow
    return x + _RADIUS + grow


def _edge_x(geometry: str, offset: float, grow: float = 0.0) -> float:
    """Return the x of the outline ``offset`` above or below the centre.

    Arrows drawn inside the particle end here rather than at its widest point,
    so that none of them pokes through a curved outline.
    """
    x = _CENTRE[0]
    radius = _RADIUS + grow
    height = abs(offset)

    if geometry == "Slab":
        return x + _SLAB * _RADIUS + grow
    if geometry == "Cylinder":
        cap_ry, half_length = _CAP * radius, radius - _CAP * radius
        if height <= half_length:
            return x + radius
        rise = min((height - half_length) / cap_ry, 1.0)
        return x + radius * math.sqrt(max(1.0 - rise * rise, 0.0))
    return x + math.sqrt(max(radius * radius - height * height, 0.0))


def _core_x(geometry: str) -> float:
    """Return the x where the impermeable core ends."""
    x = _CENTRE[0]
    if geometry == "Slab":
        return x + _CORE * _SLAB * _RADIUS
    return x + _CORE * _RADIUS


def _core_shape(scene: Scene, geometry: str) -> None:
    """Draw the impermeable core, which is the particle thinned radially.

    Only the resolved direction shrinks: a slab keeps its extent and a cylinder
    its length, because the core of either is a thinner body of the same kind.
    """
    x, y = _CENTRE
    radius = _CORE * _RADIUS

    if geometry == "Slab":
        half_width = _CORE * _SLAB * _RADIUS
        scene.add(Rect((x - half_width, y - _RADIUS), 2 * half_width, 2 * _RADIUS, Style.CORE))
    elif geometry == "Cylinder":
        cap_ry, half_length = _CAP * radius, _RADIUS - _CAP * _RADIUS
        scene.add(
            Rect((x - radius, y - half_length), 2 * radius, 2 * half_length, Style.CORE, stroked=False),
            Ellipse((x, y - half_length), radius, cap_ry, Style.CORE, stroked=False),
            EllipseArc((x, y - half_length), radius, cap_ry, 180, 360, Style.CORE),
            Polyline(((x - radius, y - half_length), (x - radius, y + half_length)), Style.CORE),
            Polyline(((x + radius, y - half_length), (x + radius, y + half_length)), Style.CORE),
            Ellipse((x, y + half_length), radius, cap_ry, Style.CORE),
        )
    else:  # Sphere
        scene.add(Ellipse((x, y), radius, radius, Style.CORE))


#: Where the domain label sits, as an offset from the centre in scene units.
_DOMAIN_AT = (-0.52 * _RADIUS, -0.62 * _RADIUS)


def _channel(geometry: str, has_core: bool) -> tuple[float, float, float]:
    """Return (inner x, lower y, upper y) of the pore the transport runs in."""
    x, y = _CENTRE
    inner = _core_x(geometry) if has_core else x
    return inner, y + _CHANNEL_RISE - _CHANNEL_HALF, y + _CHANNEL_RISE + _CHANNEL_HALF


def _pores(scene: Scene, geometry: str, has_core: bool, reserved: tuple) -> None:
    """Fill the solid matrix with pores, around whatever else is drawn in it.

    The pore space is what ``c^{\\p}`` lives in and what ``D^{\\p}`` runs
    through, so the particle has to show it; without it a bound phase on the
    pore wall, and surface diffusion along that wall, have nowhere to be.
    """
    x, y = _CENTRE
    core = _CORE * _RADIUS if has_core else 0.0
    reach = 6

    for row in range(-reach, reach + 1):
        for column in range(-reach, reach + 1):
            radius = _PORE_SIZE + 0.035 * ((row * 2 + column) % 3)
            centre = (
                x + (column + 0.5 * (row % 2)) * _PORE_STEP + _JITTER[(row * 5 + column) % 7],
                y + row * _PORE_STEP + _JITTER[(row * 3 + column * 2) % 7],
            )
            if not _inside(geometry, centre, radius + 0.12):
                continue
            if math.hypot(centre[0] - x, centre[1] - y) < core + radius + 0.1:
                continue
            if any(_overlaps(centre, radius + 0.08, box) for box in reserved):
                continue
            scene.add(Ellipse(centre, radius, radius, Style.PORE))


def _inside(geometry: str, point: Point, margin: float) -> bool:
    """Is a disc of radius ``margin`` around ``point`` within the particle?"""
    offset = point[1] - _CENTRE[1]
    if abs(offset) + margin > _RADIUS:
        return False
    return abs(point[0] - _CENTRE[0]) + margin <= _edge_x(geometry, offset) - _CENTRE[0]


def _overlaps(centre: Point, radius: float, box: tuple[float, float, float, float]) -> bool:
    """Does a disc meet an axis-aligned box given as (x0, y0, x1, y1)?"""
    x0, y0, x1, y1 = box
    return x0 - radius < centre[0] < x1 + radius and y0 - radius < centre[1] < y1 + radius


def _particle_body(scene: Scene, particle) -> None:
    """Draw the particle, its porous matrix, its film and its core."""
    geometry = getattr(particle, "geometry", "Sphere")
    limiting = not getattr(particle, "nonlimiting_filmDiff", False)
    has_core = getattr(particle, "has_core", False)

    if limiting:
        _shape(scene, geometry, _FILM, Style.FILM)
    _shape(scene, geometry, 0.0, Style.MATRIX)
    if has_core:
        _core_shape(scene, geometry)

    x, y = _CENTRE
    inner, low, high = _channel(geometry, has_core)
    domain = (x + _DOMAIN_AT[0], y + _DOMAIN_AT[1])
    _pores(
        scene,
        geometry,
        has_core,
        reserved=(
            # The drawn pore, the domain label, and the caption under the pore.
            (inner - _CHANNEL_HALF - 0.1, low - 0.1, _edge_x(geometry, 0.0), high + 0.1),
            (domain[0] - 0.5, domain[1] - 0.32, domain[0] + 0.5, domain[1] + 0.32),
            (inner + 0.3, low - 0.62, inner + 1.5, low),
        ),
    )

    grown = _surface_x(geometry, _FILM if limiting else 0.0)
    scene.add(Label(domain, r"\Omega^{\p}", Anchor.CENTER, size=1.05))
    # A leader makes clear that the label names the surface, not the film.
    scene.add(
        Polyline((_surface_point(geometry), (grown + 0.45, y + 1.25 * _RADIUS)), Style.ANNOTATION),
        Label((grown + 0.5, y + 1.25 * _RADIUS), r"\Gamma^{\p}", Anchor.LEFT, size=0.95),
    )


def _transport(scene: Scene, particle) -> None:
    """Draw one pore enlarged, and the transport that runs along it.

    Everything the particle model resolves happens here: the flux arrives
    through the film at the pore mouth, moves along the pore in the liquid,
    binds to the pore wall, and -- where the model allows it -- moves along
    that wall in the bound phase.  Drawing it in a pore rather than across a
    blank disc is what keeps ``D^{\\p}`` and ``D^{\\s}`` apart.
    """
    geometry = getattr(particle, "geometry", "Sphere")
    resolved = getattr(particle, "resolution", "1D") == "1D"
    y = _CENTRE[1]
    inner, low, high = _channel(geometry, getattr(particle, "has_core", False))

    # The pore opens at the surface, so its outer end follows the outline; its
    # inner end is rounded, because it carries on into the rest of the matrix.
    steps = 12
    mouth = tuple(
        (_edge_x(geometry, high - y - (high - low) * step / steps), high - (high - low) * step / steps)
        for step in range(steps + 1)
    )
    if getattr(particle, "has_core", False):
        # An impermeable core closes the pore off with a flat wall.
        blind = ((inner, low), (inner, high))
    else:
        blind = tuple(
            (
                inner - _CHANNEL_HALF * math.sin(math.pi * step / steps),
                y + _CHANNEL_RISE - _CHANNEL_HALF * math.cos(math.pi * step / steps),
            )
            for step in range(steps + 1)
        )
    scene.add(
        Polyline(mouth + blind, Style.PORE, closed=True, filled=True, stroked=False),
        Polyline(blind + mouth[:1], Style.WALL),
        Polyline((blind[0], mouth[-1]), Style.WALL),
    )

    if not getattr(particle, "nonlimiting_filmDiff", False):
        # Film diffusion: the flux crosses the film into the mouth of the pore.
        offset = _CHANNEL_RISE
        outside = _surface_x(geometry, _FILM) + 0.8
        lands = _edge_x(geometry, offset) + 0.04
        scene.add(
            Arrow((outside, y + offset), (lands, y + offset), Style.FLOW),
            Label((0.5 * (outside + lands), y + offset + 0.12), r"k^{\mathrm{f}}", Anchor.ABOVE, size=0.95),
        )

    wall_y = low + _BOUND
    if resolved:
        # Pore diffusion runs along the liquid in the pore; surface diffusion
        # runs along the wall, through the bound molecules sitting on it.
        pore_y = y + _CHANNEL_RISE + 0.06
        outer = _edge_x(geometry, pore_y - y) - 0.22
        scene.add(
            Arrow((outer, pore_y), (inner + 0.2, pore_y), Style.ANNOTATION, double=True),
            Label((0.5 * (inner + outer), pore_y + 0.1), r"D^{\p}", Anchor.ABOVE, size=0.9),
        )
        if getattr(particle, "has_surfDiff", False):
            outer = _edge_x(geometry, low - y) - 0.26
            scene.add(
                Arrow((outer, wall_y), (inner + 0.2, wall_y), Style.ANNOTATION, double=True),
                Label((0.5 * (inner + outer), low - 0.12), r"D^{\s}", Anchor.BELOW, size=0.9),
            )

    if getattr(particle, "has_binding", False):
        # The bound phase sits on the pore wall, which is the only place it
        # can be; drawn last so the molecules ride on the arrow along it.
        place = inner + 0.42
        while place < _edge_x(geometry, low - y) - 0.3:
            scene.add(Ellipse((place, wall_y), _BOUND, _BOUND, Style.SOLID))
            place += 0.44


def _radial_coordinate(scene: Scene, particle) -> None:
    """Draw the coordinate a 1D particle resolves, below the particle."""
    if getattr(particle, "resolution", "1D") != "1D":
        return

    geometry = getattr(particle, "geometry", "Sphere")
    x, y = _CENTRE
    baseline = y - _RADIUS - (_FILM if not getattr(particle, "nonlimiting_filmDiff", False) else 0.0) - 0.85
    surface = _surface_x(geometry)
    tick = 0.13

    scene.add(
        Arrow((x, baseline), (surface + 0.7, baseline), Style.AXIS),
        Polyline(((x, baseline - tick), (x, baseline + tick)), Style.AXIS),
        Polyline(((surface, baseline - tick), (surface, baseline + tick)), Style.AXIS),
        Label((surface, baseline - tick - 0.05), r"R^{\mathrm{p}}", Anchor.BELOW, size=0.9),
        Label((surface + 0.7, baseline + 0.1), "r", Anchor.ABOVE_RIGHT),
    )
    if getattr(particle, "has_core", False):
        # Written above the axis: on a slab the two ticks sit close together,
        # and two labels below would run into each other.
        core_x = _core_x(geometry)
        scene.add(
            Polyline(((core_x, baseline - tick), (core_x, baseline + tick)), Style.AXIS),
            Label((core_x, baseline + tick + 0.02), r"R^{\mathrm{c}}", Anchor.ABOVE, size=0.9),
        )
    else:
        scene.add(Label((x, baseline - tick - 0.05), "0", Anchor.BELOW, size=0.9, italic=False))


def _phase_box(scene: Scene, row: float, symbol: str, caption: str, style: Style) -> None:
    """Draw one phase of the chain: its concentration and what it is."""
    corner, width, height = (_CHAIN_X - _CHAIN_HALF, row - _BOX_HALF), 2 * _CHAIN_HALF, 2 * _BOX_HALF
    scene.add(
        Rect(corner, width, height, style, stroked=False, radius=0.16),
        Rect(corner, width, height, Style.WALL, filled=False, radius=0.16),
        # The two sit on opposite edges, so the padding is also what keeps them
        # apart: the widest caption, "bound phase", all but touches the symbol.
        Label((_CHAIN_X - _CHAIN_HALF + 0.3, row), symbol, Anchor.LEFT, size=1.0),
        Label((_CHAIN_X + _CHAIN_HALF - 0.3, row), caption, Anchor.RIGHT, size=0.78, italic=False),
    )


def _phase_chain(scene: Scene, particle) -> None:
    """Draw the phases of the model and how they are coupled."""
    bulk, pore, solid = _ROWS
    binds = getattr(particle, "has_binding", False)

    _phase_box(scene, bulk, r"c^{\b}", "bulk liquid", Style.BULK)
    _phase_box(scene, pore, r"c^{\p}", "pore liquid", Style.PORE)
    if binds:
        _phase_box(scene, solid, r"c^{\s}", "bound phase", Style.SOLID)

    if getattr(particle, "nonlimiting_filmDiff", False):
        # No film resistance: the two concentrations are in equilibrium, which
        # is drawn as an identity rather than as a transport step.
        middle = 0.5 * (bulk + pore)
        for offset in (-0.09, 0.09):
            scene.add(
                Polyline(
                    ((_CHAIN_X - 0.22, middle + offset), (_CHAIN_X + 0.22, middle + offset)),
                    Style.ANNOTATION,
                )
            )
        scene.add(Label((_CHAIN_X - 0.4, middle), "no film resistance", Anchor.RIGHT, size=0.78, italic=False))
    else:
        scene.add(
            Arrow((_CHAIN_X, bulk - _BOX_HALF - 0.1), (_CHAIN_X, pore + _BOX_HALF + 0.1), Style.FLOW, double=True),
            Label((_CHAIN_X + 0.25, 0.5 * (bulk + pore)), r"k^{\mathrm{f}}", Anchor.LEFT, size=0.9),
            Label((_CHAIN_X - 0.25, 0.5 * (bulk + pore)), "film diffusion", Anchor.RIGHT, size=0.78, italic=False),
        )

    if binds:
        scene.add(
            Arrow((_CHAIN_X, pore - _BOX_HALF - 0.1), (_CHAIN_X, solid + _BOX_HALF + 0.1), Style.FLOW, double=True),
            Label((_CHAIN_X + 0.25, 0.5 * (pore + solid)), _binding_symbol(particle), Anchor.LEFT, size=0.9),
            Label(
                (_CHAIN_X - 0.25, 0.5 * (pore + solid)),
                _binding_caption(particle),
                Anchor.RIGHT,
                size=0.78,
                italic=False,
            ),
        )


def _binding_symbol(particle) -> str:
    """Return the symbol written next to the binding arrow."""
    if getattr(particle, "req_binding", False):
        return r"f^{\mathrm{bind}} = 0"
    return r"f^{\mathrm{bind}}"


def _binding_caption(particle) -> str:
    """Name the binding regime, which decides whether the step is resolved."""
    return "rapid equilibrium" if getattr(particle, "req_binding", False) else "binding"


def build_particle_scene(particle, name: str = "particle") -> Scene:
    """Return the figure of one particle type.

    ``particle`` is a :class:`~src.model_particle.Particle`; it is typed loosely
    so that this module does not depend on the Streamlit-backed dataclass.
    """
    scene = Scene(name=name, xmin=0.0, ymin=0.25, xmax=13.0, ymax=7.8, title="Particle model")

    _particle_body(scene, particle)
    _transport(scene, particle)
    _radial_coordinate(scene, particle)
    _phase_chain(scene, particle)

    if getattr(particle, "resolution", "1D") != "1D":
        # A homogeneous particle has no interior coordinate to resolve.
        scene.add(
            Label(
                (_CENTRE[0], _CENTRE[1] - _RADIUS - 0.75), "homogeneous particle", Anchor.BELOW, size=0.85, italic=False
            )
        )

    return scene


__all__ = ["build_particle_scene"]
