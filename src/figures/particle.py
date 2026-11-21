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


def _particle_body(scene: Scene, particle) -> None:
    """Draw the particle, its film and its core."""
    geometry = getattr(particle, "geometry", "Sphere")
    limiting = not getattr(particle, "nonlimiting_filmDiff", False)

    if limiting:
        _shape(scene, geometry, _FILM, Style.FILM)
    _shape(scene, geometry, 0.0, Style.PORE)
    if getattr(particle, "has_core", False):
        _core_shape(scene, geometry)

    x, y = _CENTRE
    grown = _surface_x(geometry, _FILM if limiting else 0.0)
    scene.add(Label((x - 0.52 * _RADIUS, y - 0.58 * _RADIUS), r"\Omega^{\p}", Anchor.CENTER, size=1.05))
    # A leader makes clear that the label names the surface, not the film.
    scene.add(
        Polyline((_surface_point(geometry), (grown + 0.45, y + 1.25 * _RADIUS)), Style.ANNOTATION),
        Label((grown + 0.5, y + 1.25 * _RADIUS), r"\Gamma^{\p}", Anchor.LEFT, size=0.95),
    )


def _transport(scene: Scene, particle) -> None:
    """Draw the transport that runs along the particle's radial coordinate."""
    geometry = getattr(particle, "geometry", "Sphere")
    resolved = getattr(particle, "resolution", "1D") == "1D"
    surface = _surface_x(geometry)
    x, y = _CENTRE

    if not getattr(particle, "nonlimiting_filmDiff", False):
        # Film diffusion: the flux crosses the film and stops at the surface.
        # The label rides above the arrow so that wide geometries stay clear
        # of the phase chain to the right.
        outside = _surface_x(geometry, _FILM) + 0.75
        lands = _edge_x(geometry, 0.95) + 0.04
        scene.add(
            Arrow((outside, y + 0.95), (lands, y + 0.95), Style.FLOW),
            Label((0.5 * (outside + surface), y + 1.05), r"k^{\mathrm{f}}", Anchor.ABOVE, size=0.95),
        )

    if resolved:
        # Pore diffusion runs inwards along r; surface diffusion runs beside it
        # in the bound phase, which is why the two arrows are drawn in parallel.
        inner = _core_x(geometry) + 0.1 if getattr(particle, "has_core", False) else x + 0.15
        for offset, symbol, shown in (
            (-0.15, r"D^{\p}", True),
            (-1.0, r"D^{\s}", getattr(particle, "has_surfDiff", False)),
        ):
            if not shown:
                continue
            edge = _edge_x(geometry, offset) - 0.18
            scene.add(
                Arrow((edge, y + offset), (inner, y + offset), Style.ANNOTATION, double=True),
                Label((0.5 * (edge + inner), y + offset + 0.1), symbol, Anchor.ABOVE, size=0.9),
            )


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
        Label((surface, baseline - tick - 0.05), r"P^{\p}", Anchor.BELOW, size=0.9),
        Label((surface + 0.7, baseline + 0.1), "r", Anchor.ABOVE_RIGHT),
    )
    if getattr(particle, "has_core", False):
        # Written above the axis: on a slab the two ticks sit close together,
        # and two labels below would run into each other.
        core_x = _core_x(geometry)
        scene.add(
            Polyline(((core_x, baseline - tick), (core_x, baseline + tick)), Style.AXIS),
            Label((core_x, baseline + tick + 0.02), r"P^{\mathrm{core}}", Anchor.ABOVE, size=0.9),
        )
    else:
        scene.add(Label((x, baseline - tick - 0.05), "0", Anchor.BELOW, size=0.9, italic=False))


def _phase_box(scene: Scene, row: float, symbol: str, caption: str, style: Style) -> None:
    """Draw one phase of the chain: its concentration and what it is."""
    corner, width, height = (_CHAIN_X - _CHAIN_HALF, row - _BOX_HALF), 2 * _CHAIN_HALF, 2 * _BOX_HALF
    scene.add(
        Rect(corner, width, height, style, stroked=False, radius=0.16),
        Rect(corner, width, height, Style.WALL, filled=False, radius=0.16),
        Label((_CHAIN_X - _CHAIN_HALF + 0.42, row), symbol, Anchor.LEFT, size=1.0),
        Label((_CHAIN_X + _CHAIN_HALF - 0.42, row), caption, Anchor.RIGHT, size=0.78, italic=False),
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
