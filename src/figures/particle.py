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
_RADIUS = 2.05  # its outer radius (half-thickness for a slab)
_FILM = 0.42  # thickness of the liquid film around it
_CORE = 0.42  # core radius, as a fraction of the particle radius

_CHAIN_X = 9.9  # centre of the phase chain
_CHAIN_HALF = 1.95  # half-width of a phase box
_BOX_HALF = 0.46  # half-height of a phase box
_ROWS = (7.0, 4.5, 2.0)  # bulk, pore liquid, bound phase


def _shape(scene: Scene, geometry: str, grow: float, style: Style, silhouette_only: bool = False) -> None:
    """Add the particle's outline, grown by ``grow``, in ``style``.

    One function for all three geometries keeps the film, the body and the core
    in step: each is the same shape at a different size.  ``silhouette_only``
    suppresses the interior edges, which is what the film layer wants: its far
    cap would otherwise show through the particle drawn on top of it.
    """
    x, y = _CENTRE
    radius = _RADIUS + grow

    if geometry == "Slab":
        # Seen edge-on: resolved across its thickness, which is the half-width
        # the radial coordinate runs along, and extended in the other direction.
        half_height = radius * 1.3
        scene.add(Rect((x - radius, y - half_height), 2 * radius, 2 * half_height, style))
    elif geometry == "Cylinder":
        # Drawn with a short extrusion, so it is not taken for a sphere.
        depth, cap_rx = 0.42 * radius, 0.34 * radius
        scene.add(
            Rect((x - depth, y - radius), 2 * depth, 2 * radius, style, stroked=False),
            Ellipse((x - depth, y), cap_rx, radius, style, stroked=not silhouette_only),
            Ellipse((x + depth, y), cap_rx, radius, style, stroked=not silhouette_only),
            Polyline(((x - depth, y + radius), (x + depth, y + radius)), style),
            Polyline(((x - depth, y - radius), (x + depth, y - radius)), style),
        )
        if silhouette_only:
            scene.add(
                EllipseArc((x - depth, y), cap_rx, radius, 90, 270, style),
                EllipseArc((x + depth, y), cap_rx, radius, -90, 90, style),
            )
    else:  # Sphere
        scene.add(Ellipse((x, y), radius, radius, style))


def _surface_point(geometry: str) -> tuple[float, float]:
    """Return a point on the particle surface for the boundary label's leader."""
    x, y = _CENTRE
    if geometry == "Slab":
        return x + 0.5 * _RADIUS, y + 1.3 * _RADIUS
    if geometry == "Cylinder":
        return x + 0.3 * _RADIUS, y + _RADIUS
    return x + 0.73 * _RADIUS, y + 0.68 * _RADIUS


def _surface_x(geometry: str, grow: float = 0.0) -> float:
    """Return the x where the radial coordinate leaves the particle."""
    x = _CENTRE[0]
    radius = _RADIUS + grow
    if geometry == "Slab":
        return x + radius
    if geometry == "Cylinder":
        return x + 0.42 * radius + 0.34 * radius
    return x + radius


def _particle_body(scene: Scene, particle) -> None:
    """Draw the particle, its film and its core."""
    geometry = getattr(particle, "geometry", "Sphere")
    limiting = not getattr(particle, "nonlimiting_filmDiff", False)

    if limiting:
        _shape(scene, geometry, _FILM, Style.FILM, silhouette_only=True)
    _shape(scene, geometry, 0.0, Style.PORE)
    if getattr(particle, "has_core", False):
        _shape(scene, geometry, -(1.0 - _CORE) * _RADIUS, Style.CORE)

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
        # Film diffusion: the flux crosses the film into the particle. The
        # label rides above the arrow so that wide geometries stay clear of
        # the phase chain to the right.
        outside = _surface_x(geometry, _FILM) + 0.75
        scene.add(
            Arrow((outside, y + 0.95), (surface - 0.15, y + 0.95), Style.FLOW),
            Label((0.5 * (outside + surface), y + 1.05), r"k^{\mathrm{f}}", Anchor.ABOVE, size=0.95),
        )

    if resolved:
        # Pore diffusion runs inwards along r; surface diffusion runs beside it
        # in the bound phase, which is why the two arrows are drawn in parallel.
        inner = x + _CORE * _RADIUS + 0.1 if getattr(particle, "has_core", False) else x + 0.15
        scene.add(
            Arrow((surface - 0.25, y - 0.15), (inner, y - 0.15), Style.ANNOTATION, double=True),
            Label(((surface + inner) / 2, y - 0.05), r"D^{\p}", Anchor.ABOVE, size=0.9),
        )
        if getattr(particle, "has_surfDiff", False):
            scene.add(
                Arrow((surface - 0.25, y - 1.0), (inner, y - 1.0), Style.ANNOTATION, double=True),
                Label(((surface + inner) / 2, y - 0.9), r"D^{\s}", Anchor.ABOVE, size=0.9),
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
        core_x = x + _CORE * _RADIUS
        scene.add(
            Polyline(((core_x, baseline - tick), (core_x, baseline + tick)), Style.AXIS),
            Label((core_x, baseline - tick - 0.05), r"P^{\mathrm{core}}", Anchor.BELOW, size=0.9),
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
        for offset in (-0.07, 0.07):
            scene.add(
                Polyline(
                    ((_CHAIN_X - 0.3, middle + offset), (_CHAIN_X + 0.3, middle + offset)),
                    Style.ANNOTATION,
                )
            )
        scene.add(Label((_CHAIN_X + 0.5, middle), r"c^{\b} = c^{\p}", Anchor.LEFT, size=0.85))
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
