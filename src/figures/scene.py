"""Backend-independent description of a model figure.

A :class:`Scene` is a flat list of geometric primitives placed in a
right-handed coordinate frame (``x`` to the right, ``y`` **up**) and carrying a
*semantic* style name rather than a literal colour.  Both properties exist so
that a second renderer can be added without touching the figure builders:

* the y-flip that SVG needs happens once, in :mod:`src.figures.svg`, and a TikZ
  renderer (whose frame already points up) simply omits it;
* colours are resolved per backend -- to CSS custom properties in the app, to
  TikZ style names in a LaTeX export.

Scene units are roughly centimetres: a column figure is about 13 units wide.
Builders live in :mod:`src.figures.column` and :mod:`src.figures.particle`;
nothing in this module imports Streamlit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

Point = tuple[float, float]


class Style(str, Enum):
    """Semantic role of a primitive, resolved to colours by the renderer."""

    BULK = "bulk"  # interstitial (bulk) liquid volume
    BULK_FACE = "bulk-face"  # the lit end cap of a body, a shade lighter
    VOID = "void"  # space outside the bed, e.g. a radial column's channel
    WALL = "wall"  # column wall and other hard outlines
    HIDDEN = "hidden"  # edges behind the body, drawn dashed
    PARTICLE = "particle"  # particle (bead) body
    MATRIX = "matrix"  # the solid skeleton a porous particle is made of
    PORE = "pore"  # pore liquid inside a particle
    SOLID = "solid"  # bound (solid) phase
    CORE = "core"  # impermeable particle core
    FILM = "film"  # liquid film around a particle
    FLOW = "flow"  # convective flow arrows
    AXIS = "axis"  # resolved coordinate axes
    ANNOTATION = "annotation"  # labels, leader lines, brackets


class Anchor(str, Enum):
    """Placement of a label relative to its anchor point."""

    CENTER = "center"
    LEFT = "left"  # text starts at the point
    RIGHT = "right"  # text ends at the point
    ABOVE = "above"
    BELOW = "below"
    ABOVE_LEFT = "above left"
    ABOVE_RIGHT = "above right"
    BELOW_LEFT = "below left"
    BELOW_RIGHT = "below right"


@dataclass(frozen=True)
class Rect:
    """Axis-aligned rectangle given by its lower-left corner, width and height."""

    corner: Point
    width: float
    height: float
    style: Style
    filled: bool = True
    stroked: bool = True
    radius: float = 0.0


@dataclass(frozen=True)
class Ellipse:
    """Full ellipse. A circle is an ellipse with ``rx == ry``."""

    center: Point
    rx: float
    ry: float
    style: Style
    filled: bool = True
    stroked: bool = True


@dataclass(frozen=True)
class EllipseArc:
    """Part of an ellipse, from ``start_angle`` to ``end_angle`` in degrees.

    Angles are measured counter-clockwise from the positive x-axis, as in TikZ's
    ``arc[start angle=..., end angle=...]``, so the same numbers serve both
    backends.  Arcs are never filled.  ``arrow`` puts a head at the end angle,
    which is how the angular coordinate of a 3D model is drawn.
    """

    center: Point
    rx: float
    ry: float
    start_angle: float
    end_angle: float
    style: Style
    arrow: bool = False


@dataclass(frozen=True)
class Polyline:
    """Straight-line path through ``points``; closed paths may be filled.

    ``arrow`` puts a head on the last point, which is how a flow path that
    meanders through the packing is drawn.
    """

    points: tuple[Point, ...]
    style: Style
    closed: bool = False
    filled: bool = False
    stroked: bool = True
    arrow: bool = False


@dataclass(frozen=True)
class Arrow:
    """Straight arrow from ``start`` to ``end``.

    ``double`` puts a head on both ends, which is how dispersion (as opposed to
    convection) and reversible binding are drawn.
    """

    start: Point
    end: Point
    style: Style
    double: bool = False


@dataclass(frozen=True)
class Label:
    """Text placed at ``position``.

    ``text`` is written with the same ``\\b \\p \\s \\l`` shorthands as the
    equations and is expanded by :mod:`src.figures.labels`, so a figure follows
    the selected parameter format just like the equations do.
    """

    position: Point
    text: str
    anchor: Anchor = Anchor.CENTER
    style: Style = Style.ANNOTATION
    size: float = 1.0  # multiple of the scene's base font size
    italic: bool = True


Primitive = Rect | Ellipse | EllipseArc | Polyline | Arrow | Label


@dataclass
class Scene:
    """A figure: its primitives, its extent and the id prefix of its SVG defs.

    ``name`` must be unique among the figures shown on one page because the
    renderer derives arrow-head marker ids from it.  Deriving them from the name
    rather than from a counter keeps the rendered output deterministic, which is
    what lets the reference tests compare it byte for byte.
    """

    name: str
    xmin: float
    ymin: float
    xmax: float
    ymax: float
    items: list[Primitive] = field(default_factory=list)
    title: str = ""
    font_size: float = 0.52

    def add(self, *primitives: Primitive) -> None:
        """Append primitives in drawing order (later ones paint on top)."""
        self.items.extend(primitives)

    @property
    def width(self) -> float:
        return self.xmax - self.xmin

    @property
    def height(self) -> float:
        return self.ymax - self.ymin


__all__ = [
    "Anchor",
    "Arrow",
    "Ellipse",
    "EllipseArc",
    "Label",
    "Point",
    "Polyline",
    "Primitive",
    "Rect",
    "Scene",
    "Style",
]
