"""The annotation conventions every figure follows.

The builders draw different things -- a column, a particle, an isotherm -- but
they annotate them the same way, and the three rules below are the ones a
reader relies on to tell the parts of a figure apart.  They live here so that
"the figures are consistent" is something the code says once.

* **An extent is a dimension line**: a plain line with a tick at each end, the
  bounds written beyond the ticks and the symbol in the middle.  It carries no
  arrowhead, because a head reads as a direction and a coordinate is not one --
  ``z`` runs from 0 to ``L^{b}``, it does not point at anything.
* **A boundary names itself from outside** the body it belongs to, with a short
  arrow from the label onto the boundary, so the name cannot be mistaken for a
  label of the volume it encloses.
* **A transfer that runs both ways is a double-headed arrow** (``double=True``).
  Film diffusion, pore and surface diffusion and binding are exchanges, not
  one-way fluxes, and only the convective flow path carries a single head.
"""

from __future__ import annotations

from src.figures.scene import Anchor, Arrow, Label, Point, Polyline, Scene, Style

#: Half the length of the tick that closes a dimension line, in scene units.
TICK = 0.13

#: Distance between a tick and the bound written beyond it.
_GAP = 0.05


def dimension(
    scene: Scene,
    start: Point,
    end: Point,
    symbol: str,
    start_label: str | None = None,
    end_label: str | None = None,
    flip: bool = False,
    size: float = 0.9,
) -> None:
    """Draw the extent from ``start`` to ``end`` as a dimension line.

    The line is horizontal or vertical; the ticks stand perpendicular to it and
    the labels go below a horizontal line and to the right of a vertical one,
    or on the other side when ``flip`` is set.
    """
    vertical = abs(end[1] - start[1]) > abs(end[0] - start[0])
    side = -1.0 if flip else 1.0
    # Where a label sits relative to the line, and which way the ticks run.
    if vertical:
        out = (side * (TICK + _GAP), 0.0)
        tick = (TICK, 0.0)
        ends = Anchor.LEFT if side > 0 else Anchor.RIGHT
    else:
        out = (0.0, -side * (TICK + _GAP))
        tick = (0.0, TICK)
        ends = Anchor.BELOW if side > 0 else Anchor.ABOVE

    scene.add(Polyline((start, end), Style.AXIS))
    for point, text in ((start, start_label), (end, end_label)):
        scene.add(
            Polyline(((point[0] - tick[0], point[1] - tick[1]), (point[0] + tick[0], point[1] + tick[1])), Style.AXIS)
        )
        if text is not None:
            scene.add(Label((point[0] + out[0], point[1] + out[1]), text, ends, size=size, italic=not text.isdigit()))

    # The symbol goes on the other side of the line from the bounds: on a short
    # extent -- an annulus between two radii, a thin cap -- the three labels
    # would otherwise be written over one another.
    middle = (0.5 * (start[0] + end[0]) - out[0], 0.5 * (start[1] + end[1]) - out[1])
    scene.add(Label(middle, symbol, _opposite(ends), size=1.0))


def boundary(scene: Scene, point: Point, direction: Point, text: str, reach: float = 0.7, size: float = 1.0) -> None:
    """Name the boundary at ``point``, from ``direction`` outside the body.

    ``direction`` is the unit vector pointing away from the body; the label
    sits at the far end of a short arrow that runs back onto the boundary.
    """
    tail = (point[0] + direction[0] * reach, point[1] + direction[1] * reach)
    anchor = (point[0] + direction[0] * (reach + 0.07), point[1] + direction[1] * (reach + 0.07))
    scene.add(
        Arrow(tail, (point[0] + direction[0] * 0.1, point[1] + direction[1] * 0.1), Style.ANNOTATION),
        Label(anchor, text, _away(direction), size=size),
    )


def _opposite(anchor: Anchor) -> Anchor:
    """Return the anchor that places a label on the other side of a line."""
    return {
        Anchor.BELOW: Anchor.ABOVE,
        Anchor.ABOVE: Anchor.BELOW,
        Anchor.LEFT: Anchor.RIGHT,
        Anchor.RIGHT: Anchor.LEFT,
    }[anchor]


def _away(direction: Point) -> Anchor:
    """Return the anchor that keeps a label clear of an arrow out of a body."""
    if abs(direction[1]) >= abs(direction[0]):
        return Anchor.ABOVE if direction[1] > 0 else Anchor.BELOW
    return Anchor.LEFT if direction[0] > 0 else Anchor.RIGHT
