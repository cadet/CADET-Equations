"""Render a :class:`~src.figures.scene.Scene` as a scalable vector graphic.

The renderer is the only place that knows about SVG.  It performs the two
backend-specific steps the scene deliberately leaves open:

* the y-flip, because scene coordinates point up and SVG's point down;
* the colour lookup, which resolves a :class:`~src.figures.scene.Style` to a
  CSS custom property defined in :mod:`src.theme`.

Going through CSS variables rather than literal colours is what lets a figure
follow the light/dark appearance without a rerun, and it is why figures are
rendered inline with ``st.markdown`` instead of ``st.image``: an ``<img>``
element cannot inherit the page's custom properties.

Every coordinate is rounded to two decimals so that cosmetic changes to the
geometry maths do not ripple through the reference figures in
``test/data/ref_svg``.
"""

from __future__ import annotations

import math
from xml.sax.saxutils import escape

from src.figures.labels import label_runs, label_text
from src.figures.palette import colour
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

#: Palette entry to fill and to stroke each semantic style with, and the stroke
#: width. ``None`` means "no paint". An element carries the literal colour in
#: its presentation attribute and the palette name in its class, which is how
#: the same figure can follow the page's appearance and still open with its
#: colours in a drawing program -- see src/figures/palette.py.
_PAINT: dict[Style, tuple[str | None, str | None, float]] = {
    Style.BULK: ("bulk", None, 0.0),
    Style.BULK_FACE: ("bulk-face", None, 0.0),
    Style.VOID: ("void", None, 0.0),
    Style.WALL: (None, "wall", 0.055),
    Style.HIDDEN: (None, "wall", 0.04),
    Style.PARTICLE: ("particle", "wall", 0.04),
    Style.MATRIX: ("matrix", "wall", 0.04),
    Style.PORE: ("pore", "wall", 0.04),
    Style.SOLID: ("solid", "wall", 0.04),
    Style.CORE: ("core", "wall", 0.04),
    Style.FILM: ("film", "film-edge", 0.035),
    Style.FLOW: (None, "flow", 0.075),
    Style.FLUX: (None, "flow", 0.045),
    Style.AXIS: (None, "axis", 0.05),
    Style.ANNOTATION: (None, "annotation", 0.03),
}

#: The typeface of every label, declared on the root element.
#:
#: A figure that names no font is set in whatever the renderer defaults to --
#: Streamlit's own face in the app, a serif in a browser opening the file, a
#: system sans in Inkscape -- and the labels are then a different width in each.
#: That moves them: the name in a phase box is anchored on the right edge, so a
#: wider face grows it towards the symbol anchored on the left until the two
#: collide. The stack starts with Streamlit's face, so the app keeps the
#: appearance it has, and falls back to faces of the same proportions.
_FONT_STACK = "'Source Sans 3', 'Source Sans Pro', 'Segoe UI', Helvetica, Arial, sans-serif"

#: Styles drawn with a dashed stroke.
_DASHED = {Style.HIDDEN: "0.18 0.14"}

#: Size of a style's arrow head, relative to the default. A head is measured in
#: stroke widths, so a short arrow needs a smaller one: at full size two heads
#: are longer than the arrow between two touching phases and all that is left
#: of it is a blot.
_HEAD_SCALE = {Style.FLUX: 0.4}

#: Vertical offset of a super- or subscript, as a multiple of the font size.
_SUPERSCRIPT_RISE = 0.36
_SUBSCRIPT_DROP = 0.20
_SCRIPT_SCALE = 0.72

#: Distance between a label's anchor point and the text, in scene units.
_LABEL_PAD = 0.18


def _number(value: float) -> str:
    """Format a coordinate with at most two decimals and no trailing zeros."""
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if text in ("", "-0") else text


def _style_attributes(style: Style, filled: bool, stroked: bool) -> str:
    """Return the fill/stroke attributes and theme classes for ``style``."""
    fill, stroke, width = _PAINT[style]
    classes = []

    if filled and fill:
        parts = [f'fill="{colour(fill)}"']
        classes.append(f"cadet-fill-{fill}")
    else:
        parts = ['fill="none"']

    if stroked and stroke:
        parts.append(f'stroke="{colour(stroke)}"')
        classes.append(f"cadet-stroke-{stroke}")
        parts.append(f'stroke-width="{_number(width)}"')
        parts.append('stroke-linecap="round"')
        parts.append('stroke-linejoin="round"')
        if style in _DASHED:
            parts.append(f'stroke-dasharray="{_DASHED[style]}"')

    if classes:
        parts.append(f'class="{" ".join(classes)}"')
    return " ".join(parts)


class _Renderer:
    """Converts one scene into SVG elements, flipping y as it goes."""

    def __init__(self, scene: Scene, var_format: str):
        self.scene = scene
        self.var_format = var_format
        self.arrow_styles: set[Style] = set()

    # -- coordinate handling ------------------------------------------------

    def point(self, point: Point) -> tuple[float, float]:
        """Map a scene point into the SVG frame."""
        return point[0], self.scene.ymax + self.scene.ymin - point[1]

    def coordinates(self, point: Point) -> str:
        x, y = self.point(point)
        return f"{_number(x)},{_number(y)}"

    # -- primitives ---------------------------------------------------------

    def rect(self, item: Rect) -> str:
        x, y = self.point((item.corner[0], item.corner[1] + item.height))
        radius = f' rx="{_number(item.radius)}"' if item.radius else ""
        return (
            f'<rect x="{_number(x)}" y="{_number(y)}" width="{_number(item.width)}"'
            f' height="{_number(item.height)}"{radius} '
            f"{_style_attributes(item.style, item.filled, item.stroked)}/>"
        )

    def ellipse(self, item: Ellipse) -> str:
        x, y = self.point(item.center)
        return (
            f'<ellipse cx="{_number(x)}" cy="{_number(y)}" rx="{_number(item.rx)}" ry="{_number(item.ry)}" '
            f"{_style_attributes(item.style, item.filled, item.stroked)}/>"
        )

    def ellipse_arc(self, item: EllipseArc) -> str:
        """Draw an elliptical arc given by its start and end angle."""
        sweep = item.end_angle - item.start_angle
        # A full turn would start and end on the same point, where SVG's arc
        # command draws nothing at all; it is split so that it draws a circle.
        if abs(sweep) >= 360.0:
            sweep = 180.0 if sweep > 0 else -180.0
            halves = (
                EllipseArc(item.center, item.rx, item.ry, item.start_angle, item.start_angle + sweep, item.style),
                EllipseArc(
                    item.center, item.rx, item.ry, item.start_angle + sweep, item.start_angle + 2 * sweep, item.style
                ),
            )
            return "".join(self.ellipse_arc(half) for half in halves)

        start = self._on_ellipse(item, item.start_angle)
        end = self._on_ellipse(item, item.end_angle)
        large_arc = 1 if abs(sweep) > 180 else 0
        # Scene angles run counter-clockwise; the y-flip reverses that in SVG.
        sweep_flag = 0 if sweep > 0 else 1
        head = ""
        if item.arrow:
            self.arrow_styles.add(item.style)
            head = f' marker-end="url(#{self.scene.name}-head-{item.style.value})"'
        return (
            f'<path d="M {self.coordinates(start)} A {_number(item.rx)} {_number(item.ry)} '
            f'0 {large_arc} {sweep_flag} {self.coordinates(end)}" '
            f"{_style_attributes(item.style, False, True)}{head}/>"
        )

    @staticmethod
    def _on_ellipse(item: EllipseArc, angle: float) -> Point:
        radians = math.radians(angle)
        return (
            item.center[0] + item.rx * math.cos(radians),
            item.center[1] + item.ry * math.sin(radians),
        )

    def polyline(self, item: Polyline) -> str:
        commands = ["M " + self.coordinates(item.points[0])]
        commands += ["L " + self.coordinates(point) for point in item.points[1:]]
        if item.closed:
            commands.append("Z")
        head = ""
        if item.arrow:
            self.arrow_styles.add(item.style)
            head = f' marker-end="url(#{self.scene.name}-head-{item.style.value})"'
        return f'<path d="{" ".join(commands)}" {_style_attributes(item.style, item.filled, item.stroked)}{head}/>'

    def arrow(self, item: Arrow) -> str:
        """Draw an arrow; ``double`` adds a second head at the start."""
        self.arrow_styles.add(item.style)
        marker = f"{self.scene.name}-head-{item.style.value}"
        head = f' marker-end="url(#{marker})"'
        if item.double:
            head += f' marker-start="url(#{marker}-start)"'
        x1, y1 = self.point(item.start)
        x2, y2 = self.point(item.end)
        return (
            f'<line x1="{_number(x1)}" y1="{_number(y1)}" x2="{_number(x2)}" y2="{_number(y2)}" '
            f"{_style_attributes(item.style, False, True)}{head}/>"
        )

    def label(self, item: Label) -> str:
        """Draw a text label, expanding its LaTeX shorthand into tspans."""
        size = self.scene.font_size * item.size
        x, y = self.point(item.position)
        anchor, baseline, dx, dy = _label_placement(item.anchor)
        ink = _PAINT[item.style][1] or "annotation"
        style = ' font-style="italic"' if item.italic else ""

        spans, shift = [], 0.0
        for run in label_runs(item.text, self.var_format):
            # A nested script, as in P^{L^{b}}, rides one step further out.
            level = abs(run.shift)
            step = _SUPERSCRIPT_RISE if run.shift > 0 else _SUBSCRIPT_DROP
            direction = -1.0 if run.shift > 0 else 1.0
            target = direction * size * step * (1.0 + 0.7 * (level - 1)) if level else 0.0
            run_size = size * _SCRIPT_SCALE**level
            attributes = f' font-size="{_number(run_size)}"'
            if not run.italic:
                attributes += ' font-style="normal"'
            if abs(target - shift) > 1e-9:
                attributes += f' dy="{_number(target - shift)}"'
                shift = target
            spans.append(f"<tspan{attributes}>{escape(run.text)}</tspan>")

        return (
            # The space between two words is a run, and so a tspan, of its own.
            # Without xml:space a renderer may trim it as insignificant
            # whitespace, which is how "film diffusion" becomes "filmdiffusion"
            # in Inkscape; no label starts or ends with a space, so preserving
            # it changes nothing else.
            f'<text xml:space="preserve" x="{_number(x + dx)}" y="{_number(y + dy)}" '
            f'font-size="{_number(size)}" fill="{colour(ink)}" class="cadet-fill-{ink}" '
            f'text-anchor="{anchor}" dominant-baseline="{baseline}"{style}>'
            f"{''.join(spans)}</text>"
        )

    # -- document -----------------------------------------------------------

    def markers(self) -> str:
        """Return the arrow-head definitions used by this scene."""
        definitions = []
        for style in sorted(self.arrow_styles, key=lambda item: item.value):
            ink = _PAINT[style][1] or "annotation"
            scale = _HEAD_SCALE.get(style, 1.0)
            length, width = 6 * scale, 4.8 * scale
            path = (
                f"M 0 0 L {_number(length)} {_number(0.5 * width)} "
                f"L 0 {_number(width)} L {_number(0.2 * length)} {_number(0.5 * width)} Z"
            )
            for suffix, orientation in (("", "auto"), ("-start", "auto-start-reverse")):
                definitions.append(
                    f'<marker id="{self.scene.name}-head-{style.value}{suffix}"'
                    f' markerWidth="{_number(length)}" markerHeight="{_number(width)}"'
                    f' refX="{_number(0.9 * length)}" refY="{_number(0.5 * width)}" orient="{orientation}"'
                    f' markerUnits="strokeWidth"><path d="{path}" fill="{colour(ink)}"'
                    f' class="cadet-fill-{ink}"/></marker>'
                )
        return "".join(definitions)


def _label_placement(anchor: Anchor) -> tuple[str, str, float, float]:
    """Return (text-anchor, dominant-baseline, dx, dy) for an anchor."""
    horizontal = {
        Anchor.CENTER: ("middle", 0.0),
        Anchor.ABOVE: ("middle", 0.0),
        Anchor.BELOW: ("middle", 0.0),
        Anchor.LEFT: ("start", _LABEL_PAD),
        Anchor.RIGHT: ("end", -_LABEL_PAD),
        Anchor.ABOVE_LEFT: ("end", -_LABEL_PAD),
        Anchor.ABOVE_RIGHT: ("start", _LABEL_PAD),
        Anchor.BELOW_LEFT: ("end", -_LABEL_PAD),
        Anchor.BELOW_RIGHT: ("start", _LABEL_PAD),
    }[anchor]
    # SVG y grows downwards, so "above" is a negative offset.
    vertical = {
        Anchor.CENTER: ("middle", 0.0),
        Anchor.LEFT: ("middle", 0.0),
        Anchor.RIGHT: ("middle", 0.0),
        Anchor.ABOVE: ("baseline", -_LABEL_PAD),
        Anchor.ABOVE_LEFT: ("baseline", -_LABEL_PAD),
        Anchor.ABOVE_RIGHT: ("baseline", -_LABEL_PAD),
        Anchor.BELOW: ("hanging", _LABEL_PAD),
        Anchor.BELOW_LEFT: ("hanging", _LABEL_PAD),
        Anchor.BELOW_RIGHT: ("hanging", _LABEL_PAD),
    }[anchor]
    return horizontal[0], vertical[0], horizontal[1], vertical[1]


#: A scene this many units wide is displayed at _REFERENCE_MAX. Narrower
#: scenes are displayed proportionally smaller so that text is the same size
#: in every figure on the page, rather than growing as the drawing shrinks.
_REFERENCE_WIDTH, _REFERENCE_MAX = 13.0, 46.0


def render_svg(scene: Scene, var_format: str = "CADET", max_width: str | None = None) -> str:
    """Return ``scene`` as a standalone, responsive SVG element."""
    if max_width is None:
        max_width = _number(_REFERENCE_MAX * scene.width / _REFERENCE_WIDTH) + "rem"
    renderer = _Renderer(scene, var_format)

    body = []
    for item in scene.items:
        if isinstance(item, Rect):
            body.append(renderer.rect(item))
        elif isinstance(item, Ellipse):
            body.append(renderer.ellipse(item))
        elif isinstance(item, EllipseArc):
            body.append(renderer.ellipse_arc(item))
        elif isinstance(item, Polyline):
            body.append(renderer.polyline(item))
        elif isinstance(item, Arrow):
            body.append(renderer.arrow(item))
        elif isinstance(item, Label):
            body.append(renderer.label(item))
        else:  # pragma: no cover - guards against an unhandled primitive
            raise TypeError(f"Cannot render primitive of type {type(item).__name__}")

    defs = renderer.markers()
    defs_element = f"<defs>{defs}</defs>" if defs else ""
    title = escape(label_text(scene.title, var_format)) if scene.title else ""
    title_element = f"<title>{title}</title>" if title else ""

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" role="img" class="cadet-figure" '
        f'data-figure="{escape(scene.name)}" font-family="{_FONT_STACK}" '
        f'viewBox="{_number(scene.xmin)} {_number(scene.ymin)} '
        f'{_number(scene.width)} {_number(scene.height)}" '
        f'style="width:100%;max-width:{max_width};height:auto;display:block;margin:0.5em auto;">'
        f"{title_element}{defs_element}{''.join(body)}</svg>"
    )


__all__ = ["render_svg"]
