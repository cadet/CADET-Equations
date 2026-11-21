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

#: Stroke and fill for each semantic style, as (fill, stroke, stroke width).
#: ``None`` means "no paint"; names resolve against the variables in src.theme.
_PAINT: dict[Style, tuple[str | None, str | None, float]] = {
    Style.BULK: ("var(--cadet-fig-bulk)", None, 0.0),
    Style.BULK_FACE: ("var(--cadet-fig-bulk-face)", None, 0.0),
    Style.VOID: ("var(--cadet-fig-void)", None, 0.0),
    Style.WALL: (None, "var(--cadet-fig-wall)", 0.055),
    Style.HIDDEN: (None, "var(--cadet-fig-wall)", 0.04),
    Style.PARTICLE: ("var(--cadet-fig-particle)", "var(--cadet-fig-wall)", 0.04),
    Style.PORE: ("var(--cadet-fig-pore)", "var(--cadet-fig-wall)", 0.04),
    Style.SOLID: ("var(--cadet-fig-solid)", "var(--cadet-fig-wall)", 0.04),
    Style.CORE: ("var(--cadet-fig-core)", "var(--cadet-fig-wall)", 0.04),
    Style.FILM: ("var(--cadet-fig-film)", "var(--cadet-fig-film-edge)", 0.035),
    Style.FLOW: (None, "var(--cadet-fig-flow)", 0.075),
    Style.AXIS: (None, "var(--cadet-fig-axis)", 0.05),
    Style.ANNOTATION: (None, "var(--cadet-fig-annotation)", 0.03),
}

#: Styles drawn with a dashed stroke.
_DASHED = {Style.HIDDEN: "0.18 0.14"}

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
    """Return the fill/stroke attributes for ``style``."""
    fill, stroke, width = _PAINT[style]
    parts = [f'fill="{fill}"' if filled and fill else 'fill="none"']
    if stroked and stroke:
        parts.append(f'stroke="{stroke}"')
        parts.append(f'stroke-width="{_number(width)}"')
        parts.append('stroke-linecap="round"')
        parts.append('stroke-linejoin="round"')
        if style in _DASHED:
            parts.append(f'stroke-dasharray="{_DASHED[style]}"')
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
        return f'<path d="{" ".join(commands)}" {_style_attributes(item.style, item.filled, item.stroked)}/>'

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
        fill = _PAINT[item.style][1] or "var(--cadet-fig-annotation)"
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
            f'<text x="{_number(x + dx)}" y="{_number(y + dy)}" font-size="{_number(size)}" '
            f'fill="{fill}" text-anchor="{anchor}" dominant-baseline="{baseline}"{style}>'
            f"{''.join(spans)}</text>"
        )

    # -- document -----------------------------------------------------------

    def markers(self) -> str:
        """Return the arrow-head definitions used by this scene."""
        definitions = []
        for style in sorted(self.arrow_styles, key=lambda item: item.value):
            colour = _PAINT[style][1] or "var(--cadet-fig-annotation)"
            for suffix, orientation, path in (
                ("", "auto", "M 0 0 L 6 2.4 L 0 4.8 L 1.2 2.4 Z"),
                ("-start", "auto-start-reverse", "M 0 0 L 6 2.4 L 0 4.8 L 1.2 2.4 Z"),
            ):
                definitions.append(
                    f'<marker id="{self.scene.name}-head-{style.value}{suffix}" markerWidth="6"'
                    f' markerHeight="4.8" refX="5.4" refY="2.4" orient="{orientation}"'
                    f' markerUnits="strokeWidth"><path d="{path}" fill="{colour}"/></marker>'
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


def render_svg(scene: Scene, var_format: str = "CADET", max_width: str = "46rem") -> str:
    """Return ``scene`` as a standalone, responsive SVG element."""
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
        f'viewBox="{_number(scene.xmin)} {_number(scene.ymin)} '
        f'{_number(scene.width)} {_number(scene.height)}" '
        f'style="width:100%;max-width:{max_width};height:auto;display:block;margin:0.5em auto;">'
        f"{title_element}{defs_element}{''.join(body)}</svg>"
    )


__all__ = ["render_svg"]
