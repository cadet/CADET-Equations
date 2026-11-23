"""The colours a figure is drawn in, as light and dark appearance pairs.

The palette lives here rather than in :mod:`src.theme` so that it has a single
definition and no Streamlit import.  Two consumers read it, and between them
they let one figure serve both the app and a drawing program.

:mod:`src.figures.svg` paints every element with the *light* colour, written
out as a literal in the ``fill`` and ``stroke`` presentation attributes, and
adds a class naming the style it came from.  :mod:`src.theme` turns the palette
into a rule per class, which the browser applies on top because a presentation
attribute loses against any stylesheet, so the page gets the appearance it asks
for while the file keeps a colour of its own.

A figure may therefore not be painted through a CSS custom property.  Inkscape
and most other SVG editors implement neither ``var()`` nor its fallback, and an
unreadable paint leaves ``fill`` at its initial value -- black -- which is what
a figure saved from a theme-only palette looks like once it is opened outside a
browser.
"""

from __future__ import annotations

#: Semantic name -> (light appearance, dark appearance).
PALETTE: dict[str, tuple[str, str]] = {
    "bulk": ("#dce7f0", "#16324a"),
    "bulk-face": ("#eef4f9", "#1e3f5c"),
    "void": ("#ffffff", "#0e1117"),
    "wall": ("#0f172a", "#d7dee8"),
    "particle": ("#ffffff", "#2b3545"),
    "matrix": ("#d4e0ea", "#2b3545"),
    "pore": ("#f7fafd", "#19212e"),
    "solid": ("#93b9d8", "#44688f"),
    "core": ("#9fb0c2", "#5c6b80"),
    "film": ("#fff4e5", "#2e2718"),
    "film-edge": ("#b26a00", "#f0c070"),
    "flow": ("#023d6b", "#6fb3e0"),
    "axis": ("#475569", "#94a3b8"),
    "annotation": ("#0f172a", "#d7dee8"),
}


def colour(name: str) -> str:
    """Return the literal colour an element of style ``name`` is painted with."""
    return PALETTE[name][0]


def css_rules(scope: str = "") -> str:
    """Return one rule per painted style, for the page to repaint a figure.

    Written as ``light-dark()`` so that the rules follow the appearance without
    a rerun, and scoped by ``scope`` -- a selector the figures sit inside.
    """
    rules = []
    for name, (light, dark) in PALETTE.items():
        for property_name in ("fill", "stroke"):
            rules.append(f"{scope}.cadet-{property_name}-{name} {{ {property_name}: light-dark({light}, {dark}); }}")
    return "\n".join(rules)
