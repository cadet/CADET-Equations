"""The colours a figure is drawn in, as light and dark appearance pairs.

The palette lives here rather than in :mod:`src.theme` so that it has a single
definition and no Streamlit import.  Two consumers read it: the theme turns it
into the ``--cadet-fig-*`` custom properties of the page, and
:mod:`src.figures.svg` writes the light appearance into every paint as the
fallback of its ``var()``.

That fallback is what makes a saved figure usable on its own.  A file opened
outside the app has no page to inherit the custom properties from, and CSS
treats a reference to an undefined one as invalid, so ``fill`` falls back to
its initial value -- black.  Without the fallback every standalone figure is a
silhouette.
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


def paint(name: str) -> str:
    """Return the paint for ``name``: its variable, light appearance as fallback."""
    return f"var(--cadet-fig-{name}, {PALETTE[name][0]})"


def css_variables(indent: str = "    ") -> str:
    """Return one ``--cadet-fig-*`` declaration per line, for a CSS rule."""
    return "\n".join(
        f"{indent}--cadet-fig-{name}: light-dark({light}, {dark});" for name, (light, dark) in PALETTE.items()
    )
