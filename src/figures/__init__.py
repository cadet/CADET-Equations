"""Model figures: scalable drawings of the configured model.

The package mirrors :mod:`src.equations`: pure functions that turn a model into
text, here SVG rather than LaTeX, with nothing importing Streamlit.  A figure is
built in two steps -- a builder assembles a backend-independent
:class:`~src.figures.scene.Scene`, and a renderer turns it into output -- so
that a TikZ renderer can be added later without rewriting the figures.
"""

from src.figures.column import build_column_scene
from src.figures.particle import build_particle_scene
from src.figures.scene import Scene
from src.figures.svg import render_svg


def column_figure(column_model, var_format: str = "CADET") -> str:
    """Return the SVG of ``column_model``'s interstitial volume."""
    return render_svg(build_column_scene(column_model), var_format)


def particle_figures(column_model, var_format: str = "CADET") -> list[str]:
    """Return one SVG per distinct particle type of ``column_model``.

    Each figure gets its own scene name because the renderer derives its arrow
    head ids from it and several figures share one page.
    """
    types = list(getattr(column_model, "par_type_counts", None) or getattr(column_model, "particle_models", []) or [])
    return [
        render_svg(build_particle_scene(particle, f"particle-{index + 1}"), var_format)
        for index, particle in enumerate(types)
    ]


__all__ = [
    "Scene",
    "build_column_scene",
    "build_particle_scene",
    "column_figure",
    "particle_figures",
    "render_svg",
]
