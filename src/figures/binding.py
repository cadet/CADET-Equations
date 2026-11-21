"""Build the scene of the binding model figure.

Only binding models with a characteristic picture get one:

============  =========================================================
model         figure
============  =========================================================
``Linear``    the isotherm, a straight line of slope k^a / k^d
``Langmuir``  the isotherm, saturating at the capacity q^max
``SMA``       the exchange itself, after Brooks and Cramer (1992): a
              protein binds to nu sites and displaces that many
              counter-ions, while shielding sigma further sites
``Arbitrary`` none -- an unspecified isotherm has nothing to draw, and
              the particle figure already shows the exchange
============  =========================================================

The symbols match :mod:`src.equations`, so the figure and the binding term
below it name the same quantities.
"""

from __future__ import annotations

from src.figures.scene import (
    Anchor,
    Arrow,
    Ellipse,
    Label,
    Polyline,
    Rect,
    Scene,
    Style,
)

#: Binding models this module can draw.
DRAWABLE = ("Linear", "Langmuir", "SMA")

# %% Isotherm panel

_ORIGIN = (1.75, 1.15)
_WIDTH, _HEIGHT = 5.6, 4.0
_SAMPLES = 48


def _isotherm_axes(scene: Scene) -> None:
    """Draw the concentration axes an isotherm is plotted against."""
    x, y = _ORIGIN
    scene.add(
        Arrow((x, y), (x + _WIDTH + 0.55, y), Style.AXIS),
        Arrow((x, y), (x, y + _HEIGHT + 0.55), Style.AXIS),
        Label((x + _WIDTH + 0.6, y + 0.05), r"c^{\p}", Anchor.ABOVE_RIGHT, size=0.95),
        Label((x - 0.15, y + _HEIGHT + 0.6), r"c^{\s}", Anchor.ABOVE_LEFT, size=0.95),
    )


def _linear(scene: Scene) -> None:
    """A linear isotherm: the bound phase follows the pore liquid."""
    x, y = _ORIGIN
    scene.add(
        Polyline(((x, y), (x + _WIDTH, y + _HEIGHT * 0.92)), Style.FLOW),
        Label(
            (x + _WIDTH * 0.62, y + _HEIGHT * 0.57),
            r"k^{\mathrm{a}} / k^{\mathrm{d}}",
            Anchor.ABOVE_LEFT,
            size=0.95,
        ),
    )


def _langmuir(scene: Scene) -> None:
    """A Langmuir isotherm: the bound phase saturates at the capacity."""
    x, y = _ORIGIN
    plateau = _HEIGHT * 0.78
    steepness = 5.5  # fixes how quickly the curve approaches the plateau

    curve = tuple(
        (
            x + fraction * _WIDTH,
            y + plateau * steepness * fraction / (1.0 + steepness * fraction),
        )
        for fraction in (index / _SAMPLES for index in range(_SAMPLES + 1))
    )
    scene.add(
        Polyline(((x, y + plateau), (x + _WIDTH + 0.3, y + plateau)), Style.HIDDEN),
        Label((x - 0.15, y + plateau), r"q^{\mathrm{max}}", Anchor.RIGHT, size=0.9),
        Polyline(curve, Style.FLOW),
    )


# %% Steric mass action

_SURFACE_Y = 1.9
_TICK = 0.38
_SITE_SPACING = 0.62
_FIRST_SITE = 1.3
_SITES = 15
_PROTEIN_SITES = (6, 7, 8)  # the nu sites the protein occupies
_SHIELDED_SITES = (5, 9)  # the sigma sites it sterically blocks


def _site_x(index: int) -> float:
    return _FIRST_SITE + index * _SITE_SPACING


def _bracket(scene: Scene, first: int, last: int, depth: float, text: str) -> None:
    """Draw a brace under a run of binding sites and name it."""
    left, right = _site_x(first) - 0.2, _site_x(last) + 0.2
    top, bottom = _SURFACE_Y - depth, _SURFACE_Y - depth - 0.22
    scene.add(
        Polyline(((left, top), (left, bottom), (right, bottom), (right, top)), Style.ANNOTATION),
        Label((0.5 * (left + right), bottom - 0.08), text, Anchor.BELOW, size=0.9),
    )


def _sma(scene: Scene) -> None:
    """The steric mass action exchange, after Brooks and Cramer (1992)."""
    # The resin surface, carrying Lambda charged ligands.
    scene.add(Polyline(((0.9, _SURFACE_Y), (_site_x(_SITES - 1) + 0.6, _SURFACE_Y)), Style.WALL))
    for index in range(_SITES):
        x = _site_x(index)
        scene.add(Polyline(((x, _SURFACE_Y), (x, _SURFACE_Y + _TICK)), Style.WALL))
        if index not in _PROTEIN_SITES and index not in _SHIELDED_SITES:
            # A free ligand holds a counter-ion.
            scene.add(Ellipse((x, _SURFACE_Y + _TICK + 0.2), 0.19, 0.19, Style.SOLID))

    # The protein, bound across nu ligands.
    left, right = _site_x(_PROTEIN_SITES[0]) - 0.42, _site_x(_PROTEIN_SITES[-1]) + 0.42
    bottom = _SURFACE_Y + _TICK + 0.06
    scene.add(
        Rect((left, bottom), right - left, 1.25, Style.PARTICLE, radius=0.5),
        Label((0.5 * (left + right), bottom + 0.62), r"c^{\s}_{i}", Anchor.CENTER, size=0.95),
    )

    _bracket(scene, _PROTEIN_SITES[0], _PROTEIN_SITES[-1], 0.15, r"\nu_{i}")
    _bracket(scene, _SHIELDED_SITES[0], _SHIELDED_SITES[0], 0.85, r"\sigma_{i}")
    _bracket(scene, _SHIELDED_SITES[1], _SHIELDED_SITES[1], 0.85, r"\sigma_{i}")
    _bracket(scene, 0, _SITES - 1, 1.65, r"\Lambda")

    # The protein arrives from the liquid; nu counter-ions are released.
    scene.add(
        Arrow((2.1, 4.9), (2.1, bottom + 1.45), Style.FLOW),
        Label((2.25, 4.6), r"c^{\p}_{i}", Anchor.LEFT, size=0.95),
        Arrow((_site_x(_SITES - 1) - 0.3, bottom + 0.4), (_site_x(_SITES - 1) - 0.3, 4.9), Style.FLOW),
        Label((_site_x(_SITES - 1) - 0.15, 4.6), r"c^{\s}_{0}", Anchor.LEFT, size=0.95),
    )

    # Named above the surface, so the swatch is not taken for a bound ion.
    scene.add(
        Ellipse((4.4, 4.55), 0.19, 0.19, Style.SOLID),
        Label((4.75, 4.55), "counter-ion", Anchor.LEFT, size=0.8, italic=False),
    )


#: Isotherm plots share their axes; the SMA cartoon brings its own frame.
_ISOTHERMS = {"Linear": _linear, "Langmuir": _langmuir}


def build_binding_scene(binding_model: str, name: str = "binding") -> Scene | None:
    """Return the figure of ``binding_model``, or ``None`` if it has none."""
    if binding_model == "SMA":
        scene = Scene(name=name, xmin=0.4, ymin=-0.7, xmax=11.4, ymax=5.6, title="Steric mass action binding")
        _sma(scene)
        return scene

    draw = _ISOTHERMS.get(binding_model)
    if draw is None:
        return None

    scene = Scene(name=name, xmin=0.0, ymin=0.3, xmax=9.3, ymax=6.4, title=binding_model + " binding isotherm")
    _isotherm_axes(scene)
    draw(scene)
    return scene


__all__ = ["DRAWABLE", "build_binding_scene"]
