"""Build the scene of the interstitial volume (column) figure.

The body is chosen from the column geometry, the coordinate arrows from the
resolution, and the annotations from the transport settings.  Only the axial
flow cylinder exists so far; the remaining geometries follow the same pattern
and plug into :func:`build_column_scene`.

The drawings are adapted from the TikZ sources of the CADET modelling papers,
so a figure in the app and a figure in a paper show the same object.
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

#: Layout of the axial flow cylinder, in scene units.
_LEFT, _RIGHT = 2.1, 10.3  # the flat part of the body
_AXIS_Y = 4.3  # height of the cylinder axis
_BODY_RY = 1.75  # half-height of the body
_CAP_RX = 0.62  # how far an end cap bulges out


def _cylinder_body(scene: Scene) -> None:
    """Draw a horizontal cylinder with the near end cap facing the outlet."""
    left_centre = (_LEFT, _AXIS_Y)
    right_centre = (_RIGHT, _AXIS_Y)
    top, bottom = _AXIS_Y + _BODY_RY, _AXIS_Y - _BODY_RY

    # Fill first: the barrel and both caps share one silhouette.
    scene.add(
        Ellipse(left_centre, _CAP_RX, _BODY_RY, Style.BULK, stroked=False),
        Rect((_LEFT, bottom), _RIGHT - _LEFT, 2 * _BODY_RY, Style.BULK, stroked=False),
        Ellipse(right_centre, _CAP_RX, _BODY_RY, Style.BULK_FACE, stroked=False),
    )

    # Then the outline. The barrel contributes only its two straight edges, and
    # the far half of the inlet cap is hidden behind the body.
    scene.add(
        EllipseArc(left_centre, _CAP_RX, _BODY_RY, 90, 270, Style.WALL),
        EllipseArc(left_centre, _CAP_RX, _BODY_RY, -90, 90, Style.HIDDEN),
        Polyline(((_LEFT, top), (_RIGHT, top)), Style.WALL),
        Polyline(((_LEFT, bottom), (_RIGHT, bottom)), Style.WALL),
        Ellipse(right_centre, _CAP_RX, _BODY_RY, Style.WALL, filled=False),
    )


def _flow_arrows(scene: Scene) -> None:
    """Draw the inlet and outlet flow arrows and the boundary labels."""
    inlet_x, outlet_x = _LEFT - _CAP_RX, _RIGHT + _CAP_RX
    scene.add(
        Arrow((inlet_x - 1.35, _AXIS_Y), (inlet_x - 0.2, _AXIS_Y), Style.FLOW),
        Arrow((outlet_x + 0.2, _AXIS_Y), (outlet_x + 1.35, _AXIS_Y), Style.FLOW),
        Label((inlet_x - 0.75, _AXIS_Y + 0.22), r"\Gamma^{\mathrm{in}}", Anchor.ABOVE),
        Label((outlet_x + 0.75, _AXIS_Y + 0.22), r"\Gamma^{\mathrm{out}}", Anchor.ABOVE),
        Label((inlet_x - 0.75, _AXIS_Y - 0.28), "Flow", Anchor.BELOW, size=0.85, italic=False),
    )


def _wall_label(scene: Scene) -> None:
    """Point out the column wall from above the barrel."""
    wall_x = _LEFT + 0.35 * (_RIGHT - _LEFT)
    top = _AXIS_Y + _BODY_RY
    scene.add(
        Arrow((wall_x, top + 0.75), (wall_x, top + 0.08), Style.ANNOTATION),
        Label((wall_x, top + 0.82), r"\Gamma^{\mathrm{wall}}", Anchor.ABOVE),
    )


def _axial_coordinate(scene: Scene) -> None:
    """Draw the axial coordinate the 1D, 2D and 3D models resolve."""
    axis_y = _AXIS_Y - _BODY_RY - 0.95
    outlet = _RIGHT + _CAP_RX
    tick = 0.13
    scene.add(
        Arrow((_LEFT, axis_y), (outlet + 0.75, axis_y), Style.AXIS),
        Polyline(((_LEFT, axis_y - tick), (_LEFT, axis_y + tick)), Style.AXIS),
        Polyline(((outlet, axis_y - tick), (outlet, axis_y + tick)), Style.AXIS),
        Label((_LEFT, axis_y - tick - 0.04), "0", Anchor.BELOW, italic=False),
        Label((outlet, axis_y - tick - 0.04), r"L^{\b}", Anchor.BELOW),
        Label((outlet + 0.75, axis_y + 0.1), "z", Anchor.ABOVE_RIGHT),
    )


def build_column_scene(column_model) -> Scene:
    """Return the figure of ``column_model``'s interstitial volume.

    ``column_model`` is a :class:`~src.model_column.Column`; it is typed loosely
    to keep this module free of a dependency on the Streamlit-backed dataclass.
    """
    scene = Scene(
        name="column",
        xmin=0.0,
        ymin=0.6,
        xmax=12.6,
        ymax=7.6,
        title="Interstitial volume of the " + column_model.model_name(),
    )

    _cylinder_body(scene)
    _flow_arrows(scene)
    _wall_label(scene)
    if column_model.has_axial_coordinate:
        _axial_coordinate(scene)

    scene.add(Label(((_LEFT + _RIGHT) / 2, _AXIS_Y), r"\Omega^{\b}", Anchor.CENTER, size=1.15))

    return scene


__all__ = ["build_column_scene"]
