"""Small rendering utilities for the Streamlit equation generator.

This module provides helper functions used by the UI to display
availability badges and formatted equation text. Documentation is
kept concise and focused on public behaviour.
"""

import streamlit as st
from typing import List

from src.utils import format_variables


_BADGE_TOOLTIP_CSS = """<style>
.badge-container {
    position: relative;
    display: inline-block;
    margin-right: 6px;
}
.badge-container .badge-tooltip {
    display: none;
    position: absolute;
    top: 110%;
    left: 0;
    z-index: 1000;
    background: white;
    border: 1px solid #ddd;
    border-radius: 8px;
    padding: 8px 12px;
    font-size: 0.82em;
    min-width: 260px;
    max-width: 400px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.15);
    white-space: normal;
    color: #333;
}
.badge-container:hover .badge-tooltip {
    display: block;
}
</style>"""


def availability_badge_html(name: str, available: int, details=None) -> str:
    """Return an HTML badge indicating a model's availability.

    Parameters
    ----------
    name:
        Label shown on the badge.
    available:
        Availability flag: -1 (not present), 0 (approx.), 1 (present).
    details:
        Optional list of (label, value) tuples shown on hover.

    Returns
    -------
    str
        HTML snippet suitable for inline display in Streamlit.
    """
    # model not available default
    color_bg = "#fdecea"
    color_fg = "#b71c1c"
    icon = "not supported"

    if available == 1:
        color_bg = "#dce7f0"
        color_fg = "#023d6b"
        icon = "supported"
    elif available == 0:
        color_bg = "#fff4e5"
        color_fg = "#b26a00"
        icon = "approximation"

    badge_style = (
        f'background-color:{color_bg};'
        f'color:{color_fg};'
        f'padding:4px 10px;'
        f'border-radius:12px;'
        f'font-size:0.85em;'
        f'display:inline-block;'
        f'cursor:default;'
    )

    tooltip_html = ""
    if details:
        rows = "".join(
            f'<div style="margin:2px 0;"><b>{label}:</b> {value}</div>'
            for label, value in details
        )
        tooltip_html = f'<span class="badge-tooltip">{rows}</span>'

    return (
        f'<span class="badge-container">'
        f'<span style="{badge_style}">{name}: {icon}</span>'
        f'{tooltip_html}'
        f'</span>'
    )


def render_availability_badges(model) -> None:
    """Render availability badges with hover tooltips for solver details."""
    solver_info = model.solver_details()
    badges = availability_badge_html("CADET-Core", model.available_CADET_Core(),
                                     solver_info.get("CADET-Core"))
    badges += availability_badge_html("CADET-Process", model.available_CADET_Process(),
                                      solver_info.get("CADET-Process"))
    badges += availability_badge_html("CADET-Semi-Analytic", model.available_CADET_SemiAnalytic(),
                                      solver_info.get("CADET-Semi-Analytic"))
    html = _BADGE_TOOLTIP_CSS + f'<div style="margin-bottom:1em;">{badges}</div>'
    st.markdown(html, unsafe_allow_html=True)


def write_and_save(output: str, var_format: str, file_content: List[str], as_latex: bool = False):
    """Format and render `output` while collecting it for export.

    The function applies variable formatting, appends the formatted
    string to ``file_content`` and writes it to the Streamlit app
    as LaTeX or plain text.
    """
    output = format_variables(output, var_format)

    if output is not None:
        file_content.append(output)
        if as_latex:
            st.latex(output)
        else:
            st.write(output)


__all__ = ["availability_badge_html", "render_availability_badges", "write_and_save"]
