"""Small rendering utilities for the Streamlit equation generator.

This module provides helper functions used by the UI to display
availability badges and formatted equation text. Documentation is
kept concise and focused on public behaviour.
"""

import streamlit as st
from typing import List

from src.utils import format_variables


def availability_badge_html(name: str, available: int) -> str:
    """Return an HTML badge indicating a model's availability.

    Parameters
    ----------
    name:
        Label shown on the badge.
    available:
        Availability flag: -1 (not present), 0 (approx.), 1 (present).

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
        color_bg = "#e6f4ea"
        color_fg = "#137333"
        icon = "supported"
    elif available == 0:
        color_bg = "#fff4e5"
        color_fg = "#b26a00"
        icon = "approximation"

    return (
        f'<span style="'
        f'background-color:{color_bg};'
        f'color:{color_fg};'
        f'padding:4px 10px;'
        f'border-radius:12px;'
        f'font-size:0.85em;'
        f'margin-right:6px;'
        f'display:inline-block;">'
        f'{name}: {icon}'
        f'</span>'
    )


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


__all__ = ["availability_badge_html", "write_and_save"]
