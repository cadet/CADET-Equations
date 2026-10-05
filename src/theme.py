"""Theme-aware colors for the custom HTML/CSS injected by the app.

Streamlit themes its own widgets from the ``[theme.light]`` / ``[theme.dark]``
sections in ``.streamlit/config.toml``, but the raw HTML this app renders needs
its own palette. The variables below use the CSS ``light-dark()`` function and
therefore resolve against the ``color-scheme`` that Streamlit sets on the app
container, so they follow the appearance chosen under Settings > Appearance
without the app having to rerun.
"""

import streamlit as st

_THEME_CSS = """
<style>
.stApp {
    /* CADET blue, lightened for dark mode where the brand navy is unreadable */
    --cadet-link: light-dark(#023d6b, #6fb3e0);
    --cadet-link-hover: light-dark(#145a86, #9ecdf0);
    --cadet-accent: light-dark(#023d6b, #6fb3e0);

    --cadet-button-bg: light-dark(#023d6b, #0b5d94);
    --cadet-button-bg-hover: light-dark(#145a86, #1a7ab8);
    --cadet-button-fg: #ffffff;

    --cadet-panel-bg: light-dark(#f0f4f8, #1a2230);
    --cadet-panel-fg: light-dark(#0f172a, #e5e7eb);

    --cadet-tooltip-bg: light-dark(#ffffff, #1a2230);
    --cadet-tooltip-fg: light-dark(#333333, #e5e7eb);
    --cadet-tooltip-border: light-dark(#dddddd, #3a4556);
    --cadet-tooltip-shadow: light-dark(rgba(0, 0, 0, 0.15), rgba(0, 0, 0, 0.5));

    --cadet-badge-supported-bg: light-dark(#dce7f0, #16324a);
    --cadet-badge-supported-fg: light-dark(#023d6b, #9ecdf0);
    --cadet-badge-approx-bg: light-dark(#fff4e5, #3d2f12);
    --cadet-badge-approx-fg: light-dark(#b26a00, #f0c070);
    --cadet-badge-unsupported-bg: light-dark(#fdecea, #45201d);
    --cadet-badge-unsupported-fg: light-dark(#b71c1c, #f3a8a3);
}
</style>
"""


def apply_theme_css() -> None:
    """Define the CADET color variables for the current page."""
    st.markdown(_THEME_CSS, unsafe_allow_html=True)
