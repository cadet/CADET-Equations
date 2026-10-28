"""
@author: jmbr
"""

import bibtexparser
import streamlit as st
from bibtexparser.bparser import BibTexParser

from src.theme import apply_theme_css

st.logo("images/logo_CADET.png", size="large", link=None, icon_image=None)

st.set_page_config(
    page_title="Documentation",
    page_icon=":material/docs:",
)

apply_theme_css()

# Helper funciton to use citations in standard markdown (streamlit)
with open("CITATION.bib", encoding="utf-8") as bibtex_file:
    parser = BibTexParser(common_strings=True)
    bib_database = bibtexparser.load(bibtex_file, parser=parser)

citation_links = {}
formatted_refs = []


def strip_braces(value: str) -> str:
    """Drop BibTeX grouping braces, e.g. in "Eric {von Lieres}", from a field value."""
    return value.replace("{", "").replace("}", "")


for entry in bib_database.entries:
    key = entry.get("ID", "unknown").lower()
    author = strip_braces(entry.get("author", "Unknown Author"))
    title = strip_braces(entry.get("title", "Untitled"))
    chapter = entry.get("chapter", "")
    if chapter != "":
        chapter = f", chapter {chapter}"
    year = entry.get("year", "n.d.")
    journal = entry.get("journal", "")
    if journal == "":
        journal = entry.get("publisher", "")
    journal = strip_braces(journal)
    doi = entry.get("doi", "").removeprefix("https://doi.org/")

    main_author = author.split(" and ")[0].split(",")[0]
    label = f"{main_author} et al., {year}" if "and" in author else f"{main_author}, {year}"
    anchor_id = f"ref-{key}"

    citation_links[key] = f"[{label}](#{anchor_id})"

    doi_html = f' <a href="https://doi.org/{doi}" target="_blank">DOI: {doi}</a>' if doi else ""
    formatted_refs.append(f'<a id="{anchor_id}"></a>{author} ({year}). <i>{title}</i>{chapter}. {journal}.{doi_html}')


st.markdown(
    """
<style>
.cadet-hero {
    background-color: var(--cadet-panel-bg);
    color: var(--cadet-panel-fg);
    padding: 1.5rem;
    border-left: 5px solid var(--cadet-accent);
    border-radius: 8px;
    line-height: 1.6;
}

.cadet-hero-title {
    font-size: 2rem;
    font-weight: bold;
    color: var(--cadet-accent);
    font-family: serif;
}

.cadet-hero a {
    color: var(--cadet-link);
}

.cadet-hero a:hover {
    color: var(--cadet-link-hover);
}
</style>
<div class="cadet-hero">
<span class="cadet-hero-title">CADET-Equations</span><br>
is a Python tool for <b>generating model equations</b> based on a <b>unified and composable modeling framework</b> for packed-bed chromatography (<a href="#ref-breuer2026unified">Breuer et al., 2026</a>).

This framework embodies the CADET group's modeling philosophy, and its principles are reflected throughout our software stack.
[CADET-Core](https://cadet.github.io/master/index.html#cadet-core), for instance, incorporates much of the framework's modular design, enabling it to simulate the majority of models that can be configured here.

CADET-Equations provides a simplistic user interface that allows users to define models, and outputs the corresponding equations in <b>LaTeX format</b> and as <b>PDF</b>.
</div>
""",
    unsafe_allow_html=True,
)

st.markdown(
    r"""
    ## Features

    - Generation of comprehensive model equations used in packed-bed liquid chromatography
    - Configuration via model upload:
      - [CADET model file](https://cadet.github.io/master/interface/index.html) (HDF5)
      - simplistic configuration file: download and upload your session state as json file
    - Downloadable models as PDF and LaTeX files. We thankfully ask to be cited: <a href="#ref-breuer2026unified">Breuer et al., 2026</a>
    - Availability indicator for CADET-Core and CADET-Process. Note that some models can be approximated even if the corresponding equations are not implemented exactly

    ## Contribute

    We welcome contributions - such as [reporting bugs](https://github.com/cadet/CADET-Equations/issues/new?template=bug_report.yml) or [requesting features](https://github.com/cadet/CADET-Equations/issues/new?template=feature_request.yml) on GitHub, or by reaching out with modeling questions on the [CADET forum](https://forum.cadet-web.de/).

""",
    unsafe_allow_html=True,
)

st.markdown("---")
st.markdown(
    """
    ## Citation

    If you find the tool useful and use it for your own work, please cite the corresponding publication of the
    underlying unified and composable modeling framework:
    Breuer, J. M., Leweke, S., Schmölder, J., Berger, A., Gassner, G. J., von Lieres, E. (2026).
    *A unified and composable framework for modeling packed-bed liquid chromatography*.
    [https://doi.org/10.48550/arXiv.2610.33269](https://doi.org/10.48550/arXiv.2610.33269)
""",
    unsafe_allow_html=True,
)

st.markdown("---")
st.markdown("## References", unsafe_allow_html=True)
for ref in formatted_refs:
    st.markdown(ref, unsafe_allow_html=True)
