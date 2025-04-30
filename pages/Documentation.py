# -*- coding: utf-8 -*-
"""
@author: jmbr
"""

import streamlit as st
import bibtexparser
from bibtexparser.bparser import BibTexParser

st.logo("images/logo_CADET.png", size="medium", link=None, icon_image=None)

st.set_page_config(
    page_title="Documentation",
    page_icon=":material/docs:",
)

# Helper funciton to use citations in standard markdown (streamlit)
with open("CITATION.bib", "r", encoding="latin1") as bibtex_file:
    parser = BibTexParser(common_strings=True)
    bib_database = bibtexparser.load(bibtex_file, parser=parser)

citation_links = {}
formatted_refs = []

for entry in bib_database.entries:
    key = entry.get("ID", "unknown").lower()
    author = entry.get("author", "Unknown Author")
    title = entry.get("title", "Untitled")
    chapter = entry.get("chapter", "")
    if not chapter == "":
        chapter = f", chapter {chapter}"
    year = entry.get("year", "n.d.")
    journal = entry.get("journal", "")
    if journal == "":
        journal = entry.get("publisher", "")
    doi = entry.get("doi", "")

    main_author = author.split(" and ")[0].split(",")[0]
    label = f"{main_author} et al., {year}" if "and" in author else f"{main_author}, {year}"
    anchor_id = f"ref-{key}"

    citation_links[key] = f"[{label}](#{anchor_id})"

    doi_html = f' <a href="https://doi.org/{doi}" target="_blank">DOI: {doi}</a>' if doi else ""
    formatted_refs.append(f'<a id="{anchor_id}"></a>{author} ({year}). <i>{title}</i>{chapter}. {journal}.{doi_html}')


st.markdown("""
<div style="background-color: #f0f4f8; padding: 1.5rem; border-left: 5px solid #1f77b4; border-radius: 8px; line-height: 1.6;">
<span style="font-size: 2rem; font-weight: bold; color: #1f77b4; font-family: serif;">CADET-Equations</span><br>
is a Python tool for <b>generating model equations</b> based on a <b>unified and holistic modeling framework</b> for packed-bed chromatography (<a href="#ref-leweke2021">Leweke, 2021</a>).

This framework embodies the CADET group's modeling philosophy, and its principles are reflected throughout our software stack.
[CADET-Core](https://cadet.github.io/master/index.html#cadet-core), for instance, incorporates much of the framework's modular design, enabling it to simulate the majority of models that can be configured here.

CADET-Equations provides a simplistic user interface that allows users to define models, and outputs the corresponding equations in <b>LaTeX format</b> and as <b>PDF</b>.
</div>
""", unsafe_allow_html=True)


st.markdown(
"""
    ## Features

    - Generation of comprehensive model equations used in packed-bed liquid chromatography
    - Configuration via model upload:
      - [CADET model file](https://cadet.github.io/master/interface/index.html) (HDF5)
      - simplistic configuration file: download and upload your session state as json file
    - Downloadable models as PDF and LaTeX files. We thankfully ask to be cited

    ## Limitations

    - [Specific Binding Models](https://github.com/cadet/CADET-Equations/issues/15), binding is currently defined by an arbitrary function $f^\mathrm{bind}$
    - [Per-component parameterization](https://github.com/cadet/CADET-Equations/issues/7), e.g. of film diffusion kinetics
    - [Per-particle-type parameterization](https://github.com/cadet/CADET-Equations/issues/21), e.g. of film diffusion kinetics
    - [Reactions](https://github.com/cadet/CADET-Equations/issues/6)

    ## Contribute

    We welcome contributions - such as [reporting bugs on GitHub](https://github.com/cadet/CADET-Equations/issues), or by reaching out with modeling questions on the [CADET forum](https://forum.cadet-web.de/).

"""
)

st.markdown("---")
st.markdown("## References", unsafe_allow_html=True)
for ref in formatted_refs:
    st.markdown(ref, unsafe_allow_html=True)
