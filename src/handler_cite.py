"""
Bibliography utilities for parsing and formatting BibTeX entries.

This module provides functions to load a BibTeX file, extract author names,
and format in-text citations and references for both plain text and LaTeX output.
"""

import re
from pathlib import Path

import streamlit as st


def _strip_bibtex_value(value: str) -> str:
    """Remove outer braces/quotes and collapse whitespace for BibTeX values."""
    if value is None:
        return ""
    value = value.strip().rstrip(",")
    if (value.startswith("{") and value.endswith("}")) or (value.startswith('"') and value.endswith('"')):
        value = value[1:-1]
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def load_bibliography(bib_path: str) -> dict:
    """Parse a simple BibTeX file into a dict keyed by citation key."""
    bib_text = Path(bib_path).read_text(encoding="utf-8")
    entries = {}

    for match in re.finditer(r"@(\w+)\s*\{\s*([^,]+),", bib_text):
        entry_type = match.group(1).strip().lower()
        key = match.group(2).strip()
        body_start = match.end()
        next_match = re.search(r"\n@", bib_text[body_start:])
        body_end = body_start + next_match.start() if next_match else len(bib_text)
        body = bib_text[body_start:body_end]

        fields = {}
        for field_match in re.finditer(r"(?m)^\s*([A-Za-z]+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|\"[^\"]*\")\s*,?", body):
            field = field_match.group(1).lower()
            value = _strip_bibtex_value(field_match.group(2))
            fields[field] = value

        fields["entry_type"] = entry_type
        fields["key"] = key
        entries[key] = fields

    return entries


def _author_last_names(author_field: str) -> list[str]:
    """Extract author last names from a BibTeX author list."""
    if not author_field:
        return []
    names = []
    for raw_author in [part.strip() for part in author_field.split(" and ") if part.strip()]:
        cleaned = raw_author.replace("{", "").replace("}", "")
        if "," in cleaned:
            last_name = cleaned.split(",", 1)[0].strip()
        else:
            tokens = cleaned.split()
            last_name = tokens[-1] if tokens else cleaned
        names.append(last_name)
    return names


def format_intext_citation(entry: dict) -> str:
    """Return an author-year in-text citation like 'Leweke et al. (2021)'."""
    authors = _author_last_names(entry.get("author", ""))
    year = entry.get("year", "n.d.")

    if not authors:
        return f"({year})"
    if len(authors) == 1:
        return f"{authors[0]} ({year})"
    if len(authors) == 2:
        return f"{authors[0]} and {authors[1]} ({year})"
    return f"{authors[0]} et al. ({year})"


def format_reference_plain(entry: dict) -> str:
    """Return a compact plain-text reference string for Streamlit output."""
    authors = entry.get("author", "")
    year = entry.get("year", "n.d.")
    title = entry.get("title", "")
    venue = entry.get("journal") or entry.get("school") or ""
    doi = entry.get("doi", "")

    ref = f"{authors} ({year}). {title}."
    if venue:
        ref += f" {venue}."
    if doi:
        ref += f" DOI: {doi}."
    return ref


def _escape_latex(text: str) -> str:
    """Escape a plain string for safe inclusion in LaTeX text mode."""
    replacements = {
        "\\": r"\textbackslash{}",
        "{": r"\{",
        "}": r"\}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def format_reference_latex(entry: dict) -> str:
    """Return a BibTeX-entry-like line for LaTeX thebibliography."""
    return _escape_latex(format_reference_plain(entry))


def cite(citation_key: str, bibliography_entries: dict, used_citation_keys: list) -> str:
    """Register citation key and return in-text citation text."""
    if citation_key not in bibliography_entries:
        return f"[{citation_key}]"
    if citation_key not in used_citation_keys:
        used_citation_keys.append(citation_key)
    return format_intext_citation(bibliography_entries[citation_key])


def cite_html(citation_key: str, bibliography_entries: dict, used_citation_keys: list) -> str:
    """Register citation key and return linked HTML in-text citation."""
    citation_text = cite(citation_key, bibliography_entries, used_citation_keys)
    if citation_key not in bibliography_entries:
        return citation_text
    return f"<a href='#ref-{citation_key}' style='color:#023d6b;'>{citation_text}</a>"


def render_references(bibliography_entries: dict, used_citation_keys: list, file_content: list):
    """Render bibliography to Streamlit output and LaTeX export content."""
    if not used_citation_keys:
        return

    st.write("### References")

    file_content.append(r"\begin{thebibliography}{99}")

    for key in used_citation_keys:
        entry = bibliography_entries.get(key)
        if entry is None:
            continue
        st.markdown(f"- <a id='ref-{key}'></a>{format_reference_plain(entry)}", unsafe_allow_html=True)
        file_content.append(rf"\bibitem{{{key}}} " + format_reference_latex(entry))

    file_content.append(r"\end{thebibliography}")
