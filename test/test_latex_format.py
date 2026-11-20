r"""
Unit tests for the canonical LaTeX layout produced by src/latex_format.py.

The fixpoint test over test/data/ref_latex is the actual regression guard: if
every reference file is already normalized, then no generated document can
drift out of the canonical layout without a reference file changing too.
"""

import glob
import os
import re

import pytest

from src.latex_format import MATH_ENVIRONMENTS, normalize, normalize_document

REF_LATEX_DIR = os.path.join(os.path.dirname(__file__), "data", "ref_latex")
REF_LATEX_FILES = sorted(glob.glob(os.path.join(REF_LATEX_DIR, "*.tex")))

_BEGIN_RE = re.compile(r"\\begin\{([A-Za-z]+\*?)\}")
_END_RE = re.compile(r"\\end\{([A-Za-z]+\*?)\}")


# %% The reference files are the contract


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("path", REF_LATEX_FILES, ids=lambda p: os.path.basename(p)[:-4])
def test_reference_files_are_normalized(path):
    """Every reference .tex must be a fixpoint of the formatter."""
    with open(path, encoding="utf-8") as handle:
        content = handle.read()
    assert normalize_document(content) == content, (
        f"{os.path.basename(path)} is not in canonical layout; regenerate the reference files."
    )


@pytest.mark.ci
@pytest.mark.reference
@pytest.mark.parametrize("path", REF_LATEX_FILES, ids=lambda p: os.path.basename(p)[:-4])
def test_reference_files_obey_the_layout_rules(path):
    """Spot-check the rules directly, independently of the formatter."""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().split("\n")

    depth = 0
    for number, line in enumerate(lines, start=1):
        where = f"{os.path.basename(path)}:{number}"
        assert line == line.rstrip(), f"{where}: trailing whitespace"
        assert "\t" not in line, f"{where}: tab character"

        stripped = line.strip()
        begin, end = _BEGIN_RE.match(stripped), _END_RE.match(stripped)
        if begin or end:
            assert stripped.startswith("\\"), f"{where}: content before \\begin/\\end"
            assert _BEGIN_RE.fullmatch(stripped) or _END_RE.fullmatch(stripped) or stripped.endswith("}"), (
                f"{where}: content after \\begin/\\end"
            )
        if end:
            depth -= 0 if end.group(1) == "document" else 1
        assert not stripped or line.startswith("    " * depth), f"{where}: wrong indentation"
        if begin:
            depth += 0 if begin.group(1) == "document" else 1
            assert "*" not in begin.group(1) or begin.group(1) not in MATH_ENVIRONMENTS, f"{where}: unnumbered equation"
    assert depth == 0, f"{os.path.basename(path)}: unbalanced environments"


# %% Layout rules


@pytest.mark.ci
@pytest.mark.unit_test
def test_begin_and_end_get_their_own_line():
    formatted = normalize(r"text \begin{align} a &= b. \end{align} more text")
    assert formatted.split("\n") == ["text", r"\begin{align}", "    a", "    &= b.", r"\end{align}", "more text"]


@pytest.mark.ci
@pytest.mark.unit_test
def test_rows_and_alignment_markers_open_their_own_line():
    formatted = normalize(r"\begin{align} a &= b, \\ c &= d. \end{align}")
    assert formatted.split("\n") == [
        r"\begin{align}",
        "    a",
        "    &= b,",
        r"    \\",
        "    c",
        "    &= d.",
        r"\end{align}",
    ]


@pytest.mark.ci
@pytest.mark.unit_test
def test_terms_of_a_sum_are_broken_at_top_level_operators():
    formatted = normalize(r"\begin{align} a &= - b + c - d. \end{align}")
    assert formatted.split("\n") == [
        r"\begin{align}",
        "    a",
        "    &= - b",
        "        + c",
        "        - d.",
        r"\end{align}",
    ]


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize(
    "body",
    [
        r"a &= \frac{1 - b}{c}.",  # inside braces
        r"a &= \left( 1 - b \right).",  # inside \left ... \right
        r"a &= (1 - b).",  # inside parentheses
        r"a &= b^{-1}.",  # an exponent sign
        r"a &= b, \quad k = 1, \ldots, N - M.",  # in the annotation after the equation
    ],
)
def test_operators_that_are_not_top_level_terms_are_left_alone(body):
    formatted = normalize(r"\begin{align} " + body + r" \end{align}")
    assert formatted.split("\n") == [r"\begin{align}", "    a", "    &" + body.split("&", 1)[1], r"\end{align}"]


@pytest.mark.ci
@pytest.mark.unit_test
def test_second_alignment_column_keeps_its_markers():
    formatted = normalize(r"\begin{alignat}{2} a &= b & &\qquad\text{on }(0, T). \end{alignat}")
    assert formatted.split("\n") == [
        r"\begin{alignat}{2}",
        "    a",
        "    &= b",
        r"    & & \qquad\text{on }(0, T).",
        r"\end{alignat}",
    ]


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_row_without_a_marker_gets_one():
    assert "    &= b." in normalize(r"\begin{align} a = b. \end{align}").split("\n")


@pytest.mark.ci
@pytest.mark.unit_test
def test_environments_without_alignment_keep_their_rows_unmarked():
    """``equation`` takes no alignment tab, so none may be inserted."""
    assert "&" not in normalize(r"\begin{equation} a = b. \end{equation}")


@pytest.mark.ci
@pytest.mark.unit_test
def test_every_equation_is_numbered():
    formatted = normalize(r"\begin{align*} a &= b. \end{align*}")
    assert formatted.startswith(r"\begin{align}") and formatted.endswith(r"\end{align}")


@pytest.mark.ci
@pytest.mark.unit_test
def test_nested_environments_are_indented_once_per_level():
    raw = "\n".join([r"\begin{document}", r"\begin{itemize}", r"\item one", r"\end{itemize}", r"\end{document}"])
    formatted = normalize(raw)
    assert formatted.split("\n") == [
        r"\begin{document}",
        r"\begin{itemize}",
        r"    \item one",
        r"\end{itemize}",
        r"\end{document}",
    ]


# %% Whitespace


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize(
    "raw, expected",
    [
        ("a\n\n\n\nb", "a\n\nb"),  # at most one blank line
        ("text\n\n\\begin{align}\na &= b.\n\\end{align}", "text\n\\begin{align}\n    a\n    &= b.\n\\end{align}"),
        ("\\begin{align}\na &= b.\n\\end{align}\n\ntext", "\\begin{align}\n    a\n    &= b.\n\\end{align}\ntext"),
        ("a  b   c", "a b c"),  # no double spaces in text
        ("one. Two.", "one.\nTwo."),  # one sentence per line
        ("see Breuer et al. (2026) here.", "see Breuer et al. (2026) here."),  # not at an abbreviation
        ("e.g. Something happens.", "e.g. Something happens."),
        ("trailing   \nspace", "trailing\nspace"),
    ],
)
def test_whitespace_rules(raw, expected):
    assert normalize(raw) == expected


@pytest.mark.ci
@pytest.mark.unit_test
def test_blank_lines_inside_environments_are_dropped():
    assert "\n\n" not in normalize("\\begin{itemize}\n\\item one\n\n\\item two\n\\end{itemize}")


@pytest.mark.ci
@pytest.mark.unit_test
def test_space_before_punctuation_is_removed():
    assert normalize(r"\begin{align} a &= b , \end{align}").count("&= b,") == 1


@pytest.mark.ci
@pytest.mark.unit_test
def test_a_document_ends_with_exactly_one_newline():
    assert normalize_document("a\n\n\n") == "a\n"


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("path", REF_LATEX_FILES[:5], ids=lambda p: os.path.basename(p)[:-4])
def test_normalize_is_idempotent(path):
    with open(path, encoding="utf-8") as handle:
        content = handle.read()
    once = normalize_document(content)
    assert normalize_document(once) == once
