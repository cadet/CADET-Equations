r"""Canonical formatting for the generated LaTeX output.

The document is assembled from many small snippets, so a plain concatenation
is only ever as consistent as the snippets are.  :func:`normalize` re-formats
an assembled document (or any self-contained snippet) into one canonical
shape.  That keeps the generated ``.tex`` readable and diffable, and it makes
the reference files in ``test/data/ref_latex`` a meaningful contract: a diff
there is a change of content, never a change of layout.

Rules
-----
* ``\begin{...}`` and ``\end{...}`` always start a line and nothing else
  follows them on that line.
* The body between ``\begin`` and ``\end`` is indented one level.  The
  ``document`` environment does not add a level.
* Inside a math environment, every row break ``\\`` sits on a line of its own,
  every alignment marker opens a line, and the terms of a sum are broken at
  top-level ``+`` / ``-`` and indented one further level.
* Every row carries an alignment marker: a lone top-level ``=`` becomes ``&=``.
* Every equation is numbered, so starred math environments are unstarred.
* One sentence per line in running text; lines are never wrapped at a fixed
  column.
* At most one blank line anywhere; never one inside an environment, and never
  one directly before or after a display equation.
* No tabs, no trailing whitespace, no space before ``,``, ``.`` or ``;``.
"""

import re

INDENT = "    "

#: Environments whose body is laid out as equation rows.
MATH_ENVIRONMENTS = frozenset({"align", "alignat", "eqnarray", "equation", "flalign", "gather", "multline"})

#: Math environments that align their rows, i.e. that accept ``&`` markers.
ALIGNMENT_ENVIRONMENTS = frozenset({"align", "alignat", "eqnarray", "flalign"})

#: Environments that do not indent their body.
NO_INDENT_ENVIRONMENTS = frozenset({"document"})

_ENV_ARGS = r"(?:\{[^{}]*\}|\[[^\[\]]*\])*"
_BEGIN_RE = re.compile(r"\\begin\{([A-Za-z]+\*?)\}" + _ENV_ARGS)
_END_RE = re.compile(r"\\end\{([A-Za-z]+\*?)\}")
_COMMAND_RE = re.compile(r"\\(?:[A-Za-z]+|.)")

# A segment body starting with one of these is glued to its alignment marker,
# so that the marker and the operator read as one token (``&=``, ``&+``).
_OPERATOR_RE = re.compile(r"^(?:[=+\-<>]|\\(?:leq|geq|le|ge|neq|approx|equiv|sim|propto)\b)")

# Words ending in a period that do not end a sentence.
_ABBREVIATIONS = frozenset(
    {"e.g.", "i.e.", "cf.", "etc.", "al.", "vs.", "resp.", "approx.", "eq.", "eqs.", "fig.", "ref.", "no."}
)


# --------------------------------------------------------------------------- #
# depth-aware scanning
# --------------------------------------------------------------------------- #
def _top_level_tokens(text: str, wanted: frozenset) -> list[tuple[int, str]]:
    r"""Return ``(position, token)`` for every wanted token at top level.

    Top level means outside ``{}`` groups, outside ``()``/``[]`` groups and
    outside ``\left ... \right`` pairs, so that a minus sign inside
    ``\frac{1 - \varepsilon}{...}`` is never mistaken for a term separator.
    The two-character row break ``\\`` may be requested as a token as well.
    """
    tokens: list[tuple[int, str]] = []
    brace = group = left = 0
    i, n = 0, len(text)
    while i < n:
        char = text[i]
        if char == "\\":
            if text.startswith(r"\\", i):
                if brace == group == left == 0 and r"\\" in wanted:
                    tokens.append((i, r"\\"))
                i += 2
                continue
            match = _COMMAND_RE.match(text, i)
            command = match.group(0)
            if command == r"\left":
                left += 1
            elif command == r"\right":
                left -= 1
            i = match.end()
            continue
        if char == "{":
            brace += 1
        elif char == "}":
            brace -= 1
        elif char in "([":
            group += 1
        elif char in ")]":
            group -= 1
        elif brace == group == left == 0 and char in wanted:
            tokens.append((i, char))
        i += 1
    return tokens


def _split_at_top_level(text: str, token: str) -> list[str]:
    """Split ``text`` at every top-level occurrence of ``token``."""
    parts, previous = [], 0
    for position, _ in _top_level_tokens(text, frozenset({token})):
        parts.append(text[previous:position])
        previous = position + len(token)
    parts.append(text[previous:])
    return parts


# --------------------------------------------------------------------------- #
# math environments
# --------------------------------------------------------------------------- #
def _tidy_math_line(line: str) -> str:
    """Drop the space in front of trailing punctuation."""
    return re.sub(r"\s+([,.;])", r"\1", line)


def _ensure_alignment_marker(row: str) -> str:
    """Give a row without any alignment marker one in front of its relation."""
    if _top_level_tokens(row, frozenset({"&"})):
        return row
    relations = _top_level_tokens(row, frozenset({"="}))
    if not relations:
        return row
    position = relations[0][0]
    return (row[:position].rstrip() + " &" + row[position:]).strip()


def _split_terms(body: str) -> list[str]:
    r"""Split the terms of a sum at top-level ``+`` / ``-``.

    A leading sign stays with its term, and everything from the first
    top-level comma onwards (``, \quad k = 1, \ldots, N - M``) is left alone,
    since it annotates the equation rather than contributing a term.
    """
    body = body.strip()
    if not body:
        return [""]

    commas = _top_level_tokens(body, frozenset({","}))
    annotation_start = commas[0][0] if commas else len(body)

    cuts = []
    for position, _ in _top_level_tokens(body, frozenset({"+", "-"})):
        if position == 0 or position >= annotation_start:
            continue
        preceding = body[:position].rstrip()
        if preceding and (preceding[-1].isalnum() or preceding[-1] in "}])|"):
            cuts.append(position)

    terms, previous = [], 0
    for cut in cuts:
        terms.append(body[previous:cut].strip())
        previous = cut
    terms.append(body[previous:].strip())
    return [term for term in terms if term]


def _format_row(row: str, level: int, aligned: bool = True) -> list[str]:
    """Lay out a single equation row as one line per alignment segment.

    ``aligned`` is false for environments such as ``equation`` that take no
    alignment markers; their rows are only broken into terms.
    """
    row = row.strip()
    if not aligned:
        terms = _split_terms(row)
        return [INDENT * level + _tidy_math_line(terms[0])] + [
            INDENT * (level + 1) + _tidy_math_line(term) for term in terms[1:]
        ]

    row = _ensure_alignment_marker(row)
    columns = _split_at_top_level(row, "&")

    # Pair each column with the marker that introduces it, merging the markers
    # of empty columns (``& &`` in an alignat annotation) into the next one.
    segments: list[tuple[str, str]] = []
    marker = ""
    for index, column in enumerate(columns):
        if index:
            marker = (marker + " &").strip()
        body = column.strip()
        if not body:
            continue
        segments.append((marker, body))
        marker = ""
    if marker:
        segments.append((marker, ""))

    lines = []
    for marker, body in segments:
        terms = _split_terms(body)
        if not marker:
            head = terms[0]
        elif _OPERATOR_RE.match(terms[0]):
            head = marker + terms[0]
        else:
            head = (marker + " " + terms[0]).strip()
        lines.append(INDENT * level + _tidy_math_line(head))
        for term in terms[1:]:
            lines.append(INDENT * (level + 1) + _tidy_math_line(term))
    return lines


def _format_math_body(body: str, level: int, environment: str) -> list[str]:
    """Lay out the body of a math environment."""
    body = re.sub(r"\s+", " ", body).strip()
    if not body:
        return []
    aligned = environment in ALIGNMENT_ENVIRONMENTS
    lines = []
    for index, row in enumerate(_split_at_top_level(body, r"\\")):
        if index:
            lines.append(INDENT * level + r"\\")
        lines.extend(_format_row(row, level, aligned))
    return lines


# --------------------------------------------------------------------------- #
# running text
# --------------------------------------------------------------------------- #
def _ends_with_abbreviation(text: str) -> bool:
    word = text.split()[-1] if text.split() else ""
    return word.lower() in _ABBREVIATIONS or bool(re.fullmatch(r"[A-Za-z]\.", word))


def _split_sentences(line: str) -> list[str]:
    """Split running text after sentence-ending punctuation."""
    sentences, start, in_math, i = [], 0, False, 0
    while i < len(line):
        char = line[i]
        if char == "\\":
            i += 2
            continue
        if char == "$":
            in_math = not in_math
        elif not in_math and char in ".!?" and line[i + 1 : i + 2] == " ":
            following = line[i + 2 :].lstrip()
            if following[:1].isupper() and not _ends_with_abbreviation(line[: i + 1]):
                sentences.append(line[start : i + 1].strip())
                start = i + 2
        i += 1
    sentences.append(line[start:].strip())
    return [sentence for sentence in sentences if sentence]


def _format_text_line(line: str, level: int) -> list[str]:
    """Lay out one line of running text, one sentence per line.

    A list entry stays on a single line: ``\\item`` and ``\\bibitem`` delimit
    the entries, so breaking inside one would make the list unreadable.
    """
    line = re.sub(r" {2,}", " ", line.strip())
    line = _tidy_math_line(line)
    if re.match(r"\\(?:bib)?item\b", line):
        return [INDENT * level + line]
    return [INDENT * level + sentence for sentence in _split_sentences(line)]


# --------------------------------------------------------------------------- #
# document structure
# --------------------------------------------------------------------------- #
def _number_all_equations(text: str) -> str:
    """Unstar math environments so that every equation carries a number."""
    starred = "|".join(sorted(MATH_ENVIRONMENTS))
    return re.sub(r"(\\(?:begin|end)\{(?:" + starred + r"))\*\}", r"\1}", text)


def _explode_environments(text: str) -> str:
    r"""Put every ``\begin`` and ``\end`` on a line of its own."""

    def put_on_own_line(match: re.Match) -> str:
        before = "" if re.search(r"(?:^|\n)[ \t]*$", text[: match.start()]) else "\n"
        after = "" if re.match(r"[ \t]*(?:\n|$)", text[match.end() :]) else "\n"
        return before + match.group(0) + after

    text = _BEGIN_RE.sub(put_on_own_line, text)
    return _END_RE.sub(put_on_own_line, text)


def _is_display_begin(line: str) -> bool:
    match = _BEGIN_RE.match(line.strip())
    return bool(match) and match.group(1) in MATH_ENVIRONMENTS


def _is_display_end(line: str) -> bool:
    match = _END_RE.match(line.strip())
    return bool(match) and match.group(1) in MATH_ENVIRONMENTS


def _tidy_blank_lines(lines: list[str]) -> list[str]:
    """Collapse blank runs and drop blank lines around display equations."""
    collapsed: list[str] = []
    for line in lines:
        if not line.strip():
            if not collapsed or not collapsed[-1].strip():
                continue
            collapsed.append("")
        else:
            collapsed.append(line)

    result: list[str] = []
    for index, line in enumerate(collapsed):
        if not line.strip():
            following = collapsed[index + 1] if index + 1 < len(collapsed) else ""
            if _is_display_begin(following) or _is_display_end(result[-1] if result else ""):
                continue
        result.append(line)
    return result


def normalize(tex: str) -> str:
    """Return ``tex`` in the canonical layout described in the module docstring.

    The function is idempotent: normalizing already normalized text is a no-op.
    It only moves whitespace around, with two deliberate exceptions: starred
    math environments are unstarred, and a row without an alignment marker
    gets one in front of its relation.
    """
    text = tex.replace("\r\n", "\n").replace("\r", "\n").replace("\t", " ")
    text = _number_all_equations(text)
    text = _explode_environments(text)

    lines: list[str] = []
    stack: list[str] = []
    math_body: list[str] | None = None

    for raw_line in text.split("\n"):
        line = raw_line.strip()
        level = sum(1 for env in stack if env not in NO_INDENT_ENVIRONMENTS)

        begin = _BEGIN_RE.fullmatch(line)
        end = _END_RE.fullmatch(line)

        if begin and math_body is None:
            lines.append(INDENT * level + line)
            stack.append(begin.group(1))
            if begin.group(1) in MATH_ENVIRONMENTS:
                math_body = []
            continue

        if end and stack and stack[-1] == end.group(1):
            if math_body is not None:
                lines.extend(_format_math_body(" ".join(math_body), level, end.group(1)))
                math_body = None
            stack.pop()
            closing_level = sum(1 for env in stack if env not in NO_INDENT_ENVIRONMENTS)
            lines.append(INDENT * closing_level + line)
            continue

        if math_body is not None:
            if line:
                math_body.append(line)
            continue

        if not line:
            if not any(env not in NO_INDENT_ENVIRONMENTS for env in stack):
                lines.append("")
            continue

        lines.extend(_format_text_line(line, level))

    return "\n".join(_tidy_blank_lines(lines))


def normalize_document(tex: str) -> str:
    """Normalize a whole document: no leading blank lines, one trailing newline."""
    return normalize(tex).strip("\n") + "\n"


# --------------------------------------------------------------------------- #
# building display equations
# --------------------------------------------------------------------------- #
def display_equation(body: str, environment: str = "align", columns: int | None = None) -> str:
    r"""Wrap an equation body in a numbered display environment.

    Call sites only supply the equation itself; the layout inside the
    environment is imposed by :func:`normalize` further down the line.
    """
    arguments = "" if columns is None else "{" + str(columns) + "}"
    return "\\begin{" + environment + "}" + arguments + "\n" + body.strip() + "\n\\end{" + environment + "}"


def display_rows(rows, environment: str = "align", columns: int | None = None, terminator: str = ".") -> str:
    r"""Wrap several equation rows in one display environment.

    The rows are separated by ``,\\`` and the last one is closed with
    ``terminator``.
    """
    body = ",\n" + "\\\\" + "\n"
    return display_equation(body.join(row.strip() for row in rows) + terminator, environment, columns)


__all__ = [
    "INDENT",
    "MATH_ENVIRONMENTS",
    "display_equation",
    "display_rows",
    "normalize",
    "normalize_document",
]
