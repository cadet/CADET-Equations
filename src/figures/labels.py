"""Turn the LaTeX shorthand used in figure labels into plain text runs.

Figure labels are written with the same ``\\b \\p \\s \\l`` shorthands as the
equations (see :mod:`src.utils`), so a figure honours the selected parameter
format for free: in ``Legacy`` format a label reading ``c^{\\s}`` comes out as
``q``, exactly as the equations do.

SVG has no math renderer, so a label is reduced to a list of :class:`Run`
objects -- a piece of text plus the vertical shift and the upright/italic flag
it should be drawn with.  The SVG backend maps those to ``<tspan>`` elements.  A
TikZ backend would ignore this module and emit the LaTeX source unchanged, which
is why the conversion lives here and not in :mod:`src.figures.svg`.
"""

from __future__ import annotations

from dataclasses import dataclass

from src.utils import format_variables

#: Vertical position of a run, in units of the base font size.
NORMAL, SUPERSCRIPT, SUBSCRIPT = 0, 1, -1

_SYMBOLS = {
    r"\alpha": "α",
    r"\beta": "β",
    r"\gamma": "γ",
    r"\Delta": "Δ",
    r"\delta": "δ",
    r"\varepsilon": "ε",
    r"\epsilon": "ε",
    r"\theta": "θ",
    r"\lambda": "λ",
    r"\mu": "μ",
    r"\nu": "ν",
    r"\rho": "ρ",
    r"\sigma": "σ",
    r"\phi": "φ",
    r"\varphi": "φ",
    r"\Gamma": "Γ",
    r"\Lambda": "Λ",
    r"\Omega": "Ω",
    r"\ell": "ℓ",
    r"\partial": "∂",
    r"\nabla": "∇",
    r"\infty": "∞",
    r"\cdot": "·",
    r"\times": "×",
    r"\to": "→",
    r"\rightleftharpoons": "⇌",
    r"\in": "∈",
    r"\dots": "…",
    r"\ldots": "…",
    r"\approx": "≈",
    r"\leq": "≤",
    r"\geq": "≥",
}

#: Symbols that are upright in mathematical typesetting: capital Greek letters,
#: relations and operators. Lowercase Greek letters stay italic.
_UPRIGHT_SYMBOLS = {
    r"\Delta",
    r"\Gamma",
    r"\Lambda",
    r"\Omega",
    r"\nabla",
    r"\infty",
    r"\cdot",
    r"\times",
    r"\to",
    r"\rightleftharpoons",
    r"\in",
    r"\dots",
    r"\ldots",
    r"\approx",
    r"\leq",
    r"\geq",
}

#: Commands that only affect spacing or sizing and carry no glyph.
_SPACING = {r"\,": " ", r"\;": " ", r"\:": " ", r"\ ": " ", r"\quad": "  ", r"\qquad": "    "}
_IGNORED = (r"\left", r"\right", r"\big", r"\Big", r"\bigg", r"\Bigg", r"\displaystyle")


@dataclass(frozen=True)
class Run:
    """A stretch of label text sharing one vertical shift and one font style."""

    text: str
    shift: int = NORMAL
    italic: bool = True


def _balanced_group(source: str, index: int) -> tuple[str, int]:
    """Return the contents of the ``{...}`` group starting at ``index``.

    ``index`` must point at the opening brace.  The returned index points just
    past the matching closing brace.
    """
    depth = 0
    for position in range(index, len(source)):
        if source[position] == "{":
            depth += 1
        elif source[position] == "}":
            depth -= 1
            if depth == 0:
                return source[index + 1 : position], position + 1
    # Unbalanced input: treat the remainder as the group.
    return source[index + 1 :], len(source)


def _argument(source: str, index: int) -> tuple[str, int]:
    """Return the argument of a command or script at ``index``.

    The argument is either a braced group or, in LaTeX fashion, the single
    character that follows.
    """
    if index >= len(source):
        return "", index
    if source[index] == "{":
        return _balanced_group(source, index)
    return source[index], index + 1


def _command(source: str, index: int) -> tuple[str, int]:
    """Return the control sequence starting at the backslash at ``index``."""
    position = index + 1
    while position < len(source) and source[position].isalpha():
        position += 1
    if position == index + 1:  # a single non-letter command such as "\,"
        position += 1
    return source[index:position], position


def _convert(source: str, shift: int, italic: bool, runs: list[Run]) -> None:
    """Append the runs of ``source`` to ``runs``, inheriting shift and style."""
    index = 0
    while index < len(source):
        char = source[index]

        if char == "\\":
            command, index = _command(source, index)
            if command == r"\mathrm" or command == r"\text" or command == r"\operatorname":
                argument, index = _argument(source, index)
                _convert(argument, shift, False, runs)
            elif command == r"\vec" or command == r"\bar" or command == r"\hat":
                argument, index = _argument(source, index)
                _convert(argument, shift, italic, runs)
            elif command in _SYMBOLS:
                runs.append(Run(_SYMBOLS[command], shift, italic and command not in _UPRIGHT_SYMBOLS))
            elif command in _SPACING:
                runs.append(Run(_SPACING[command], shift, italic))
            elif command in _IGNORED:
                continue
            else:  # unknown command: drop the backslash, keep the name
                runs.append(Run(command[1:], shift, italic))

        elif char in "^_":
            argument, index = _argument(source, index + 1)
            # A script inside a script, as in P^{L^{b}}, nests one level
            # further out rather than sitting at the same height.
            step = SUPERSCRIPT if char == "^" else SUBSCRIPT
            _convert(argument, shift + step if shift * step > 0 else step, italic, runs)

        elif char == "{":
            group, index = _balanced_group(source, index)
            _convert(group, shift, italic, runs)

        elif char == "$":  # labels are plain math; stray delimiters are noise
            index += 1

        else:
            # Digits, operators and punctuation are upright by math convention.
            runs.append(Run(char, shift, italic and char.isalpha()))
            index += 1


def _merge(runs: list[Run]) -> list[Run]:
    """Join neighbouring runs that share a shift and a font style."""
    merged: list[Run] = []
    for run in runs:
        if not run.text:
            continue
        if merged and merged[-1].shift == run.shift and merged[-1].italic == run.italic:
            merged[-1] = Run(merged[-1].text + run.text, run.shift, run.italic)
        else:
            merged.append(run)
    return merged


def label_runs(text: str, var_format: str = "CADET") -> list[Run]:
    """Expand ``text`` for ``var_format`` and split it into drawable runs."""
    runs: list[Run] = []
    _convert(format_variables(text, var_format), NORMAL, True, runs)
    return _merge(runs)


def label_text(text: str, var_format: str = "CADET") -> str:
    """Return the label as flat text, used for accessible descriptions."""
    return "".join(run.text for run in label_runs(text, var_format))


__all__ = ["NORMAL", "SUBSCRIPT", "SUPERSCRIPT", "Run", "label_runs", "label_text"]
