# -*- coding: utf-8 -*-
"""

This file implements utility helpers

"""

import re
from typing import List

_CADET_PATTERNS = [
    (re.compile(r"\\l(?![a-zA-Z])"), r"\\mathrm{\\ell}"),
    (re.compile(r"\\b(?![a-zA-Z])"), r"\\mathrm{b}"),
    (re.compile(r"\\p(?![a-zA-Z])"), r"\\mathrm{p}"),
    (re.compile(r"\\s(?![a-zA-Z])"), r"\\mathrm{s}"),
]

_LEGACY_PATTERNS = [
    (re.compile(r"c\^\{\\l\}(?![a-zA-Z])"), r"c"),
    (re.compile(r"c\^\{\\b\}(?![a-zA-Z])"), r"c"),
    (re.compile(r"V\^\{\\l\}(?![a-zA-Z])"), r"V^{\\mathrm{\\ell}}"),
    (re.compile(r"V\^\{\\b\}(?![a-zA-Z])"), r"V^{\\mathrm{\\ell}}"),
    (re.compile(r"c\^\{\\p\}_\\{i\\}(?![a-zA-Z])"), r"c_{\\mathrm{p},i}"),
    (re.compile(r"c\^\{\\p\}_i(?![a-zA-Z])"), r"c_{\\mathrm{p},i}"),
    (re.compile(r"c\^\{\\p\}_\\{j,i\\}(?![a-zA-Z])"), r"c_{\\mathrm{p},j,i}"),
    (re.compile(r"c\^\{\\p\}(?![a-zA-Z])"), r"c_{\\mathrm{p}}"),
    (re.compile(r"c\}\\^\{\\p\\}(?![a-zA-Z])"), r"c}_{\\mathrm{p}}"),
    (re.compile(r"c\}\\^\{\\l\\}(?![a-zA-Z])"), r"c}"),
    (re.compile(r"c\^\{\\s\}(?![a-zA-Z])"), r"q"),
    (re.compile(r"c\}\\^\{\\s\\}(?![a-zA-Z])"), r"q}"),
    (re.compile(r"\\p(?![a-zA-Z])"), r"\\mathrm{p}"),
    (re.compile(r"\\s(?![a-zA-Z])"), r"\\mathrm{s}"),
]

def format_variables(input_str: str, var_format: str, mul_pars: bool = False) -> str:
    """format variable symbols according to the chosen format.

    var_format must be either "CADET" or "Legacy".
    """
    if not isinstance(var_format, str):
        raise ValueError(f"Format {var_format} is not supported. (expected 'CADET' or 'Legacy')")

    if var_format == "CADET":
        for patt, repl in _CADET_PATTERNS:
            input_str = patt.sub(repl, input_str)
    elif var_format == "Legacy":
        for patt, repl in _LEGACY_PATTERNS:
            input_str = patt.sub(repl, input_str)
        # Merge various types of adjacent subscripts produced by sequential replacements
        input_str = re.sub(r"c_\{\\mathrm\{([^}]*)\}\}_\{([^}]*)\}", r"c_{\\mathrm{\1},\2}", input_str)
        input_str = re.sub(r"c_\{\\mathrm\{([^}]*)\}\}_([A-Za-z0-9]+)", r"c_{\\mathrm{\1},\2}", input_str)
        input_str = re.sub(r"c_\{([^}]*)\}_\{([^}]*)\}", r"c_{\1,\2}", input_str)
        input_str = re.sub(r"c_\{([^}]*)\}_([A-Za-z0-9]+)", r"c_{\1,\2}", input_str)
    else:
        raise ValueError(f"Format {var_format} is not supported. (expected 'CADET' or 'Legacy')")

    return input_str


__all__ = ["format_variables"]
