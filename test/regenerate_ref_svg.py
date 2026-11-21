"""Regenerate the reference figures in ``test/data/ref_svg``.

Run from the repository root with the repository on the path::

    PYTHONPATH=. python test/regenerate_ref_svg.py            # all models
    PYTHONPATH=. python test/regenerate_ref_svg.py GRM CSTR   # just these

Then check ``git status``: the files that change are exactly the models the
edit affected, which is the cheapest way to see its blast radius.
"""

import os
import sys

sys.path[:0] = [os.path.abspath("."), os.path.abspath("test")]

from test_refFigures import FIGURE_MODELS, render_figures  # noqa: E402

TEST_DIR = os.path.abspath("test")
REFERENCE_DIR = os.path.join(TEST_DIR, "data", "ref_svg")


def main(names: list[str]) -> None:
    os.makedirs(REFERENCE_DIR, exist_ok=True)
    for name in names or FIGURE_MODELS:
        document = render_figures(name, TEST_DIR)
        path = os.path.join(REFERENCE_DIR, name + ".svg")
        # Written with the platform line ending, as the reference LaTeX is;
        # the test reads both back in text mode, so either one compares equal.
        with open(path, "w", encoding="utf-8") as file:
            file.write(document)
        print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1:])
