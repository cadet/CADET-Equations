"""Regenerate the reference figures in ``test/data/ref_svg/<model>/``.

Run from the repository root with the repository on the path::

    PYTHONPATH=. python test/regenerate_ref_svg.py            # all models
    PYTHONPATH=. python test/regenerate_ref_svg.py GRM CSTR   # just these

Then check ``git status``: the files that change are exactly the models the
edit affected, which is the cheapest way to see its blast radius.
"""

import os
import sys

sys.path[:0] = [os.path.abspath("."), os.path.abspath("test")]

from test_refFigures import FIGURE_MODELS, reference_dir, render_figures  # noqa: E402

TEST_DIR = os.path.abspath("test")


def main(names: list[str]) -> None:
    for model_name in names or FIGURE_MODELS:
        figures = render_figures(model_name, TEST_DIR)
        directory = reference_dir(model_name, TEST_DIR)
        os.makedirs(directory, exist_ok=True)

        # A figure the model no longer draws would otherwise stay behind and
        # fail the test, which compares the directory listing as well.
        for stale in os.listdir(directory):
            if stale.endswith(".svg") and stale[: -len(".svg")] not in figures:
                os.remove(os.path.join(directory, stale))
                print("removed", os.path.join(directory, stale))

        for figure_name, figure in figures.items():
            path = os.path.join(directory, figure_name + ".svg")
            # Written with the platform line ending, as the reference LaTeX is;
            # the test reads both back in text mode, so either one compares equal.
            with open(path, "w", encoding="utf-8") as file:
                file.write(figure)
            print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1:])
