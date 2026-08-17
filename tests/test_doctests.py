"""
Runs the ``Examples`` blocks in the docstrings as tests.

The documentation executes those blocks at build time with jupyter-sphinx, so a stale
example becomes a visibly wrong number on a published page. Running them here as well
means the suite catches it first, on every supported Python, without needing Sphinx
installed.

This is not hypothetical: five of the first batch of examples written for
:mod:`physics` quoted values that were close to, but not, what the code returns. Each
was a plausible hand-computed figure. The suite caught all five.
"""

import doctest
import importlib
import unittest

import _support  # noqa: F401  (path setup)

# Modules whose docstring examples are meant to run. Deliberately a list rather than a
# directory walk: the CLI tools take a DEM and write files, so their examples are
# illustrative rather than executable, and figures builds plots.
#
# site_searcher is included for the routines that take plain arrays -- the terrain
# derivatives, the slope band, the capacity lattice. Its pipeline stages need a DEM on
# disk and are exercised by the regression tests instead.
#
# **No example may open a path relative to the working directory.** CI runs this job
# with `working-directory: tests`, so `data/whatever.csv` resolves under tests/ there
# and under the repository root when run by hand -- passing locally and failing in CI,
# which is exactly what happened. Build the path from the module's own __file__, as the
# aperture examples do, or write an example that needs no file at all.
# combine_experiments, figures and fetch_roads were added after 41 of their examples
# turned out never to have been executed, and two of them were wrong: geographic_extent
# documented a bottom edge of -16.005 where the code returns -15.995, which would make
# a 100-row raster of 0.01 degree pixels 1.01 degrees tall, and canyon_geometry expected
# `7.6` from a value NumPy 2 reprs as `np.float64(7.6)`. Both are published in the API
# reference, and both had been read many times without being run.
#
# figures is here despite building plots: its examples assert on figure *structure* --
# axis counts, figure widths -- which needs no display, and `_support` has already
# forced MPLBACKEND to Agg by the time these import. Only the modules whose examples
# genuinely cannot run headless stay out.
DOCTESTED = ("oroscope.physics", "oroscope.scoring", "oroscope.arrival_scan",
             "oroscope.aperture", "oroscope.explain", "oroscope.site_searcher",
             "oroscope.combine_experiments", "oroscope.figures",
             "oroscope.fetch_roads")

OPTIONS = doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS


def load_tests(loader, tests, ignore):
    """unittest hook: adds one DocTestSuite per module listed above."""
    for name in DOCTESTED:
        module = importlib.import_module(name)
        tests.addTests(doctest.DocTestSuite(module, optionflags=OPTIONS))
    return tests


class TestEveryModuleIsCovered(unittest.TestCase):
    """A module gaining examples should gain them here too."""

    def test_listed_modules_all_import(self):
        for name in DOCTESTED:
            self.assertTrue(importlib.import_module(name), name)

    def test_physics_examples_exist(self):
        """Guards against the list silently covering nothing."""
        from oroscope import physics
        found = doctest.DocTestFinder().find(physics)
        with_examples = [t for t in found if t.examples]
        self.assertGreater(len(with_examples), 5,
                           "physics should carry runnable examples")


if __name__ == "__main__":
    unittest.main()
