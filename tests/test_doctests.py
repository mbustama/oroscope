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
DOCTESTED = ("oroscope.physics", "oroscope.scoring", "oroscope.arrival_scan",
             "oroscope.aperture", "oroscope.explain", "oroscope.site_searcher")

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
