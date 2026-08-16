"""
Documentation that can go stale silently.

A wrong number in a docstring is caught by ``test_doctests``; a *missing* option in the
README is caught by nothing, and the README had drifted to documenting 34 of 83 flags
plus one, ``--fresnel_buffer``, that no longer existed. Prose is not testable, but its
coverage of a generated surface is, so that is what these pin.
"""

import argparse
import io
import os
import re
import sys
import unittest
from contextlib import redirect_stdout

from _support import REPO_ROOT  # noqa: F401  (also sets up sys.path)

from oroscope import site_searcher as ss


def cli_flags():
    """Every long option the tool accepts, read from the parser rather than a list."""
    captured = {}
    real = argparse.ArgumentParser.parse_args

    def intercept(self, *args, **kwargs):
        captured["parser"] = self
        raise SystemExit(0)

    argparse.ArgumentParser.parse_args = intercept
    argv = sys.argv
    sys.argv = ["site_searcher.py"]
    try:
        with redirect_stdout(io.StringIO()):
            ss.main()
    except SystemExit:
        pass
    finally:
        argparse.ArgumentParser.parse_args = real
        sys.argv = argv

    return {a.option_strings[0] for a in captured["parser"]._actions
            if a.option_strings and a.dest != "help"}


def read(*parts):
    with open(os.path.join(REPO_ROOT, *parts), encoding="utf-8") as f:
        return f.read()


class TestTheCliPageDocumentsTheCli(unittest.TestCase):
    """
    ``docs/source/cli.rst`` is the canonical option reference; drift in either
    direction is a defect. It moved there from the README, which is now
    library-first — but the reference still has to be complete wherever it lives.
    """

    @classmethod
    def setUpClass(cls):
        cls.flags = cli_flags()
        cls.page = read("docs", "source", "cli.rst")
        cls.documented = set(re.findall(r"``(--[a-z_0-9]+)``", cls.page))

    def test_every_option_is_documented(self):
        missing = sorted(self.flags - self.documented)
        self.assertEqual(missing, [], f"cli.rst does not document: {missing}")

    # Options belonging to the other console scripts and to the repository's own
    # full-DEM runner, which this page documents too.
    OTHER_TOOLS = {"--only", "--dry-run", "--north", "--south", "--west", "--east",
                   "--labels", "--out", "--mode", "--require", "--no_image", "--sweep",
                   "--keep_runs", "--open_topography_api_key",
                   # oroscope-fetch-roads
                   "--places", "--places_only", "--bbox", "--classes", "--step_deg",
                   # oroscope-combine
                   "--roads"}

    def test_no_option_is_documented_that_does_not_exist(self):
        """``--fresnel_buffer`` outlived the code by some months."""
        phantom = sorted(self.documented - self.flags - self.OTHER_TOOLS)
        self.assertEqual(phantom, [], f"cli.rst documents non-existent flags: {phantom}")

    def test_the_option_table_carries_types_and_defaults(self):
        """A reference without defaults sends the reader to the source anyway."""
        self.assertIn("- Default", self.page)
        self.assertIn("- What it does", self.page)

    def test_the_precedence_order_is_stated_the_way_the_code_resolves_it(self):
        """
        The command line wins. It used to lose to the config file, silently, and the
        documentation described the losing behaviour long after it was fixed.
        """
        section = self.page[self.page.index("Where things are resolved from"):]
        cli = section.lower().index("option you actually typed")
        config = section.lower().index("configuration file")
        self.assertLess(cli, config,
                        "the page must state that a typed option beats the config")

    def test_the_readme_leads_with_code_rather_than_the_command_line(self):
        """
        Most people use this as a library, so the README's first example should be a
        call rather than a shell line.
        """
        readme = read("README.md")
        body = readme[readme.index("# Oroscope"):]
        first_python = body.index("```python")
        first_usage = body.index("oroscope --config_path")
        self.assertLess(first_python, first_usage,
                        "the README should show a call before a command line "
                        "(a `pip install` line before either is fine)")


class TestEveryModuleIntroducesItself(unittest.TestCase):
    """
    ``functions.rst`` documents each module with ``automodule``, which renders the
    module docstring as the section's introduction. A module without one contributes a
    bare list of functions to the published API reference.
    """

    MODULES = ("oroscope", "oroscope.site_searcher", "oroscope.arrival_scan",
               "oroscope.physics", "oroscope.scoring", "oroscope.aperture",
               "oroscope.explain", "oroscope.combine_experiments",
               "oroscope.crop_dem", "oroscope.sensitivity", "oroscope.figures",
               "oroscope.fetch_dem", "oroscope.generate_env")

    def test_all_modules_have_a_docstring(self):
        import importlib
        bare = []
        for name in self.MODULES:
            module = importlib.import_module(name)
            if not (module.__doc__ or "").strip():
                bare.append(name)
        self.assertEqual(bare, [], f"modules with no docstring: {bare}")

    def test_the_docstrings_say_something(self):
        """A one-word docstring satisfies the check above and helps nobody."""
        import importlib
        thin = []
        for name in self.MODULES:
            module = importlib.import_module(name)
            if len((module.__doc__ or "").split()) < 8:
                thin.append(name)
        self.assertEqual(thin, [], f"modules with a near-empty docstring: {thin}")


class TestTheNotebooksCallRealNames(unittest.TestCase):
    """
    Notebooks are executed in CI, which is what makes their claim to work checkable —
    except notebooks 7 and 8, which drive whole searches and read stored full-DEM
    results, and are excluded because executing them on every push costs far more than
    it checks.

    So the drift they are exposed to is checked here instead, statically: every
    ``ss.<name>`` and ``explain.<name>`` the generator writes must still exist. That is
    the failure mode a rename produces, and it is the one the excluded execution would
    otherwise have caught.
    """

    @classmethod
    def setUpClass(cls):
        cls.source = read("tools", "make_notebooks.py")

    def attributes(self, alias):
        return set(re.findall(rf"\b{alias}\.([a-zA-Z_][a-zA-Z_0-9]*)", self.source))

    def test_every_site_searcher_name_exists(self):
        missing = sorted(n for n in self.attributes("ss") if not hasattr(ss, n))
        self.assertEqual(missing, [],
                         f"notebooks call site_searcher names that do not exist: {missing}")

    def test_every_explain_name_exists(self):
        from oroscope import explain
        missing = sorted(n for n in self.attributes("explain") if not hasattr(explain, n))
        self.assertEqual(missing, [],
                         f"notebooks call explain names that do not exist: {missing}")

    def test_every_generated_notebook_is_committed(self):
        """A generator entry with no notebook beside it means someone forgot to run it."""
        names = set(re.findall(r'"(\d\d_[a-z_]+\.ipynb)":', self.source))
        self.assertTrue(names, "found no notebook names in the generator")
        for name in sorted(names):
            path = os.path.join(REPO_ROOT, "notebooks", name)
            self.assertTrue(os.path.exists(path), f"{name} is generated but not committed")

    def test_the_excluded_notebooks_are_still_excluded_on_purpose(self):
        """
        If either ever becomes cheap enough to execute, these are the lines to delete —
        but silently dropping the exclusion while they are still expensive would slow
        every push by an hour and a half.
        """
        workflow = read(".github", "workflows", "lint.yml")
        for name in ("07_explaining_a_run.ipynb", "08_the_full_dem.ipynb"):
            self.assertIn(name, workflow,
                          f"{name} must be named in the CI workflow, excluded or not")


class TestTheDocumentedPublicSurfaceIsReal(unittest.TestCase):
    """``__all__`` drives autodoc, so a name that has gone stale is a broken page."""

    def test_everything_in_all_exists(self):
        for name in ss.__all__:
            self.assertTrue(hasattr(ss, name), f"site_searcher.__all__ names {name}, "
                                               f"which does not exist")

    def test_explain_all_exists(self):
        from oroscope import explain
        for name in explain.__all__:
            self.assertTrue(hasattr(explain, name), f"explain.__all__ names {name}")


if __name__ == "__main__":
    unittest.main()
