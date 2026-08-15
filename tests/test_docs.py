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

import site_searcher as ss


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


class TestReadmeDocumentsTheCli(unittest.TestCase):
    """The option tables are the reference; drift in either direction is a defect."""

    @classmethod
    def setUpClass(cls):
        cls.flags = cli_flags()
        cls.readme = read("README.md")
        cls.documented = set(re.findall(r"`(--[a-z_0-9]+)`", cls.readme))

    def test_every_option_is_documented(self):
        missing = sorted(self.flags - self.documented)
        self.assertEqual(missing, [], f"README does not document: {missing}")

    def test_no_option_is_documented_that_does_not_exist(self):
        """``--fresnel_buffer`` outlived the code by some months."""
        phantom = sorted(self.documented - self.flags)
        self.assertEqual(phantom, [], f"README documents non-existent flags: {phantom}")

    def test_the_precedence_order_is_stated_the_way_the_code_resolves_it(self):
        """
        The command line wins. It used to lose to the config file, silently, and the
        README described the losing behaviour long after it was fixed.
        """
        section = self.readme[self.readme.index("## 4."):]
        section = section[:section.index("### Complete List")]
        cli = section.lower().index("command line")
        config = section.lower().index("config file")
        self.assertLess(cli, config,
                        "the README must state that the command line beats the config")


class TestEveryModuleIntroducesItself(unittest.TestCase):
    """
    ``functions.rst`` documents each module with ``automodule``, which renders the
    module docstring as the section's introduction. A module without one contributes a
    bare list of functions to the published API reference.
    """

    MODULES = ("site_searcher", "arrival_scan", "physics", "scoring", "aperture",
               "explain", "combine_experiments", "crop_dem", "sensitivity",
               "figures", "fetch_dem", "generate_env")

    def test_all_modules_have_a_docstring(self):
        bare = []
        for name in self.MODULES:
            module = __import__(name)
            if not (module.__doc__ or "").strip():
                bare.append(name)
        self.assertEqual(bare, [], f"modules with no docstring: {bare}")

    def test_the_docstrings_say_something(self):
        """A one-word docstring satisfies the check above and helps nobody."""
        thin = []
        for name in self.MODULES:
            module = __import__(name)
            if len((module.__doc__ or "").split()) < 8:
                thin.append(name)
        self.assertEqual(thin, [], f"modules with a near-empty docstring: {thin}")


class TestTheNotebooksCallRealNames(unittest.TestCase):
    """
    Notebooks are executed in CI, which is what makes their claim to work checkable —
    except notebook 7, which drives whole searches and reads stored full-DEM results,
    and is excluded because executing it on every push costs more than it checks.

    So the drift it is exposed to is checked here instead, statically: every
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
        import explain
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

    def test_the_excluded_notebook_is_still_excluded_for_a_reason(self):
        """
        If 07 ever becomes cheap enough to execute, this is the line to delete — but
        silently dropping the exclusion while it is still expensive would slow every
        push by an hour and a half.
        """
        workflow = read(".github", "workflows", "lint.yml")
        self.assertIn("07_running_a_search.ipynb", workflow,
                      "notebook 07 must be named in the CI workflow, excluded or not")


class TestTheDocumentedPublicSurfaceIsReal(unittest.TestCase):
    """``__all__`` drives autodoc, so a name that has gone stale is a broken page."""

    def test_everything_in_all_exists(self):
        for name in ss.__all__:
            self.assertTrue(hasattr(ss, name), f"site_searcher.__all__ names {name}, "
                                               f"which does not exist")

    def test_explain_all_exists(self):
        import explain
        for name in explain.__all__:
            self.assertTrue(hasattr(explain, name), f"explain.__all__ names {name}")


if __name__ == "__main__":
    unittest.main()
