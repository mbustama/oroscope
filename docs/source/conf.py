# -*- coding: utf-8 -*-
r"""Sphinx configuration for the Oroscope documentation.

Build the HTML documentation with::

    pip install -r docs/requirements.txt
    cd docs && make html

The result lands in ``docs/build/html``.

The library modules live in ``src/`` as top-level modules rather than inside a
package, so that directory is put on ``sys.path`` below and ``autodoc`` imports
them by their bare names, exactly as user code does.
"""

import os
import pathlib
import re
import sys

sys.path.insert(0, os.path.abspath('../../src'))   # so `import oroscope` works uninstalled

# Importing site_searcher pulls in matplotlib.pyplot, so Sphinx's own process needs a
# non-interactive backend before autodoc imports anything.
#
# Set on the matplotlib module rather than through the MPLBACKEND environment variable,
# and the distinction is the whole reason the figures on physics.rst were blank. An
# environment variable is inherited by child processes, and jupyter-sphinx runs the
# `jupyter-execute` blocks in a *separate kernel*: Agg reached that kernel, overrode the
# inline backend that captures figures as PNGs, and every diagram silently produced no
# output. The pages built clean and showed nothing.
import matplotlib                                  # noqa: E402

matplotlib.use('Agg')

# -- Project information -----------------------------------------------------

project = 'Oroscope'
copyright = '2026, Mauricio Bustamante'
author = 'Mauricio Bustamante'


# The version lives in pyproject.toml, which is what `python -m build` reads, so that
# file is the only one entitled to state it. Keeping a second copy here would mean two
# hand-edits per release with nothing checking they agree.
#
# Read from the file rather than from package metadata so the documentation builds from
# a source tree that has not been installed, which is what `cd docs && make html` does;
# metadata is the fallback for a build running from an installed distribution.
def _project_version():
    """Returns the version, from pyproject.toml or from installed metadata."""
    pyproject = pathlib.Path(__file__).resolve().parents[2] / 'pyproject.toml'
    try:
        found = re.search(r'^version = "([^"]+)"', pyproject.read_text(), re.M)
    except OSError:
        found = None
    if found is not None:
        return found.group(1)

    from importlib.metadata import version as installed_version

    return installed_version('oroscope')


release = _project_version()
version = '.'.join(release.split('.')[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',       # API reference, generated from the docstrings
    'sphinx.ext.mathjax',       # Renders the LaTeX math in the docstrings
    'sphinx.ext.viewcode',      # Links API entries to highlighted source
    'sphinx.ext.intersphinx',   # Cross-links to the Python, NumPy and SciPy docs
    'numpydoc',                 # Parses the numpydoc-style docstrings
    'sphinx_copybutton',        # Copy-to-clipboard button on code blocks
    'sphinxcontrib.bibtex',     # References page (refs.bib)
    'myst_parser',              # Lets pages .. include:: the Markdown docs
    'jupyter_sphinx',           # Runs the narrative examples at build time
]

# numpydoc would otherwise try to document the members of every class it meets; the
# library exposes module-level routines and a couple of small record types.
numpydoc_show_class_members = False

# numpydoc validates that every parameter in a signature appears in the docstring.
# Left on deliberately: the whole point of this pass was that the parameters were not
# documented, and a check nobody runs is a check that stops being true.
numpydoc_validation_checks = {'PR01', 'PR02', 'PR03'}

bibtex_bibfiles = ['refs.bib']

# Only the inventories actually cross-referenced. Each one costs a network round trip
# per build, and a fetch that fails is a warning, which the -W build treats as fatal.
intersphinx_timeout = 10
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
}

# Keep the source order of the routines, which groups them the way the pipeline runs,
# rather than sorting them alphabetically.
autodoc_member_order = 'bysource'

# The scan kernels are compiled by Numba at import. autodoc only needs their
# signatures and docstrings, both of which survive the decorator.
autodoc_mock_imports = []

# Type annotations are rendered in the signature line, where they read as part of the
# call. The docstrings state units and meaning, which the annotation cannot; the two
# are complementary rather than duplicates -- `float` and "elevation angle in degrees,
# positive upward" answer different questions.
autodoc_typehints = 'signature'
autodoc_typehints_format = 'short'

# `from __future__ import annotations` makes every annotation a string, so autodoc
# needs the module's globals to resolve them back into links.
autodoc_type_aliases = {}

master_doc = 'index'
templates_path = ['_templates']
exclude_patterns = ['_build']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_title = 'Oroscope %s' % version

# `logo_only = False` keeps html_title -- and so the version -- visible beneath the
# image rather than letting it stand alone. The logo is an SVG, which sphinx-rtd-theme
# serves as-is: it stays sharp at any sidebar width, and it is transparent outside its
# disc, so it sits correctly on the theme's background without a matte.
html_logo = '_static/oroscope_logo.svg'

html_theme_options = {
    'logo_only': False,
}
