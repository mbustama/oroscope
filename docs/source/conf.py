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

sys.path.insert(0, os.path.abspath('../../src'))

# Importing site_searcher pulls in matplotlib.pyplot, which will try to reach for a
# display if one is not ruled out first. autodoc imports every documented module, so
# this has to happen before Sphinx does.
os.environ.setdefault('MPLBACKEND', 'Agg')

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

master_doc = 'index'
templates_path = ['_templates']
exclude_patterns = ['_build']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_title = 'Oroscope %s' % version

html_theme_options = {
    'logo_only': False,
}
