Oroscope: terrain site search for particle-astrophysics observatories
=====================================================================

.. image:: https://github.com/mbustama/oroscope/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/mbustama/oroscope/actions/workflows/tests.yml
   :alt: tests

.. image:: https://github.com/mbustama/oroscope/actions/workflows/lint.yml/badge.svg
   :target: https://github.com/mbustama/oroscope/actions/workflows/lint.yml
   :alt: Code Quality

.. image:: https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg
   :target: https://mbustama.github.io/oroscope/
   :alt: Documentation

.. image:: https://img.shields.io/pypi/v/oroscope.svg
   :target: https://pypi.org/project/oroscope/
   :alt: PyPI

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: License: GPL v3

.. image:: https://img.shields.io/badge/python-3.9+-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python 3.9+

.. image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Code style: ruff

**Oroscope** searches digital elevation models for ground that can host a
particle-astrophysics observatory. Greek *oros*, mountain, and *skopein*, to look at.

It answers one structural question, and that is what lets a single engine serve
experiments that look nothing alike:

   From this patch of ground, is there a target surface at the right range, in the
   right direction, at the right relative orientation, with the right matter behind it?

GRAND wants terrain a few degrees below the horizon and tens of kilometres away, to
catch radio from air showers started by Earth-skimming tau neutrinos. TAMBO wants a
canyon wall two to five kilometres across, to catch the particles themselves. They
differ in their numbers, not in their structure — so both are configurations of the
same scan, and a search can report where each is viable and where the two coincide.

.. jupyter-execute::
   :hide-code:
   :hide-output:

   # As physics.rst: jupyter-sphinx runs these blocks in a plain kernel, which renders
   # a Figure as its text repr rather than as a PNG unless the inline backend is
   # switched on. Without this the diagram below was the literal text
   # "<Figure size 1020x415 with 2 Axes>" -- on the front page, in a build that
   # reported no error.
   %matplotlib inline
   %config InlineBackend.figure_formats = ['svg']

.. jupyter-execute::
   :hide-code:

   from oroscope import figures
   # Assigned, not left bare: the inline backend flushes the figure at the end of the
   # cell *and* the returned Figure renders its own repr, so calling this bare put the
   # same diagram on the page twice. physics.rst uses the same form for the same reason.
   _ = figures.walk_mechanism()

What it does
------------

#. **Screens terrain** by slope, aspect, altitude and exclusion zones.
#. **Scans arrival directions** from every surviving pixel, finding where a tau could
   exit and how much rock lies behind it — one profile walk per (candidate, azimuth),
   serving every elevation bin at once.
#. **Scores** the result against per-experiment criteria, each a named component.
#. **Cleans up** morphologically, labels sites, and places detectors on a lattice in
   real ground coordinates.
#. **Writes** GeoTIFF, KML, PNG and JSON, plus a selection funnel, a provenance record,
   and a plain-language account of what was found and why.

What it is not
--------------

* **Not a shower simulation.** Development comes from a Gaisser–Hillas profile and an
  isothermal atmosphere. No Monte Carlo, no detector response, no trigger model.
* **Not a flux or exposure calculation.** It reports geometry and acceptance-shaped
  quantities; turning those into an event rate needs a flux and a response table.
* **Not a substitute for a site visit.** There is no slope stability, land access,
  power or cost.
* **Not externally validated.** The physics is checked against closed-form synthetic
  terrain and against itself. It has not been compared with any collaboration's own
  simulated acceptance. :doc:`assumptions` is blunt about what that means.

Getting started
---------------

.. code-block:: shell

   pip install oroscope

or, for a clone with the notebooks, the configurations and the test suite:

.. code-block:: shell

   git clone https://github.com/mbustama/oroscope.git
   cd oroscope
   pip install -e .

.. code-block:: shell

   oroscope --config_path config/grand_colca_config.json

See :doc:`installation` and :doc:`quickstart` for the longer version, and
:doc:`physics` for what is actually being computed.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   installation
   quickstart
   howitworks
   cli
   data
   notebooks
   physics
   assumptions
   glossary
   implementation
   functions
   references

Author
------

Written by Mauricio Bustamante (mbustamante@gmail.com). Bug reports and questions are
best raised as `GitHub issues <https://github.com/mbustama/oroscope/issues>`_.

License
-------

Released under the `GNU General Public License v3
<https://www.gnu.org/licenses/gpl-3.0>`_. The full text ships with the source, as
``LICENSE`` in the repository root.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
