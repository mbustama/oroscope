API reference
=============

Generated from the docstrings. Every ``Examples`` block is **executed when this page is
built**, so what is shown is what the code returns rather than numbers written beside
it.

Physics
-------

Closed-form physics with no terrain in it: atmosphere, Earth chord, tau range and exit
probability, shower development, geomagnetic field, Cherenkov footprint. Self-contained
and usable on its own.

.. automodule:: oroscope.physics
   :members:
   :undoc-members:

The arrival scan
----------------

The scan kernel: profile walking, column depth, Fresnel clearance and radio-noise
line-of-sight. Compiled with Numba; see :doc:`physics` for what it computes.

.. automodule:: oroscope.arrival_scan
   :members:
   :undoc-members:

Scoring
-------

Score shapes — band, saturating, ramp — and their composition into a single figure of
merit. Read the warning in :doc:`assumptions` about thresholding a product before
choosing ``min_score``.

.. automodule:: oroscope.scoring
   :members:
   :undoc-members:

Aperture
--------

Aperture estimate, tabulated response, and inference of a response curve from a
published one.

.. automodule:: oroscope.aperture
   :members:
   :undoc-members:

The pipeline
------------

Screening, morphology, capacity and outputs, plus the command-line entry point.

.. automodule:: oroscope.site_searcher
   :members:
   :undoc-members:

Explaining a run
----------------

Turning a results dictionary into an account of what was found and why: which
constraint set the size of the answer, what held the surviving sites back, and which
of the numbers are assumptions. Pure, so an old results file can be explained months
later with no DEM and no pipeline.

.. automodule:: oroscope.explain
   :members:
   :undoc-members:

Tools
-----

.. automodule:: oroscope.crop_dem
   :members:

.. automodule:: oroscope.combine_experiments
   :members:

.. automodule:: oroscope.sensitivity
   :members:

Figures
-------

The schematics used throughout this documentation, as importable functions so they can
be restyled and reused.

.. automodule:: oroscope.figures
   :members:


Fetching data
-------------

The one-shot tools that bring a region's inputs onto disk are setup rather than library,
which is why ``import oroscope`` does not re-export them and why their command-line side
is documented on :doc:`cli` instead. Only the names other pages link into are listed
here: :data:`oroscope.fetch_dem.REGIONS` is cross-referenced from :doc:`data`, and a
cross-reference with no target renders as plain text rather than failing the build, so it
sat there unresolved and unnoticed.

.. automodule:: oroscope.fetch_dem
   :members: REGIONS, download_dem, generate_and_patch_config
