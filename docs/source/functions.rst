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

.. automodule:: physics
   :members:
   :undoc-members:

The arrival scan
----------------

The scan kernel: profile walking, column depth, Fresnel clearance and radio-noise
line-of-sight. Compiled with Numba; see :doc:`physics` for what it computes.

.. automodule:: arrival_scan
   :members:
   :undoc-members:

Scoring
-------

Score shapes — band, saturating, ramp — and their composition into a single figure of
merit. Read the warning in :doc:`assumptions` about thresholding a product before
choosing ``min_score``.

.. automodule:: scoring
   :members:
   :undoc-members:

Aperture
--------

Aperture estimate, tabulated response, and inference of a response curve from a
published one.

.. automodule:: aperture
   :members:
   :undoc-members:

The pipeline
------------

Screening, morphology, capacity and outputs, plus the command-line entry point.

.. automodule:: site_searcher
   :members:
   :undoc-members:

Tools
-----

.. automodule:: crop_dem
   :members:

.. automodule:: combine_experiments
   :members:

.. automodule:: sensitivity
   :members:

Figures
-------

The schematics used throughout this documentation, as importable functions so they can
be restyled and reused.

.. automodule:: figures
   :members:
