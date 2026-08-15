Tutorial notebooks
==================

Six worked notebooks live in `notebooks/
<https://github.com/mbustama/oroscope/tree/main/notebooks>`_, numbered in reading
order. Each carries its figures inline, so they can be read on GitHub without being
run, and each ends with links to the previous and next.

**None of them needs a real DEM.** Every notebook builds its own terrain, because the
DEMs are a quarter of a gigabyte and are not in the repository — a tutorial only its
author can run is not a tutorial.

To run them rather than read them:

.. code-block:: shell

   pip install "oroscope[notebooks]"
   jupyter lab notebooks/

.. note::

   The notebooks are not built into this documentation. Executing them on every docs
   build would be slow, and they are more useful where their outputs are already
   stored. They are executed in CI on every push, which is what makes their claim to
   work checkable rather than merely plausible.

Start here
----------

`1. Getting started <https://github.com/mbustama/oroscope/blob/main/notebooks/01_getting_started.ipynb>`_
   The question the tool answers, a piece of terrain, one scan, and what the
   observables mean.

`2. The arrival scan <https://github.com/mbustama/oroscope/blob/main/notebooks/02_the_arrival_scan.ipynb>`_
   How one profile walk fills every elevation bin at once, why azimuths rather than
   bins set the cost, and why the ground at your own feet blocks every steep
   downward ray.

The physics
-----------

`3. The physics toolkit <https://github.com/mbustama/oroscope/blob/main/notebooks/03_physics_toolkit.ipynb>`_
   Tau decay lengths, shower development in grammage rather than metres, the Earth
   chord, and geomagnetic emission — all usable with no terrain at all.

`4. Criteria and scoring <https://github.com/mbustama/oroscope/blob/main/notebooks/04_criteria_and_scoring.ipynb>`_
   The three score shapes and their composition, plus two traps: a saturating score
   whose scale is wrong stops discriminating, and a product score has no safe
   threshold.

Two experiments
---------------

`5. GRAND and TAMBO <https://github.com/mbustama/oroscope/blob/main/notebooks/05_grand_and_tambo.ipynb>`_
   The same engine, two experiments, and exactly which numbers differ. Recovers a
   canyon's wall slope from the terrain it was built with.

`6. Combining and sensitivity <https://github.com/mbustama/oroscope/blob/main/notebooks/06_combining_and_sensitivity.ipynb>`_
   Joint, union and co-location, then the more important question: how much a result
   depends on its assumptions. Read this one before quoting a capacity.
