Tutorial notebooks
==================

Twelve worked notebooks live in `notebooks/
<https://github.com/mbustama/oroscope/tree/main/notebooks>`_, numbered in reading
order. Each carries its figures inline, so they can be read on GitHub without being
run, and each ends with links to the previous and next.

**None of them needs a real DEM.** Every notebook builds its own terrain, because the
DEMs are a quarter of a gigabyte and are not in the repository — a tutorial only its
author can run is not a tutorial. Notebook 7 will *use* the Colca DEM for three of its
eight animations when it is there, because those three are about what a criterion does
to real ground, but it falls back to synthetic terrain and says on each figure which it
used.

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

Seeing it move
--------------

`7. Animating the mechanism <https://github.com/mbustama/oroscope/blob/main/notebooks/07_animating_the_mechanism.ipynb>`_
   The eight animations in ``tools/make_animations.py``, built and explained — what
   each one argues, and why the six rejected candidates are better served by a static
   figure. Ends with how to turn the MP4s into animated GIFs, both through the tool and
   with plain ``ffmpeg``.

Running one
-----------

`8. Explaining a run <https://github.com/mbustama/oroscope/blob/main/notebooks/08_explaining_a_run.ipynb>`_
   The pipeline as an ordinary Python call, and how to read what it says: configuration
   as data, the memory pre-flight, the results dictionary it hands back, the funnel and
   its binding constraint, which sites are actually in the result, why each one is good
   and what held it back. Two searches — one that finds ground and one that finds none,
   because the empty result is the case a bare results file serves worst.

Real ground, region by region
-----------------------------

These share a numbering block so a new region slots in without renumbering its
neighbours. Each **reads** results produced locally by ``tools/run_full_dem.py`` rather
than running the searches itself.

`9. Arequipa, the full DEM <https://github.com/mbustama/oroscope/blob/main/notebooks/09_arequipa_dem.ipynb>`_
   The run that had never been done: GRAND alone, TAMBO alone, and the combination,
   over the whole DEM rather than a crop, at 1 arc-second.

`10. Ancash <https://github.com/mbustama/oroscope/blob/main/notebooks/10_ancash_dem.ipynb>`_
   The same three searches over the Cordillera Blanca and the Callejón de Huaylas, at
   the same resolution and with every transferable criterion unchanged — so a
   difference in the answer is a difference in the ground rather than in the question.

*11. Lima* — not written yet; the number is reserved.

`12. Peru, all of it <https://github.com/mbustama/oroscope/blob/main/notebooks/12_peru_dem.ipynb>`_
   The search run over a whole country — 339 million pixels at 3 arc-seconds, in four
   minutes. As much about reading a coarse answer honestly as about Peru: why the
   resolution is forced rather than chosen, why the area is a bracket and the site
   count is worse than the area, and one worry that was checked and came back fine.

Going deeper
------------

`13. Turning the knobs <https://github.com/mbustama/oroscope/blob/main/notebooks/13_turning_the_knobs.ipynb>`_
   Each parameter turned, and what moves. Two of the five sit on cliffs and a third
   changes one reported number while leaving another alone. The striding penalty is
   shown to depend on **how thin the accepted ground is**, which is why the same knob
   costs 1.51× at Colca and 23.0× on the Callejón de Huaylas. Runs in seconds on terrain
   it builds itself.

.. note::

   Notebooks 7 through 10 are not executed in CI. Seven renders eight animations and
   wants an ``ffmpeg`` the runner does not have; eight drives whole searches; nine and
   ten read stored full-DEM results that take tens of minutes to produce. They are run
   locally when a configuration changes; ``tests/test_docs.py`` checks statically that
   the API names they call still exist, which is the drift the execution would have
   caught.

   Notebook 12 *is* executed, despite reading a store it will not find there. Every
   cell that needs the store is guarded and every other cell computes from the library,
   so on a bare runner it degrades to the prose and the arithmetic rather than failing —
   which is the behaviour worth having checked on every push.
