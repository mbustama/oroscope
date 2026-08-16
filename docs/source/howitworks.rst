How the search works
====================

This page introduces the vocabulary. **Screening**, **striding**, **the arrival scan**,
**scoring**, **closing**, **pruning** — these words appear throughout the documentation,
the code and the run summaries, and none of them is guessable from ordinary English.
Read this before anything else that quotes a number.

.. contents::
   :local:
   :depth: 2


The one question
----------------

Everything the tool does reduces to a single question, asked once per pixel of an
elevation model:

   From this patch of ground, is there a **target surface** at the **right range**, in
   the **right direction**, at the **right relative orientation**, with the **right
   matter behind it**?

That is what makes one engine serve two experiments that look nothing alike. GRAND wants
terrain a few degrees below the horizon and tens of kilometres away, to catch radio from
air showers started by Earth-skimming tau neutrinos. TAMBO wants a canyon wall two to
five kilometres across, to catch the particles themselves. **They differ in their
numbers, not in their structure.**


The seven stages
----------------

.. jupyter-execute::
   :hide-code:

   from oroscope import figures
   figures.pipeline_stages()

The figure is a real run — TAMBO over the full Ancash DEM, 68.6 million pixels — and the
important thing about it is that it has **two halves that work in opposite directions**.

Everything down to *scoring* **removes** candidates. Everything below it **rebuilds a map
from what survived**, which is why the count rises again at closing. Reading a funnel
table as though every row were a filter is the commonest way to misread one.

Stage by stage:

**1. Screening.** A cheap per-pixel test: slope inside the experiment's deployable band,
optionally altitude, aspect, distance to a road, and distance from a radio-noise zone.
Slope is the one that matters — GRAND needs 3–25°, TAMBO 20–60°, and a pixel has only one
slope, which is why co-location is decided here more than anywhere else.

**2. Striding.** Of the pixels that survive screening, one in ``candidate_stride`` is
kept. This is **cost control, not a criterion**: the arrival scan is the expensive stage
and its cost is linear in the number of candidates. Striding is unbiased in *acceptance*
— measured identical at strides 1 and 5 — but it has a consequence at the far end of the
pipeline that is easy to miss. See :ref:`striding-and-closing`.

**3. The arrival scan.** The heart of the tool. From each surviving candidate, walk
outward along several bearings and, for each, find where the ray first meets terrain. A
direction is *accepted* when that intersection lies in the configured distance window, the
elevation window, and (optionally) strikes ground steep enough. One walk fills every
elevation bin at once, so the **azimuth count sets the cost** and the elevation binning is
nearly free.

**4. Scoring.** Each accepted candidate is given named component scores in [0, 1] —
column depth, exit distance, accepted solid angle, shower development, tau decay, and
others as configured — which are combined into one number and cut at ``min_score``. See
:ref:`why-a-product-is-a-cliff`.

**5. Closing.** A morphological operation that fills holes. Its job is to undo striding:
turn a lattice of isolated accepted marks back into a region a detector array could
actually occupy.

**6. Pruning and selection.** Connected regions are labelled, then dropped if they are
too narrow (``min_width_km``), too small, or hold too few detector positions
(``min_sub_array_size``). What remains is the answer.


.. _striding-and-closing:

Striding and closing are one decision, not two
----------------------------------------------

.. jupyter-execute::
   :hide-code:

   from oroscope import figures
   figures.striding_and_closing()

Striding leaves a gap of ``candidate_stride`` pixels between kept candidates. Closing
repairs it — but only if its structuring element is **larger than that gap**. Below the
gap the marks never touch and the mask stays a scatter of isolated pixels; above it the
region reappears almost intact.

**The transition is at the gap and it is abrupt.** In the figure a 3-pixel element
recovers 0.04× of the accepted set and a 5-pixel one recovers 0.68× — seventeen times
more, for two pixels of element.

This is not a hypothetical. It is the mechanism behind a real **4.75× under-report** of
TAMBO's area at Colca, and a **291× one** on the steeper ground of the Callejón de
Huaylas, where the accepted strips are narrower still. GRAND never suffers it: its 1 km
closing element bridges a 154 m stride gap without noticing.

.. warning::

   The damage is not done to the area measurement. It is done in **pruning**: a
   fragmented mask becomes thousands of tiny regions, and almost none of them clears
   ``min_sub_array_size``. Every run warns when the element cannot bridge the gap
   (:func:`~oroscope.site_searcher.warn_stride_outruns_closing`). **Heed it, raise
   ``gap_close_km``, or read the area as a lower bound and say so.**


.. _why-a-product-is-a-cliff:

Why a product score has no safe threshold
------------------------------------------

.. jupyter-execute::
   :hide-code:

   from oroscope import figures
   figures.score_composition()

Components are combined by **multiplication**, so a candidate must be good at everything
and one bad component sinks it. That is defensible — these are requirements, not
preferences — but it has a consequence for any threshold placed on the result.

A product of several numbers in [0, 1] piles up near zero however good the terrain is.
So ``min_score`` does not express a mild preference: it sits on a cliff, and **where it
lands depends on how many components happen to be enabled**. Adding a component moves
every score down and silently tightens the cut.

Measured on one real search, cuts of 0.0, 0.35 and 0.5 gave 45,928, 2,056 and **zero**
detector positions. ``score_percentile`` is the scale-free alternative and is
scored against the same ranking; see :doc:`assumptions`.


What the answer is, and is not
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - It is
     - It is not
   * - A **ranking** of ground by how well it satisfies named geometric and physical
       criteria.
     - An **event rate**. No flux, no cross-section and no detector response are folded
       in. The scores are not apertures.
   * - A **capacity estimate** for an arbitrarily placed array.
     - An optimised layout. Detectors are placed from each region's bounding-box corner,
       not fitted — see :doc:`implementation`.
   * - Reproducible from a configuration file and a DEM.
     - Validated against an external simulation. Nothing here has been.

:doc:`assumptions` is the blunt version of this table, and every run prints its own.


Where to go next
----------------

- :doc:`glossary` — every term on this page, defined, with links to where it is used.
- :doc:`notebooks` — notebook 7 animates the mechanisms above; notebook 8 drives the
  pipeline and reads its output; notebooks 9 to 12 are real regions.
- :doc:`physics` — what the geometry gets multiplied by.
- :doc:`implementation` — how it is actually built, for anyone changing it.
