Implementation notes
====================

For someone changing the code, or trying to work out why it is shaped the way it is.
Nothing here is needed to *use* the tool — that is :doc:`quickstart` and :doc:`cli` —
and nothing here is physics, which is :doc:`physics`. This is the rest: the decisions
that are invisible from outside and expensive to rediscover.

Almost every entry below exists because something went wrong once. The measurements are
in ``docs/ROADMAP.md`` §6, which is the long form of this page.

.. contents::
   :local:
   :depth: 2


The shape of a run
------------------

One process, six phases, and the expensive one is third:

.. list-table::
   :header-rows: 1
   :widths: 18 20 62

   * - Phase
     - Cost
     - Note
   * - Load DEM
     - seconds
     - The GeoTIFF is converted once to a ``.npy`` beside it and thereafter
       **memory-mapped**, so the kernel can evict it. That is why the memory estimate
       excludes it.
   * - Screening
     - seconds
     - Tiled, ``tile_size`` at a time. Tests the *squared* gradient, so neither a square
       root nor an arctangent is formed over the tile.
   * - Arrival scan
     - **minutes to tens of minutes**
     - Numba, ``parallel=True``. Linear in candidates × azimuths. Elevation bins are
       nearly free.
   * - Morphology
     - seconds
     - Closing and pruning, tiled with a halo.
   * - Capacity
     - seconds
     - Lattice packing per region.
   * - Save and draw
     - seconds
     - **A separate memory peak.** See below.


Memory, which is the thing that actually breaks
------------------------------------------------

Three separate quantities get confused, and confusing them has cost this project a
machine.

**1. The search's anonymous memory.** What
:func:`~oroscope.site_searcher.estimate_peak_memory_gb` models: the per-candidate
observable arrays and, dominating them, the ~36 arrays live at once inside
:func:`~oroscope.scoring.compose`. It deliberately excludes the memory-mapped DEM.

**2. The map's memory.** A *separate* peak landing on top of the first at the very end.
The map renders at ``downsample_factor * 2``; measured at **~190 bytes per viz pixel**
plus ~130 MB of matplotlib fixed cost, by
:func:`~oroscope.site_searcher.estimate_visualisation_memory_gb`. Three runs in one
session finished their searches and then died drawing the picture, with the JSON and the
GeoTIFF already written. If a run reports its numbers but has no PNG, this is why.

**3. The address-space cap.** ``--max_memory_gb`` is applied as ``RLIMIT_AS``, which caps
**virtual** address space and therefore counts every mapping — the DEM cache, the
ping-pong buffers, the arenas numba and BLAS reserve. So it is not comparable with either
estimate above.

.. warning::

   **A cap above available memory is not a cap.** ``RLIMIT_AS`` protects the machine only
   if the process reaches it *before* the kernel runs out; above that line the OOM killer
   gets there first and picks its victim by size rather than by fault.
   :func:`~oroscope.site_searcher.preflight_memory` warns, returns
   ``cap_exceeds_available``, and with ``refuse=True`` raises instead. ``0`` disables
   capping entirely and ``tools/run_full_dem.py`` refuses it outright.

Two calibration points, both measured rather than assumed:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Run
     - Estimate
     - Peak RSS
     - Peak virtual
   * - Arequipa, 129 Mpx, ``4 / 5``
     - 5.08 GiB
     - **6.59 GiB**
     - **7.80 GiB**
   * - Huaylas crop, 11.4 Mpx, ``1 / 1``
     - 3.27 GiB
     - **4.66 GiB**
     - **5.86 GiB**

Both pairs were re-measured on the post-audit code, sampling ``VmHWM`` and ``VmPeak``
from ``/proc`` once a second across the process tree, so the table is one vintage.
Arequipa had read 5.68 GiB RSS with the virtual column blank, and a blank there is the
dangerous kind of gap: a cap sized from the RSS figure is sized from the wrong quantity.
Three attempts settled it — a 6.5 GiB cap died in the final analysis refusing 69 MiB, a
7.5 GiB cap carried GRAND but lost TAMBO's map to ``std::bad_alloc`` in the AGG backend,
and 8.0 GiB completed. **The gap between the two columns is 1.2 GiB, and that gap is the
whole reason the estimate must not be used to set the cap.**

Huaylas moved the other way, 5.31 → 4.66 GiB resident, which is the combination's
modelled memory getting 0.65 GiB cheaper. Measurements age in both directions, and
neither direction is safe to assume.

The estimator is calibrated on a *strided* run and is roughly **2× optimistic at
candidate_stride 1**. Treat it as a lower bound and size a cap from measurement when the
sampling is unusual.

**``candidate_stride`` is the memory lever, not ``downsample_factor``.** Candidates are
taken on the native grid, so downsampling scales the labelling arrays as its inverse
square and barely touches the dominant term. This is the single most useful fact for
planning a large run and it is the opposite of what most people reach for.


Numerics and geometry
---------------------

**Two Earth radii, deliberately.** The particle trajectory is not refracted, so it uses
the true radius; the radio path is, so Fresnel clearance uses the inflated *radio* radius
(``k`` ≈ 4/3). Mixing them is a real error of a few hundred metres of apparent drop over
80 km.

**Geographic pixels are not square in metres.** A degree of longitude shrinks with the
cosine of the latitude, so a 1 arc-second pixel spans ~30.7 m north–south and ~29.5 m
east–west at 17° S. The longitude scale is evaluated at the DEM's *centre* latitude so
the residual error is spread evenly rather than accumulating toward one edge. Every
metric quantity goes through :class:`~oroscope.site_searcher.MapGrid`.

**The elevation binning is free; the azimuths are not.** The walk tracks the running
maximum of the apparent terrain angle, which only increases, so each new value claims a
contiguous band of bins and one pass fills them all. Doubling the bins costs almost
nothing; doubling the azimuths doubles the run.

**Detector placement is anchored, not fitted.** ``count_grid_capacity`` lays the lattice
from each region's bounding-box corner and counts positions on usable ground. Capacity is
therefore an estimate for an *arbitrarily placed* array, not the best achievable packing.
Positions are laid out in metres and only then looked up in the pixel grid — the earlier
version converted the spacing to an integer pixel stride and truncated three separate
times, overcounting by 7.4% at GRAND's 1 km and **58% at TAMBO's 100 m**.


Things that were silently wrong once
-------------------------------------

Each of these produced plausible output while being incorrect, which is what made them
expensive.

- **A preset name iterated character by character.** ``rfi_zones`` reached the pipeline
  unresolved, and the pipeline iterates whatever it is given — so a search that believed
  it was excluding five zones excluded none, and printed ``RFI Zones: 8 active``, one per
  letter. Fixed by making :func:`~oroscope.site_searcher.resolve_rfi_zones` the single
  translation, called by all three callers.
- **A component that appeared when switched off.** Whether the geomagnetic weighting had
  been applied was judged by comparing weighted and plain solid angles — but a candidate
  that accepted no direction has a ratio of zero by construction, and those zeros stood in
  as evidence. Judged on viable candidates only, now.
- **A results prefix that named one experiment.** Outputs were
  ``grand_search_results_*`` even for TAMBO runs. Renamed to ``oroscope_results_*``;
  readers still accept the old prefix, so older runs load.
- **Funnel rows read by position.** A run with :term:`RFI zones` carries an extra stage,
  so index 4 means different things in two regions. Reading positionally made GRAND's
  acceptance at Arequipa look like 20% when it is 60.1%. Read by **name**.


Conventions worth keeping
--------------------------

- **A library must not choose the matplotlib backend.** CI asserts that importing
  ``oroscope`` leaves it untouched. Tools may set ``Agg``; the package may not.
- **Examples in docstrings are executed.** ``tests/test_doctests.py`` runs every
  ``Examples`` block, so the values in them must be *computed*, not predicted.
- **Notebooks are generated**, from ``tools/make_notebooks.py``. Edit the generator,
  never the ``.ipynb``. They are executed in CI except where the cost is prohibitive, and
  what replaces that execution is a static check that every API name they call exists.
- **Negative results are recorded** in ``docs/ROADMAP.md`` so they are not retried.
- **Lint from the repository root**: ``ruff check .`` also lints the notebooks.


Where the numbers live
----------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Path
     - What it holds
   * - ``output/<run>/``
     - Everything a run writes: GeoTIFF, world file, KML, PNG, results JSON, provenance,
       explanation. Gitignored.
   * - ``results/<region>_full/``
     - The small, readable artefacts a notebook reads: results JSON, provenance,
       explanation, combined report. Committed.
   * - ``results/region_comparison.md``
     - Every region against every other, regenerated by ``tools/compare_regions.py``.
   * - ``docs/ROADMAP.md`` §6
     - The measurements behind everything on this page.
