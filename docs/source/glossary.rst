Glossary
========

Every term this project uses in a sense that is not the ordinary English one, with a
link to where it is explained and where it bites.

.. glossary::
   :sorted:

   acceptance
      The fraction of scanned candidates that survived to the morphology. **Two funnel
      rows can answer that and they mean different things.** ``directions accepted`` is
      the *geometry* alone — candidates for which the :term:`arrival scan` accepted at
      least one direction. Where a score cut is in force a second row follows it
      (``score >= 0.35``, or ``score in top 25%``), and *that* is the count the reported
      sites, area and capacity were built from. Quoting the geometric row beside them
      overstates acceptance badly: at TAMBO Arequipa it reads 63% where the run kept
      9.7%.

      **Read it by stage name, never by position** — a run with :term:`RFI zones`
      carries an extra funnel row, so the same index means different things in two
      regions. Unbiased under :term:`striding`, which is why striding is described as
      cost control. See :doc:`howitworks`.

   arrival scan
      The central operation: from a candidate pixel, walk outward along each bearing of
      the :term:`azimuth fan` and find where the ray first meets terrain. A direction is
      accepted when that intersection falls inside the distance window, the elevation
      window and any target-slope band. One walk fills every elevation bin at once, so
      the azimuth count sets the cost. See :doc:`howitworks` and
      :func:`~oroscope.arrival_scan.scan`.

   aspect
      The compass direction a slope faces. The :term:`azimuth fan` is centred on it, so
      a candidate looks out across the ground that falls away in front of it.

   azimuth fan
      The bearings the :term:`arrival scan` tests at one candidate: ``n_azimuths``
      directions within ``azimuth_half_width_deg`` of the :term:`aspect`. Cost is linear
      in this number. :func:`~oroscope.arrival_scan.azimuth_fan` returns the offsets;
      notebook 7 animates the sweep.

   candidate
      A pixel that survived :term:`screening` and :term:`striding`, and is therefore
      given to the :term:`arrival scan`. Candidates are taken on the **native** grid, not
      the downsampled one, which is why :term:`downsample factor` is the weaker of the
      two memory levers.

   closing
      The morphological operation that fills holes in the accepted mask, turning the
      lattice :term:`striding` leaves back into a region an array could occupy. Its
      structuring element is ``gap_close_km``, defaulting to the antenna spacing.
      **It must be larger than the stride gap**; see :doc:`howitworks` and
      :func:`~oroscope.site_searcher.warn_stride_outruns_closing`.

   co-location
      Siting two experiments so that one deployment serves both. **Not the same as the**
      :term:`joint region`: a partner array does not have to stand on ground both
      experiments accept, because what couples them is a shared line of sight to the same
      massif rather than a shared footprint. Measured across three regions, a GRAND array
      of 100 antennas fits within ~10 km of the best TAMBO site, 1,000 within 20–30 km and
      5,000 within 40–60 km. :func:`~oroscope.combine_experiments.colocation_capacity`
      and :func:`~oroscope.combine_experiments.smallest_radius_for`.

   column depth
      Rock traversed along an arrival direction, in g/cm². Accumulated over the whole
      profile walk, so it is bounded by where the walk stopped rather than by the
      target's true thickness — read it as a lower bound. See :doc:`physics`.

   composition
      How component scores are combined: ``product`` (unforgiving, the default),
      ``mean`` (compensating) or ``min`` (weakest link). Under ``product`` the weights
      act as exponents. See :func:`~oroscope.scoring.compose`.

   downsample factor
      ``downsample_factor``. Reduces the grid on which the mask is labelled and **area
      is measured**, but not the grid candidates are taken on. So it saves memory as its
      inverse square on the labelling arrays only, and it costs area for thin features —
      about 30% for a canyon strip. Capacity is still counted at full resolution, which
      is why area and capacity can disagree. See :doc:`assumptions`.

   funnel
      The per-stage record of how many pixels survived: DEM pixels, screening,
      striding, directions accepted, the score cut where one applies, closing, pruning,
      selection. The geometry and the score cut are separate rows, because they are
      separate questions and a reader needs to know which one bound the run. Printed by
      every run and stored in the results JSON. **Its two halves run in opposite
      directions** — see :doc:`howitworks`. :class:`~oroscope.site_searcher.Funnel`.

   joint region
      Ground accepted by two experiments at once, so a single site could host both. Its
      size relative to the smaller experiment's mask is near-constant across regions —
      about 73% at unbiased sampling, about 44% at ``4 / 5``, and **the two must never be
      mixed**. See :doc:`notebooks` (notebook 11).

   min score
      ``min_score``. The cut applied to the composed score. The dominant assumption in
      the project: a product piles up near zero, so any threshold in the middle sits on a
      cliff whose position moves when a component is added.
      :doc:`howitworks` shows why; ``score_percentile`` is the scale-free alternative.

   pruning
      Dropping labelled regions that are too narrow (``min_width_km``), too small, or
      hold too few detector positions (``min_sub_array_size``). This is where a
      fragmented mask is actually lost — not in the area measurement. See
      :doc:`howitworks`.

   RFI zones
      Radio-frequency-interference exclusion zones: circles or polygons around
      settlements and industry that candidates must avoid. Supplied as a preset name, a
      JSON string or an explicit list. A preset name reaching the pipeline unresolved
      used to be iterated character by character, silently excluding nothing; see
      :func:`~oroscope.site_searcher.resolve_rfi_zones`.

   screening
      The cheap per-pixel test that runs before anything expensive: slope band, and
      optionally altitude, aspect, road distance and :term:`RFI zones`. Slope is the
      criterion that decides co-location, because a pixel has one slope and both
      experiments must accept it.

   scoring
      Turning the :term:`arrival scan`'s observables into named component scores in
      [0, 1] and combining them by :term:`composition`. The components are retained
      per site so a result can be *described* rather than merely ranked. See
      :func:`~oroscope.scoring.score_candidates`.

   striding
      ``candidate_stride``. Keeping one surviving pixel in N as a :term:`candidate`.
      Cost control rather than a criterion, and unbiased in :term:`acceptance` — but it
      leaves a gap that :term:`closing` must bridge, and where it does not the reported
      area collapses. **The memory lever**, more than :term:`downsample factor`.
      See :doc:`howitworks`.

   target
      The terrain a ray strikes: the far canyon wall for TAMBO, a distant massif for
      GRAND. Distinct from the ground the array stands on, which is why there are two
      slope criteria — ``min_slope_deg``/``max_slope_deg`` for the near ground and
      ``min_target_slope_deg`` for what is struck.

   viz downsample
      The map renders at ``downsample_factor * 2``, so its raster is a quarter the area
      of the labelling grid. It is a **separate memory peak** landing on top of the
      search's, and it is the stage that most often fails; see :doc:`implementation` and
      :func:`~oroscope.site_searcher.estimate_visualisation_memory_gb`.
