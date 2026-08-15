Assumptions and limitations
===========================

What this code takes for granted, what it does not model, and which numbers are
choices rather than measurements. :doc:`physics` derives the criteria; this page says
where they stop being trustworthy.

It is deliberately blunt. A site search produces authoritative-looking areas and
detector counts, and the only defence against those being over-read is writing down
what is behind them.

.. contents::
   :local:
   :depth: 2


Numbers that are assumptions, not measurements
----------------------------------------------

These are the ones to check before quoting a result.

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - Parameter
     - Default
     - What it rests on
   * - ``beta``, the tau energy-loss constant
     - :math:`0.6\times10^{-6}\,(E/{\rm EeV})^{0.20}` cm²/g
     - Estimated from mass scaling, in the range (0.4–1.0)×10⁻⁶. Moves the
       production-and-escape optimum in proportion, though not the siting conclusion.
       Worth pinning to whatever the collaboration uses.
   * - Geomagnetic declination
     - −6.9° (Arequipa IGRF)
     - **Does not follow the site.** Inclination is derived from the DEM's own
       coordinates through a centred dipole, but the dipole is unreliable for
       declination (−0.2° against a measured −6.9°) and is deliberately not used for
       it. Supply the IGRF value per site.
   * - ``grammage_band_fraction``
     - 0.1
     - How far down the Gaisser–Hillas profile still counts as a usable shower. A
       choice about detector capability, not a property of the shower.
   * - ``min_target_slope_deg``
     - unset
     - When set for a canyon search, the floor separating a wall from a hillside.
   * - ``decay_spectral_index``
     - 2.0
     - The flux slope the decay term is folded against. May be pinned to a value or
       given as a ``(low, high)`` pair to marginalise over. Folding replaced a single
       representative energy, which was far worse: see below.
   * - ``min_score``
     - 0.0
     - A cut on a product of six components. See below.


The two that dominate a canyon result
--------------------------------------

Measured by ``oroscope-sensitivity`` against a Colca TAMBO baseline, varying one
parameter at a time:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - Parameter
     - Low
     - Baseline
     - High
   * - ``decay_spectral_index`` (folded)
     - 1.5 → 7205
     - 2.0 → **9717**
     - 2.7 → 10 495
   * - ``min_score``
     - 0.0 → 45 928
     - 0.35 → **2056**
     - 0.5 → **0**
   * - ``min_target_slope_deg``
     - 0° → 7442
     - 25° → **2056**
     - 35° → **0**

**The decay term is now folded over the spectrum, and that fixed it.** Evaluated at a
single representative energy it *was* the answer rather than an approximation to one:
across TAMBO's own 3 PeV – 1 EeV reach the reported capacity ran from 10 878 to zero,
because the decay length runs from 147 m to 49 km against a ~3 km crossing. Folded
against the flux, the same result varies by 1.46× across a plausible range of spectral
index. The remaining exposure is the index itself, and it can be marginalised rather
than chosen.

What has *not* been done is the full job: the number of detected events is
:math:`\int \Phi(E)\,A(E)\,P_{\rm decay}(E)\,{\rm d}E`, and the acceptance
:math:`A(E)` is exactly what no available table supplies. The weight here is the flux
alone.

**A product score has no safe threshold.** Six components each in :math:`[0,1]`
multiply to a distribution piled up near zero, so ``min_score`` anywhere in the middle
sits on a cliff — 22× the baseline at 0.0, zero at 0.5. Prefer ranking sites and taking
the best :math:`N`; a weighted geometric mean would also spread the distribution.


Reported area is not physics-accepted area
-------------------------------------------

Three things stand between the accepted pixels and the number in the results file, and
only the first is physics.

**Morphological closing inflates.** Measured with a stride-1 control run at Colca,
closing with a 1 km element more than doubles the area — **2.29×**. So a reported
4580 km² corresponds to about 2120 km² the physics actually accepted. Closing is not
wrong: a site has to be a deployable region rather than a scatter of pixels. But the
reported figure is an upper bound, and the two should not be conflated.

**Candidate striding is unbiased, which was worth checking.** ``candidate_stride: 5``
samples one pixel in five. The same control run gives acceptance of 60.1% at both
strides, and the stride-corrected area matches the stride-1 truth to 0.05%, so this
costs nothing and saves five sixths of the scan.

**Area and capacity are measured on different grids.** Per-site ``area_km2`` comes from
the downsampled map, capacity from the full-resolution mask. At ``downsample_factor``
greater than 1 a feature only a few pixels wide loses area it keeps detectors on — for
a canyon strip that is a ~30% discrepancy. Use ``downsample_factor: 1`` for thin
features, or read the two as measuring different things.


Physics that is not modelled
----------------------------

* **Neutral-current regeneration.** Only charged-current attenuation is applied to the
  Earth chord, so the suppression of downgoing directions is somewhat overstated.
* **No shower simulation.** Everything about shower development comes from a
  Gaisser–Hillas profile and an isothermal atmosphere. There is no Monte Carlo, no
  detector response, and no trigger model.
* **No flux, no cross-section spectrum, no exposure.** The search reports geometry and
  acceptance-shaped quantities. Turning those into an event rate needs a flux and a
  response table, and :mod:`aperture` only estimates the folding.
* **Column depth is bounded by the walk.** The profile is walked out to
  ``max_dist_km`` and the depth histogram accumulates over that path, so a short-range
  search reports a depth set by where the walk stopped rather than by the target's
  thickness. At TAMBO's 5 km this makes the depth band score ~1 for everything; it
  degrades gracefully, but it is not measuring the far wall.
* **Isothermal atmosphere**, single scale height. Adequate for grammage over a few
  kilometres; not a substitute for a real profile at large zenith angles.
* **Standard rock everywhere.** One density for the whole crust. Real geology varies,
  and column depth scales with it directly.


Terrain and data limitations
----------------------------

* **The DEM is the world.** A 30 m SRTM or AW3D30 tile cannot resolve a cliff face,
  and the search inherits every artefact the DEM has — voids, stripes, vegetation
  canopy read as ground. Nodata is carried as NaN and excluded rather than
  interpolated.
* **Sub-pixel detector spacing is permitted and is a continuum limit.** At spacings
  finer than the DEM's pixels, capacity is area divided by area-per-detector; the
  terrain mask cannot say whether those positions are individually usable.
* **The layout is anchored, not fitted.** Detectors are placed from each site's
  bounding-box corner rather than optimised over placement, so capacity is an estimate
  for an arbitrarily-placed array, not the best achievable packing.
* **No logistics beyond distance-to-road.** No slope stability, land access, power,
  or cost.
* **Geographic coordinates only.** Pixel sizes come from a local flat-Earth conversion
  at the DEM's centre latitude; the tool does not reproject, and a projected DEM is not
  handled.


Where the results have been checked
-----------------------------------

Worth stating, because it bounds how much the rest should be trusted.

* Terrain fixtures are **synthetic with closed-form answers** — slope, aspect, target
  distance and canyon geometry are known analytically, so tests assert against
  arithmetic rather than against a previous run. The fixtures themselves are verified
  before the code that uses them.
* The far-wall slope measurement **recovers a canyon fixture's own wall slope
  exactly** at 15°, 25°, 35° and 45°, and at Colca it recovers 34.7–44.3° against a
  published ~40°.
* Capacity **matches analytic lattice density to 2%** from 1000 m down to 60 m
  spacing, on both square and triangular grids.
* Golden-file regression pins whole-pipeline output on synthetic terrain and a real
  DEM crop.

What has *not* happened: no comparison against a collaboration's own simulated
acceptance, and no end-to-end validation against a site chosen by other means. Until
one of those exists, treat absolute numbers as internally consistent rather than
externally verified.


A prediction worth falsifying
-----------------------------

The Earth-absorption treatment implies the effective arrival window **narrows with
energy**: its lower edge climbs from −4.4° at 100 PeV to −0.9° at 10 EeV. If a
collaboration's simulated window does not narrow that way, one of the two treatments
has the absorption wrong. That is the most useful single check anyone with an
independent simulation could run against this code.
