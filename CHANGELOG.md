# Changelog

Notable changes, newest first. Measured deltas are quoted where a change moved a number.

## 0.5.0 — 2026-08-17

First published release. The version jumps from the unreleased 0.1.0 because what is
being released is the audited pipeline and a set of results measured on it, not the
prototype that number described.

### Added
- **One import: `import oroscope`.** The modules moved from flat top-level names into a
  real package, whose `__init__` re-exports the whole public surface — 151 names — while
  the submodules stay importable when a narrower namespace reads better. The notebooks'
  `sys.path` insert is gone; `pip install -e .` is the only setup step.
- **A dedicated CLI page** in the documentation, with the complete option reference —
  all 87, with types and defaults, generated from the parser so it cannot drift. The
  README and the quickstart now lead with code, since that is how most people use this.
- **Every run explains itself.** A plain-language summary — what was found, which
  funnel stage set the size of the answer and the parameter behind it, which named
  score component held each site back, and which numbers are assumptions with their
  measured sensitivity. Printed last and saved as `explanation.txt`; on by default,
  `--no_explain` suppresses it. `explain.explain_results(results)` is pure, so an old
  results file can be re-explained with no DEM and nothing re-run.
- **Named score components are now stored per site**, so a site can be attributed
  rather than only ranked — both ways round. The summary reports **why each site is
  good**: the criteria it satisfies, each with the measurement that earned it ("1.08 sr
  of accepted sky, targets at 3,137 m, striking 39° terrain"), and the one that held it
  back. On the Colca configs `solid_angle` is the weakest component at 16 of 16 TAMBO
  sites — and at every selected site in every region searched: 62 of 62 in Ancash, 84 of
  84 in Lima, 85 of 85 in Arequipa.
- **Site records carry coordinates**: centre latitude/longitude and a bounding box, so
  a reader can find the ground without opening the raster in a GIS.
- **The combination explains itself too.** `oroscope-combine` prints and saves an
  account of the overlay: what each experiment brings, how much ground they can share,
  and **which screening band decides that** — a pixel has one slope and both
  experiments must accept it, so co-location is settled there. What each asks of the
  *view* may differ freely, and is reported as no obstacle.
- The summary also reports the energy the geometry favours, and **what to try next**,
  as concrete commands chosen from what the run actually did.
- **CLI/library parity.** `max_memory_gb` is a pipeline parameter and the memory
  estimate, warning and cap are `preflight_memory()`; `load_config`,
  `generate_config` and `default_config` are ordinary functions rather than a literal
  inside `main()`; and the pipeline **returns its results dictionary** instead of
  `None`, so callers no longer re-read the JSON it just wrote.
- **Two experiments from one engine.** TAMBO is now a configuration rather than a code
  path, alongside GRAND, and `oroscope-combine` overlays two runs to report joint,
  union and co-location.
- **Per-role slope criteria.** The scan records how steep the terrain it strikes is,
  along the arrival azimuth, so a canyon search can require a *far wall* as well as a
  deployable near one. Recovers a fixture's wall slope exactly at 15°, 25°, 35°, 45°.
- **Tau decay folded over a spectrum**, with the index pinned or marginalised over a
  range. Replaces a single representative energy, which *was* the answer rather than an
  approximation to it: 10 878 → 0 detector positions across TAMBO's own energy reach,
  against 1.46× across a plausible range of spectral index.
- **Rank-based candidate selection** (`--score_percentile`) and best-N site selection
  (`--stop_at_target`), because the default score is a product whose distribution piles
  up near zero and so has no safe absolute threshold.
- **The DEM's origin is read from its own tiepoint**, and a supplied origin that
  disagrees is reported rather than silently honoured.
- `oroscope-sensitivity`, which varies one parameter at a time and tabulates what moves.
- Memory safeguards: a pre-flight estimate, an address-space cap, and subprocess
  isolation per sweep point.
- Documentation: physics and assumptions pages, eight tutorial notebooks, and
  reproducible figures in `src/oroscope/figures.py`.
- **Two more notebooks.** **7** drives the pipeline from Python and reads what it
  says, about a run that finds ground and one that finds none — the empty result being
  the case a bare results file serves worst. **8** is the full Arequipa DEM: it *reads*
  stored results rather than producing them, since three searches at half an hour each,
  against CI that executes every notebook on every push, is not a tutorial.
  `tools/run_arequipa_full.py` produces that store locally into `results/arequipa_full/`;
  regenerate it when a configuration changes. Configurations for both experiments over
  the full DEM are in `config/`.

### Fixed
- **The funnel blamed the arrival window for cuts `min_score` made.** It recorded the
  post-cut count under the name `directions accepted`, so geometry and the score cut
  were one row. Separated on the Ancash TAMBO run, the scan keeps **82%** of what
  striding hands it and `min_score` takes **8.4×**. Every ratio derived from that row
  was measuring two things at once — including the "acceptance is identical at 14.0%"
  that demonstrated striding was unbiased. It still is, at both stages independently:
  geometry 81.640% against 81.612%, score cut 8.59% against 8.60%.
- **A store could be half-refreshed and look current.** `results/huaylas_full/` held
  three committed control files no run produces, four months stale at the old spacing,
  and `results/arequipa_full/manifest.json` listed three files while the store held
  nine — because `--only` wrote a manifest of what *that* invocation stored. A reader
  following the manifest's own "read this first" would have concluded GRAND had never
  run. The manifest now separates what the store holds from what the run wrote, and
  dates everything else under `also_present`.
- **An animation was broken by a guard added in the same release, and CI could not see
  it.** Teaching `scoring.compose` to refuse an all-zero weighting broke
  `product_collapse`, whose first frame *is* an all-zero weighting — the empty product
  before any criterion applies. Notebook 07 is excluded from CI because the runner has
  no ffmpeg. Fixed in the animation, where the empty product is meaningful, with a test
  that builds every animation's first frame and needs no ffmpeg.
- **A figure was published twice.** `howitworks.rst` was the only page without a
  `%matplotlib inline` setup block, so its first block rendered once as the cell's
  result and once more as the backend flushed — the seven-stage funnel appeared twice
  in a row on the live site.
- **A parameter's default depended on which door you came in by.** The pipeline's
  signature, the argparse parser and `default_config()` disagreed on **ten**
  parameters — `search_mode` was `single` from Python and `distributed` from a shell,
  `min_dist_km` was 30 km against 10. The signature was the odd one out in every case
  and now matches the documented values; a test compares all three pairwise. Two were
  in the template, and one of those was a hazard: `origin_lat`/`origin_lon` defaulted
  to `0.0`, a *valid* coordinate in the Gulf of Guinea, so a forgotten placeholder
  georeferenced a run to the wrong continent instead of failing. They are `null` now,
  meaning "read the DEM's own tiepoint". The bundled configurations set every
  parameter explicitly and are unaffected: Colca is still 15 sites and 9717 detectors.
- **`import oroscope` forced the matplotlib backend.** `combine_experiments` called
  `matplotlib.use("Agg")` at module level — harmless for a standalone module, not for a
  package front door: it reached into every caller's session and overrode the inline
  backend, so notebooks captured no figures at all. Chosen in `main()` now, where the
  command line actually needs it. A library must not decide how its user's figures are
  rendered.
- **The map title said "GRAND site search" on every run**, including TAMBO's, as did the
  console banner and the KML placemark names.
- **A `geomagnetic` score component appeared in runs that had switched the weighting
  off.** Whether it was applied is judged by comparing the weighted solid angle with
  the plain one, and a candidate that accepted no directions has a ratio of zero by
  construction — so those zeros stood in as evidence of weighting. Harmless under a
  product composition (it was multiplying by one, and TAMBO's numbers are unchanged),
  wrong under `mean`, and it listed a disabled criterion among the reasons a site
  was good.
- **`oroscope-combine` crashed on any single-mode run.** Two faults in the same path,
  both found by new tests: `capacity_of` did not catch the `ValueError` from
  `int('N/A')` — which is what `search_mode: single` writes — and the console summary
  applied thousands-grouping to a string when a run had no capacity to report. Either
  one took the whole combination down.
- **`main()` leaked its log file and never restored `sys.stdout`.** It tees both
  streams into the run's log; the swap and the open handle outlived the call, so a
  process running it twice stacked a `TeeLogger` on the previous one and leaked a
  handle each time.
- **Documentation had drifted from the code.** The README documented 34 of 83 CLI
  options and one, `--fresnel_buffer`, that no longer existed; it described the old
  precedence rule, config-over-command-line, which was fixed long ago; the startup
  banner still described a single ray cast to a target mountain; and three modules,
  including `site_searcher` itself, had no docstring at all. `tests/test_docs.py` now
  pins each of those.
- **`--output_directory_base_with_given_json` lost to the config file.** It resolved
  before the merge loop and so kept the precedence every other flag had shed, meaning
  it was silently ignored whenever a config set it.
- **The library did not create its own output directory**, so a caller who passed a
  path that did not exist got `FileNotFoundError` from inside numpy's `open_memmap`,
  naming a scratch buffer rather than the directory.
- **The results file listed more sites than were selected, with nothing saying which.**
  With `--stop_at_target`, `sites` holds everything that cleared the thresholds while
  `total_sites`, `total_capacity` and the exported raster cover only the selection — so
  anything totalling the list over-reported. Measured on a synthetic run: 2 sites and
  243.9 km² against a mask holding 1 site and 215.7 km². Each record now carries
  `selected`, and the run summary counts and sums the selection. Its area now equals
  the exported mask's exactly, in plain, truncated, downsampled and single modes.
- **`oroscope-combine` read a stale mask.** It took the alphabetically first `.tif` in
  a run directory, and the pre-rename `grand_search_results_*` prefix sorts before
  `oroscope_results_*` — so a re-run directory was overlaid using its superseded mask,
  silently. Corrected Colca figures: TAMBO 44.5 → **83.6 km²**, joint 26.4 → **50.1
  km²**, union 4598.3 → **4613.7 km²**. GRAND's own numbers are unchanged. *Those
  absolute values are the size of this fix, not the figures this release ships: the
  150 m re-run below moved Colca's TAMBO area to 203.0 km² and its joint to 123.3.*
- **Capacity was over-counted twice.** Integer pixel stamping inflated it by 7.4% at
  1 km spacing and 58% at 100 m; and capacity was counted over each site's *bounding
  box*, which also contains other sites — 38% on a canyon network, 2.07× on a synthetic
  case.
- **A 256th site crashed the run** (`uint8` overflow) after the physics had been paid for.
- **A matplotlib figure leaked per run**, which took a ten-point sweep to 6.9 GB and
  into the OOM killer.
- **The command line lost to the config file**, silently, making every flag a no-op
  whenever a generated config was used.
- Per-site work was O(sites × pixels); now flat in the site count.
- The declared Python 3.9 floor was false — six modules used `X | None` without
  `from __future__ import annotations`.
- The digitised published curves were never committed, so a clone failed eight tests.
- Documentation figures rendered as text, not images.

### Changed
- **Every published number comes from the audited pipeline.** All six regions were
  re-searched and every derived figure regenerated from them, so `results/`, the
  notebooks and the documentation agree. GRAND barely moved; TAMBO moved everywhere.
  Arequipa's TAMBO area went 111.9 → **1,036.9 km²** and its capacity 9,024 → **49,271**;
  GRAND's area 88,527.5 → 88,208.2 km².

  | region | sampling | GRAND km² | TAMBO sites / km² | joint km² | share of TAMBO |
  | --- | --- | --- | --- | --- | --- |
  | colca | 1 / 5 | 4,569.4 | 16 / 203.0 | 123.3 | 60.7% |
  | huaylas | 1 / 1 | 8,249.5 | 32 / 291.3 | 228.7 | **78.5%** |
  | cajatambo | 1 / 1 | 5,573.8 | 44 / 774.5 | 591.7 | **76.4%** |
  | ancash | 4 / 5 | 42,791.9 | 62 / 740.0 | 411.1 | 55.6% |
  | lima | 4 / 5 | 51,209.0 | 84 / 915.4 | 509.8 | 55.7% |
  | arequipa | 4 / 5 | 88,208.2 | 85 / 1,036.9 | 619.1 | 59.7% |

  **Quote the crops, never the departments, for TAMBO** — department figures are
  striding artefacts. The joint share is two rows and never one: ~56–60% strided,
  ~76–79% unbiased.
- **TAMBO runs at the published 150 m detector spacing**, not 100 m, and
  `gap_close_km` defaults to that spacing. Both striding penalties collapse as a result:
  Colca **4.75× → 1.51×**, the Callejón de Huaylas **291× → 23.0×**. A 3-pixel closing
  element recovers 0.04× of the accepted set and a 5-pixel one 0.68×, so 100 m sat on
  the wrong side of that cliff and 150 m does not.
- **The solid-angle score is scale-free.** `solid_angle_half_fraction` is a fraction of
  the sky the fan and the arrival window actually make available; the steradian form
  silently encoded the fan width, so a value tuned for one fan was wrong for any other.
  `solid_angle_half_sr` still works and still means steradians.
- **Closing inflation at Colca re-measured, 2.29× → 2.35×**, from the same stride-1
  control that produced the original. `explain.AREA_INFLATION_AT_COLCA` is the only
  place the figure lives, and a test asserts against the constant rather than a literal.
- **A conclusion was withdrawn.** The joint no longer "barely moves with scale": at
  100 m the Colca crop gave 50.1 km² and the whole Arequipa DEM 50.2; at 150 m it is
  123.3 against 619.1. The old invariance was an artefact of an element too small to
  reconnect what striding cut apart.
- **Documentation figures render as SVG and carry a zoom control** — a magnifier on
  each figure and an overlay, as static assets rather than an extension.
- **Figure labels capitalise their first word** — axis labels, titles, legend entries
  and annotations, everywhere. Stated at the top of `figures.py` so it holds.
- Renamed to **oroscope**; outputs are `oroscope_results_*` (the old prefix still reads).
- `src/setup.py` → `src/oroscope/fetch_dem.py`, so `pip install` no longer runs the downloader.
- Packaged: `pip install -e .`, five console scripts, CI on Python 3.9–3.13.
- Every criterion is a configuration knob; nothing that shapes a result is hard-coded.

### Known limitations
See [assumptions and limitations](https://mbustama.github.io/oroscope/assumptions.html).
The short version: reported area is ~2.35× the physics-accepted area because of
morphological closing, the detector acceptance *A(E)* is not modelled, Askaryan
emission is not modelled, and nothing has been checked against an external simulation.
