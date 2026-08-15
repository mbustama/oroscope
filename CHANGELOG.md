# Changelog

Notable changes, newest first. Measured deltas are quoted where a change moved a number.

## Unreleased

### Added
- **Every run explains itself.** A plain-language summary — what was found, which
  funnel stage set the size of the answer and the parameter behind it, which named
  score component held each site back, and which numbers are assumptions with their
  measured sensitivity. Printed last and saved as `explanation.txt`; on by default,
  `--no_explain` suppresses it. `explain.explain_results(results)` is pure, so an old
  results file can be re-explained with no DEM and nothing re-run.
- **Named score components are now stored per site**, so a weak site can be attributed
  rather than only reported. On the Colca configs the attribution is unambiguous:
  `solid_angle` is the weakest component at 15 of 15 TAMBO sites.
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
- Documentation: physics and assumptions pages, six tutorial notebooks, and reproducible
  figures in `src/figures.py`.

### Fixed
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
- Renamed to **oroscope**; outputs are `oroscope_results_*` (the old prefix still reads).
- `src/setup.py` → `src/fetch_dem.py`, so `pip install` no longer runs the downloader.
- Packaged: `pip install -e .`, five console scripts, CI on Python 3.9–3.13.
- Every criterion is a configuration knob; nothing that shapes a result is hard-coded.

### Known limitations
See [assumptions and limitations](https://mbustama.github.io/oroscope/assumptions.html).
The short version: reported area is ~2.3× the physics-accepted area because of
morphological closing, the detector acceptance *A(E)* is not modelled, and nothing has
been checked against an external simulation.
