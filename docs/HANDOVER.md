# Handover brief — Oroscope

Written to be fed into a fresh session. It assumes no memory of the previous one.

**Repository:** `mbustama/oroscope` (renamed from `site-search`; GitHub redirects the old
name). Local path `~/Research/GRAND/oroscope`, with a `site_search` symlink left beside
it for anything that still points at the old path — **delete it when convenient.**

**Branch:** `dev`, 39 commits ahead of `main`, all pushed. **Head at handover: `34887d9`.**
**Tests:** 370, stdlib `unittest`, ~18 s. **CI:** 8 jobs, all green.

`main` is protected by a repository ruleset: no direct push, no force-push, no deletion,
pull request required, seven status checks. **You cannot push to `main`** — work on `dev`
and open a PR.

---

> **Status update.** Both immediate tasks in §1 are **delivered** — see `docs/ROADMAP.md`
> §6.23 (`--explain`) and §6.24 (CLI/library parity). Doing so turned up a bug in
> `oroscope-combine` that makes **the §5.1 table below wrong for TAMBO**; the corrected
> numbers are in that table and in ROADMAP §6.23b. §1 is kept as the record of what was
> asked for.

## 1. The immediate tasks

### 1.1 `--explain`, on by default

The owner asked for a human-readable summary of a run: *"here is what was found and
why"*. Everything needed is already in the results JSON, but a reader has to assemble
the story themselves, and these runs are meant to be handed to other people.

**It must default to on.** The owner was explicit. Add `--no_explain` to suppress.

What it should draw on, all of which already exists:

- the **funnel** (`results["funnel"]`) — survivors after each filter. The line where the
  count collapses is the constraint responsible, and that is the single most useful
  thing to say when a search returns little or nothing.
- `results["regions"]` — labelled regions → past area threshold → past capacity
  threshold → selected.
- per-site records, already sorted by capacity, each with 34 `arrival_scan` fields
  (mean/median/p90 of solid angle, exit distance, column depth, horizon, grammage,
  Earth chord, altitude, far-wall slope, score and its named components).
- `results["parameters"]` — every resolved knob, so the summary can name the ones that
  did the work.
- `provenance.json` — git commit, DEM sha256, package versions.

Worth saying in the output, because they are the things a reader will otherwise get
wrong: that **reported area is ~2.3× the physics-accepted area** (morphological
closing, §5.2), that the score components are named so a weak site can be *attributed*,
and which parameters are assumptions rather than measurements.

Write it as a function that takes the results dict and returns a string, called by the
CLI — not as printing scattered through the pipeline. That way the library gets it too
(see §1.2) and it is testable without running a search.

### 1.2 Everything the CLI can do, the library must do too

The owner asked for full parity. Measured, the CLI flags with no
`find_grand_regions_interactive()` equivalent are:

| flag | status |
| --- | --- |
| `--max_memory_gb` | **A real gap.** Applied only in `main()`; a library user must call `apply_memory_cap()` themselves. |
| `--nearest_sampling`, `--no_geomagnetic`, `--require_sky`, `--include_near_field`, `--no_print_info` | Fine — negative-form aliases whose positive form (`bilinear_sampling`, `use_geomagnetic`, `require_terrain`, `exclude_near_field`, `print_info`) is a library parameter. |
| `--config_path`, `--config_preset`, `--generate_config`, `--output_directory_base_with_given_json` | CLI-level concerns, but a library user would reasonably want `load_config(path)` and `generate_config(path, preset)` as functions. Worth adding. |

Also only in `main()` and worth moving into the library: the **pre-flight memory
estimate and warning**, and the **origin resolution print**. Both are useful to anyone
driving the pipeline in a loop.

The pipeline function returns `None`. It writes files and the caller reads the JSON
back. **Returning the results dict** would make the library genuinely usable and costs
nothing — `tests/_support.py` and `sensitivity.py` both currently re-read the file it
just wrote.

---

## 2. What does NOT need to be run again

Read this before measuring anything. Several of these cost real time or crashed the
machine.

**Do not re-run:**

- **The benchmark baseline.** `bench/baseline.json` was refreshed at head `34887d9` on a
  quiet machine (1-minute load 0.89) and matches the current code. Only re-run after a
  change that should move a stage timing, with `python bench/benchmark.py --update`.
- **The Colca searches.** Both configs were run at this head; the outputs in
  `output/grand_colca_config/`, `output/tambo_colca_config/` and
  `output/combined_colca/` are current. Numbers in §5.
- **The sensitivity sweeps.** Both the single-energy and the spectrum-folded sweeps are
  recorded in `docs/ROADMAP.md` §6.20–6.21 with their tables. Re-running costs ~10
  minutes and will reproduce them.
- **The stride-1 control run.** `config/grand_colca_stride1.json` exists and its result
  is recorded (§5.2): striding is unbiased, closing inflates 2.29×.
- **The notebooks.** All six are committed *with their outputs*. Regenerate only if you
  change `tools/make_notebooks.py`, and then re-execute — but note the trap in §3.
- **DEM downloads.** `input/dem/` holds `arequipa_SRTMGL1.tif`, `lima_AW3D30.tif` and the
  derived `colca.tif` crop. **The Arequipa DEM already covers Colca Canyon** — verified,
  1673 m of incision against the published ~1.5 km. No new download is needed for either
  experiment.

**Do not re-derive:** everything in §6. Those are measured facts, several of which
contradicted a confident prior.

**Do not retry:** everything in §7. Those were tried and rejected with numbers.

---

## 3. Environment, machine, and three traps

- **Use the conda env `sssearch`**: `conda activate sssearch`. It has numpy, scipy,
  numba, tifffile, imagecodecs, matplotlib, tqdm, sphinx and the docs stack, coverage,
  ruff, jupyter. It does **not** have pytest, which is why the suite is stdlib
  `unittest`.
- The package is installed editable (`pip install -e .`), so `import physics` works from
  anywhere and the five console scripts (`oroscope`, `oroscope-combine`,
  `oroscope-crop`, `oroscope-sensitivity`, `oroscope-fetch-dem`) are on `PATH`.
- **Cap parallelism at 8 cores.** The machine has 12 but is shared.

**Trap 1 — memory.** A ten-point sensitivity sweep once reached 6.9 GB and was killed by
the OOM killer, taking other work with it. The cause was ours (a leaked matplotlib
figure per run) and is fixed, but the safeguards matter: every run now prints an
estimate against available memory, caps its own address space at 80% of available, and
`sensitivity.py` runs each point in a subprocess. The machine has 15 GB with typically
~6 free. **For the full DEM use `downsample_factor: 4`** — the estimator says 2.3 GiB
against 4.5 GiB at `downsample_factor: 1`.

**Trap 2 — timings are unreliable here.** It is a hybrid CPU: 2 P-cores (4.6 GHz, CPUs
0–3) and 8 E-cores (3.4 GHz), running at about a third of rated clock under load. The
same unchanged code measured 43.6 s and 39.8 s in consecutive runs. **A/B alternating
inside one process, single-threaded, on a subsample.** See ROADMAP §6.12.

**Trap 3 — forcing a matplotlib backend breaks image capture, silently.** This bit twice.
Setting `MPLBACKEND=Agg` as an environment variable propagates into child kernels and
overrides the inline backend, so notebooks stored no figures and documentation pages
rendered figures as the text `<Figure size ...>`. Both built clean and reported no
error. `conf.py` now uses `matplotlib.use('Agg')` — module-level, not inherited — and
the docs carry a `%matplotlib inline` setup cell. **If you touch either, check that
images are actually produced**, don't trust a green build.

---

## 4. Repo map

| path | what |
| --- | --- |
| `src/site_searcher.py` | 3517 lines. Pipeline, CLI, screening, morphology, capacity, outputs, memory guards. |
| `src/arrival_scan.py` | The scan kernel: profile walking, column depth, Fresnel, RFI line-of-sight. Numba. |
| `src/physics.py` | Closed-form physics, no terrain: atmosphere, shower profile, Earth chord, tau range and decay, geomagnetic, Cherenkov. |
| `src/scoring.py` | Score shapes (band, saturating, ramp) and composition. |
| `src/aperture.py` | Aperture estimate, tabulated response, `infer_response()`. |
| `src/combine_experiments.py` | Overlays two or more runs: joint, union, co-location. |
| `src/crop_dem.py` | Cuts a lat/lon window out of a DEM. |
| `src/sensitivity.py` | One-at-a-time parameter sweeps, each point in a subprocess. |
| `src/figures.py` | The publication figures, as functions returning `Figure`. |
| `src/fetch_dem.py` | Downloads DEMs. Was `setup.py`, whose name hijacked `pip install`. |
| `tests/` | 370 tests. `synthetic.py` builds terrain with closed-form answers. |
| `tools/make_notebooks.py` | Generates the six tutorials. Edit here, not the `.ipynb`. |
| `docs/source/` | Sphinx. `physics.rst` derives the criteria; `assumptions.rst` is the blunt list of what the numbers rest on. |
| `docs/ROADMAP.md` | ~1500 lines. The durable record: every phase, every measurement, every negative result. **Read §6.11, §6.12, §6.20–6.22.** |
| `bench/benchmark.py` | Per-stage timings and peak RSS, gated at 30% regression. |
| `config/` | `grand_colca_config.json`, `tambo_colca_config.json` (same crop, combinable), plus arequipa/lima and the stride-1 diagnostic. |

## 5. What the tool does now, and its current numbers

It screens a DEM by slope/aspect/altitude/exclusion zones, scans arrival directions from
each survivor, scores against per-experiment criteria, cleans up morphologically, labels
sites and places detectors on a lattice, then writes GeoTIFF/world file/KML/PNG/JSON plus
a funnel and a provenance record.

**GRAND and TAMBO are configurations, not code paths.** That was the phase 2 claim and it
held: adding an experiment means writing a JSON file.

### 5.1 Colca, at this head

**Corrected.** The figures first written here came from `oroscope-combine` reading a
stale mask — it took the alphabetically first `.tif`, and the pre-rename
`grand_search_results_*` prefix sorts before `oroscope_results_*`. Fixed; see ROADMAP
§6.23b. GRAND's own numbers were unaffected.

| | area | sites | capacity | of its own area in the joint |
| --- | --- | --- | --- | --- |
| GRAND | 4580.2 km² | 1 | 5317 | 1.1% *(was 0.6%)* |
| TAMBO | **83.6 km²** *(was 44.5)* | 15 | 9717 | 59.9% *(was 59.3%)* |
| **joint** | **50.1 km²** *(was 26.4)* | | | |
| **union** | **4613.7 km²** *(was 4598.3)* | | | |

**Co-location is decided by slope, not arrival geometry.** GRAND's 3–25° deployable band
against Colca's ~40° walls leaves only a 20–25° sliver.

### 5.2 Numbers to quote carefully

- **Reported area is ~2.29× the physics-accepted area.** Morphological closing, measured
  with a stride-1 control. GRAND's 4580 km² corresponds to ~2120 km² actually accepted.
- **Candidate striding is unbiased** — acceptance identical at strides 1 and 5, and the
  stride-corrected area matches the stride-1 truth to 0.05%.
- **TAMBO's capacity now varies by 1.46× across a plausible spectral index**, having
  varied without bound (10878 → 0) when the decay was evaluated at a single energy.
  `min_score` is now the dominant assumption at 2.38× to 0.20×.

---

## 6. Key measured facts — do not re-derive

| finding | value |
| --- | --- |
| Slope depends on measurement baseline | median 17.8° at ~61 m, 10.8° at 1 km |
| Morphological closing inflates area | **2.29×**, measured at stride 1 |
| Candidate striding | unbiased; stride-5 matches stride-1 area to 0.05% |
| Capacity over-count, integer stamping | +7.4% at 1 km, **+58% at 100 m** — fixed |
| Capacity over-count, bounding box vs region | **+38%** on a canyon network, 2.07× synthetic — fixed |
| The ±3° window sits below the horizon almost everywhere | median horizon 7.3° |
| Geomagnetic asymmetry | east-facing targets worth 3.7× north-facing |
| Earth absorption narrows the window | −4.4° at 100 PeV, −2.0° at 1 EeV, −0.9° at 10 EeV |
| Colca crossing supplies | ~170 g/cm² across 2 km, ~390 across the full 4.5 km |
| Shower maximum needs | 561 g/cm² at 3 PeV, 700 at 1 EeV |
| Tau decay length | 147 m at 3 PeV, 49 km at 1 EeV |
| Bilinear vs nearest sampling | acceptance +13.4%, 1.44× cost |
| Azimuth locality penalty | 1.14× (the sweep premise, measured and rejected) |
| Scan thread scaling | 1.85× at 2 threads, **3.70× at 8** — hardware, not scheduling |
| Far-wall slope recovered at Colca | 34.7–44.3°, median 38.6°, against a published ~40° |

## 7. Tried and REJECTED — do not repeat

| attempt | result |
| --- | --- |
| Whole-raster azimuthal sweep | Premise measured wrong; locality penalty only 1.14×. |
| Hoisting the per-sample division | **Slower**, 24.5 → 26.1 s. |
| Solving the ray's exit distance up front | ~5% but broke 405 of 40,000 candidates. |
| Tabulating the scan's per-bin transcendentals | Bit-identical, **0.975×**. Reverted. |
| Short-circuiting the histogram's bin search | Bit-identical, **0.992×**. Reverted. |
| Tuning `numba.set_parallel_chunksize` | ~3%. Block-dealing already handles it. |

**The lesson:** per (candidate, azimuth) the bin loop runs 12 times and the profile walk
~2700, a ratio of 225:1. **The scan is not compute-bound.** Four separate arithmetic
optimisations returned nothing. Anything further must reduce *samples* or *memory
traffic*, not flops.

## 8. Working conventions — please keep

1. **Measure before optimising, and after.** Several confident hypotheses here were wrong.
2. **Run examples, don't read them.** Eight docstring examples were plausible, close, and
   wrong; every one was caught by executing it. `tests/test_doctests.py` runs them all.
3. **Negative results go in `docs/ROADMAP.md`** so they are not retried.
4. **Fixtures are verified before the code that uses them** (`tests/test_fixtures.py`).
5. Regenerate goldens deliberately: `cd tests && UPDATE_GOLDEN=1 python -m unittest test_regression`.
6. Commit messages explain *why* and state measured deltas.
7. The roadmap is updated in the same commit as the code.
8. **A failing test is more often a wrong test than wrong code** — it has been, six times
   here. The recurring cause is forgetting that a detector on the ground has every steep
   downward direction blocked by the ground at its feet.

## 9. Remaining work, ranked

**Physics**
1. **The detector acceptance `A(E)` is not modelled.** An event rate is
   ∫Φ(E)·A(E)·P(E)dE; the weight used is the flux alone. *Partially doable now*: both
   published curves are in `data/` and `aperture.infer_response()` divides one by our
   geometric model to recover everything else. A better weight than flat, but it
   inherits the published site's geometry.
2. **`min_score` is the dominant assumption.** `--score_percentile` and
   `--stop_at_target` now exist as rank-based alternatives; the configs still use the
   absolute cut. Consider switching them.
3. Column depth is still bounded by the walk unless `max_range_km` is set.
4. Neutral-current regeneration not modelled — Earth-chord suppression overstated.
5. β, the tau energy-loss constant, is an estimate (0.4–1.0×10⁻⁶). Needs a collaboration
   value.
6. Geomagnetic **declination** does not follow the site (inclination does). Needs IGRF.

**Verification**
7. **Nothing has been checked against an external simulation.** The Earth-absorption
   prediction in §6 is the cheapest such test and is ready for someone to run.
8. **The full Arequipa DEM has never been run** — every number is from crops. The owner
   wants this *after* `--explain`. Use `downsample_factor: 4`; expect ~25–30 min.

**Software**
9. ~~`--explain`~~ — done, ROADMAP §6.23.
10. ~~CLI/library parity~~ — done, ROADMAP §6.24.
11. The pipeline still requires `cd src` for relative config paths.
12. No release: nothing on PyPI, Pages not yet deployed (the workflow exists and fires
    from `main`), so several README badges are dark by design.

## 10. Open questions for the owner

1. **IGRF declination per site.** Inclination follows the DEM's coordinates via a dipole;
   declination falls back to Arequipa's −6.9°. The dipole is unreliable for declination
   (−0.2° against a measured −6.9°) and is deliberately not used for it.
2. **β**, as above.
3. **The TAMBO assumptions**, all flagged in `config/tambo_colca_config.json` and
   `docs/source/assumptions.rst`: the 20–60° near-wall band, the 25° far-wall floor, the
   ±20° arrival window, the 0.1 shower-content fraction, γ = 2.0, and `min_score` 0.35.
4. **A prediction worth checking against the collaboration's own acceptance:** the
   effective arrival window should narrow with energy, its lower edge climbing from
   −4.4° at 100 PeV to −0.9° at 10 EeV. If their simulated window does not narrow that
   way, one of the two treatments has the absorption wrong.
5. **A logo.** The design brief is `docs/LOGO.md` -- feed its paragraph to an
   image-generation session. It lands in `docs/source/_static/`, and `conf.py` has the
   `html_logo` line commented out ready to enable.
