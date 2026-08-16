# Handover brief — Oroscope

Written to be fed into a fresh session. It assumes no memory of the previous one.

**Repository:** `mbustama/oroscope`, **public**. Local path `~/Research/GRAND/oroscope`,
with a `site_search` symlink beside it for anything still pointing at the old path —
**delete it when convenient.**

**Branch:** `dev`. **Head at handover: `8c108b1`.** `main` contains everything on `dev`
(merged as PR #2), so the two are level apart from merge commits.
**Tests:** 541, stdlib `unittest`, ~30 s. **CI:** 8 jobs, all green.
**Documentation:** live at <https://mbustama.github.io/oroscope/>, deployed from `main`.

`main` is protected by a repository ruleset: no direct push, no force-push, no deletion,
pull request required, seven status checks, **zero approvals**. So you *can* land work
yourself — open a PR from `dev`, wait for the checks, merge — but you cannot push to
`main` directly.

---

## 1. The full Arequipa DEM ✅ done

**Run on 2026-08-16.** The store is populated, notebook 8 shows the real numbers, and the
findings are in `docs/ROADMAP.md` §6.26. Headline: GRAND 88,527.5 km² in one region with
101,948 detectors; TAMBO 111.9 km² across 26 sites spread over 310 km; joint 50.2 km²,
which is the same 50 km² the Colca crop already had. Both experiments are bound by the
same funnel stage the crops were, so the crops were representative — but TAMBO's
acceptance halves at full scale (9.7% against 17.5%), which says Colca is exceptional
canyon rather than typical.

Read §6.26 before re-running anything here. What follows is how to reproduce it.

**It cost 26 minutes, not ninety:** GRAND 24.2 min, TAMBO 1.2 min. TAMBO is cheap because
its targets are 2-5 km away against GRAND's 10-40.

**It needs `--max-memory-gb`.** The first attempt died 23 minutes in against the default
cap, because the memory estimator under-predicted by 2.5×. That is fixed (§6.26a), but
the default ceiling is still 80% of available, which is not enough on this machine when
the desktop holds half of RAM. Pass the flag.

### 1.1 How to run it

```bash
conda activate sssearch
cd ~/Research/GRAND/oroscope
python tools/run_arequipa_full.py --dry-run                  # costs nothing, starts nothing
python tools/run_arequipa_full.py --max-memory-gb 7.0        # GRAND, TAMBO, then the combination
```

The dry run printed this on 2026-08-16, and is the check to repeat first:

```text
DEM:       input/dem/arequipa_SRTMGL1.tif
estimate:  5.08 GiB at downsample_factor 4
available: 6.4 GiB
would run: grand, tambo, then combine
expected:  ~25 min for grand, ~1 min for tambo
store:     results/arequipa_full
```

**Set the cap deliberately.** 7.0 GiB worked against ~7.9 GiB free and a measured
5.68 GiB peak. Do not pass 0 — that disables the cap and invites the OOM killer, which
is trap 1. And note that `conda activate` may need `conda init` first; calling
`~/anaconda3/envs/sssearch/bin/python` directly is equivalent and avoids it.

`--only grand` or `--only tambo` runs one search and skips the combination, if you would
rather do them separately.

**It writes two places.** The full outputs — GeoTIFF, world file, KML, PNG, log — land in
`output/arequipa_full_{grand,tambo,combined}/`, which is gitignored. The small artefacts
— results JSON, provenance, explanation, a few hundred KB — are copied into
`results/arequipa_full/`, **which is committed**, along with a `manifest.json` recording
when and from what. That store is **populated** as of 2026-08-16.

### 1.2 Why the store exists

[Notebook 8](../notebooks/08_the_full_dem.ipynb) *reads* those results rather than
producing them. GRAND alone is 25 minutes, against CI that executes every notebook on
every push, for something that changes only when a configuration does. (TAMBO is ~1
minute; the whole store rebuilds in 26, not the ninety once assumed.) This works only
because
`explain.explain_results()` is a pure function of the results dictionary — no DEM,
nothing re-run.

Notebooks 7 and 8 are therefore **excluded from the CI execution job**, and
`tests/test_docs.py` checks statically that every API name they call still exists.

After the run: re-execute notebook 8 so its stored outputs show the real numbers.

```bash
cd notebooks && env -u MPLBACKEND jupyter nbconvert --execute --inplace 08_the_full_dem.ipynb
```

### 1.3 The configurations

`config/grand_arequipa_full.json` and `config/tambo_arequipa_full.json`. They are the
Colca crop configs with exactly three changes: the full DEM, `origin_lat`/`origin_lon`
null so the corner is read from the file's own tiepoint, and `downsample_factor: 4`
instead of 1. Every criterion is otherwise unchanged, **deliberately** — the point of
this run is scale, not a different question.

`downsample_factor: 4` has a price worth stating wherever the numbers are read: area is
measured on the downsampled mask while capacity is measured at full resolution, so a
feature a few pixels wide keeps its detectors and loses area. That matters more for
TAMBO's canyon strips than for GRAND's blobs. **Read these areas as lower bounds.**

### 1.4 What was looked at, in this order — all four answered in roadmap §6.26

1. **Is the binding constraint the same one the crops found?** This is the most
   consequential question in the whole run. If a full DEM is bound by a different funnel
   stage than its crops were, the crops were not representative and every number derived
   from them needs re-reading. For comparison, at Colca both were bound by
   `directions accepted` — GRAND keeping 60.1%, TAMBO 17.5%.
2. **The area**, against the crop scaled up, and against the closing factor *this run
   reports for itself* rather than the 2.29× quoted from Colca.
3. **The site count and their spread.** A crop cannot say whether the good ground is one
   region or fifty scattered ones. That is a deployment question, not a physics one.
4. **The weakest score component.** On the crops it is `solid_angle` at 15 of 15 TAMBO
   sites and at GRAND's single site. If that holds at full scale it is a statement about
   the criterion — about `solid_angle_half_sr` — rather than about Peru.

Record the answers in `docs/ROADMAP.md` §6.26, which is where the expectation is written
down, and correct §6.26 if the run contradicts it.

---

## 2. What does NOT need to be run again

Read this before measuring anything. Several of these cost real time or crashed the
machine.

**Do not re-run:**

- **The Colca searches.** Both configs were re-run at this head; the outputs in
  `output/grand_colca_config/`, `output/tambo_colca_config/` and
  `output/combined_colca/` are current, and carry the site coordinates and named score
  components added this session. Numbers in §5.
- **The sensitivity sweeps.** Both the single-energy and spectrum-folded sweeps are in
  `docs/ROADMAP.md` §6.20–6.21 with their tables. Re-running costs ~10 minutes and will
  reproduce them.
- **The stride-1 control runs**, at both GRAND and TAMBO settings.
  `config/grand_colca_stride1.json` and `config/tambo_colca_stride1.json` exist and their
  results are recorded. Striding is unbiased in acceptance at both. At TAMBO's 100 m
  element it costs **4.75× of the reported area** — roadmap §6.34.
- **The notebooks.** All eight are committed *with their outputs* and were executed
  against the installed package at this head.
- **DEM downloads.** `input/dem/` holds `arequipa_SRTMGL1.tif`, `lima_AW3D30.tif` and the
  derived `colca.tif` crop. **The Arequipa DEM already covers Colca Canyon** — verified,
  1673 m of incision against the published ~1.5 km.

**The benchmark baseline is refreshed** (2026-08-16), and the way it was done is the
part to carry forward. This machine still cannot resolve a 30% change on a short stage —
two identical passes gave `arequipa_900/ray_tracing` 1.228 s then 1.820 s, 48% apart — so
rather than wait for a quiet machine the harness was taught to cope:
`--repeat N` keeps the per-stage **minimum** (timing noise is one-sided: nothing runs
faster than its true cost), `spread_pct` records what the machine could resolve, and the
gate skips any stage whose spread exceeds half of it. **Re-measure with
`python bench/benchmark.py --update --repeat 5`**, never `--update` alone. Only
`arequipa_2500` is genuinely gateable here. ROADMAP §6.37; §6.25 for the history.

**Do not re-derive:** everything in §6. **Do not retry:** everything in §7.

---

## 3. Environment, machine, and the traps

- **Use the conda env `sssearch`**: `conda activate sssearch`. It has numpy, scipy,
  numba, tifffile, imagecodecs, matplotlib, tqdm, sphinx and the docs stack, coverage,
  ruff, jupyter. It does **not** have pytest, which is why the suite is stdlib
  `unittest`.
- The package is installed editable (`pip install -e .`), so `import oroscope` works from
  anywhere and the five console scripts (`oroscope`, `oroscope-combine`, `oroscope-crop`,
  `oroscope-sensitivity`, `oroscope-fetch-dem`) are on `PATH`.
- **Cap parallelism at 8 cores.** The machine has 12 but is shared.
- **The bundled configs resolve `dem_path` relative to `src/`.** Run them from there, or
  pass an absolute path. This is a known wart, §9.11.

**Trap 1 — memory.** A ten-point sensitivity sweep once reached 6.9 GB and was killed by
the OOM killer, taking other work with it. The cause was ours (a leaked matplotlib figure
per run) and is fixed, but the safeguards matter: every run prints an estimate against
available memory, caps its own address space at 80% of available, and `sensitivity.py`
runs each point in a subprocess. The machine has 15 GB with 6–7 typically free.

**Trap 2 — timings are unreliable here.** A hybrid CPU: 2 P-cores (4.6 GHz) and 8
E-cores (3.4 GHz), running at about a third of rated clock under load. **A/B alternating
inside one process, single-threaded, on a subsample.** Never trust a before/after taken
in consecutive whole-suite runs. See ROADMAP §6.12 and §6.25.

**Trap 3 — the matplotlib backend, which has now bitten three times.**

1. `MPLBACKEND=Agg` as an *environment variable* propagates into child kernels and
   overrides the inline backend, so notebooks stored no figures and documentation pages
   rendered figures as the text `<Figure size ...>`. Both built clean and reported no
   error.
2. `index.rst` lacked the `%matplotlib inline` setup cell the other pages carry, so the
   front-page diagram published as that same literal text.
3. **`import oroscope` forced the backend.** `combine_experiments` called
   `matplotlib.use("Agg")` at module level — harmless for a standalone module, fatal for
   a package front door, where it reached into every caller's session and killed inline
   figure capture in all eight notebooks.

The rule: **a library must not decide how its user's figures are rendered.** The CI
packaging job now asserts that importing `oroscope` leaves the backend untouched. And
whenever you touch a figure path, **check that images are actually produced** — count
`image/png` outputs in the notebooks, or grep the built HTML for `Figure size`. Do not
trust a green build.

**Trap 4 — lint the way CI does.** `ruff check .` from the repository root. It lints the
notebooks too, and `ruff check src/ tests/` does not — an unused `plt` import in a
notebook that draws nothing failed CI after a local check had passed.

**Trap 5 — `pgrep -f "some/script.py"` matches its own command line.** A wait loop built
on it never exits, because it is waiting for itself. Two zombie tasks this session.

**Trap 6 — notebooks 7 and 8 cannot be executed locally with the repo's kernelspec.**
The `python3` kernel registered on this machine points at `~/anaconda3/bin/python3`,
which is base and has no `oroscope`, so `jupyter nbconvert --execute` fails with
`ModuleNotFoundError` on the import cell. CI does not hit this because it installs the
package into the runner's default python — and notebooks 7 and 8 are the two CI never
executes, so nothing catches it. This appeared when the flat modules became a package
and the notebooks' `sys.path` insert was removed. Either register a kernelspec for
`sssearch`, or point `JUPYTER_PATH` at one:

```bash
mkdir -p /tmp/k/kernels/python3 && cat > /tmp/k/kernels/python3/kernel.json <<'JSON'
{"argv": ["/home/mbustamante/anaconda3/envs/sssearch/bin/python", "-m",
          "ipykernel_launcher", "-f", "{connection_file}"],
 "display_name": "Python 3", "language": "python"}
JSON
cd notebooks && env -u MPLBACKEND JUPYTER_PATH=/tmp/k jupyter nbconvert --execute --inplace 08_the_full_dem.ipynb
```

---

## 4. Repo map

| path | what |
| --- | --- |
| `src/oroscope/__init__.py` | The package front door: re-exports 131 names. `import oroscope` is the whole setup. |
| `src/oroscope/site_searcher.py` | 3945 lines. Pipeline, CLI, config files, screening, morphology, capacity, outputs, memory guards. |
| `src/oroscope/arrival_scan.py` | The scan kernel: profile walking, column depth, Fresnel, RFI line-of-sight. Numba. |
| `src/oroscope/physics.py` | Closed-form physics, no terrain: atmosphere, shower profile, Earth chord, tau range and decay, geomagnetic, Cherenkov. |
| `src/oroscope/explain.py` | The run summary and the combination summary. Pure functions of the results dict. |
| `src/oroscope/scoring.py` | Score shapes (band, saturating, ramp) and composition. |
| `src/oroscope/aperture.py` | Aperture estimate, tabulated response, `infer_response()`. |
| `src/oroscope/combine_experiments.py` | Overlays two or more runs: joint, union, co-location. |
| `src/oroscope/crop_dem.py` | Cuts a lat/lon window out of a DEM. |
| `src/oroscope/sensitivity.py` | One-at-a-time parameter sweeps, each point in a subprocess. |
| `src/oroscope/figures.py` | The publication figures. **States the label convention** (§8.9). |
| `src/oroscope/fetch_dem.py` | Downloads DEMs. Was `setup.py`, whose name hijacked `pip install`. |
| `tests/` | 541 tests across 14 files. `synthetic.py` builds terrain with closed-form answers. |
| `tools/make_notebooks.py` | Generates the eight tutorials. **Edit here, not the `.ipynb`.** Only rewrites what changed. |
| `tools/run_arequipa_full.py` | The full-DEM runner. §1. |
| `results/arequipa_full/` | The committed store notebook 8 reads. Currently empty but for its README. |
| `docs/source/cli.rst` | The command line, with the complete 82-option reference generated from the parser. |
| `docs/source/assumptions.rst` | The blunt list of what the numbers rest on. |
| `docs/ROADMAP.md` | ~2000 lines. The durable record. **Read §6.11, §6.12, §6.20–6.33.** |
| `bench/benchmark.py` | Per-stage timings and peak RSS, gated at 30% regression. |

## 5. Current numbers, at this head

| | area | sites | capacity |
| --- | --- | --- | --- |
| GRAND, Colca crop | 4580.2 km² | 1 | 5317 |
| TAMBO, Colca crop | 83.6 km² | 15 | 9717 |
| **joint** | 50.1 km² | | Jaccard 0.0109 |
| **union** | 4613.7 km² | | |
| GRAND, full DEM | 88,527.5 km² | 1 | 101,948 |
| TAMBO, full DEM | 111.9 km² | 26 | 9024 |
| **joint, full DEM** | 50.2 km² | | Jaccard 0.0006 |
| **union, full DEM** | 88,589.2 km² | | |

The crop rows are measured at `downsample_factor: 1` and the full-DEM rows at 4, so
**TAMBO's two areas are not commensurable** — a canyon strip loses ~30% of its area to
downsampling while keeping its detectors. Compare the acceptance rates instead, which are
measured on the same grid in both. §6.26 of the roadmap does this properly.

**Co-location is decided by slope, not arrival geometry.** GRAND's 3–25° deployable band
against Colca's ~40° walls leaves only a 20–25° sliver — 23% of the narrower band. What
the two ask of the *view* (distance window, arrival elevations) differs freely and is no
obstacle: a pixel has one slope and both must accept it.

### 5.1 Numbers to quote carefully

- **Reported area is not physics-accepted area.** Morphological closing inflated it
  **2.29×** at Colca, measured against a stride-1 control. Each run now reports the
  factor for itself: **2.19× for GRAND** — an independent check on that number, agreeing
  to 4% — but **0.53× for TAMBO**, whose 100 m element cannot bridge the gaps
  `candidate_stride: 5` leaves.
- **TAMBO's area is a lower bound by 4.75×, now measured** (roadmap §6.34). The stride-1
  control at TAMBO settings separates the two effects the 0.53× conflated: acceptance is
  unbiased (17.494% against 17.491%), closing alone inflates 1.17×, and the fragmentation
  of a mask marked one pixel in five and closed with an element too small to reconnect it
  costs **4.75×**. Read TAMBO's Colca area as **~397 km², not 83.6**, and its capacity as
  **~45,856, not 9,717**. The full-DEM 111.9 km² is under-reported for the same reason
  plus downsampling, and cannot be measured directly here (26.8M candidates, ~10 GiB).
  GRAND is unaffected at either scale. Every run now warns when the element cannot outrun
  the stride gap.
- **`solid_angle` is the weakest score component at 15 of 15 TAMBO sites** and at GRAND's
  single site. The Colca result is set almost entirely by `solid_angle_half_sr`.
- **TAMBO's capacity varies by 1.46× across a plausible spectral index**, having varied
  without bound (10878 → 0) when the decay was evaluated at a single energy. `min_score`
  is the dominant assumption at 2.38× to 0.20×.

## 6. Key measured facts — do not re-derive

| finding | value |
| --- | --- |
| Slope depends on measurement baseline | median 17.8° at ~61 m, 10.8° at 1 km |
| Morphological closing inflates area | **2.29×** at a 1 km element (stride-1 control) |
| The same, from a run's own funnel | 2.19× GRAND, **0.53× TAMBO** (100 m element) |
| Candidate striding, acceptance | unbiased at both elements: 60.1% vs 60.1% (GRAND), 17.491% vs 17.494% (TAMBO) |
| Candidate striding, **area** | costs **4.75×** at a 100 m element, nothing at 1 km — the element must outrun the stride gap (154 m at stride 5) |
| Capacity over-count, integer stamping | +7.4% at 1 km, **+58% at 100 m** — fixed |
| Capacity over-count, bounding box vs region | **+38%** on a canyon network — fixed |
| The ±3° window sits below the horizon almost everywhere | median horizon 7.3° |
| Geomagnetic asymmetry | east-facing targets worth 3.7× north-facing |
| Earth absorption narrows the window | −4.4° at 100 PeV, −2.0° at 1 EeV, −0.9° at 10 EeV |
| Colca crossing supplies | ~170 g/cm² across 2 km, ~390 across the full 4.5 km |
| Shower maximum needs | 561 g/cm² at 3 PeV, 700 at 1 EeV |
| Tau decay length | 147 m at 3 PeV, 49 km at 1 EeV |
| Bilinear vs nearest sampling | acceptance +13.4%, 1.44× cost |
| Scan thread scaling | 1.85× at 2 threads, **3.70× at 8** — hardware, not scheduling |
| Far-wall slope recovered at Colca | 34.7–44.3°, median 38.6°, against a published ~40° |
| Storing the named score components | 1.62× on the per-site aggregation, +48 ms on a 19 s search |

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
~2700, a ratio of 225:1. **The scan is not compute-bound.** Anything further must reduce
*samples* or *memory traffic*, not flops.

## 8. Working conventions — please keep

1. **Measure before optimising, and after.** Several confident hypotheses here were wrong.
2. **Run examples, don't read them.** `tests/test_doctests.py` executes every `Examples`
   block. Eight were plausible, close, and wrong.
3. **Negative results go in `docs/ROADMAP.md`** so they are not retried.
4. **Fixtures are verified before the code that uses them** (`tests/test_fixtures.py`).
5. Regenerate goldens deliberately:
   `cd tests && UPDATE_GOLDEN=1 python -m unittest test_regression`.
6. **Lint as CI does: `ruff check .` from the root.** Trap 4.
7. Commit messages explain *why* and state measured deltas.
8. The roadmap is updated in the same commit as the code.
9. **Figure labels capitalise their first word** — axis labels, titles, legend entries and
   annotations, in the notebooks and in the maps the pipeline writes. Labels beginning
   with a function name or a symbol are left alone. Stated at the top of `figures.py`.
10. **Documentation is library-first.** Show `import oroscope` and a call before showing a
    shell line. The command line has its own page.
11. **A failing test is more often a wrong test than wrong code** — it has been, seven
    times here. The recurring cause is forgetting that a detector on the ground has every
    steep downward direction blocked by the ground at its feet.

## 9. Remaining work, ranked

**Physics**
1. **The detector acceptance `A(E)` is not modelled.** An event rate is
   ∫Φ(E)·A(E)·P(E)dE; the weight used is the flux alone. *Partially doable now*: both
   published curves are in `data/` and `aperture.infer_response()` divides one by our
   geometric model to recover everything else.
2. **`min_score` is the dominant assumption.** `--score_percentile` exists as the
   scale-free alternative; the configs still use the absolute cut. Consider switching.
3. ~~Column depth is bounded by the walk unless `max_range_km` is set.~~ ✅ measured,
   2026-08-16. It is a **6.4× under-report** at the default, and walking 4× the distance
   window fixes it with an *identical* selection. But at 12× the same run keeps 6.0% of
   directions against 17.5% and collapses — the knob is not monotone. Configs unchanged;
   the run now reports the factor and the cliff. Roadmap §6.39.
4. Neutral-current regeneration not modelled — Earth-chord suppression overstated.
5. ~~β, the tau energy-loss constant, is an estimate.~~ ✅ configurable, 2026-08-16 —
   `physics.set_tau_energy_loss()`. Still wants a collaboration value, but adopting one
   no longer means editing source. **Note it does not affect a search**: β enters tau
   range and survival through rock, which the search does not model; the search uses
   the decay length E/m·cτ, which carries no β. The run's summary had been listing it
   as an assumption behind its numbers, and no longer does. Roadmap §6.38.
6. ~~Geomagnetic **declination** does not follow the site.~~ ✅ socket added,
   2026-08-16. `physics.set_declination_model(fn)` and `declination_from_grid()`; the
   pipeline consults them before the constant fallback. **No IGRF table is shipped on
   purpose** — its coefficients are somebody else's published work and typing them from
   memory would give plausible wrong declinations. Install `ppigrf`/`pyIGRF` and pass
   its function, or supply a NOAA grid. Roadmap §6.40.

**Verification**
7. **Nothing has been checked against an external simulation.** The Earth-absorption
   prediction in §6 is the cheapest such test and is ready for someone to run.
8. ~~**The full Arequipa DEM.**~~ ✅ done, 2026-08-16. §1 and roadmap §6.26.

**Software**
9. ~~**A stride-1 control at TAMBO settings.**~~ ✅ done, 2026-08-16. It is a lower bound
   by **4.75×**. Roadmap §6.34; `config/tambo_colca_stride1.json`.
10. ~~**The config→pipeline translation is duplicated three times.**~~ ✅ done,
    2026-08-16. `config_to_pipeline_kwargs()` / `run_from_config()`; all three callers
    use them. It had already cost a bug: the sweep child never resolved `rfi_zones`, so
    a preset name was iterated character by character and a sweep on a GRAND config ran
    with **no exclusion zones while printing `RFI Zones: 8 active`**. The recorded
    sweeps are unaffected (they used TAMBO's `"none"`). Roadmap §6.35.
11. ~~**The pipeline resolves paths relative to the working directory.**~~ ✅ done,
    2026-08-16. `load_config()` resolves `dem_path`, `road_map_path` and `resume_dir`
    against the configuration's own directory, and the output base follows the same
    rule, so a search runs identically from anywhere. No shipped config changed:
    `config/` and `src/` are both one level below the root, so `../input/...` names the
    same file either way. A cwd-relative path still works, with a warning. **Note
    `oroscope-fetch-dem` is not covered** — it writes `../input/dem/` and `../config/`
    relative to the working directory and has no configuration file to be relative to,
    so it still wants running from `src/`.
12. ~~**Refresh `bench/baseline.json` on a quiet machine.**~~ ✅ done, 2026-08-16 —
    without a quiet machine, by measuring what this one can resolve first. Two
    identical passes disagreed by **48%** on a short stage, so `--repeat N` now keeps
    the per-stage **minimum** (noise is one-sided), `spread_pct` records the resolution
    alongside the measurement, and the gate ignores any stage whose spread exceeds half
    of it. Refreshed with `--repeat 5`. **Only `arequipa_2500` is really gateable here**
    (ray tracing 8.4% spread); `synthetic_900/ray_tracing` spreads 149.6% and is now
    excluded rather than failing builds at random. Roadmap §6.37.
13. **No release.** Nothing on PyPI. The logo question is settled: it is a **PNG
    everywhere** (1024×1024 RGBA), because PyPI does not render SVG and a project
    carrying both formats eventually ships two different logos.

## 10. Open questions for the owner

1. **IGRF declination per site.** The mechanism now exists (§6.40) and wants feeding:
   either add `ppigrf` as a dependency and pass its function to
   `physics.set_declination_model()`, or export a NOAA declination grid covering the
   Arequipa DEM into `data/` and load it with `declination_from_grid()`. Which of those
   the collaboration prefers is the open question. Until then the constant −6.9°
   fallback stands, which is right for southern Peru.
2. **β**, as above.
3. **The TAMBO assumptions**, all flagged in `config/tambo_colca_config.json` and
   `docs/source/assumptions.rst`: the 20–60° near-wall band, the 25° far-wall floor, the
   ±20° arrival window, the 0.1 shower-content fraction, γ = 2.0, and `min_score` 0.35.
   Add to that **`solid_angle_half_sr`**, which §5.1 shows is what the TAMBO result
   actually turns on.
4. **A prediction worth checking against the collaboration's own acceptance:** the
   effective arrival window should narrow with energy, its lower edge climbing from
   −4.4° at 100 PeV to −0.9° at 10 EeV. If their simulated window does not narrow that
   way, one of the two treatments has the absorption wrong.

---

## 11. What the previous session did

Context on why things are as they are. All of it is in `docs/ROADMAP.md` §6.23–6.33 with
the measurements.

**Delivered:**

- **`--explain`, on by default.** Every run prints a plain-language account of itself and
  saves it as `explanation.txt`: what was found, which funnel stage set the size of the
  answer and the parameter behind it, **why each site is good** (the criteria it
  satisfies, each with the measurement that earned it) and what held it back, the closing
  factor measured from the run itself, which numbers are assumptions, what energy the
  geometry favours, and what to try next. `oroscope-combine` gained the same treatment,
  including *which screening band decides co-location*.
- **CLI/library parity.** `max_memory_gb`, `preflight_memory()`, `load_config()`,
  `generate_config()`, `default_config()`, and the pipeline **returns its results dict**.
- **`import oroscope`.** The flat modules became a real package; the notebooks' `sys.path`
  insert is gone.
- **Site records carry coordinates** and their named score components.
- **Documentation:** a CLI page with the full option reference, library-first README and
  quickstart, the logo, and eight notebooks (7 explains a run — including one that finds
  nothing — and 8 is the full DEM).
- **Coverage 65% → 74%**; `combine_experiments` 30% → 88%, `crop_dem` 18% → 90%.

**Bugs found, several of which had been producing plausible wrong numbers:**

| what | effect |
| --- | --- |
| `oroscope-combine` read a stale `.tif` | TAMBO 44.5 → **83.6 km²**, joint 26.4 → **50.1** |
| The results file listed more sites than were selected | area over-reported wherever `stop_at_target` was used |
| `oroscope-combine` crashed on any single-mode run | two faults in one path |
| A `geomagnetic` component appeared with the weighting *off* | a disabled criterion listed among a site's strengths |
| `main()` leaked its log file and never restored `sys.stdout` | |
| Ten defaults disagreed across three sources | `origin` defaulted to 0.0 — a *valid* coordinate in the Gulf of Guinea |
| `import oroscope` forced the matplotlib backend | all notebook figures silently lost |
| The map title said "GRAND site search" on every run | including TAMBO's |
| The README documented 34 of 83 options, one imaginary | plus the precedence rule backwards |

**One worth reading in full (ROADMAP §6.30):** the first version of the combination
summary reasoned about *why* two experiments could share ground, and was confidently
wrong — it concluded GRAND and TAMBO "cannot share ground at all", printed directly above
the 50.1 km² they demonstrably share. It had treated the viewing windows as shared
constraints when only properties of the ground itself are shared. **A summary that
reasons can be wrong in ways one that merely restates cannot, and it will be wrong
persuasively.** Check any such feature against a case whose answer you already know.
