# Handover brief — GRAND / TAMBO site-search tool

Written to be fed into a fresh session. It assumes no memory of the previous one.

**Branch:** `dev`, 19 commits ahead of `main`, all pushed to `origin/dev`.
**Head at handover:** `a7b4ec1`.
**Tests:** 250, `unittest`, ~25 s, no dependencies beyond what the tool already needs.

**Phase 3 has since been taken a second pass** — see §6.4 for which leads are done,
which were measured and closed, and which remain. The headline conclusion is that the
arrival scan is **not compute-bound**: four separate arithmetic optimisations of the
inner loop have now failed to help, because the profile walk outweighs everything else
in the kernel by 225:1. Further scan work has to reduce samples or memory traffic.

**Before trusting any timing on this machine, read ROADMAP §6.12.** It is a hybrid
CPU (2 P-cores + 8 E-cores) that runs at about a third of its rated clock under load,
and the same unchanged code measured 43.6 s and 39.8 s in consecutive runs. A/B
alternating inside one process, single-threaded, on a subsample.

---

## 1. What the tool does

It searches Digital Elevation Models of the Peruvian Andes for sites suitable for
neutrino observatories. Originally GRAND-only (radio detection of air showers from
Earth-skimming tau neutrinos); being generalised to cover TAMBO (particle detection
across a deep canyon) and, eventually, other experiments.

The pipeline: screen terrain by slope/aspect/altitude/exclusion zones → scan arrival
directions from each surviving pixel to find where a tau could exit and how much rock
lies behind it → score the results → clean up morphologically → label sites and
compute capacity → write GeoTIFF/KML/PNG/JSON.

```bash
cd src && python site_searcher.py --config_path ../config/arequipa_config.json
```

## 2. Environment and machine

- **Use the conda env `sssearch`**: `conda activate sssearch`. It has numba, scipy,
  tifffile, imagecodecs, matplotlib, joblib, tqdm, psutil. It does **not** have pytest,
  which is why the suite is stdlib `unittest`.
- **Cap parallelism at 8 cores.** The machine has 12 but is shared; the owner asked
  explicitly not to saturate it. `bench/benchmark.py` is pinned to 8 via `BENCH_CORES`.
- A long-running job of the owner's (`m1_2nu.py`) is often at ~100% of one core.
  **Check `uptime` before trusting any timing.** The benchmark harness prints a warning
  when the 1-minute load average exceeds 1.0, and records it in `baseline.json`.
- Real DEMs live in `input/dem/` (gitignored, ~250 MB). Tests that need them skip
  automatically when absent.

## 3. Repo map

| path | what |
| --- | --- |
| `src/site_searcher.py` | 2033 lines. Pipeline, CLI, screening, morphology, capacity, outputs. |
| `src/arrival_scan.py` | The scan kernel: profile walking, column depth, Fresnel, RFI line-of-sight. Numba. |
| `src/physics.py` | Closed-form physics: atmosphere, Earth chord, tau range and exit probability, geomagnetic field, Cherenkov footprint. Pure Python/NumPy, no terrain. |
| `src/scoring.py` | Score shapes (band, saturating, ramp) and composition. |
| `src/aperture.py` | Aperture estimate, tabulated response, `infer_response()`. |
| `src/fetch_dem.py` | Downloads DEMs. Was `setup.py`, whose name hijacked `pip install`; renamed with packaging. |
| `src/generate_env.py` | AST-based conda env generator. |
| `tests/` | 250 tests. `synthetic.py` builds terrain with closed-form answers. |
| `bench/benchmark.py` | Per-stage timings and peak RSS on fixed cases, gated at 30% regression. |
| `data/` | Hand-digitised published curves (GRAND Fig. 25, TAMBO Fig. 3). |
| `docs/ROADMAP.md` | 1159 lines. The full record: physics review, every phase, every measurement. **Read it.** |

## 4. Working conventions

These have worked well; please keep them.

1. **Measure before optimising, and measure after.** Several confident hypotheses in
   this project turned out to be wrong when measured (§7).
2. **Every change that alters results gets quantified**, via the golden files, and the
   delta goes in the commit message.
3. **Negative results are recorded**, in `docs/ROADMAP.md`, so they are not retried.
4. **Fixtures are verified before the code that uses them.** `tests/test_fixtures.py`
   checks that the synthetic terrain has the geometry it claims.
5. Regenerate goldens deliberately: `cd tests && UPDATE_GOLDEN=1 python -m unittest test_regression`.
6. Commit messages explain *why*, and state measured deltas.
7. The roadmap is the durable record; update it in the same commit as the code.

## 5. State of each phase

**Phase 0 — foundations. Done.** Test harness, benchmarks, selection funnel (per-filter
survivor counts, in the log and results JSON), provenance (git commit, DEM sha256,
resolved parameters, package versions).

**Phase 1 — physics core. Done, and beyond its original scope.** The engine scans
(azimuth, elevation) arrival directions, finds the first terrain intersection, and
integrates the sub-surface chord for column depth. Scores are composable with named
components. Added along the way: geomagnetic weighting, atmospheric grammage, Earth
chord, Cherenkov footprint, RFI line-of-sight shielding, muon shielding, an estimated
production-and-escape optimum, and both published aperture curves digitised.

**Phase 2 — generalisation. Not started.** See §8.

**Phase 3 — performance. In progress, and the immediate task.** See §6.

**Phase 4 — usability. Not started.** See §9.

---

## 6. Phase 3 — the immediate task

### 6.1 Where the time goes now

End-to-end on a 2500² Arequipa crop, 8 cores, full physics:

| stage | seconds | share |
| --- | --- | --- |
| arrival scan | 43.6 | **91%** |
| morphology | 1.6 | 3% |
| topographic screen | 0.7 | 1.5% |
| capacity analysis | 0.8 | 1.5% |

Phase 3 so far took the whole run from 82.0 s to 47.9 s (30.3 s if bilinear sampling
is disabled). The full Arequipa DEM is ~20× the pixels, so roughly 15 minutes.

### 6.2 Already tried and REJECTED — do not repeat

| attempt | result |
| --- | --- |
| Whole-raster azimuthal sweep | Premise measured wrong. Locality penalty across azimuths is only **1.14×**, not the large effect assumed; the prefetcher handles a constant row stride and 8 threads hide the latency. Would buy ≤14% for a large rewrite that does not extend to a (θ,φ) scan with sub-surface chords. |
| Hoisting the per-sample division | **Slower**, 24.5 s → 26.1 s. The division was already pipelined; the replacement needed two multiplies per sample for the horizon. |
| Solving the ray's exit distance up front to drop the bounds test | ~5%, but broke 405 of 40,000 candidates. The bound is `(cols − c0)/dc_per_m`, not `(cols − 1 − c0)/dc_per_m`; reproducing `int()` truncation at every edge is error-prone. Not worth it. |
| Shuffling candidates for load balance | Works (scaling 3.3× → 5.0×) but destroys locality; net ~10%. Superseded by block-dealing, which is what is in the code. |
| Tabulating the scan's per-bin transcendentals (was lead (f), and more) | Bit-identical, measured **0.975×**. Reverted. |
| Short-circuiting the histogram's bin search | Bit-identical, measured **0.992×**. Reverted. |
| Tuning `numba.set_parallel_chunksize` | ~3% at 8 threads (8.44 s → 8.17 s). Block-dealing already handles the imbalance. |

**The lesson from the last two is the important one:** per (candidate, azimuth) the bin
loop runs 12 times and the profile walk runs ~2700 samples, a ratio of 225:1. Anything
outside the walk is under 1% of the kernel. **The scan is not compute-bound** — do not
spend more effort on flops per sample. See ROADMAP §6.11–6.12.

### 6.3 Already delivered in Phase 3

- **Separable morphology** — a rectangle of ones factorises into a column and a row,
  so O(N(h+w)) not O(Nhw). Bit-identical. 10.5 s → 1.6 s.
- **Slope-space scan loop** — all comparisons are monotonic in elevation angle, so
  comparing against pre-computed tangents removes an arctangent per sample.
- **Block-dealt candidate ordering** — `prange` schedules statically and candidates
  arrive in spatial order, so threads got uniformly-short or uniformly-long work.
  Dealing blocks of neighbours round-robin balances without losing locality. 31.3 → 23.3 s.
- **Bilinear profile sampling** — accuracy fix that costs 1.44×. Changes acceptance by
  +13.4%. `--nearest_sampling` opts out.
- **Streaming DEM cache** — non-evictable memory 623 → 132 MiB, and now bounded by
  block size rather than DEM size.
- **Topographic screen in gradient space** — slope band tested as `tan²(min) ≤ dx²+dy²
  ≤ tan²(max)`, and every filter after slope/altitude works on the surviving subset, so
  `arctan2` is evaluated only where the aspect is read. Byte-identical candidates across
  six filter configurations; 1.25× to 4.47× depending on how much the filters reject.
- **Per-site passes made O(N)** — one `ndimage.mean`, one lookup-table recolouring and
  one stable `argsort` grouping replace three O(sites × pixels) loops. Flat in the site
  count: 0.78 s → 0.22 s at 100 sites, and 900 sites now costs the same as 16.
- **A crash fixed on the way**: `labeled_viz` was `uint8`, so selecting a 256th site
  raised `OverflowError` after the physics had already been paid for. Relevant to
  TAMBO, whose layout is many small sub-arrays. Pinned by `TestManySites`.

### 6.4 Leads — status

**Done:** (a) screening transcendentals, (b) per-site `labeled == site_id`,
(c) `summarize_observables_by_site`, (d) the unused `joblib` import, (f) azimuth
sin/cos, (g) RFI zone masking — (f) measured neutral and was reverted; the rest are in.

**Measured and closed without a change:** (e) the Fresnel pass is 6.5% of the scan
(43.55 s → 46.57 s), not the large share suspected, so fusing it has a 6.5% ceiling.
(i) multiprocessing — screening and morphology together are now ~3% of runtime.

**Still open:** (h) fusing `uniform_filter` with `np.gradient` — note this only bites
when `slope_baseline_m` is set, which the shipped configs do not set. (j) the two-pass
DEM cache build, which is first-run latency rather than throughput. (k) an early exit
in the walk, which is the one lead with real headroom left and is described below.

The original text of the remaining leads follows.

**(a) The screening stage computes `arctan`, `sqrt` and `arctan2` over every pixel.**
`src/site_searcher.py:382-383`:
```python
slope  = np.degrees(np.arctan(np.sqrt(dx**2 + dy**2)))
aspect = np.degrees(np.arctan2(-dx, dy)) % 360
```
The slope is used only for a band comparison. `min ≤ atan(√g) ≤ max` is equivalent to
`tan²(min) ≤ dx²+dy² ≤ tan²(max)`, which needs neither `sqrt` nor `arctan`. And the
aspect is only ever read at surviving pixels (`aspect[cr, cc]`), so `arctan2` could be
evaluated on the masked subset instead of the whole tile. On a 2500² map that is tens
of millions of transcendental calls avoided. This is the same trick already applied
inside the scan kernel, not yet applied here.

**(b) Per-site work is O(sites × pixels).** `src/site_searcher.py:962`
`site_mask_ds = (labeled == site_id)` scans the whole downsampled array once per site,
and `:992` `np.isin(labeled, final_selection_ids)` is similar. `find_objects` slices
are already computed a few lines above — use them to restrict to each site's bounding
box, or build a lookup table indexed by label.

**(c) `summarize_observables_by_site` in `site_searcher.py` masks the full accepted
array once per site.** Same O(S·N) pattern; a single pass with `np.bincount` over site
ids would do it in one sweep.

**(d) `joblib` is imported but unused.** `src/site_searcher.py:10` — `Parallel` and
`delayed` were only needed by the legacy ray-caster, which was removed. Dropping the
import (and the dependency, if nothing else needs it) is free.

**(e) The Fresnel clearance pass re-walks the path for every accepted cell.**
`arrival_scan._min_clearance_ratio`. It runs only for accepted directions, but with
high acceptance that is most of them, and it duplicates the main walk's coordinate
arithmetic. Fusing it into the main walk, or capping the number of directions it
evaluates, is worth measuring — it is a plausible chunk of the 1.44× bilinear cost.

**(f) `sin`/`cos` of the azimuth are recomputed per candidate per azimuth.** With
`use_aspect=True` the azimuth is `aspect_i + offset_a`, so angle-addition from
per-candidate `sin/cos(aspect)` and precomputed `sin/cos(offset)` would replace two
transcendentals with four multiplies and two adds. Minor but cheap.

**(g) RFI zone masking loops over zones with full-array operations per zone**
(`get_candidates_chunked`). With a dozen zones that is a dozen passes; could be one.

**(h) `terrain_derivatives` runs `uniform_filter` then `np.gradient`** — two passes
over the block. A single convolution with the derivative-of-box kernel would do both.

**(i) Multiprocessing.** Everything currently parallel uses Numba threads. The
screening and morphology stages are tile-parallel and embarrassingly so, but are now a
small share of runtime — measure before bothering. If the scan ever needs to scale past
one machine's cores, splitting candidates across processes is straightforward since
they share only the read-only DEM.

**(j) `build_elevation_cache` writes a temporary raw file then converts it** — two
passes over the disk. Decoding segment-by-segment and converting in one pass would
halve the I/O, at the cost of dealing with tifffile's segment API.

**(k) The scan walks to `max_range` always.** There is no early exit, and adding one
is subtle: the depth histogram accumulates over the whole path, and the reported
horizon would be truncated if the walk stopped early. If an early exit is attempted,
be explicit about what it does to `horizon_deg`, which is a reported observable.

### 6.5 Method that worked

- Establish a baseline on a fixed input, then change one thing.
- **Warm the JIT before timing.** A measurement that appeared to show a large memory
  locality effect was JIT compilation in the first case. This cost real effort.
- Check correctness against the previous implementation on a large sample, not just
  the unit tests — several of these changes are numerically subtle.
- Beware the machine's load average.

---

## 7. Key measured facts

Do not re-derive these; they are in `docs/ROADMAP.md` with context.

| finding | value |
| --- | --- |
| Slope depends on measurement baseline | median 17.8° at ~61 m, 10.8° at 1 km |
| Single-azimuth ray under-accepted | rejected 2.2× more than it accepted |
| Morphological closing inflates area | reported area was ~18× the physics-validated area |
| Capacity over-count from integer stamping | +7.4% at 1 km spacing, **+58% at 100 m** |
| The ±3° window sits below the horizon almost everywhere | median horizon 7.3° |
| Geomagnetic asymmetry | east-facing targets worth 3.7× north-facing |
| Shower maturity at Arequipa | 1.07 × X_max |
| Antennas across the radio footprint at 1 km spacing | 0.38 |
| Production-and-escape optimum | 12 km of rock at 100 PeV rising to 23 km at 10 EeV |
| Earth absorption cuts the window | −4.4° at 100 PeV, −2.0° at 1 EeV, −0.9° at 10 EeV |
| Azimuth locality penalty | 1.14× (the sweep premise, measured and rejected) |
| Bilinear vs nearest sampling | acceptance +13.4%, 9.8% of candidates flip |

---

## 8. Phase 2 — generalisation (after Phase 3)

The goal the owner set at the outset: make the tool serve more than GRAND, and combine
criteria across experiments.

**The unifying observation:** GRAND and TAMBO ask the same structural question — *from
this patch of ground, is there a target surface at the right range, in the right
direction, at the right relative orientation, with the right matter behind it?* They
differ in numbers, not structure. One scan engine already answers both.

**Proposed layering:** terrain → criteria (local and view) → experiment spec (YAML) →
region/layout → combination (`all` / `any` / weighted, plus *joint* vs *union* and a
co-location report).

**Known blocker, must be fixed first.** `count_grid_capacity` stamps integer strides,
which over-counts capacity by **+58% at TAMBO's 100 m spacing** — only ~3 pixels span
one detector separation on a 30 m DEM, and three separate `int()` truncations compound.
TAMBO capacity must come from usable area and wall geometry analytically, or from a
resampled grid. Characterisation tests already pin the current behaviour in
`tests/test_capacity.py`; they will fail loudly when it is fixed, which is intended.

**Channels to support**, from the owner's guidance:

| channel | accepted arrival directions |
| --- | --- |
| GRAND neutrinos | −3° to +3° about the horizon |
| GRAND cosmic rays | above the horizon, unless a nearby mountain blocks |
| TAMBO | facing a canyon |

The cosmic-ray channel **inverts** the test — terrain is an obstruction, not a target —
and `require_terrain=False` already expresses that in the same kernel. Any framework
that cannot express both GRAND channels from one scan is not general enough.

**TAMBO parameters** (ref. [2] and owner guidance): 5,000 detection units, 100 m
spacing to start (published nominal is 150 m), triangular grid, Colca Canyon ~1.5 km
deep with ~4.5 km between valley sides, τ range in the valley 50 m – 5 km, shower
3–10 km long and 200 m across, energy reach ~3 PeV – 1 EeV.

Two consequences already established: Colca's walls are ~40°, well outside GRAND's
3–25° deployable band, so **the slope criterion must be per-experiment and probably
per-role** (the far wall wants to be steep, the near wall deployable); and 5,000 units
at 150 m need ~97 km² of wall, i.e. a long strip rather than a compact blob, which the
current `min_width_km` opening would destroy. **Layout models must be per-experiment.**

Also for TAMBO: no Fresnel term, no geomagnetic dependence, footprint set by lateral
particle spread, and `--grammage_mode particle` (the band, not the radio threshold),
because particle content dies after shower maximum whereas radio simply propagates.

---

## 9. Phase 4 — usability (sketch)

- **Auto-detect `origin_lat`/`origin_lon` from the GeoTIFF tiepoint.** Verified present
  and matching the current configs to ~1e-4°. Removes the most error-prone input.
- ~~Rename `src/setup.py`~~ — done; it is now `src/fetch_dem.py`.
- Real packaging: `pyproject.toml`, console entry points, pinned env including `imagecodecs`.
- rasterio/pyproj for CRS and outputs, retiring the hand-written `.tfw`.
- `--explain` funnel report; parameter sweeps with a sensitivity table.
- Revisit the config > fallback > **CLI** precedence — CLI losing to a config file is
  the reverse of what users expect.

---

## 10. Open questions for the owner

1. **IGRF field values per site.** Inclination now follows the DEM's own coordinates via
   a dipole model; declination still falls back to the Arequipa IGRF value (−6.9°) and
   should be supplied per site. The dipole is unreliable for declination (−0.2° vs a
   measured −6.9° at Arequipa) and is deliberately not used for it.
2. **β, the tau energy-loss constant.** Estimated at (0.4–1.0)×10⁻⁶ cm²/g from mass
   scaling, with an assumed energy dependence `0.6×10⁻⁶ (E/1 EeV)^0.20`. Worth pinning
   to whatever the collaboration uses; it moves the optimum in proportion, though not
   the siting conclusion.
3. **Neutral-current regeneration** is not modelled — only charged-current attenuation —
   so the Earth-chord suppression is somewhat overstated.
4. **A prediction worth checking against the collaboration's own acceptance:** the
   effective arrival window should narrow with energy, its lower edge climbing from
   −4.4° at 100 PeV to −0.9° at 10 EeV. If their simulated window does not narrow that
   way, one of the two treatments has the absorption wrong.

---

## 11. Lessons worth carrying forward

- **My confident hypotheses were wrong more than once.** The memory-locality premise
  for the sweep, the flat-with-energy optimum, the harmonic tau range, the claim that
  grammage should be scored as a band for radio. Each survived until measured or
  re-derived. Measure or derive; do not assert.
- **Failing tests were more often wrong tests than wrong code** — four separate times.
  The recurring cause was forgetting that a detector on the ground has every downward
  direction blocked by the ground at its own feet. When a test fails, check the test's
  premise before the code.
- **A rewrite fixed a bug nobody was looking for**, twice: the slope-space rewrite
  exposed an `int()` truncation that counted terrain below the acceptance window, and
  the physics re-derivation found the τ range formula was wrong. Rewrites are worth
  doing partly for this.
- **The benchmark harness's load warning earns its keep.** A "regression" in morphology
  turned out to be the owner's own job competing for CPU.
