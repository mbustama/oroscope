# Handover brief — Oroscope

Written to be fed to a fresh session. It assumes no memory of the previous one.

**Repository:** `mbustama/oroscope`, **public**. Local path `~/Research/GRAND/oroscope`.
A `site_search` symlink sits beside it for anything pointing at the old path; the owner
knows about it and has chosen to keep it.

**Branch:** `dev`, head **`559dc4f`**, pushed. `main` is 7 commits behind — the last
merge was PR #4. **Tests:** 598, stdlib `unittest`, ~30 s. **CI:** 8 jobs, all green.
**Documentation:** <https://mbustama.github.io/oroscope/>, deployed from `main`.

`main` is protected: no direct push, PR required, seven status checks, **zero
approvals**. You can open a PR and wait for checks, but **you cannot merge** — the
permission classifier blocks `gh pr merge`. Open the PR, report the green checks, and
hand the merge command to the owner. Do not try to work around it.

---

## 0. Start here — the four things to do next

The previous session ended mid-request. **Nothing below has been started.** Do these,
roughly in this order; 1 and 2 are one piece of work.

### 1. Four new animations

`tools/make_animations.py` already exists and produces four (`the_walk`, `the_funnel`,
`stride_and_closing`, `energy_window`) as MP4 + GIF. Add these four, agreed after an
explicit "truly useful, not trivially so" filter — the test is whether the intermediate
states carry the argument, not whether it moves:

1. **The azimuth fan at one candidate.** `the_walk` animates *elevation*; nothing
   animates *azimuth*. The scan tests `n_azimuths` bearings within `azimuth_half_width_deg`
   of aspect (`arrival_scan.azimuth_fan(n, half_width)` returns the offsets). Sweep the
   bearing and show which find a wall at the right range and which find sky or the
   candidate's own hillside. Completes the mechanism pair, and acceptance-versus-bearing
   is a polar quantity static figures render badly.
2. **The product score collapsing as components are added.** The strongest of the four.
   `min_score` is the dominant assumption, and the reason a product threshold is
   treacherous is *dynamic*: each component multiplied in drags the whole distribution
   toward zero while the cut stays put. Animate the histogram with the cut fixed and
   watch six components walk the population under the line. Ties directly to the
   measured sweep in notebook 8 (100% → 46% across `min_score` 0.00–0.15).
3. **A criterion biting on real ground.** The notebook sweeps give *how much*; this
   gives *where*. Morph the accepted mask over Colca as `min_target_slope_deg` crosses
   the wall-slope distribution and watch the rims let go. Cheap now — a small run is
   seconds.
4. **The tau's fate through rock.** Propagate a tau with energy and survival falling
   against the production-and-escape optimum (`physics.production_escape_optimum_gcm2`,
   5.7×10⁶ g/cm² at 1 EeV; `physics.tau_survival`). Corrects the commonest
   misconception in the problem — *more rock is not better* — which the code already
   gets right via the depth **band** but which currently lives only in prose.

**Rejected, with reasons — do not build these:** a continuous closing-element sweep
(`stride_and_closing` already covers it), sin α over azimuth (a polar plot does it
statically and better), decay probability vs distance × energy (a contour plot is the
right tool), Earth chord vs elevation (one static diagram), Cherenkov footprint growth
(a simple monotone relation), detector packing filling a site (decorative).

### 2. A new notebook that generates the animations

Notebook **09**, generated from `tools/make_notebooks.py` — **edit the generator, never
the `.ipynb`**. It should build all eight animations and explain each one's point.

The owner asked specifically: **assume the notebook produces MP4 only, and include
instructions for converting to animated GIF.** `tools/make_animations.py` already does
both (`--format mp4,gif`), so the notebook should show the MP4 path and then document
the conversion — both the built-in `--format gif` route and a plain `ffmpeg` recipe, so
a reader with an MP4 from anywhere can convert it.

Notebooks 7 and 8 are excluded from the CI execution job; 09 will need the same
treatment if it runs searches, and `tests/test_docs.py` must still see every API name it
calls. Check `.github/workflows/lint.yml` — the exclusion is a literal `rm` of the two
filenames.

### 3. A Peru-wide search

The owner wants "areas all across Peru, to see the full power of the code". **This needs
sizing before it is attempted** — I did the arithmetic, do not redo it:

Peru is roughly lat −18.4…0.0, lon −81.4…−68.6, so ~18.4° × 12.8°.

| DEM | pixels | `candidate_stride` 5 | 15 | 30 |
| --- | --- | --- | --- | --- |
| SRTM 1-arcsec (30 m) | 3,052 Mpx | 110 GiB | 39 GiB | 22 GiB |
| SRTM 3-arcsec (90 m) | 339 Mpx | 12.7 GiB | **4.8 GiB** | 2.8 GiB |

The machine has ~7–8 GiB free. **So: 3-arcsec at stride 15 is the only single-run option**,
or tile the country at 1-arcsec and combine. Two things that matter:

- **`candidate_stride` is the memory lever, not `downsample_factor`.** Candidates are
  taken on the native grid, so downsampling barely touches the dominant term (roadmap
  §6.26a). This is the single most useful fact for planning the run.
- **Raising the stride costs area unless `gap_close_km` rises with it.** At stride 15 the
  gap is 15 px; the closing element must outrun it or the mask fragments and the area
  collapses — measured at **4.75×** for TAMBO (§6.34). Every run now warns
  (`warn_stride_outruns_closing`), so heed it. Read a Peru-wide area as a coarse survey
  number and say so.

Also needed: a DEM. `oroscope-fetch-dem` downloads from OpenTopography with a free API
key; `oroscope-fetch-roads --dem <path> --places` then gets roads and towns for the map.
Ask the owner for the key rather than guessing at one.

### 4. A question to think about, not to code yet

> After identifying a combined site — say one in the Colca — can the code establish the
> winning configuration for a *specific realization*: ~100 TAMBO units sprinkled with a
> few GRAND antennas?

**Short answer: not yet, and the gap is specific.** My analysis, for the next session to
build on:

*What the code does now.* `count_grid_capacity` lays a regular hex or square lattice
from the site's **bounding-box corner** and counts positions landing on usable ground.
It is **anchored, not fitted** — the docstring says so. That answers *how many fit*, not
*where they should go*.

*What is missing, in order of difficulty.*

1. **The per-pixel score is computed and then thrown away.** `scoring.score_candidates`
   returns a score per candidate; only per-site aggregates (mean/p50/p90) reach the
   results file. Retaining it as a raster is a small change and unlocks everything else.
2. **There is no placement routine.** With a score raster, a greedy or blue-noise
   placement of N units maximising summed score subject to a minimum spacing is
   straightforward and would genuinely answer "where".
3. **"Winning" is only as meaningful as the objective.** The score is a *ranking proxy*,
   not an event rate. Optimising it is defensible but it is not maximising detected
   neutrinos. That needs the differential acceptance `A(E)` — the outstanding physics
   ask (§9.1). `--decay_weight_by` now takes `flux`, `acceptance` or
   `flux_times_acceptance`, but no real `A(E)` table exists, and one *inferred* from a
   published integral curve is demonstrably unsafe (§6.42: it returns zero sites).

*One structural point worth making to the owner.* For a joint array the two experiments
barely compete for ground — a 100 m TAMBO strip sits inside a single 1 km GRAND cell.
The real constraint is the **per-role slope band**: GRAND needs 3–25° and TAMBO 20–60°,
so within one joint site the GRAND antennas belong on the gentler shoulders and the
TAMBO units on the wall. A layout tool should respect that split, and the pipeline
already computes both bands. That is probably the most useful thing to say before any
optimisation is written.

---

## 1. What the previous session did

23 commits, `1fa8810` → `559dc4f`. All of it is in `docs/ROADMAP.md` §6.26 and §6.34–6.43
with the measurements. Condensed:

**The full Arequipa DEM was run** — the first search over all 128.6 Mpx, and the answer
to the question the whole project had been building toward.

| | area km² | sites | capacity | bound by | kept |
| --- | --- | --- | --- | --- | --- |
| GRAND, full DEM | 88,527.5 | 1 | 101,948 | directions accepted | 61.6% |
| TAMBO, full DEM | 111.9 | 26 | 9,024 | directions accepted | 9.7% |
| joint | 50.2 | | | | Jaccard 0.0006 |

Both bound by the same stage the crops were, so the crops were representative about
*what* limits the answer. The rate differs: GRAND 61.6% vs the crop's 60.1%
(indistinguishable), TAMBO 9.7% vs 17.5% (half). **Colca is ordinary ground for GRAND
and exceptional for TAMBO.** GRAND is one contiguous region holding 10× its target;
TAMBO is 26 sites over 310 km. The joint barely moved from the crop's 50.1 km² — 21×
more ground yielded no more co-locatable ground.

**Ten roadmap items closed** (1–6, 9–12). The three that changed a published number:

- **TAMBO's area is low by 4.75×** (§6.34). The stride-1 control separated two effects
  the funnel's 0.53× conflated: acceptance is unbiased (17.494% vs 17.491%) but the mask
  is closed *before* area is measured, and a 100 m element cannot bridge the 154 m gap
  stride 5 leaves. Read Colca as **~397 km² and ~45,856 detectors**. Every run now warns.
- **A silent bug in the sweep path** (§6.35). It passed `rfi_zones` through as the preset
  *name*, which the pipeline iterated character by character — printing `RFI Zones: 8
  active`, one per letter, while excluding nothing. Recorded sweeps are unaffected (they
  used TAMBO's `"none"`). Fixed by `config_to_pipeline_kwargs()` / `run_from_config()`,
  now the single translation for all three callers.
- **Column depth is truncated 6.4× by the walk** (§6.39), and walking 4× further fixes it
  with a byte-identical selection — but 12× collapses the run to 6.0% acceptance.

**The memory estimator was wrong by 2.5×** (§6.26a) and killed the first full-DEM attempt
23 minutes in. It counted only the arrays the scan returns; the peak is inside
`scoring.compose`, ~36 per-candidate arrays. Now 7.21 GiB at `downsample_factor` 1 and
5.08 at 4.

**Maps were substantially rebuilt.** Geographic axes, scale bar, north arrow, altitude
colorbar matched to the panel height, greyscale relief on the overlay so the categories
own the colour, water and nodata distinguished, roads from OpenStreetMap, real town
markers, `--reveal` frames for talks. Details in §2 below.

**Docs and notebooks.** Physics page rewritten with GRAND/TAMBO signal physics treated
separately, a joint particle+radio section, and what each signal *tells you about the
shower*. Notebook 8 documents every assumption, shows its maps inline, and sweeps four
parameters. References are their own page with INSPIRE entries.

---

## 2. Owner preferences — follow these without being asked

Learned by correction over this session. They are consistent and worth honouring.

**Figures**

- **No titles.** The caption carries it. This applies to every map.
- **Legend outside the axes, at the top.** Four columns unless told otherwise.
- **Legend text minimal**: no counts, no parentheticals. "Roads", not "Roads (230, OSM)".
  "Both — 50.1 km²", not with the percentages appended.
- **Colorbar height must match the plot panel.** Use `ss.attach_colorbar`, which takes
  the space from the axes' own divider. `fig.colorbar(fraction=...)` sizes against the
  *figure* and overshoots.
- **Scale bar and north arrow** on every map. `ss.add_scale_bar`, `ss.add_north_arrow`.
- **Grey base with a colorbar** beats a colourful base. When the overlay needs colour,
  the terrain gives it up — greyscale altitude with the bar reads better and leaves the
  hues to the categories.
- **Roads green**, not neutral dark: a thin dark line vanishes into hillshade exactly
  where the ground is steep, which is where the sites are.
- **Only a few labels.** 6 labels, 25 markers is the current setting and was arrived at
  by being told twice that more was too many.
- **Attribution goes in the caption, not on the figure.**
- **Colour scales track the data**, not a fixed range (`ss.altitude_limits`).
- **Progressive-reveal frames** are wanted for talks, and the layout must be pixel-identical
  between frames. Verify it, do not assume it.

**Working style**

- **Proceed without asking** when the path is clear. They have said so explicitly, twice.
  Reserve questions for genuine forks.
- **Do not do trivial work.** When asked for options, filter hard and say what you
  rejected and why.
- **Source data, never invent it.** Town coordinates came from OpenStreetMap rather than
  memory; bibliography entries came verbatim from INSPIRE. This matters to them.
- **Notebooks are educational.** Figures must be shown inline, not merely saved.
- They will feed a brief to a fresh session rather than let context run out — write
  handovers accordingly.

**Project conventions (from the repo, unchanged)**

- Measure before optimising, and after. Several confident hypotheses here were wrong.
- Run examples, don't read them — `tests/test_doctests.py` executes every `Examples`
  block. **Doctest values must be computed, not predicted**; I got them wrong three
  separate times this session by writing plausible numbers.
- Negative results go in `docs/ROADMAP.md` so they are not retried.
- Lint as CI does: `ruff check .` **from the root** — it lints the notebooks too.
- Commit messages explain *why* and state measured deltas.
- The roadmap is updated in the same commit as the code.
- Figure labels capitalise their first word.

---

## 3. Environment and traps

- **`conda activate sssearch` fails** (`conda init` not run). Call the interpreter
  directly: `/home/mbustamante/anaconda3/envs/sssearch/bin/python`. Same environment.
- **`gh` needs `GIT_CONFIG_NOSYSTEM=1`** — the sandbox blocks `/etc/gitconfig` and every
  `gh` call fails without it. `gh pr create --body-file` cannot read from the scratchpad;
  pipe on stdin with `--body-file -`.
- **Jupyter's `python3` kernelspec points at base anaconda**, which has no `oroscope`, so
  `nbconvert --execute` and the docs' `jupyter_sphinx` blocks both fail with
  `ModuleNotFoundError`. Write a kernelspec into a scratch dir and pass `JUPYTER_PATH`:

  ```bash
  mkdir -p /tmp/k/kernels/python3 && cat > /tmp/k/kernels/python3/kernel.json <<'JSON'
  {"argv": ["/home/mbustamante/anaconda3/envs/sssearch/bin/python", "-m",
            "ipykernel_launcher", "-f", "{connection_file}"],
   "display_name": "Python 3", "language": "python"}
  JSON
  cd notebooks && env -u MPLBACKEND JUPYTER_PATH=/tmp/k jupyter nbconvert --execute --inplace 08_the_full_dem.ipynb
  ```

- **`ffmpeg` is the snap build** and cannot write outside its confinement — it will not
  write into `/tmp/claude-*`. Write into the repo's gitignored `output/`, or extract
  frames with PIL.
- **Memory.** The desktop holds ~8 GB of 15. The full DEM needs `--max-memory-gb 7.0`
  passed explicitly; the default 80%-of-available cap is not enough. Never pass 0 — that
  disables the cap and invites the OOM killer.
- **`pgrep -f "script.py"` matches its own command line**, and so does an `echo` of the
  same string in the same command. Use `ps` and check the output.
- **Matplotlib backend.** A library must not choose it. CI asserts that importing
  `oroscope` leaves it untouched. When touching a figure path, **check images are
  actually produced** — count `image/png` outputs, or grep built HTML for `Figure size`.
  Do not trust a green build.
- **Notebook size.** `show_figure` in the generator picks JPEG over PNG above 250 KiB;
  terrain maps are photographic and PNG made notebook 8 4.2 MB.

---

## 4. Repo map (what changed)

| path | what |
| --- | --- |
| `src/oroscope/site_searcher.py` | Pipeline, CLI, config, screening, morphology, capacity, outputs. Now also `run_from_config`, `config_to_pipeline_kwargs`, the map furniture (`add_scale_bar`, `add_north_arrow`, `attach_colorbar`, `add_roads`, `add_settlements`, `altitude_limits`), and `warn_stride_outruns_closing`. |
| `src/oroscope/fetch_roads.py` | **New.** Downloads roads and populated places from OpenStreetMap via Overpass. `oroscope-fetch-roads --dem <path> --places`. |
| `src/oroscope/combine_experiments.py` | The overlay. Now geographic, greyscale-with-colorbar, road- and town-aware, and `--reveal`. |
| `src/oroscope/physics.py` | Adds `set_tau_energy_loss`, `set_declination_model`, `declination_from_grid`, `nc_regeneration_factor`, and `weight_by` on the spectrum fold. |
| `tools/make_animations.py` | **New.** Four animations, MP4 + GIF. Extend for §0.1. |
| `tools/make_notebooks.py` | Generates the notebooks. **Edit here, not the `.ipynb`.** |
| `tools/run_arequipa_full.py` | Full-DEM runner. Finds roads/places automatically, passes `--reveal`. |
| `input/roads/` | Gitignored. 8,780 roads and 1,268 places for the Arequipa DEM. |
| `docs/source/refs.bib` | INSPIRE entries verbatim. Do not retype a record from memory. |

---

## 5. Current numbers

| | area km² | sites | capacity |
| --- | --- | --- | --- |
| GRAND, Colca crop | 4,580.2 | 1 | 5,317 |
| TAMBO, Colca crop | 83.6 → **read as ~397** | 15 → 29 | 9,717 → **~45,856** |
| GRAND, full DEM | 88,527.5 | 1 | 101,948 |
| TAMBO, full DEM | 111.9 (low, see §1) | 26 | 9,024 |
| joint, full DEM | 50.2 | | Jaccard 0.0006 |

**Quote these with their caveats.** TAMBO's areas are lower bounds by ~4.75× from
striding and ~30% again from downsampling at full scale. The joint is limited by TAMBO's
mask, so it is a floor. GRAND is unaffected.

---

## 6. Do not repeat

- **The full Arequipa DEM run.** Done, 26 minutes (GRAND 26.8, TAMBO 1.7), store
  populated, notebook 8 executed against it.
- **The stride-1 controls.** Both GRAND and TAMBO. §6.34.
- **The sensitivity sweeps.** §6.20–6.21.
- **The benchmark baseline.** Refreshed with `--repeat 5`; only `arequipa_2500` is
  gateable on this host. Re-measure with `--update --repeat 5`, never `--update` alone.
- **The memory estimator investigation.** §6.26a.
- **The road and place downloads** for Arequipa — they are in `input/roads/`.
- **Everything in ROADMAP §6 and §7.** Read §6.26, §6.26a, §6.34–6.43 before measuring
  anything.

## 7. Still open

1. **`A(E)` remains the outstanding physics ask.** The selector exists; no real
   differential table does. An inferred one is unsafe (§6.42).
2. **Askaryan (charge-excess) emission is not modelled** — found this session, not yet
   acted on. Only the geomagnetic term is implemented, so `sin α → 0` scores exactly
   zero, and under a *product* composition that zero rejects the site outright when
   charge-excess would leave it ~10–20% efficient. Peru is near the magnetic equator, so
   it is precisely north–south geometries that get zeroed. With a 0.14 floor the
   published 3.7× east-over-north ratio compresses to ~3.4× and nothing lands on zero.
   Lesser omissions found alongside it: tau decay branching (~17% to muons, a flat
   normalisation), |B| magnitude (only sin α is used), and the galactic radio background.
3. **`min_score` → `score_percentile`.** 0.35 ≡ percentile 22.8 on Colca, and a scan
   shows **no knee**. Switching restates published numbers: the owner's call, not yours.
4. **Nothing has been checked against an external simulation.** The Earth-absorption
   prediction (window edge −4.4° at 100 PeV to −0.9° at 10 EeV) is the cheapest such
   test and is ready to run.
5. **No release.** Nothing on PyPI. The logo is a PNG everywhere, which settles the one
   packaging question that was open.
