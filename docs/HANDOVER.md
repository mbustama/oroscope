# Handover brief — Oroscope

Written to be fed to a fresh session. It assumes no memory of the previous one.

**The next session's job is to audit the code.** §0 says what to audit and where the
bodies are likely buried. Everything after that is context.

**Repository:** `mbustama/oroscope`, **public**. Local path `~/Research/GRAND/oroscope`.
A `site_search` symlink sits beside it for anything pointing at the old path; the owner
knows about it and has chosen to keep it.

**Branch:** `dev`, head **`c946960`**, pushed, **CI all green** — 8 jobs: `docs`, `ruff`,
`Notebooks execute`, and tests on Python 3.9 through 3.13. `main` is **21 commits
behind**; the last merge was PR #4 and **there is no open PR**. **Tests:** 609, stdlib
`unittest`, ~30 s. **Documentation:** <https://mbustama.github.io/oroscope/>, from `main`.

`main` is protected: no direct push, PR required, seven status checks, **zero approvals**.
You can open a PR and wait for checks, but **you cannot merge** — the permission
classifier blocks `gh pr merge`. Open the PR, report the checks, hand the merge command to
the owner. Do not try to work around it.

---

## 0. The audit — start here

Nothing is half-finished. What follows is where to look, in the order I would look.

### The known defect, still unfixed on purpose

**`physics.tau_exit_probability` under-resolves its own integral** (§6.44). It spreads
`samples` uniformly over `[0, X]`, but only interactions within about one tau range of the
far surface contribute, so as `X` grows the spacing outruns the only region that matters.
At 3 PeV and `X` = 10⁹ g/cm² the default is **8× the converged value and inverts the trend
with depth**. `depth_band_from_energy` inherits it and returns a band excluding the true
1 EeV optimum by an order of magnitude.

**No published number is affected** — every config leaves `depth_band_gcm2` null. Left
alone because the fix changes a physics function's outputs and wants its own commit with a
test pinning converged values. The substitution `u = X − x` on a log grid is the obvious
form. **Owner's call.**

### Where I would audit hardest, and why

1. **The memory model.** Three quantities get confused and it has cost a machine —
   anonymous memory (`estimate_peak_memory_gb`), the map's separate peak
   (`estimate_visualisation_memory_gb`), and the `RLIMIT_AS` cap. The estimator is
   **~2× optimistic at `candidate_stride` 1** (measured: 3.27 GiB predicted, 5.31 GiB
   RSS, 6.51 GiB virtual). It is calibrated on a *strided* run and has exactly one
   data point at stride 1. `n_scoring_arrays` deserves a second calibration.
2. **The `combine` step's memory is not modelled at all.** It renders at full resolution
   where the search's map renders at `downsample_factor * 2`, so it needs ~4× the map's
   estimate. It failed twice this session and was re-run standalone both times. This is a
   real gap in `estimate_visualisation_memory_gb`.
3. **Funnel rows are read by position in places.** Reading positionally is a live trap —
   a run with RFI zones carries an extra stage, and it made GRAND's Arequipa acceptance
   look like 20% when it is 61.6%. `compare_regions.py` reads by name; check everything
   else does.
4. **`count_grid_capacity` is anchored, not fitted.** Every capacity in the project is an
   estimate for an arbitrarily placed array. Worth confirming nothing quietly treats it
   as an optimum.
5. **Numba kernels are excluded from coverage** (`^\s*@jit` in `[tool.coverage.report]`),
   so a regression there is invisible to coverage. `tests/test_arrival_scan.py` drives
   them against closed-form terrain; check that is still true after any change.
6. **The three "silently wrong once" classes** listed in `docs/source/implementation.rst`
   — a preset name iterated character by character, a component that appeared when
   switched off, a results prefix naming one experiment. Each produced plausible output
   while being incorrect. Look for more of that shape.

### Open questions, unchanged

- **`A(E)` is the outstanding physics ask.** The selector exists (`--decay_weight_by`);
  no real differential table does, and an inferred one is demonstrably unsafe (§6.42).
- **Askaryan (charge-excess) emission is not modelled.** Only the geomagnetic term is, so
  `sin α → 0` scores exactly zero and under a *product* composition that zero rejects the
  site outright. Peru is near the magnetic equator, so north–south geometries get zeroed.
- **Nothing has been checked against an external simulation.** The Earth-absorption
  prediction (window edge −4.4° at 100 PeV to −0.9° at 10 EeV) is the cheapest such test.
- **`min_score` → `score_percentile`** (§6.43). 0.35 ≡ percentile 22.8 on Colca and a scan
  shows **no knee**. Switching restates published numbers: the owner's call.

---

## 1. What the previous session did

Twenty commits, `35d7814` → `c946960`. Measurements in `docs/ROADMAP.md` §6.44–6.52.

**Animations and notebooks.** `tools/make_animations.py` went to eight animations. The
notebooks were **renumbered** at the owner's request and now run: 07 animations, 08
explaining a run, then a per-region block — 09 Arequipa, 10 Ancash, **11 reserved
originally for Lima and now Lima**, 12 Peru — then 13 turning the knobs. Generator
variables are **content-named** (`NB_ANCASH`, not `NB10`) so the next region does not force
a rename.

**Three department runs, one dataset.** Ancash and Lima were added; Lima was
**re-downloaded as SRTMGL1** to replace AW3D30, because a dataset difference would have sat
inside every comparison as a confound.

| | Arequipa | Ancash | Lima |
| --- | --- | --- | --- |
| median slope | 11.1° | 23.0° | 20.4° |
| GRAND per pixel | 1.00× | 0.91× | 0.72× |
| TAMBO per pixel | 1.00× | 2.93× | 2.09× |

The answer tracks steepness in opposite directions for the two experiments.

**The biggest finding: the striding penalty for TAMBO is not 4.75×.** On the Callejón de
Huaylas crop, run unbiased at `1 / 1` against a `4 / 5` control on identical ground, TAMBO
lost **291× in area and 386× in capacity** while GRAND moved 1.1×. Acceptance was
**identical** at 14.0% both ways. All of the loss happens between closing and selection:
the mask fragments into 7,954 regions of which one clears `min_sub_array_size`. **Every
strided TAMBO area and capacity in this project is a lower bound by a terrain-dependent
factor with no useful upper limit.**

**A correction I had to make to my own work an hour after writing it.** The joint region's
share of TAMBO's mask sat at 44.9 / 43.0 / 46.2% across the departments and looked like a
constant of the two experiments. Both **unbiased** crops give ~73%. There are two
constants, and the strided one is not real. Quote ~72–75%.

**Docs.** Four new pages — `howitworks.rst` (the vocabulary, with three new schematics in
`oroscope.figures`), `glossary.rst`, `implementation.rst`, `data.rst` — plus notebook 13
and a co-location helper.

**A memory safeguard, after the machine went down twice.** See §3, Trap 1.

---

## 2. Owner preferences — follow these without being asked

**Figures**

- **No titles.** The caption carries it. Applies to every map.
- **Legend outside the axes, at the top**, four columns unless told otherwise — but build
  it from the categories that actually occur.
- **Legend text minimal**: no counts, no parentheticals.
- **Colorbar height must match the panel.** Use `ss.attach_colorbar`.
- **Scale bar and north arrow** on every map.
- **Grey base with a colorbar** beats a colourful base. **Roads green**, not dark.
- **Only a few labels** — 6 labels, 25 markers.
- **Capitalise the first word** of every axis label, legend entry and annotation.
- **Attribution in the caption, not on the figure.**
- Publication-quality figures are wanted for talks; they are library functions so they can
  be exported at any dpi.

*Known inconsistency:* the pipeline's own search map still prints counts and parentheticals
in its legend, against the rule above. Not changed — never in scope.

**Working style**

- **Proceed without asking** when the path is clear. Said explicitly, several times.
- **Do not do trivial work.** When asked for options, filter hard and say what you
  rejected. When asked to assess, give a recommendation with reasons — the owner accepted
  "extend notebook 11 rather than add a fourteenth" on that basis.
- **Source data, never invent it.** Town coordinates from OpenStreetMap, bibliography from
  INSPIRE, API limits quoted from the vendor's page, region bounds from Nominatim.
- **Notebooks are educational**, figures shown inline, and **each region notebook carries
  the full `explanation.txt` of its runs inline** — the owner asked for this explicitly.
- **Check for zombie processes periodically.** The owner asked twice. See Trap 6.
- **Credentials:** the owner will hand over an API key when asked and expects it used
  directly. **Never commit it**; document how a reader gets their own.

**Project conventions**

- Measure before optimising, and after. Several confident hypotheses were wrong this
  session, including two of mine that I had already written down.
- **Doctest values must be computed, not predicted** — `tests/test_doctests.py` runs every
  `Examples` block.
- Negative results go in `docs/ROADMAP.md` so they are not retried.
- Lint as CI does: `ruff check .` **from the root** — it lints notebooks too.
- Commit messages explain *why* and state measured deltas.
- The roadmap is updated in the same commit as the code.

---

## 3. Environment and traps

**Trap 1 — memory, and it is the one that bites hardest.** ~8 GiB available of 15, and **a
search that reaches it kills the session, not just the run.** That happened twice. Before
launching anything:

```bash
python tools/run_full_dem.py --region <r> --dry-run
```

- `--max_memory_gb` is `RLIMIT_AS`, capping **virtual** address space, so it counts the
  memory-mapped DEM. `estimate_peak_memory_gb` estimates **anonymous** memory and excludes
  it. They are not comparable.
- **A cap above available memory is not a cap** — the OOM killer arrives first.
  `preflight_memory` now warns, returns `cap_exceeds_available`, and with `refuse=True`
  raises. `run_full_dem.py` passes `refuse=True` and **rejects `--max-memory-gb 0`**.
- **I passed `max_memory_gb=0` myself** — the exact thing my own config comment forbade.
  The advice existed and nothing enforced it. That is why the check is now a mechanism.

**Trap 2 — `conda activate sssearch` fails.** Call the interpreter directly:
`/home/mbustamante/anaconda3/envs/sssearch/bin/python`.

**Trap 3 — `gh` needs `GIT_CONFIG_NOSYSTEM=1`.** And `gh pr create --body-file` cannot read
from the scratchpad; pipe on stdin with `--body-file -`.

**Trap 4 — Jupyter's `python3` kernelspec points at base anaconda**, which has no
`oroscope`:

```bash
mkdir -p /tmp/k/kernels/python3 && cat > /tmp/k/kernels/python3/kernel.json <<'JSON'
{"argv": ["/home/mbustamante/anaconda3/envs/sssearch/bin/python", "-m",
          "ipykernel_launcher", "-f", "{connection_file}"],
 "display_name": "Python 3", "language": "python"}
JSON
cd notebooks && env -u MPLBACKEND JUPYTER_PATH=/tmp/k jupyter nbconvert --execute --inplace 13_turning_the_knobs.ipynb
```

**Trap 5 — `ffmpeg` is the snap build** and cannot write into `/tmp/claude-*`. Write into
the repo's gitignored `output/`.

**Trap 6 — a waiter that greps for a process matches itself.**
`pgrep -f "run_full_dem.py --region ancash"` finds its own command line and spins forever.
**Wait on a file or a log string, never on a process name.** Three zombies accumulated this
way before the owner noticed.

**Trap 7 — `_ = figures.foo()` in a `jupyter-execute` block renders only in the first block
of a page.** End the block with the figure as the last expression. Verify by counting
`_images/` entries in the built HTML — do not trust a green build.

**Trap 8 — escaping `\n` through the notebook generator.** The generator holds cell sources
as ordinary Python strings, so `\n` is consumed when the notebook is built. Write `\\n`, or
avoid it with a bare `print()`.

**Trap 9 — the docs build runs `sphinx -W`.** A short title underline, an undocumented
parameter, or an unreachable intersphinx inventory each fail it. **Build locally before
pushing:**

```bash
env -u MPLBACKEND JUPYTER_PATH=/tmp/k python -m sphinx -W -b html docs/source /tmp/db
```

**Trap 10 — a killed run leaves orphans.** Check `output/<run>/` for `buffer_*.npy` and
mismatched timestamps before trusting a directory.

**Trap 11 — Overpass rate-limits.** A department-sized roads+places fetch takes several
minutes and sometimes stalls. **Fetch context before starting a search** — a run resolves
its map inputs once, at the beginning.

**Disk** was at 99% at one point; now ~14 GB free. `old/` holds 3.5 GB of superseded
material if room is needed.

---

## 4. Repo map

| path | what |
| --- | --- |
| `src/oroscope/site_searcher.py` | Pipeline, CLI, config, screening, morphology, capacity, outputs, map furniture. Adds `estimate_visualisation_memory_gb`, `REFUSE_FRACTION`, and `preflight_memory(refuse=)`. |
| `src/oroscope/combine_experiments.py` | The overlay. Adds `colocation_capacity` and `smallest_radius_for`. |
| `src/oroscope/figures.py` | Physics figures plus three new schematics: `pipeline_stages`, `striding_and_closing`, `score_composition`. |
| `src/oroscope/fetch_dem.py` | Four regions, `--region`, `OPENTOPOGRAPHY_API_KEY`, no longer shells out to a `site_searcher.py` in the CWD. |
| `tools/run_full_dem.py` | Was `run_arequipa_full.py`. Region table: arequipa, ancash, lima, huaylas, cajatambo. |
| `tools/compare_regions.py` | **New.** Regenerates `results/region_comparison.md` from the stores. |
| `tools/make_notebooks.py` | Generates all 13 notebooks. **Edit here, never the `.ipynb`.** |
| `docs/source/{howitworks,glossary,implementation,data}.rst` | **New.** |
| `results/<region>_full/` | Committed small artefacts, including every `explanation.txt`. |
| `results/region_comparison.md` | Every region against every other. |

---

## 5. Current numbers

| | area km² | sites | capacity | sampling |
| --- | --- | --- | --- | --- |
| GRAND, Arequipa | 88,527.5 | 1 | 101,948 | 4 / 5 |
| TAMBO, Arequipa | 111.9 | 26 | 9,024 | 4 / 5 |
| GRAND, Ancash | 43,091.2 | 1 | 49,447 | 4 / 5 |
| TAMBO, Ancash | 174.9 | 35 | 14,290 | 4 / 5 |
| GRAND, Lima | 51,677.6 | 1 | 59,270 | 4 / 5 |
| TAMBO, Lima | 190.9 | 40 | 15,775 | 4 / 5 |
| GRAND, Peru (90 m) | 563,411 | 17 | 633,655 | 4 / 15 |
| **GRAND, Huaylas crop** | **8,294.9** | 1 | **9,609** | **1 / 1** |
| **TAMBO, Huaylas crop** | **855.1** | **109** | **98,696** | **1 / 1** |
| **GRAND, Cajatambo crop** | **5,541.1** | 1 | **6,424** | **1 / 1** |
| **TAMBO, Cajatambo crop** | **1,119.2** | **97** | **129,359** | **1 / 1** |

**Quote the crops, not the departments, for TAMBO.** The department numbers are lower
bounds by a factor between 4.75 and 291. GRAND is unaffected either way. Peru is a survey:
read its area as 4–6 × 10⁵ km², and distrust its site *count* more than its area — the
largest "site" has the whole DEM as its bounding box.

---

## 6. Do not repeat

- **Any of the five searches above.** All stored, all committed.
- **The Huaylas `4 / 5` control.** That is what measured 291×.
- **The stride-1 controls at Colca**, the sensitivity sweeps, the memory-estimator
  investigation, the 90 m slope-band check (answer: no effect), the road and place
  downloads for all four regions.
- **The animation candidate filter.** Six were rejected with reasons in §6.45.
- **Everything in ROADMAP §6 and §7.** Read §6.26a, §6.34, §6.44–6.52 before measuring.

## 7. Still open

1. **`A(E)`**, §0 above.
2. **`tau_exit_probability`'s integration**, §0 and §6.44.
3. **Askaryan emission**, §0.
4. **`min_score` → `score_percentile`**, §6.43 — the owner's call.
5. **External validation** — nothing has been checked against a simulation.
6. **A national TAMBO answer** needs 1 arc-second and tiling; reasons in
   `config/grand_peru_survey.json`.
7. **The joint-realization optimiser**, §6.52. Three steps in order: retain the per-pixel
   score as a raster, then a placement routine given a **patch-aware** feasibility test,
   then a real objective. **Optimise over the union, never the intersection** — the
   measurement says a naive formulation gives a confidently wrong answer.
8. **No release.** Nothing on PyPI.
