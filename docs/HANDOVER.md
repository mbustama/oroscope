# Handover brief — Oroscope

Written to be fed to a fresh session. It assumes no memory of the previous one.

**Repository:** `mbustama/oroscope`, **public**. Local path `~/Research/GRAND/oroscope`.
A `site_search` symlink sits beside it for anything pointing at the old path; the owner
knows about it and has chosen to keep it.

**Branch:** `dev`, head **`d333841`**, pushed. `main` is 14 commits behind — the
last merge was PR #4. **Tests:** 600, stdlib `unittest`, ~30 s. **CI:** 8 jobs.
**Documentation:** <https://mbustama.github.io/oroscope/>, deployed from `main`.

`main` is protected: no direct push, PR required, seven status checks, **zero
approvals**. You can open a PR and wait for checks, but **you cannot merge** — the
permission classifier blocks `gh pr merge`. Open the PR, report the green checks, and
hand the merge command to the owner. Do not try to work around it.

---

## 0. Start here

The previous session's four asks are all **delivered**. Nothing below is half-finished.
What follows is what is genuinely next, in the order I would take it.

### 1. Open the PR

`dev` is pushed and there is no open PR for it. Seven commits are waiting to land.

```bash
GIT_CONFIG_NOSYSTEM=1 gh pr create --base main --head dev --title "..." --body-file -
```

Then report the green checks and hand the merge command to the owner — you cannot merge
it yourself (see the note above). Watch for the docs-job flake in Trap 10 before
concluding anything is broken.

### 2. `tau_exit_probability` under-resolves its own integral — decide whether to fix it

Found this session, measured, recorded in **§6.44**, and deliberately **not fixed**. It
integrates over interaction depth with `np.linspace(0, X, samples)` — a fixed number of
points over the *whole* depth — while only interactions within about one tau range of
the far surface contribute. As `X` grows the spacing outruns the only region that
matters. At 3 PeV and `X` = 10⁹ g/cm² the default is **8× the converged value and
inverts the trend with depth**, so the curve grows a spurious maximum at the grid edge.

`depth_band_from_energy` inherits it and returns (5.6×10⁷, 2.9×10⁸) for TAMBO's
3 PeV – 1 EeV range, a band excluding the true 1 EeV optimum of 5.7×10⁶ by more than an
order of magnitude, whose low edge *rises* when the minimum energy is *lowered*.

**No published number is affected** — every config leaves `depth_band_gcm2` null, so runs
score against the default (10⁵, 10⁷) and never call it. That is why it was left alone:
the fix changes a physics function's outputs and wants its own commit with a test pinning
the converged values. The substitution `u = X − x` on a log grid is the obvious form.
**Owner's call.**

### 3. The Peru survey has an open follow-up that was started and abandoned

§6.46 reports the survey at `candidate_stride` 15 with `gap_close_km` 1.5. That closing
element is 1.5× GRAND's own 1 km antenna spacing, chosen to outrun the 1,382 m gap
stride 15 leaves — so the area carries a **declared upward bias**, bracketed at
4–6 × 10⁵ km².

The run that removes it is `candidate_stride` 10 with `gap_close_km` 1.0: the element
then equals the array's own spacing *and* still outruns the 922 m gap. **I started that
run and it killed the machine.** Estimate 6.74 GiB against ~8.7 available. Read §6.46's
second half before retrying, and see Trap 1 below. If you retry it, do it on a machine
with more headroom, or accept the stride-15 bracket and move on — the bracket is already
honest.

### 4. `A(E)` is still the outstanding physics ask

Unchanged. The selector exists (`--decay_weight_by` takes `flux`, `acceptance`,
`flux_times_acceptance`); no real differential table does. One *inferred* from a
published integral curve is demonstrably unsafe — §6.42, it returns zero sites.

### 5. If a layout tool is ever written, read §6.47 first

The joint-realization question was answered this session and the answer is not the one
the obvious approach assumes. Summary in §5 below; the whole thing is §6.47.

---

## 1. What the previous session did

Six commits, `35d7814` → `b70ef72`.

**Four more animations**, taking `tools/make_animations.py` to eight (§6.45). The filter
was whether the intermediate states carry the argument; six candidates were rejected
because a static figure does them better, and those reasons are in §6.45 so they are not
re-proposed.

- `the_azimuth_fan` — completes the pair with `the_walk`. That one sweeps elevation at
  one bearing; this sweeps the bearing.
- `product_collapse` — six real components multiplied into a real cut, one at a time.
- `slope_criterion` — `min_target_slope_deg` crossing the wall-slope distribution over
  Colca. The sweeps say how much; this says where.
- `tau_in_rock` — energy and survival falling against the production-and-escape optimum.

Three of the eight now read `input/dem/colca.tif` when present and fall back to synthetic
terrain — **saying which on the figure** — when not, because "where a criterion bites on
real ground" is not something synthetic terrain can honestly show.

**Two measurements fell out of building them,** both in §6.45:

- **The product collapse is real but not evenly shared.** 100% of viable candidates above
  `min_score` 0.35 before any component, 32.2% after six — but `solid_angle` alone takes
  it from 100% to 35.9%, and `distance` moves it by *nothing*, the scan having already
  applied that same 2–5 km window as a hard criterion before scoring saw it.
- **The wall-slope mask outlives its own median by 20°.** Half the candidates see a mean
  wall slope under 29.7°, yet a 30° floor keeps 87% of them; the half-way point is 50°.
  The criterion is per direction, the observable is a mean over accepted directions.
  **Read `target_slope_deg` as a description, never as a prediction of what a cut does.**

**Notebook 9**, generated from `tools/make_notebooks.py` as usual. Builds all eight,
explains what each argues, and documents MP4 → GIF two ways: `--format gif`, and a plain
`ffmpeg` `palettegen`/`paletteuse` recipe for an MP4 from anywhere. The palette matters —
342 KiB against pillow's 761 KiB on the same animation, and better colour.

New in the tool: `write_mp4_with_stills()`, which writes the MP4 and grabs stills in the
same pass. **One function and not two because the builders accumulate** — the ray drawn
at frame 30 is still on the axes at frame 60, which is what makes the fan fill in — so
the frames can be walked exactly once.

**The Peru-wide survey ran.** §6.46. 22,080 × 15,360 = 339 Mpx at 3 arc-seconds:
**17 sites, 563,411 km², 633,655 antenna positions, in four minutes.** Caveats in §5.

**Notebook 10** is that survey on its own, at the owner's request. It reads the stored
run rather than repeating it, but everything else computes live — the memory table, the
stride/closing table, the funnel, the area bracket, and the resolution check. **Unlike
7, 8 and 9 it is executed in CI**: every cell needing the store is guarded and nothing
else touches the filesystem outside the library, so on a bare runner it degrades to
prose and arithmetic instead of failing. Verified by running it from an empty directory
with no repo around it — the thing that would have broken there was a synthetic fallback
importing from `../tests`, so it is built inline instead. Keep that property if you edit
it.

**The joint-realization question was answered.** §6.47, and §5 below.

**A memory safeguard**, after the machine went down. §6.46 and Trap 1.

---

## 2. Owner preferences — follow these without being asked

Learned by correction across two sessions. Consistent and worth honouring.

**Figures**

- **No titles.** The caption carries it. This applies to every map.
- **Legend outside the axes, at the top.** Four columns unless told otherwise — but
  build it from the categories that actually occur; a legend naming a colour that never
  appears sends the reader hunting for nothing.
- **Legend text minimal**: no counts, no parentheticals. "Roads", not "Roads (230, OSM)".
- **Colorbar height must match the plot panel.** Use `ss.attach_colorbar`.
  `fig.colorbar(fraction=...)` sizes against the *figure* and overshoots.
- **Scale bar and north arrow** on every map. `ss.add_scale_bar`, `ss.add_north_arrow`.
- **Grey base with a colorbar** beats a colourful base.
- **Roads green**, not neutral dark: a thin dark line vanishes into hillshade exactly
  where the ground is steep, which is where the sites are.
- **Only a few labels.** 6 labels, 25 markers, arrived at by being told twice.
- **Attribution goes in the caption, not on the figure.**
- **Colour scales track the data** (`ss.altitude_limits`).
- **Progressive-reveal frames** are wanted for talks, pixel-identical between frames.
  Verify it, do not assume it.

*Known inconsistency:* the pipeline's own search map still prints counts and
parentheticals in its legend ("Site 5: 600703 DUs (533861.48 km²)"), against the rule
above. Not changed — it was not in scope — but worth raising.

**Working style**

- **Proceed without asking** when the path is clear. Said explicitly, more than twice.
  Reserve questions for genuine forks.
- **Do not do trivial work.** When asked for options, filter hard and say what you
  rejected and why.
- **Source data, never invent it.** Town coordinates from OpenStreetMap, bibliography
  verbatim from INSPIRE, API limits quoted from the vendor's own page. This matters.
- **Notebooks are educational.** Figures shown inline, not merely saved.
- They will feed a brief to a fresh session rather than let context run out — write
  handovers accordingly.
- **Credentials:** the owner will hand over an API key when asked and expects it used
  directly. **Never commit it**; document how a reader gets their own instead.

**Project conventions (from the repo, unchanged)**

- Measure before optimising, and after. Several confident hypotheses have been wrong —
  including one this session (see §5, the 90 m slope check).
- Run examples, don't read them — `tests/test_doctests.py` executes every `Examples`
  block. **Doctest values must be computed, not predicted.**
- Negative results go in `docs/ROADMAP.md` so they are not retried.
- Lint as CI does: `ruff check .` **from the root** — it lints the notebooks too.
- Commit messages explain *why* and state measured deltas.
- The roadmap is updated in the same commit as the code.
- Figure labels capitalise their first word.

---

## 3. Environment and traps

**Trap 1 — memory, and it is the one that bites hardest.** This desktop has ~8 GiB of 15
actually available, and **a search that reaches it kills the session, not just the run.**
That happened this session. Before launching any search:

```python
ss.preflight_memory(dem, downsample_factor=…, candidate_stride=…, max_memory_gb=…)
```

Read both numbers. **Do not launch when the estimate is above ~70% of available, or when
the cap is above available — ask first.** Two things to understand:

- `--max_memory_gb` is `RLIMIT_AS`, which caps **virtual** address space and so counts
  every mapping — the DEM's `.npy` cache, the ping-pong buffers. `estimate_peak_memory_gb`
  estimates **anonymous** memory and deliberately excludes them. On a 339 Mpx DEM the two
  differ by over 2 GiB.
- **A cap above available memory is not a cap.** The OOM killer arrives before `RLIMIT_AS`
  fires. The pre-flight now warns and returns `cap_exceeds_available`. When the two
  constraints cannot both be met, the configuration does not fit: raise
  `candidate_stride`, which is the memory lever (§6.26a), not the cap. Never pass 0.

**Trap 2 — `conda activate sssearch` fails** (`conda init` not run). Call the interpreter
directly: `/home/mbustamante/anaconda3/envs/sssearch/bin/python`. Same environment.

**Trap 3 — `gh` needs `GIT_CONFIG_NOSYSTEM=1`** — the sandbox blocks `/etc/gitconfig` and
every `gh` call fails without it. `gh pr create --body-file` cannot read from the
scratchpad; pipe on stdin with `--body-file -`.

**Trap 4 — Jupyter's `python3` kernelspec points at base anaconda**, which has no
`oroscope`, so `nbconvert --execute` and the docs' `jupyter_sphinx` blocks both fail with
`ModuleNotFoundError`:

```bash
mkdir -p /tmp/k/kernels/python3 && cat > /tmp/k/kernels/python3/kernel.json <<'JSON'
{"argv": ["/home/mbustamante/anaconda3/envs/sssearch/bin/python", "-m",
          "ipykernel_launcher", "-f", "{connection_file}"],
 "display_name": "Python 3", "language": "python"}
JSON
cd notebooks && env -u MPLBACKEND JUPYTER_PATH=/tmp/k jupyter nbconvert --execute --inplace 09_animating_the_mechanism.ipynb
```

**Trap 5 — `ffmpeg` is the snap build** and cannot write outside its confinement — it will
not write into `/tmp/claude-*`. Write into the repo's gitignored `output/`, or extract
frames with PIL.

**Trap 6 — `pgrep -f "script.py"` matches its own command line**, and so does an `echo` of
the same string in the same command. Use `ps` and check the output.

**Trap 7 — matplotlib backend.** A library must not choose it; CI asserts that importing
`oroscope` leaves it untouched. When touching a figure path, **check images are actually
produced** — count `image/png` outputs, or grep built HTML for `Figure size`. Do not trust
a green build.

**Trap 8 — notebook size.** `show_figure` picks JPEG over PNG above 250 KiB; terrain maps
are photographic and PNG made notebook 8 4.2 MB. Notebook 9's stills use the same trick.

**Trap 9 — a killed run leaves orphans.** The crashed search left two 339 MB
`buffer_*.npy` scratch files and an output directory mixing two runs' artefacts. Check
`output/<run>/` for `buffer_*.npy` and for mismatched timestamps before trusting it.

**Trap 10 — a red CI docs job may be a network flake, not a code bug.** Sphinx runs with
`-W`, and intersphinx fetches an inventory per entry in `intersphinx_mapping` over the
network. One unreachable host emits one warning and exits 1, with the real cause buried a
hundred lines above `build finished with problems, 1 warning`. **Before investigating,
grep the log for `failed to reach any of the inventories`.**

This bit twice on 2026-08-16, both times on `docs.scipy.org`. Fixed by removing the scipy
entry, which nothing referenced — no role in `docs/source` or `src/` resolved against it,
and the only mention of the name was a plain-text dependency row. `numpy` and `python`
stay, because docstring type fields do resolve against them; if either flakes the same
way, cache the inventory rather than removing it.

**Disk** was at 99% at one point this session and is now ~14 GB free. `old/` holds 3.5 GB
of superseded material if room is ever needed.

---

## 4. Repo map (what changed this session)

| path | what |
| --- | --- |
| `tools/make_animations.py` | Eight animations now. Adds `_colca_ground()` (real DEM or synthetic fallback), `write_mp4_with_stills()`, and the four new builders. |
| `tools/make_notebooks.py` | Adds `NB09`. **Edit here, never the `.ipynb`.** |
| `notebooks/09_animating_the_mechanism.ipynb` | Generated and executed. Excluded from CI execution. |
| `notebooks/10_the_peru_survey.ipynb` | Generated and executed. **Not** excluded — it is written to survive a runner with no store. |
| `src/oroscope/fetch_dem.py` | Adds the `peru` region (SRTMGL3), `--region`, `--output_dir`, `--config_dir`, `OPENTOPOGRAPHY_API_KEY`. No longer shells out to a `site_searcher.py` in the CWD — a pre-package leftover that broke it outside `src/`. |
| `src/oroscope/site_searcher.py` | `preflight_memory` warns when the cap exceeds available memory and returns `cap_exceeds_available`. |
| `config/grand_peru_survey.json` | **New.** The national survey. Every choice is commented with its reason. |
| `input/dem/peru_SRTMGL3.tif` | Gitignored, 302 MB, downloaded this session. |
| `tests/test_docs.py` | Checks notebook 9's `ma.<name>` calls resolve and that it builds every `BUILDERS` entry and no others. |
| `docs/ROADMAP.md` | §6.44–6.47. |

---

## 5. Current numbers

| | area km² | sites | capacity |
| --- | --- | --- | --- |
| GRAND, Colca crop | 4,580.2 | 1 | 5,317 |
| TAMBO, Colca crop | 83.6 → **read as ~397** | 15 → 29 | 9,717 → **~45,856** |
| GRAND, full Arequipa DEM | 88,527.5 | 1 | 101,948 |
| TAMBO, full Arequipa DEM | 111.9 (low) | 26 | 9,024 |
| joint, full Arequipa DEM | 50.2 | | Jaccard 0.0006 |
| **GRAND, all of Peru (90 m)** | **563,411** | **17** | **633,655** |

**Quote these with their caveats.** TAMBO's areas are lower bounds by ~4.75× from striding
and ~30% again from downsampling. The joint is limited by TAMBO's mask, so it is a floor.

**Peru specifically.** Read the area as **4–6 × 10⁵ km²**: the stride-corrected accepted
set is 407,805 km² and closing takes it to 563,411, a factor 1.38. And **"17 sites" is the
number to distrust, not the area** — the largest site's bounding box is the entire DEM,
because a 1.5 km element applied to a strided scatter over 339 Mpx merges the whole
cordillera into one component. Its accepted candidates are Andean (mean altitude 2,446 m);
the polygon enclosing them is not.

**One prediction checked and wrong**, which is the useful part: 3 arc-seconds does *not*
move the slope screen. On 20 Mpx of Arequipa the 3–25° band holds 67.6% of the map at 30 m
and 67.4% at 90 m — what is lost at the ceiling is regained at the floor. So the survey's
screening is a fact about Peru, not about the grid.

**The joint-realization answer (§6.47).** Of the 50.1 km² joint Colca mask, 92.4% is
TAMBO-band wall and only 3.63 km² is ground a GRAND antenna could stand on — in **1,702
fragments whose largest is 0.038 km²**, against the 0.866 km² one 1 km lattice cell needs.
**Not one fragment holds one antenna.** So an optimiser pointed at the intersection reports
the realization impossible. It is not: 100 TAMBO units need 0.87 km² against 49.4 km²
available, and GRAND-deployable ground sits a **median 0.92 km away** — inside a single
GRAND cell. **Optimise over the union with a per-role band constraint; the coupling is
shared line of sight, not shared footprint.**

---

## 6. Do not repeat

- **The full Arequipa DEM run.** Done, 26 minutes, store populated, notebook 8 executed.
- **The Peru survey at stride 15.** Done, four minutes, `output/grand_peru_survey/`.
- **The stride-1 controls.** Both GRAND and TAMBO. §6.34.
- **The sensitivity sweeps.** §6.20–6.21.
- **The 90 m slope-band check.** Done, §6.46 — the answer was "no effect".
- **The memory estimator investigation.** §6.26a.
- **The road and place downloads** for Arequipa — in `input/roads/`.
- **The animation candidate filter.** Six were rejected with reasons in §6.45.
- **Everything in ROADMAP §6 and §7.** Read §6.26, §6.26a, §6.34–6.47 before measuring
  anything.

## 7. Still open

1. **`A(E)`.** §0.4 above.
2. **`tau_exit_probability`'s integration.** §0.2 above, §6.44.
3. **Askaryan (charge-excess) emission is not modelled.** Only the geomagnetic term is
   implemented, so `sin α → 0` scores exactly zero, and under a *product* composition that
   zero rejects the site outright when charge-excess would leave it ~10–20% efficient.
   Peru is near the magnetic equator, so it is precisely north–south geometries that get
   zeroed. With a 0.14 floor the published 3.7× east-over-north ratio compresses to ~3.4×
   and nothing lands on zero. Lesser omissions alongside it: tau decay branching (~17% to
   muons, a flat normalisation), |B| magnitude (only sin α is used), galactic background.
4. **`min_score` → `score_percentile`.** 0.35 ≡ percentile 22.8 on Colca, and a scan shows
   **no knee**. Switching restates published numbers: the owner's call. §6.43 is the table
   to decide against, and §6.45's finding that one component does nine tenths of the
   collapse sharpens the case.
5. **Nothing has been checked against an external simulation.** The Earth-absorption
   prediction (window edge −4.4° at 100 PeV to −0.9° at 10 EeV) is the cheapest such test,
   is ready to run, and now has an animation arguing it (`energy_window`).
6. **A national TAMBO answer** needs 1 arc-second and tiling. Reasons it cannot be done at
   90 m are in `config/grand_peru_survey.json` under `_comment_no_tambo`.
7. **No release.** Nothing on PyPI.
