# Handover brief — Oroscope

Written for a fresh session with no memory of the last one. **The previous session was
an audit; this one's job is to make the stored results agree with the audited code, and
then publish.**

**Repository:** `mbustama/oroscope`, **public**. Local path `~/Research/GRAND/oroscope`.
A `site_search` symlink sits beside it for anything pointing at the old path.

**Branch:** `dev`, head **`ecdd92d`**, pushed, **CI all green** — 8 checks: `docs`,
`ruff`, `Notebooks execute`, tests on Python 3.9–3.13. **677 tests**, stdlib `unittest`,
~30 s. `main` is **119 commits behind** and there is no open PR.
**Documentation:** <https://mbustama.github.io/oroscope/>, built from `main`.

`main` is protected: PR required, seven checks, zero approvals. You can open a PR and
report its checks, but **you cannot merge** — the permission classifier blocks
`gh pr merge`. Hand the merge command to the owner. Do not work around it.

---

## 0. Start here — the one thing blocking everything else

**`results/` reflects none of the last eight commits.** The audit changed what the
pipeline computes; no search has been run since. Every stored number, every notebook
output, and 39 figures in the docs are from superseded code. Nothing is broken — the
store is simply *old*, and deliberately so: the previous session stopped rather than
ship a store it could not vouch for.

The whole to-do list is that, in order. Do not reorder it; each step consumes the one
before.

### 1. Run the searches

```bash
python tools/run_full_dem.py --region <r> --dry-run     # ALWAYS first
```

Six regions: `colca`, `huaylas`, `cajatambo`, `ancash`, `lima`, `arequipa`. Smallest
first. Both experiments and the combine per region — no `--only`. About **1.5 hours**,
almost all of it GRAND Arequipa (~25 min) and Lima (~20 min).

Expect TAMBO's numbers to move a lot and GRAND's a little. Two independent reasons, both
measured and both in the roadmap: TAMBO moved to the published 150 m spacing (§6.62), and
`gap_close_km` defaults to the detector spacing, so that also widened the closing
element. **Quote the crops, never the departments, for TAMBO** — the department figures
are striding artefacts (§6.49, §6.62).

### 2. Regenerate what derives from the runs

- `python tools/compare_regions.py` → `results/region_comparison.md`
- `explain.AREA_INFLATION_AT_COLCA` — a **module constant** (2.29) quoted in every run's
  summary, measured from a 100 m Colca run. Recompute from the fresh Colca store, do not
  hand-edit.
- `figures.pipeline_stages` — hardcodes the Ancash TAMBO funnel. It currently merges
  "Arrival scan" and "Scoring" into one bar because the stored run predates §6.53 and
  cannot separate them. A post-fix run can, so **split it back into seven stages**.

### 3. Re-execute the notebooks

`python tools/make_notebooks.py` rewrites only what changed, then re-execute those.
Notebooks 08–11 read the stores; 12 reads the Peru survey (untouched); 06, 07, 13 read
no results. See Trap 4 for the kernelspec.

### 4. Refresh the 39 documentation numbers

Inventoried by group in **roadmap §6.68**. `assumptions.rst` carries a warning at the top
admitting its figures predate the code — **delete that warning** once they are refreshed.

### 5. Then, and only then

An **ultrareview** is worth running once code and data finally agree — it is not worth
running now, because it would spend its depth re-finding drift already catalogued. Then
open the PR to `main` and hand over the merge.

---

## 1. What the last session did

Fifteen commits, `994fa62` → `ecdd92d`. **Every change is written up in `docs/ROADMAP.md`
§6.53–§6.68 with its measurement — read those rather than re-deriving.** In brief:

| § | what |
| --- | --- |
| 6.53 | The funnel now separates geometry from the score cut. `--explain` had been blaming the arrival window for cuts `min_score` made. |
| 6.54 | A mistyped `--score_weights` name was accepted and silently dropped. |
| 6.55 | The memory pre-flight was sized against the cheaper of two configs. |
| 6.56 | `tau_exit_probability` integrates in `u = X − x`; the old grid was 8× high at 3 PeV. |
| 6.57 | The combination's memory is modelled and 0.65 GiB cheaper. It had no test coverage at all. |
| 6.58 | Numba kernels are measurable via `NUMBA_DISABLE_JIT=1`. |
| 6.59–6.61 | Sweep timeout, a stale funnel reader, an all-sites area beside a selected count. |
| 6.62 | TAMBO at the published 150 m; `aperture.array_scale_factor` for the published curves. |
| **6.63, 6.67** | **The reported solid angle was the whole circle whatever the fan.** Then made scale-free: the score is now the *fraction* of available sky, so no constant needs retuning when the fan or arrival window moves. |
| 6.64 | The morphology tiling halo was half what closing needs. |
| 6.65 | Six smaller ones, including Colca finally being in the region table. |
| 6.66 | What a max-effort code review found **in the audit's own work** — including a 44× under-report and a wrong acceptance column already committed. |
| 6.68 | The documentation pass, and the inventory of what still waits on the re-run. |

The README was cut 590 → 393 lines: its "Complete List of Options" had drifted ten flags
behind `cli.rst`, which is generated against the parser and tested.

---

## 2. Owner preferences — follow without being asked

- **Proceed without asking** when the path is clear. Said repeatedly.
- **Do not do trivial work.** When asked for options, filter hard and say what you
  rejected. When asked to assess, recommend with reasons.
- **Measure, do not assert.** Several confident hypotheses were wrong last session,
  including three of mine that I had already written down. If you claim a number, compute
  it in front of the owner.
- **Say plainly when a finding is wrong**, including a reviewer's — do not change code to
  satisfy a claim you have refuted.
- **Source data, never invent it.** If a number needs a run that has not happened, say so
  rather than filling it in.
- Negative results go in `docs/ROADMAP.md` so they are not retried.
- Commit messages explain *why* and quote measured deltas. The roadmap is updated in the
  same commit as the code.
- **Figures:** no titles; legend outside the axes at the top; minimal legend text;
  `ss.attach_colorbar` so the bar matches the panel; scale bar and north arrow on every
  map; grey base with a colorbar; roads green; few labels; capitalise the first word of
  every label; attribution in the caption.
- **Notebooks are educational**, figures inline, each region notebook carrying its runs'
  full `explanation.txt`.
- **Credentials:** the owner hands over an API key when asked and expects it used
  directly. Never commit it.

---

## 3. Environment and traps

**Trap 1 — memory, and it still bites hardest.** ~8 GiB available of 15. A search that
reaches it kills the *session*, not just the run. `--dry-run` before anything. `--max_memory_gb`
is `RLIMIT_AS` and caps *virtual* address space; `estimate_peak_memory_gb` estimates
*anonymous* memory. They are not comparable. **A cap above available memory is not a cap.**

**Trap 2 — `conda activate sssearch` fails.** Call the interpreter directly:
`/home/mbustamante/anaconda3/envs/sssearch/bin/python`.

**Trap 3 — `gh` needs `GIT_CONFIG_NOSYSTEM=1`.** And `gh pr create --body-file` cannot
read the scratchpad; pipe on stdin with `--body-file -`.

**Trap 4 — Jupyter's `python3` kernelspec points at base anaconda**, which has no
`oroscope`:

```bash
mkdir -p /tmp/k/kernels/python3 && cat > /tmp/k/kernels/python3/kernel.json <<'JSON'
{"argv": ["/home/mbustamante/anaconda3/envs/sssearch/bin/python", "-m",
          "ipykernel_launcher", "-f", "{connection_file}"],
 "display_name": "Python 3", "language": "python"}
JSON
cd notebooks && env -u MPLBACKEND JUPYTER_PATH=/tmp/k jupyter nbconvert --execute --inplace 09_arequipa_dem.ipynb
```

**Trap 5 — `ffmpeg` is the snap build** and cannot write into `/tmp/claude-*`. Write into
the repo's gitignored `output/`.

**Trap 6 — a waiter that greps for a process matches itself.** Wait on a file or a log
string, never a process name.

**Trap 7 — `_ = figures.foo()` in a `jupyter-execute` block renders only in the first
block of a page.** End the block with the figure as the last expression.

**Trap 8 — escaping `\n` through the notebook generator.** Write `\\n`, or use a bare
`print()`.

**Trap 9 — the docs build runs `sphinx -W` with numpydoc validation.** An undocumented
parameter, a namedtuple field with no entry, a duplicated `Examples` heading, or
Parameters documented out of signature order each fail it:

```bash
env -u MPLBACKEND JUPYTER_PATH=/tmp/k python -m sphinx -W -b html docs/source /tmp/db
```

**Build it LAST, after the final edit.** Two commits went red last session because the
docs were built before the last change rather than after.

**Trap 10 — CI runs the tests with `working-directory: tests`.** A doctest or test that
opens a relative path passes locally and fails in CI. Build paths from `__file__`, as
`_support.REPO_ROOT` does.

**Trap 11 — `ruff check .` from the repo root lints the notebooks too.**

**Trap 12 — a killed run leaves orphans.** Check `output/<run>/` for `buffer_*.npy` and
mismatched timestamps before trusting a directory.

**Trap 13 — Overpass rate-limits.** Fetch roads and places *before* starting a search.

**Trap 14 — freeze the code before running anything long.** The previous session
invalidated an in-flight rerun three times by finding another defect while it ran. Finish
the fixes, get the suite green, *then* run.

**Disk** was at 98% with ~13 GB free. `old/` holds 3.5 GB of superseded material.

---

## 4. Still open

Unchanged from before, none of them touched last session:

1. **`A(E)`** — the outstanding physics ask. `data/` holds two *integral* published
   curves and `aperture.array_scale_factor` corrects them for array size; it cannot
   correct for the site, and `assumptions.rst` says exactly why. A real differential
   table is what is missing.
2. **Askaryan (charge-excess) emission is not modelled**, so `sin α → 0` scores exactly
   zero and a product composition rejects the site outright. Peru is near the magnetic
   equator, so north–south geometries are hit hardest.
3. **Nothing has been checked against an external simulation.** The Earth-absorption
   prediction (window edge −4.4° at 100 PeV to −0.9° at 10 EeV) is the cheapest such test.
4. **`min_score` → `score_percentile`** (§6.43) — the owner's call.
5. **A national TAMBO answer** needs 1 arc-second and tiling.
6. **The joint-realization optimiser** (§6.52). Optimise over the *union*, never the
   intersection.
7. **No release.** Nothing on PyPI.

New, from last session:

8. **§6.49 (the 291× striding penalty) and §6.51 (three departments) want re-measuring at
   150 m.** Both were measured at 100 m.
9. **The talk.** The owner is writing an abstract on co-locating a particle and a radio
   detector. The co-location result is the structural one: *a pixel has one slope*, so
   GRAND (3–25°) and TAMBO (20–60°) compete at the screening step, and the viewing
   windows are in no conflict at all. Joint share is **~73–81%** of TAMBO's mask from the
   unbiased crops — a range that depends on array design, not a constant.
