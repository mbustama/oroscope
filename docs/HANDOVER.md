# Handover brief — Oroscope

Written for a fresh session with no memory of the last one. **The previous session was
the re-run. It is done: the stores, the notebooks and the documentation now agree with
the audited code.** This one inherits a consistent repository and a physics to-do list.

**Repository:** `mbustama/oroscope`, **public**. Local path `~/Research/GRAND/oroscope`.
A `site_search` symlink sits beside it for anything pointing at the old path.

**Branch:** `dev`, level with `main`. **706 tests**, stdlib `unittest`, on Python
3.9–3.14. **`oroscope` 0.5.0 is on PyPI**, tagged `v0.5.0`, with a published GitHub
Release. **Documentation:** <https://mbustama.github.io/oroscope/>, built from `main`.

`main` is protected: PR required, seven checks, zero approvals. You can open a PR and
report its checks, but **you cannot merge** — the permission classifier blocks
`gh pr merge`. Hand the merge command to the owner. Do not work around it.

---

## 0. State of the store — good, for the first time in a while

All six regions were re-run on the audited code and every derived number was refreshed
from them. `results/` and `docs/` agree. `assumptions.rst` no longer carries a warning
admitting otherwise, because it no longer needs one.

| region | sampling | GRAND km² | TAMBO sites / km² | joint km² | share of TAMBO |
| --- | --- | --- | --- | --- | --- |
| colca | 1 / 5 | 4,569.4 | 16 / 203.0 | 123.3 | 60.7% |
| huaylas | 1 / 1 | 8,249.5 | 32 / 291.3 | 228.7 | **78.5%** |
| cajatambo | 1 / 1 | 5,573.8 | 44 / 774.5 | 591.7 | **76.4%** |
| ancash | 4 / 5 | 42,791.9 | 62 / 740.0 | 411.1 | 55.6% |
| lima | 4 / 5 | 51,209.0 | 84 / 915.4 | 509.8 | 55.7% |
| arequipa | 4 / 5 | 88,208.2 | 85 / 1,036.9 | 619.1 | 59.7% |

**Quote the crops, never the departments, for TAMBO** — the department figures are
striding artefacts. The joint share is **two rows and never one**: ~56–60% strided,
~76–79% unbiased. Roadmap §6.70.

**No ultrareview has ever been run on this work**, through the audits and the release
alike. The owner declined; `/code-review ultra` is user-triggered and billed, and you
cannot launch it.

---

## 1. The immediate job

**There isn't one.** The re-run landed, 0.5.0 is published and installable from PyPI,
the branches are level and nothing is unpushed. What follows is a backlog, not a queue.

Two things are the owner's call and neither blocks anything:

1. **Should `data/` ship in the wheel?** The two curves are hand-digitised from a Nature
   Astronomy figure and an arXiv figure, so redistributing them in a public package is a
   licensing question rather than a packaging one. Today they are repo-only, no runtime
   path breaks, and the docs say so.
2. **`bench/baseline.json` wants refreshing** with `--update --repeat 5` on an idle
   machine. Its *results* are stale and now say so on every run; its *timings* are good
   and would be spoiled by a busy machine.

**Releasing again**, when it comes to that: bump `pyproject.toml` and
`oroscope.__version__` together (a test pins them), give the changelog a `## <version>`
heading (another test pins that), merge, tag the merge commit, and publish a GitHub
Release — the `release: published` event is what triggers the upload. Rehearse first from
the Actions tab, which uploads `0.5.0.devN` to TestPyPI and exercises the OIDC exchange
that no local check reaches. Installing from TestPyPI needs
`--extra-index-url https://pypi.org/simple/`, because TestPyPI does not mirror
dependencies and will otherwise try to build a 2015 numpy from source.

---

## 2. What the last session did

One commit, `8c73445` → `c6bcf03`, written up in `docs/ROADMAP.md` **§6.69–§6.73**.
Read those rather than re-deriving. In brief:

| § | what |
| --- | --- |
| 6.69 | `AREA_INFLATION_AT_COLCA` 2.29 → **2.35**, from the same stride-1 control that produced the original. The test guarding it asserted the literal, not the constant. |
| 6.70 | The six regions. TAMBO moved everywhere, GRAND barely. `solid_angle` is the weakest component at **every** selected site in every region. |
| 6.71 | Both striding penalties collapse at 150 m: Colca **4.75× → 1.51×**, Huaylas **291× → 23.0×**. §6.49's "acceptance identical at 14.0%" was geometry × score under a geometry label. |
| 6.72 | The Arequipa cap. `--max-memory-gb` is `RLIMIT_AS`; the documented 5.68 GiB was RSS. Measured **VmHWM 6.59, VmPeak 7.80**. |
| 6.73 | Two zombie stores and the manifest that hid them; `make_notebooks.py` tracks code, not data. |

**One published conclusion did not survive.** `physics.rst` said the joint "barely moved
with scale" — Colca 50.1 km², the whole DEM 50.2. At 150 m it is 123.3 against 619.1.
The old invariance was an artefact of a closing element too small to reconnect what
striding cut apart. Rewritten, not renumbered.

---

## 3. Owner preferences — follow without being asked

- **Proceed without asking** when the path is clear. Said repeatedly.
- **Do not do trivial work.** When asked for options, filter hard and say what you
  rejected. When asked to assess, recommend with reasons.
- **Measure, do not assert.** Two confident hypotheses were wrong last session — that a
  6.5 GiB cap sized from RSS would hold, and that TAMBO's failed map was GRAND's fault.
  Both were refuted by measurement in minutes. If you claim a number, compute it.
- **Say plainly when a finding is wrong**, including a reviewer's.
- **Source data, never invent it.** If a number needs a run that has not happened, say so.
- Negative results go in `docs/ROADMAP.md` so they are not retried.
- Commit messages explain *why* and quote measured deltas. The roadmap is updated in the
  same commit as the code.
- **Figures:** no titles; legend outside the axes at the top; minimal legend text;
  `ss.attach_colorbar`; scale bar and north arrow on every map; grey base with a
  colorbar; roads green; few labels; capitalise the first word of every label;
  attribution in the caption.
- **Notebooks are educational**, figures inline, each region notebook carrying its runs'
  full `explanation.txt` **and the exact invocations that produced it** — including the
  roads fetch and any control run.
- **Credentials:** the owner hands over an API key when asked and expects it used
  directly. Never commit it.

---

## 4. Environment and traps

**Trap 1 — memory, still the one that bites.** ~8 GiB available of 15. `--dry-run`
before anything. **`--max_memory_gb` is `RLIMIT_AS` and caps *virtual address space*;
`estimate_peak_memory_gb` estimates *anonymous* memory. Never size one from the other.**
Arequipa: estimate 5.08, measured **RSS 6.59, address space 7.80**, so it needs
`--max-memory-gb 8.0`. A cap above available memory is not a cap.

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
string, never a process name. And note the searches **buffer their progress**: a log can
sit unchanged for 17 minutes while the run saturates 8 cores. Judge liveness by CPU time
in `/proc`, not by log growth.

**Trap 7 — `_ = figures.foo()` in a `jupyter-execute` block renders only in the first
block of a page.** End the block with the figure as the last expression.

**Trap 8 — escaping `\n` through the notebook generator.** Write `\\n`, or a bare `print()`.

**Trap 9 — the docs build runs `sphinx -W` with numpydoc validation:**

```bash
env -u MPLBACKEND JUPYTER_PATH=/tmp/k python -m sphinx -W -b html docs/source /tmp/db
```

**Build it LAST, after the final edit.**

**Trap 10 — CI runs the tests with `working-directory: tests`.** Build paths from
`__file__`, as `_support.REPO_ROOT` does.

**Trap 11 — `ruff check .` from the repo root lints the notebooks too.**

**Trap 12 — a killed run leaves orphans.** Check `output/<run>/` for `buffer_*.npy` and
mismatched timestamps. `output/arequipa_full/` is a known dead directory — orphan
buffers only, sparse, nothing reads it.

**Trap 13 — Overpass rate-limits.** Fetch roads and places *before* starting a search.
All six DEMs now have them, Colca included.

**Trap 14 — freeze the code before running anything long.**

**Trap 15 — CI does not execute notebooks 07–11.** 07 needs an ffmpeg the runner lacks;
08–11 need stores that take an hour to produce. So a change can break them and CI stays
green — which is exactly what happened to `product_collapse` when §6.54 added a guard.
**Run them locally whenever the pipeline or the stores change.**

**Trap 16 — `make_notebooks.py` tracks code, not data.** It reports a notebook unchanged
when only its inputs moved. Anything reading a store must be re-executed regardless.

**Disk** was at 98% with ~12 GB free. `old/` holds 3.5 GB of superseded material.

---

### 6.  Never examined, if you want depth

`old/` (3.5 GB, gitignored) · the RFI-zone and settlement coordinates · the Numba kernel
arithmetic, which coverage cannot see by construction · the changelog's older prose. And
four modules the suite barely reaches: `fetch_roads` 17%, `fetch_dem` and `generate_env`
27%, **`sensitivity` 32%** — that last one produced the sweep figures now published in
`assumptions.rst`, which is the same exposure as any trusted output from lightly-tested
code.

**The pattern worth carrying forward.** Almost everything the audits found was true
locally and false for a stranger: a config preset that worked because an obsolete DEM
sits on this machine, a release step that worked where dependencies happened to be
installed, doctests nobody ran, figure captions in code cells that a sweep over prose
could not see. None could fail here; two of them failed in the TestPyPI rehearsal, which
is where you want to find out. When checking something, prefer the question "would this
hold on a machine that has never run this project?"

## 5. Still open

The physics list, unchanged — none of it was touched by the re-run:

1. **`A(E)`** — the outstanding physics ask. `data/` holds two *integral* published
   curves and `aperture.array_scale_factor` corrects them for array size; it cannot
   correct for the site. A real differential table is what is missing.
2. **Askaryan (charge-excess) emission is not modelled**, so `sin α → 0` scores exactly
   zero and a product composition rejects the site outright. Peru is near the magnetic
   equator, so north–south geometries are hit hardest.
3. **Nothing has been checked against an external simulation.** The Earth-absorption
   prediction (window edge −4.4° at 100 PeV to −0.9° at 10 EeV) is the cheapest such test.
4. **`min_score` → `score_percentile`** (§6.43) — the owner's call. The re-run makes the
   case stronger: 0.35 is `score_percentile` **17.8**, the median candidate score is
   **0.13**, and a sweep shows no knee anywhere (6.7 → 525.0 km² across percentiles
   5 → 40).
5. **A national TAMBO answer** needs 1 arc-second and tiling.
6. **The joint-realization optimiser** (§6.52). Optimise over the *union*, never the
   intersection.
7. **~~No release.~~** 0.5.0 is on PyPI. What is left is the backlog in §1.
8. **IGRF declination per site** — inclination follows the DEM's coordinates now;
   declination still falls back to the Arequipa value.
9. **The estimator is left uncalibrated**: 5.08 GiB against a measured 6.59 on Arequipa.
   Deliberate — re-fitting `n_scoring_arrays` to chase one region is how the pre-flight
   came to be sized against the cheaper of two configs (§6.55). Fix the second column,
   not the constant.
10. **The talk.** The owner is writing an abstract on co-locating a particle and a radio
    detector. The co-location result is the structural one: *a pixel has one slope*, so
    GRAND (3–25°) and TAMBO (20–60°) compete at the screening step, while the viewing
    windows (±3° against ±20°) are in no conflict at all. Joint share is **76–79%** of
    TAMBO's mask from the unbiased crops — now measured on two of them, and a range that
    depends on array design rather than a constant.
