# Site Search — Development Roadmap

Working document. Phases 0 and 1 are scoped for implementation; phases 2–4 are
sketched and will be scoped in detail before each is started.

**Status:** planning. Nothing in phases 0–4 is implemented yet. The `dev` branch
currently contains only the georeferencing / core-scaling fixes of commit `a9843f9`.

---

## 1. Where we are

The tool searches Peruvian topography for GRAND deployment sites: it screens terrain
by slope, altitude, aspect and exclusion zones, casts one ray per candidate pixel to
find a target mountain, prunes the result morphologically, and packs an antenna grid
into whatever survives.

Commit `a9843f9` fixed the georeferencing that everything else rests on (map
resolution read from the DEM, per-axis metric pixel sizes, `--num_cores` honoured).
That work was corrective, not scientific — the criteria themselves are unchanged and
are the subject of phase 1.

### 1.1 Decisions taken (2026-08-03)

| Question | Decision |
| --- | --- |
| Energy ranges | GRAND > 100 PeV; TAMBO 1–100 PeV. Both must stay configurable. |
| "Interaction depth" | Use **column depth**, not a height proxy. |
| Output form | **Scores**, not binary masks. |
| Effective areas | Published curves exist — GRAND ref. [1] Fig. 25, TAMBO ref. [2] Fig. 3. |
| TAMBO detector separation | Start at **100 m**. (Note: ref. [2] nominal is 150 m — see §4.3.) |
| Backward compatibility | **Clean break.** Do what is best. |
| Differential acceptance table | Not available at present — see §4.10 for how phase 1 proceeds without it. |

### 1.2 Arrival-direction acceptance (guidance, 2026-08-03)

| experiment / channel | accepted arrival directions |
| --- | --- |
| GRAND, neutrinos | about **-3° to +3°** relative to the horizon |
| GRAND, cosmic rays | **above** the horizon, unless a nearby mountain blocks |
| TAMBO | must be **facing a canyon** |

Other experiments will impose different demands, so the acceptance window is a
per-experiment, per-*channel* configuration rather than a global constant. Note that
GRAND alone needs two channels with different — in fact nearly disjoint — windows,
which is a useful early stress test of the criterion framework in phase 2.

---

## 2. Physics and mathematics review

Findings are ordered by impact. Two were quantified against real Arequipa terrain
rather than inferred from the source; those measurements are reproducible from the
harness built in phase 0.

### 2.1 Measured

**Slope depends strongly on the measurement baseline, which is currently undeclared.**
Measured on a 2500 × 2500 crop of `arequipa_SRTMGL1.tif`:

| baseline | median slope | fraction passing 3–25° |
| --- | --- | --- |
| native (~61 m central difference — what the code uses) | 17.8° | 60.4% |
| ~250 m | 15.6° | 65.2% |
| ~500 m | 13.5° | 71.1% |
| ~1 km | 10.8° | 77.6% |

The antenna spacing is 1 km, so the physically meaningful slope is the one at
array-footprint scale, not at 61 m. The native baseline also absorbs DEM noise:
SRTM relative error of 2–5 m over a 61 m baseline is 2–5° of *apparent* slope,
which is the same size as the 3° minimum-slope threshold.

**The single-azimuth ray is the dominant restriction.** Over 40,000 sampled
candidates (min_dist 10 km, max_dist 40 km):

| ray directions tested | candidates accepted |
| --- | --- |
| aspect only (what the code does) | 13.4% |
| any azimuth within ±60° of aspect | 29.7% |

The single ray rejects **2.2× more sites than it accepts**. Aspect derived from a
noisy gradient is itself uncertain by tens of degrees on rough terrain, so the one
direction tested is partly arbitrary.

**Reported area is ~18× the area that passes the physics.** From the funnel on the
2500² Arequipa crop: 755,328 candidates → 104,632 ray-tracing hits (13.9%) →
2,323,104 pixels after the 1 km morphological closing (a 22× gain) → 1,922,720
pixels in the selected site. Closing, not the physics, determines most of the
reported site extent.

**Capacity is over-counted, severely at fine spacings.** `count_grid_capacity`
derives its strides with three independent `int()` truncations (`spacing_r`,
`spacing_c`, and the hex row step `v_step`), each of which shrinks the grid, so
antennas end up closer than the requested ground spacing:

| requested spacing | actual N-S | actual E-W | over-count |
| --- | --- | --- | --- |
| 1000 m (GRAND) | 829 m | 981 m | **+7.4%** |
| 200 m | 154 m | 178 m | +26.4% |
| 150 m (TAMBO, published) | 92 m | 149 m | **+42.3%** |
| 100 m (TAMBO, starting value) | 61 m | 89 m | **+58.1%** |

The error grows as the spacing approaches the pixel size. At TAMBO spacings there are
only ~3 pixels per detector spacing, so a 30 m DEM cannot represent the layout by
integer stamping at all. This is a phase 2 blocker, not a rounding nicety — see §5.1.

### 2.2 Correctness

**A. The line-of-sight test never tests line of sight.** `check_physics_chunk` walks
outward and accepts the first sample exceeding `detector + 1000 m + fresnel_buffer`.
It never checks whether nearer terrain occludes that sample. `fresnel_buffer` —
documented as clearance against intermediate terrain — is applied as a constant
offset on the target's height, not as a clearance test against anything.

The default value is a well-chosen number applied the wrong way: the first Fresnel
zone radius over a 30 km path is

```
r₁ = √(λ d₁ d₂ / (d₁ + d₂))  →  212 m at 50 MHz, 106 m at 200 MHz
```

so 200 m is the right order of magnitude for a clearance margin that is never computed.

**B. The "1 km interaction depth" is a height proxy for a column-depth requirement.**
The code substitutes "target is 1000 m taller than the detector", conflating a
vertical height difference with slant depth through rock, blind to the target's
thickness and to the orientation of its face. Now superseded by the column-depth
decision (§1.1).

**C. The target's orientation is never checked.** For the τ to exit toward the array,
the target face must point back at the detector.

**C2. The height criterion is mismatched to GRAND's acceptance window.** Requiring the
target to stand `1000 + fresnel_buffer` metres above the detector is equivalent to
demanding a minimum *elevation angle* that varies with distance:

| distance | required elevation | within ±3°? |
| --- | --- | --- |
| 10 km | +6.81° | no — above the window |
| 15 km | +4.52° | no — above the window |
| 20 km | +3.37° | no — above the window |
| 25 km | +2.66° | yes |
| 40 km | +1.58° | yes |
| 80 km | +0.59° | yes |

Two consequences, given the §1.2 window of -3° to +3°:

1. Inside ~22 km the criterion demands terrain **above** the acceptance window, so the
   targets it selects there are ones GRAND could not use anyway.
2. Being a floor, it admits nothing in the **-3° to 0°** half of the window — terrain
   at or below the detector's horizontal — even though those are valid arrival
   directions, and for Earth-skimming trajectories arguably the more important half.

A height threshold cannot express an angular acceptance band; only an explicit
elevation-angle window can. This is the concrete case for §4.11.

**D. Curvature radius 8500 km is the microwave 4/3-Earth refraction rule.** A
defensible convention at 50–200 MHz, but it should be explicit. At 80 km it changes
apparent drop from 502 m (true 6371 km geometry) to 376 m — comparable to the
Fresnel buffer itself.

**E. `target_antennas` does nothing in `distributed` mode.** The threshold applied is
`min_sub_array_size`, and the selection loop appends every qualifying site without
ever stopping at the target. The headline parameter is inert in the default mode.

**F. Reported areas are upper bounds.** Candidates are thinned 5× and then
binary-*closed* with a 1 km structuring element, so a "site" is really "everything
within ~1 km of a validated sample". Square structuring elements also impose
axis-aligned artefacts on site outlines.

### 2.3 Verdict

The tool errs in **both** directions, so the biases do not cancel predictably:

- **Too permissive:** no occlusion test (A), inflated areas (F), no target-orientation check (C)
- **Too restrictive:** single azimuth (§2.1, 2.2×), hard distance box (§3.4), binary cuts with no partial credit
- **Ill-defined:** slope scale (§2.1), depth units (B)

Current output is best read as a plausible shortlist, not a quantitative ranking.

---

## 3. Phase 0 — Foundations ✅ delivered

**Goal:** make every later change measurable. Nothing here changes results.

Phase 1 rewrites the scientific core; without a baseline we cannot tell a fix from a
regression. Phase 0 must land first.

Result-neutrality was verified explicitly: the pipeline was run on two fixed inputs
before and after the instrumentation and the outputs diffed identical.

```bash
cd tests && python -m unittest discover        # 54 tests, ~4 s, no extra dependencies
python bench/benchmark.py                      # compare against bench/baseline.json
python bench/benchmark.py --update             # rewrite the baseline
UPDATE_GOLDEN=1 python -m unittest test_regression   # after an intended change
```

The suite uses stdlib `unittest` rather than `pytest`, which is not installed in the
`sssearch` environment; `pytest` will collect it unchanged if it is ever added.

### 3.1 Test harness (`tests/`)

Synthetic terrain fixtures with analytically known answers, so correctness does not
depend on eyeballing maps. Three already exist in prototype form from the `a9843f9`
verification work and should be promoted into the suite:

- **Planar slope** at a known angle and aspect → slope/aspect recovery, and the
  scale-dependence table of §2.1 as a regression.
- **Isolated peak** at a known azimuth, distance and height → ray direction, distance
  bounds, height threshold, and an A/B against an isotropic-pixel kernel (which must
  miss it — that is the bug fixed in `a9843f9`).
- **Idealized canyon**: two opposing walls, known separation, depth and wall slope →
  the TAMBO geometry primitives of phase 2.
- **Circular exclusion zone** on a uniform plane → the zone must measure the same
  radius on the ground along both axes.

Care needed: probes must use the same rounding convention as the kernel. An early
version of the peak test used `round()` where the kernel truncates and reported
false failures against a 1-pixel target. Targets in fixtures should be several
pixels across, as real terrain features are.

### 3.2 Golden-file regression

Fixed inputs → stored summary JSON (site count, capacity, area, per-stage pixel
counts) with tolerances. Inputs: the 2500² Arequipa crop and the synthetic fixtures.
The full-region DEMs are too slow for CI but should be runnable on demand.

### 3.3 Benchmark harness (`bench/`)

Per-stage wall time and peak RSS on fixed inputs, committed as `bench/baseline.json`.
Cold runs (the cached `.npy` is removed first) so numbers are comparable between
invocations. Measured baseline:

| case | topo screen | ray tracing | morphology | capacity | peak RSS |
| --- | --- | --- | --- | --- | --- |
| synthetic_900 | 0.03 s | 0.17 s | **0.77 s** | 0.02 s | 296 MiB |
| synthetic_1800 | 0.13 s | 0.36 s | **2.54 s** | 0.02 s | 386 MiB |
| arequipa_900 | 0.04 s | 0.03 s | **0.81 s** | 0.01 s | 609 MiB |
| arequipa_2500 | 0.30 s | 0.18 s | **4.61 s** | 0.08 s | 617 MiB |

**Morphology is 77–90% of every run.** The physics is not the bottleneck; the
cleanup is. That is what phase 3's separable/running min-max morphology targets.

One further cost the harness exposed: the ray-caster is JIT-compiled inside joblib's
worker processes, so every fresh invocation pays ~0.75 s of compilation that the
parent cannot amortise. Replacing joblib with numba `prange` (phase 3) removes it
entirely. The benchmark warms the workers before timing so this does not distort the
per-case numbers.

### 3.4 Funnel instrumentation

Count and report pixels surviving each criterion, in order. For a tool that often
returns zero sites, *why* is the most valuable output it can produce, and it is
required to substantiate any "too limiting" claim.

### 3.5 Provenance

Git commit, DEM checksum, resolved config and package versions written into every run
directory.

**Exit criteria:** suite passes on `dev`; baseline timings and funnel counts committed;
`a9843f9`'s behaviour reproducible from a stored golden file.

---

## 4. Phase 1 — Physics core

**Goal:** replace the geometric core with one that computes what it claims to.
Every change lands as an *option* with the previous behaviour reproducible, so each
one's effect on site counts can be quantified rather than discovered later.

### 4.1 Terrain layer ✅ delivered

Slope and aspect are computed at an explicit `--slope_baseline_m`, on a
correspondingly smoothed DEM, so the scale dependence of §2.1 is a stated parameter
rather than an accident of `np.gradient`. Default is `None` (native resolution).

Derivatives are now taken on a **haloed block and cropped**, which removed a latent
tiling artefact: `np.gradient` falls back to one-sided differences at array edges, so
every tile boundary previously carried a wrong slope and aspect. Screening results
depended on `--tile_size`, a purely computational parameter.

This was the first deliberate result change of phase 1, and the phase 0 harness
quantified it on the Arequipa crop golden:

| quantity | before | after |
| --- | --- | --- |
| pixels passing slope | 580,154 | 580,129 |
| ray-tracing hits | 3,214 | 3,183 |
| site area | 117.55 km² | 122.42 km² |
| capacity | 153 | 151 |

Worth noting the amplification: a 25-pixel correction to the slope mask moved the
reported area by 4%, because the 1 km morphological closing magnifies every change in
the validated set. It is a further illustration of finding F.

Tiling invariance is now pinned by tests — tiled screening must reproduce the untiled
result exactly, at every tile size, with and without a slope baseline.

### 4.2 / 4.5 / 4.11 Arrival-direction scan engine ✅ delivered

`src/arrival_scan.py`, with 25 tests against terrain whose answer is known in closed
form. Not yet wired into the pipeline — that is the next step.

For each candidate and azimuth, one profile walk yields **every elevation bin at
once**. The running maximum of the terrain's elevation angle only increases, so each
new maximum claims a contiguous band of bins; and since a ray at angle θ is
underground wherever `θ_terrain(d) > θ`, binning `θ_terrain` and taking an inclusive
suffix sum gives the underground path length for all bins simultaneously. Rays
crossing several ridges accumulate all the rock they traverse, not just the first
chord.

Reported per candidate: accepted-direction count, accepted solid angle (sr), mean
distance to the exit point, maximum and mean column depth (g/cm²), and the horizon
angle. `require_terrain=False` inverts the test for the cosmic-ray channel, where
terrain is an obstruction rather than a target — the same kernel, no second code path.

**Measured on the 2500² Arequipa crop, 755,339 candidates:**

| azimuths | elevation bins | seconds | candidates accepted | median column depth |
| --- | --- | --- | --- | --- |
| 1 | 1 | 7.2 | 0.3% | 8.9×10⁶ g/cm² |
| 5 | 6 | 30.8 | 58.3% | 2.7×10⁶ g/cm² |
| 9 | 12 | 51.7 | 66.7% | 2.9×10⁶ g/cm² |
| 17 | 24 | 94.8 | 71.1% | 3.1×10⁶ g/cm² |

Two things to read from this.

*The design goal holds.* Cost scales with azimuths (~5.6 s each) and is nearly
independent of elevation sampling — quadrupling the bins from 6 to 24 is almost free,
which is exactly what the single-walk-per-azimuth construction was for.

*It is ~290× slower than the ray-caster it replaces* (51.7 s against 0.18 s), which is
the trade §4.11 anticipated. Whole-run cost on this crop goes from ~5.4 s to ~57 s.
Phase 3 targets it; note morphology is still 77–90% of the *old* run, so the profile
is now genuinely dominated by physics rather than cleanup.

*Acceptance rises from 13.9% to 66.7%.* That is finding B and finding C2 combined,
in one number: the old single ray at a distance-dependent elevation floor was
rejecting most of what a ±3° window admits.

### 4.2b A physical constraint that emerged

A detector standing on the ground has **every downward direction blocked by the ground
at its own feet**: over flat terrain a sub-horizontal ray goes underground within a
pixel or two and stays there, so its exit point is metres away and its column depth is
the entire traced path. This is correct, and it is why the decay-baseline window is
not optional — only a site whose local terrain falls away can use the lower half of an
acceptance window at all. It also says something about siting: a good GRAND site is
not merely *on* a slope, it is on a slope whose own terrain does not occlude the
sub-horizontal half of the window. Pinned by tests.

### 4.2c Deferred to phase 3: whole-raster azimuthal sweep

The original plan was to replace per-candidate scanning with a per-azimuth sweep
across the raster, where the horizon for *every* pixel comes from one running-max walk
— O(N) per azimuth with sequential access. That remains the right optimisation, but it
does not extend cleanly to a (θ, φ) scan carrying sub-surface chords, and the
curvature term is not separable in the way the convex-hull construction needs.
Correctness first: the delivered engine scans per candidate, and phase 3 optimises
against the measured 51.7 s.

Two accuracy items also remain for that work: bilinear sampling instead of `int()`
truncation, which currently biases each step by up to half a pixel, and early
termination once no further sample can qualify.

### 4.3 Visibility and Fresnel clearance

A target qualifies only if it *is* the horizon at its distance. Fresnel clearance
becomes a real test of `r₁ = √(λ d₁ d₂ / (d₁+d₂))` against the intervening profile,
with a configurable frequency band, and the clearance fraction reported as a score
rather than a pass/fail.

### 4.4 Curvature and refraction

Effective Earth radius becomes a parameter (k-factor), defaulting to the current
8500 km with the true-geometry option documented alongside.

### 4.5 Column depth

Per the §1.1 decision. Given a detector pixel P and a target surface point X, the
neutrino travels along the line P→X extended *beyond* X into the rock (the τ
continues along the neutrino direction to good approximation). Column depth is then

```
X_rock = ρ_rock · L_rock ,   L_rock = ∫ 1[z(s) < DEM(x(s), y(s))] ds
```

marching outward from X along that line and integrating the sub-surface segment.
Default ρ_rock = 2.65 g/cm³ (standard rock); reported in both g/cm² and km w.e.

Honest limitation to record in the docs: a DEM gives the surface only, so this
assumes uniform density and no sub-surface structure. It is a large improvement on a
height difference, not a geological model.

Because the τ must both be produced *and* escape, the useful column depth is a
**band with an optimum**, not a floor — which is exactly why scores (§1.1) are the
right output form.

### 4.6 Target orientation

The target face's aspect must point back at the detector within a tolerance
(criterion C).

### 4.7 Energy-derived distance windows

The current 10–80 km box is replaced by a window derived from the τ decay length:

```
L_decay = (E / m_τ) · cτ_τ ,   m_τ = 1.77686 GeV, cτ_τ = 87.03 µm
```

| E | L_decay |
| --- | --- |
| 1 PeV | 49 m |
| 10 PeV | 490 m |
| 100 PeV | 4.9 km |
| 1 EeV | 49 km |
| 10 EeV | 490 km |

Two independent consistency checks support this parameterisation:

- TAMBO's stated 1–100 PeV range gives 49 m – 4.9 km, matching the "Range 50 m–5 km"
  annotation in ref. [2] Fig. 1 exactly.
- The current hardcoded 10–80 km GRAND window corresponds to 0.2–1.6 EeV — i.e. the
  existing default silently encodes an energy assumption.

Shower development length is added on top (ref. [2] quotes 3–10 km shower length,
200 m diameter). Energy ranges stay configurable per experiment as required.

### 4.8 Scores

Every criterion returns a value in [0,1] with a documented shape, plus its diagnostic
count for the funnel. Hard physical impossibilities remain score 0. Composition rule
(product, weighted mean, or min) is a per-experiment setting — default to be chosen
once we see the score distributions on real terrain.

### 4.9 Validation against published results

The strongest verification available: **the tool should reproduce the published
apertures for the published sites.**

- TAMBO: reproduce the Colca Canyon aperture curve, ref. [2] Fig. 3, for 5,000 units
  at 150 m spacing.
- GRAND: reproduce the shape of ref. [1] Fig. 25.

Important caveat on how the published curves can be used. They are **integral**
quantities — integrated over the whole array, all geometries and the whole site — so
they cannot be applied as a per-pixel lookup. Their roles are:

1. a **validation anchor** (§4.9), and
2. an **energy-response weight** when comparing energy ranges.

Per-pixel weighting needs a *differential* acceptance in (distance, elevation angle,
column depth, energy).

### 4.10 Working without a differential acceptance table (decided 2026-08-03)

No such table is available at present. That blocks **absolute apertures**, but not
site *ranking* — which is what a site-search tool is for. Any factor that is
energy-dependent but site-independent cancels when comparing two sites for the same
experiment and energy band, so a defensible relative score is available today.

The plan adapts as follows.

1. **Compute and store the geometric observables, not just a score.** For every
   candidate: target distance, arrival elevation angle, column depth, usable azimuth
   range, and clearance. Sites carry *histograms* of these, not scalars.
2. **Score relatively**, within one experiment and energy band, from those
   observables plus the analytic factors below.
3. **Keep the physics response pluggable**: a documented default parameterisation
   with explicit assumptions, behind an interface that loads a CSV table when one
   exists.
4. **Absolute apertures come later for free.** Because the histograms are stored,
   folding them against an acceptance table is a post-processing step — no re-running
   of the expensive terrain analysis.

One factor *is* exactly analytic and needs no table, and it is the one that couples
most strongly to geometry: the probability that the tau decays inside the useful
range,

```
P_decay = exp(-d_min / L) - exp(-d_max / L),    L = (E/m_tau) * c*tau
```

The pieces that genuinely need simulation are the tau exit probability given column
depth, and the trigger probability given shower geometry. Those get parameterised
defaults with the assumptions written down, and are replaced when a table arrives.

Validation against ref. [1] Fig. 25 and ref. [2] Fig. 3 remains possible: integrating
our per-site model over the published Colca configuration should reproduce the
published aperture up to a single free normalisation. That is a real test — a model
with the wrong geometry dependence will not match the *shape* whatever the
normalisation.

### 4.11 Engine geometry: scan arrival directions, not "find a tall mountain"

Working through the column-depth calculation exposed that the current formulation is
not merely crude but structurally wrong, and that the fix reshapes the engine.

The current test asks: *is there terrain 1 km taller than me at 10-80 km?* The
physical question is: *from which arrival directions does a backward ray from this
pixel enter rock, and how much rock does it cross?*

Tracing backward from a candidate along an arrival direction (azimuth phi, elevation
angle theta):

- Rays **above** the local horizon escape to the sky and contribute nothing.
- Rays **below** the horizon strike terrain. That first intersection is the tau exit
  point; the distance to it is the decay baseline; and the chord beyond it, where the
  ray runs under the surface, is the column depth.

So the useful quantity is a scan over **(azimuth, elevation angle)** pairs, each
yielding (distance, column depth) — which is exactly the differential geometry §4.10
wants histogrammed. Note the grazing ray is the *least* useful one: at the horizon
the column depth goes to zero. Terrain that is merely tall is not the target; terrain
that subtends solid angle with rock behind it is.

The §1.2 guidance sets the scan range directly, and makes each experiment a
configuration of the same engine rather than a separate code path:

| channel | elevation scan | what the engine reports |
| --- | --- | --- |
| GRAND ν | -3° to +3° about horizontal | rock-backed solid angle, distance and column-depth histograms |
| GRAND CR | above the horizon | unobstructed sky solid angle; nearby terrain is a *penalty*, not a target |
| TAMBO | across a canyon | opposing-wall distance, depth and subtended angle |

The cosmic-ray channel is a useful check on the design: it inverts the test — terrain
in the accepted directions *reduces* the score instead of enabling it. Any framework
that cannot express both from one scan is not general enough for the experiments
after these two.

This supersedes §4.2's per-candidate ray fan as the engine's shape.

**Sequencing correction.** §4.2 claimed the whole-raster azimuthal sweep is both the
accuracy fix and a speedup. That holds for a pure horizon calculation, but the O(N)
skyline algorithm does not extend cleanly to a (theta, phi) scan with sub-surface
chords, and the curvature term is not separable in the way the convex-hull trick
needs. Correctness first: phase 1 implements a direct scan, phase 3 optimises it
against the phase 0 baseline. Expect phase 1 to be *slower* than the current code —
the harness exists to quantify exactly that.

---

## 5. Phase 2 — Generalization beyond radio detection *(sketch — to be scoped)*

The unifying observation: GRAND and TAMBO ask the same structural question — *"from
this patch of ground, is there a target surface at the right range, in the right
direction, at the right relative orientation, with the right amount of matter behind
it?"* They differ in numbers, not in structure. GRAND wants a tall mass 10–80 km away
in a forward arc; TAMBO wants a steep opposing wall a few km across a gorge. **One
sweep engine (§4.2) answers both.**

Proposed layering: terrain → criteria (local + view) → experiment spec (YAML) →
region/layout → combination (`all` / `any` / weighted, plus *joint* vs *union* and a
co-location report).

### 5.1 TAMBO parameters from ref. [2]

| quantity | value |
| --- | --- |
| detection units | 5,000 plastic scintillator |
| spacing | 150 m, triangular grid (**start at 100 m per §1.1**) |
| Colca Canyon depth | ~1.5 km |
| median distance between valley sides | 4.5 km |
| τ range in valley | 50 m – 5 km |
| air shower | 3–10 km length, 200 m diameter |
| energy reach | ~3 PeV – 1 EeV |

### 5.2 TAMBO search geometry (guidance, 2026-08-03)

| parameter | value | note |
| --- | --- | --- |
| wall-to-wall separation | **a few km**, configurable | Colca's published median is 4.5 km |
| canyon depth | **1.5 km** (Colca as the worked example) | ref. [2] Fig. 2 |
| detector spacing | 100 m to start | published nominal 150 m, §1.1 |
| layout | triangular grid | ref. [2] |

Separation is a *window*, not a single value, and the physics fixes its scale
independently: the τ must decay in the air gap, and at 1–100 PeV the decay length is
49 m – 4.9 km (§4.7). A few km is therefore exactly right, and the same reasoning
says the window should move with the configured energy range rather than being a
fixed constant. Proposed default 2–6 km, bracketing Colca.

Two things fall out of the published geometry that matter for the criterion design:

- **Colca's walls are ~40° steep.** Taking 1.5 km depth, 4.5 km rim-to-rim and a ~1 km
  floor gives `4500 = 1000 + 2*1500/tan(s)`, so `s ≈ 40.6°`. That is well outside
  GRAND's 3–25° deployable band, so the **slope criterion itself must be
  per-experiment** — and probably per-*role*: the far wall wants to be steep for τ
  exit, while the near wall must still be deployable. A single global slope band
  cannot express that.
- **Depth and separation are not independent.** Any two of {depth, separation, wall
  slope} determine the third, so the criteria should be stated in terms of the two
  that are measured (depth and separation) with wall slope derived, rather than
  filtering on all three as though they were free.

The fixtures now carry a `COLCA` preset with these numbers, verified against the
published figures, so TAMBO criteria can be developed against known-answer terrain
before touching a real DEM.

Two consequences worth flagging early:

- **Site shape differs fundamentally between experiments.** 5,000 units at 150 m
  triangular spacing need ~97 km² of usable canyon wall — a long strip along the
  canyon, not a compact blob. The current `min_width_km` opening with a square
  structuring element would destroy such a strip. Region/layout models must be
  per-experiment, not shared.
- **The 100 m vs 150 m spacing discrepancy** (§1.1 vs ref. [2]) should be confirmed;
  the roadmap assumes it is a deliberate starting point and keeps spacing configurable.
- **Integer grid stamping cannot express TAMBO's layout.** At 100–150 m spacing on a
  30 m DEM there are only 3–5 pixels per detector separation, and the truncation
  cascade of §2.1 over-counts capacity by 42–58%. Capacity for TAMBO must be computed
  analytically from usable area and wall geometry, or on a resampled grid — not by
  stamping integer strides. This is why layout models have to be per-experiment.

Ref. [2] Fig. 1 also annotates ">4 km shielding from background muons" — a possible
additional criterion on rock overburden. **To confirm** whether this is a site
requirement we should encode.

---

## 6. Phases 3–4 *(sketch — to be scoped)*

**Phase 3 — performance.** Separable/running min-max morphology (O(N) instead of
O(N·k²), targeting the 6.8 s of §3.3); numba `prange` replacing joblib in the hot
path; native-dtype DEM caching; content-addressed cache keys. Precision is improved,
not traded: bilinear sampling and DEM-resolution marching (§4.2).

**Phase 4 — usability.** Auto-detect `origin_lat`/`origin_lon` from the GeoTIFF
tiepoint (verified present, matching current configs to ~1e-4°); rename `src/setup.py`
(it is not a packaging file and that name hijacks `pip install`); real packaging;
rasterio/pyproj for CRS and outputs; `--explain` funnel report; parameter sweeps.

---

## 7. Open questions

1. **Differential acceptance (blocking for §4.9).** Is there a simulation output —
   NuTauSim, ZHAireS, or the collaborations' own chains — that can be tabulated as
   acceptance vs (distance, elevation angle, column depth, energy)? Published
   integral curves alone cannot weight individual pixels.
2. **Score composition** (§4.8): product, weighted mean, or min?
3. **TAMBO muon shielding** (§5.1): is ">4 km rock" a site-selection criterion?
4. **Spacing** (§5.1): confirm 100 m as the starting value against the published 150 m.

---

## References

[1] GRAND Collaboration (J. Alvarez-Muñiz et al.), *The Giant Radio Array for
Neutrino Detection (GRAND): Science and Design*, arXiv:1810.09994 (2018);
Sci. China Phys. Mech. Astron. (2020). Effective area: Fig. 25.
20 sub-arrays × 10,000 antennas over 10,000 km² each, E > 10⁸ GeV.

[2] TAMBO Collaboration, *Measuring the high-energy neutrino sky using the deep-valley
neutrino observatory TAMBO*, Nature Astronomy (2026),
doi:10.1038/s41550-026-02916-4. Exposure: Fig. 3.
