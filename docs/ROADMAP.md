# Site Search — Development Roadmap

Working document. Phase 0 and phase 1 are delivered; phases 2–4 are sketched and
will be scoped in detail before each is started.

**Status:** phase 0 (harness) and phase 1 (physics core) are complete on `dev`.
Sections marked ✅ are implemented and covered by tests; everything else is plan.
Section 2 is the original review of the inherited code and is kept as the record of
why the work was done — the defects it describes are fixed unless noted.

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
elevation-angle window can. This is the concrete case for §4.2.

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

## 4. Phase 1 — Physics core ✅ delivered

**Goal:** replace the geometric core with one that computes what it claims to.
Every change lands as an *option* with the previous behaviour reproducible, so each
one's effect on site counts can be quantified rather than discovered later.

All of phase 1 is implemented and covered by tests. In summary: slope is measured
over a stated baseline; the engine scans arrival directions and measures column depth
instead of looking for a tall mountain; distance windows follow from an energy range;
Fresnel clearance and the refraction k-factor are explicit; and sites are ranked by
composable scores whose components are reported separately.

Two planned items turned out not to need separate work. Target orientation (§4.6) is
subsumed by column depth. The whole-raster sweep (§4.2.2) moved to phase 3, where
optimisation belongs.

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

### 4.2 Arrival-direction scan engine ✅ delivered

*Supersedes the original §4.2 sweep plan, §4.5 column depth and §4.2 engine geometry.*

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
the trade §4.2 anticipated. Whole-run cost on this crop goes from ~5.4 s to ~57 s.
Phase 3 targets it; note morphology is still 77–90% of the *old* run, so the profile
is now genuinely dominated by physics rather than cleanup.

*Acceptance rises from 13.9% to 66.7%.* That is finding B and finding C2 combined,
in one number: the old single ray at a distance-dependent elevation floor was
rejecting most of what a ±3° window admits.

### 4.2.1 A physical constraint that emerged

A detector standing on the ground has **every downward direction blocked by the ground
at its own feet**: over flat terrain a sub-horizontal ray goes underground within a
pixel or two and stays there, so its exit point is metres away and its column depth is
the entire traced path. This is correct, and it is why the decay-baseline window is
not optional — only a site whose local terrain falls away can use the lower half of an
acceptance window at all. It also says something about siting: a good GRAND site is
not merely *on* a slope, it is on a slope whose own terrain does not occlude the
sub-horizontal half of the window. Pinned by tests.

### 4.2.2 Deferred to phase 3: whole-raster azimuthal sweep

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

### 4.2.3 Pipeline integration ✅ delivered

`--physics_mode {legacy,scan}`, default `legacy` so the previous behaviour stays
reproducible and every change remains measurable. Accepted candidates feed the same
morphology and labelling stages, and each site carries the **distribution** of its
candidates' observables (`arrival_scan` in the results JSON: mean/p50/p90 of solid
angle, distance, column depth and horizon angle) — the §4.10 structure that lets
absolute apertures be folded in later without re-running the terrain analysis.

**A/B on the 2500² Arequipa crop**, 5–25 km window, 9 azimuths, 12 bins:

| mode | wall time | sites | capacity | area | candidates passing physics |
| --- | --- | --- | --- | --- | --- |
| legacy | 12.6 s | 1 | 2,176 | 1,525 km² | 11.2% |
| scan | 77.6 s | 1 | 6,157 | 4,969 km² | 64.6% |

### 4.2.4 The window sits below the horizon almost everywhere

The scan reports a **median horizon angle of 7.3°** across accepted candidates, with
p90 at 17°. In terrain like this the whole ±3° acceptance window lies *below* the
local horizon, so essentially every direction in it strikes rock. That is why
acceptance is 64.6% rather than something selective.

The consequence matters for the design: **"does this direction hit terrain" carries
almost no discriminating power in the Andes.** The selection has to come from *how
much* rock and *at what distance* — the column-depth band and the decay probability —
not from visibility. The delivered site covers 4,969 km² of a ~5,600 km² crop, i.e.
the criterion as currently configured selects nearly the whole map.

This is not a defect of the engine, which is measuring correctly; it is the
measurement telling us where the physics has to do the work. It also raises the
priority of §4.8 (scores) and of a physically-motivated depth band over everything
else remaining in phase 1.

### 4.3 Fresnel clearance ✅ delivered

`--fresnel_frequency_mhz` enables a second pass over accepted directions, measuring
the worst clearance along the path in units of the first Fresnel radius
`r₁ = √(λ d₁ d₂ / (d₁+d₂))`, reported as a score rather than a gate. The scan already
guarantees an unobstructed line of sight — the intersection is by construction the
first terrain met — so this catches paths that *graze* an intervening ridge and suffer
diffraction loss even though the geometric path is clear.

Getting a usable measure took three corrections, each found by measuring rather than
by reasoning, and each worth recording because they are properties of the problem
rather than of the code:

1. **The far endpoint must be the shower, not the exit point.** Approaching the target
   both the clearance and r₁ go to zero, so their ratio collapses for every path
   regardless of whether anything obstructs it. The radio source is the air shower
   developing some kilometres after the τ decays, and using that as the endpoint
   removes the degeneracy.
2. **The antenna must be given a height.** A receiver at ground level always has
   terrain inside the first Fresnel zone immediately beside it.
3. **The near field must be excluded.** Even at 5 m, the measure was dominated by the
   ground next to the antenna: the median ratio swung by **28×** across mast heights
   from 0 to 100 m, a parameter with nothing to do with site quality. Skipping the
   first 500 m cuts that spread to **2.1×**, and the measure then responds to
   intervening ridges as intended. Below about 500 m it is really measuring ground
   roughness at a scale a 30 m DEM cannot resolve.

Defaults: 5 m antenna height, 500 m near-field exclusion, 3 km shower offset. On the
Arequipa crop the median best-direction clearance is then 1.7 r₁.

### 4.4 Curvature and refraction ✅ delivered

`--refraction_k` sets the effective Earth radius as a k-factor: k = 1 is true
geometry, k = 4/3 the radio convention that yields the 8500 km the tool has always
used. Over an 80 km path the apparent drop is 376 m at k = 4/3 against 502 m at
k = 1 — a difference comparable to the Fresnel clearance itself, which is why it
should be stated rather than assumed.

### 4.6 Target orientation — subsumed by column depth ✅

No separate criterion is needed. A face sloping away from the candidate is never
struck at all, and among faces that are struck, the depth measurement already
distinguishes them. Verified by tests.

Worth recording the direction of the effect, which is not the intuitive one: a
*gentler* face presents **more** rock for the same summit height, because a
near-horizontal ray runs further underground before reaching the summit. Since the τ
must escape as well as be produced, this is direct evidence that the depth criterion
wants an **optimum band rather than a floor** — and it is consistent with TAMBO
wanting steep canyon walls, where the τ has to exit through a short path.

### 4.7 Energy-derived distance windows ✅ delivered

`--energy_min_pev` / `--energy_max_pev` set the decay-baseline window from
`L = (E/m_τ)·cτ`, and the run banner reports the correspondence in both directions, so
a hand-set distance window now says out loud what energies it implies. The helpers
reproduce the published numbers: 1–100 PeV gives 49 m – 4.9 km against ref. [2]'s
quoted 50 m – 5 km, and the inherited 10–80 km GRAND default corresponds to
0.2–1.6 EeV. `decay_probability()` provides the one exactly-analytic acceptance factor.

The mapping from an energy range to a window is a stated convention, not a derivation:
it fixes the scale correctly, but the useful window also depends on acceptance details
this tool does not model.

### 4.8 Scores ✅ delivered

`src/scoring.py`. Every component returns [0, 1] with a documented shape, and the
components are reported separately so a site's weakness can be attributed rather than
disappearing into one opaque number.

| component | shape | rationale |
| --- | --- | --- |
| `depth` | band with soft flanks | the τ must be produced *and* escape, so an optimum, not a floor |
| `distance` | band over the decay-baseline window | defaults to the configured window |
| `solid_angle` | saturating, `x/(x+half)` | more is better with diminishing returns |
| `clearance` | ramp in r₁ | present only when a frequency is configured |

Composition is `product` (unforgiving), `mean` (compensating) or `min` (weakest link),
selected per experiment. A component of exactly zero sinks a product — physical
impossibilities must score zero, not merely small.

The depth band default is deliberately wide (1e5–1e7 g/cm²). The physically motivated
band for a given energy range is still open (§7), and a wide default ranks sites
without pretending to encode physics the tool has not been given.

**Measured effect** on the Arequipa crop, 5–25 km window, 9 azimuths, 12 bins:

| configuration | sites | area | capacity |
| --- | --- | --- | --- |
| legacy ray-caster | 1 | 1,525 km² | 2,176 |
| scan, no scoring floor | 1 | 4,969 km² | 6,157 |
| + depth band 3e5–5e6, `min_score` 0.35 | 1 | 3,685 km² | 4,736 |
| + Fresnel 50 MHz, k = 4/3 | 1 | 3,385 km² | 4,442 |

### 4.9 Validation ✅ delivered, within what is checkable

`src/aperture.py` separates the estimate into three parts: the **geometric aperture**
(area × solid angle, fully determined by terrain), the **analytic decay factor**
(`exp(-d_min/L) - exp(-d_max/L)`, no free parameters), and a **pluggable response**
defaulting to unity, replaceable by a table via `TabulatedResponse`. Under that split
the absolute normalisation is unknown but the energy *shape* and the site *ranking*
are not — and both are tested.

What is validated:

- aperture scales linearly with area and with solid angle;
- with a unit response the energy dependence is exactly the decay probability;
- the high-energy tail falls as 1/E, since a fixed window subtends
  `(d_max - d_min)/L`;
- a canyon baseline peaks in TAMBO's PeV band, the inherited 10–80 km GRAND window
  peaks near an EeV — so a site's geometry alone predicts the energies it suits;
- the τ decay length reproduces ref. [2]'s quoted 50 m – 5 km for 1–100 PeV;
- the Colca fixture reproduces the published 1.5 km depth and 4.5 km separation.

What is **not** validated, and why. Ref. [1] Fig. 25 and ref. [2] Fig. 3 are integral
over a whole array, all geometries and one site, so they cannot be applied per pixel;
they can anchor the normalisation once supplied *as data*. Reading numbers off a
published figure by eye is not a measurement, so the module provides the machinery to
compare against a supplied curve rather than a transcription of one. Supplying either
curve as a two-column CSV would close this.

### 4.4.1 Correction: particle geometry and radio propagation need different radii ✅

The 4/3-Earth radius was being applied to the *particle* geometry, which is wrong.
Neutrinos and taus are not refracted; they travel in straight lines, so the geometry
deciding where the tau exits uses the **true** 6371 km radius. The 4/3 convention
exists to straighten a refracted radio ray, and applies to the Fresnel clearance of
the signal path and to nothing else. The two are now separate, and `--refraction_k`
governs only the radio path.

Effect at 25 km is small (median column depth 1.989e6 → 1.981e6 g/cm², acceptance
64.6% → 64.5%) but it grows with distance, and it was free to fix.

### 4.12 Physics accounted for ✅ delivered

The sweep of §4.12 below identified seven omissions; six are now implemented in
`src/physics.py` and the scan kernel, with 44 tests. All the analytic pieces are
closed-form and checkable by hand, which is why they live apart from the terrain code.

**Measured on the Arequipa crop** (4335 m median site altitude, 10.1 km median
baseline, 9 azimuths, 12 bins):

| quantity | value | what it says |
| --- | --- | --- |
| mean `sin(α)` | **0.601** | geomagnetic weighting removes 40% of the effective acceptance |
| median score, geometry only | 0.163 | |
| median score, with geomagnetism | **0.094** | a 42% reduction: first-order, as predicted |
| shower maturity | **1.07 × X_max** | these sites sit *exactly* at shower maximum |
| footprint radius | 189 m | at 4335 m and 10.1 km |
| antennas across the footprint | **0.38** | a 1 km grid under-samples it ~2.6× |

Three of these are statements the tool could not previously make at all.

**(b) Geomagnetic weighting.** `--geomag_declination_deg` / `--geomag_inclination_deg`
weight each accepted cell by `sin(α)` to the field, giving a second solid angle
alongside the raw one. Left unweighted unless a field is supplied, since guessing a
field vector would be worse than declining to weight. Verified against the closed
form: near the magnetic equator a north-facing target retains under 5% of its raw
acceptance while an east-facing one keeps essentially all of it — two sites with
identical terrain statistics that no geometric measure can tell apart.

**(a) Atmospheric grammage.** The slant integral has a closed form, so no numerical
integration is needed, and it is checked against a 200,000-step numerical integral.
Reported per site and scored as a band around X_max. A given path at 4000 m carries
0.622 of its sea-level grammage — and the tool is choosing between sites that differ
by exactly that much altitude.

**(c) Earth chord.** `2R sin θ` for downgoing directions, with an optional
`--nu_interaction_length_gcm2` attenuation term. Reported always, weighted only when
an interaction length is supplied. Worth noting what the measurement showed: with a
5 km minimum baseline the negative half of the window is mostly rejected anyway,
because near ground blocks it (§4.2.1), so the chord term rarely engages for these
sites. It will matter more for configurations with short baselines.

**(d) Energy-dependent depth band.** `depth_band_from_energy()` derives the band from
the tau range, combining the boosted decay length and the energy-loss length
harmonically so the range grows then saturates. The energy-loss constant β is the
least certain number in the module — published values span roughly 0.4–1.0e-6 cm²/g —
so this fixes the *scale* of the useful depth rather than its precise value.

**(e) Footprint versus spacing.** The Cherenkov cone narrows with altitude, so a higher
site has a *smaller* footprint and needs a *denser* array. Scored as antennas across
the footprint diameter. The measured 0.38 quantifies GRAND's sparse-array trade.

**(f) RFI line-of-sight shielding.** Sources occluded by terrain contribute nothing;
survivors contribute as `1/d²`. This reuses the horizon machinery the scan already
needs, so the only new cost is one short walk per candidate per source.

**(g) TAMBO channel differences** remain for phase 2, where per-channel criteria live.

### 4.12b Original sweep: physics not yet accounted for

Recorded from a deliberate sweep for omissions rather than found by a failing test.
Ordered by how much each would change site *ranking*, which is what the tool exists to
produce.

**(a) Atmospheric grammage — shower development is not measured in metres.** An air
shower develops through slant depth in g/cm², not path length, and air density falls
as `exp(-h/H)` with `H ≈ 8.4 km`. At 4000 m the density is 0.62 of sea level, so a
20 km path yields ~1500 g/cm² against ~2450 at sea level. The tool currently treats a
kilometre at 2000 m and at 4500 m as the same kilometre, **while comparing candidate
sites that differ by exactly that much altitude**. Fix: integrate density along each
accepted path and express the distance criterion relative to X_max (~700 g/cm²)
instead of in metres. This makes the useful distance window altitude-dependent, as it
physically is.

**(b) Geomagnetic angle — the azimuth of a target matters, not just its existence.**
Radio emission from an air shower is dominantly geomagnetic, with amplitude
proportional to `sin(α)`, α being the angle between the shower axis and the local
field **B**. Peru sits near the magnetic equator, where **B** is close to horizontal
and roughly northward, so showers propagating north–south have strongly suppressed
geomagnetic emission while east–west ones are near maximal. The scan currently treats
all azimuths as equivalent, which is a first-order omission for GRAND: two sites with
identical terrain statistics but differently oriented targets are not equally good.
Fix: weight each accepted (azimuth, elevation) cell by `sin(α)` from an IGRF field
vector for the site. Askaryan emission is not field-dependent and keeps the
suppression finite rather than absolute.

**(c) The Earth chord dominates the negative half of the window.** Tracing backward at
angle θ below horizontal, the ray descends into the bulk Earth: the chord is
`2R sin θ`, about 220 km at −1° and 670 km at −3°. The scan truncates column depth at
`max_range` (tens of km), so it systematically under-measures downgoing directions.
The physics wants two quantities, not one: the **full chord** governs neutrino
attenuation and is analytic in θ, independent of topography; the **last tau-range
before exit** governs tau production and escape, and *is* what the DEM chord measures.
Splitting them would make the depth criterion mean something definite, and would
explain why the measured depth distribution presses against the `max_range` ceiling.

**(d) The depth optimum moves with energy.** The neutrino interaction length shortens
and the tau range lengthens as energy rises, so the optimal column depth is
energy-dependent. The band is currently a constant, which is defensible only while it
is a placeholder (§7).

**(e) Footprint versus spacing.** The radio footprint is a Cherenkov ring a few hundred
metres across, against a 1 km antenna grid, so counted antennas are a cost proxy and
not an effective area. Worse, `n − 1` falls with altitude, shrinking the Cherenkov
angle and the footprint, so a high site needs *denser* spacing than a low one for the
same trigger efficiency — an altitude–spacing coupling the layout model does not have.

**(f) RFI shielding is free and is not being used.** The scan already computes horizons.
Whether a settlement is line-of-sight visible from a candidate is therefore nearly free
to evaluate, and far more physical than a circular exclusion zone: a town behind a
ridge is not equivalent to one in plain view.

**(g) TAMBO differs in kind, not degree.** Particle detection means no Fresnel term and
no geomagnetic dependence, with the footprint set by lateral particle spread (~100–200 m,
which is what the published 150 m spacing matches). And the canyon width is limited by
grammage as much as by decay length: at 3500 m, 4.5 km of air is only ~360 g/cm²,
well short of X_max, so TAMBO observes showers still developing. This reinforces the
phase 2 conclusion that criteria are per-channel.

**Not worth modelling for ranking:** Galactic background noise sets the absolute
trigger threshold but is site-independent, so it cancels in a comparison (§4.10).

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

1. **A column-depth band for GRAND above 100 PeV.** Now the highest-value input.
   Section 4.2.4 showed the geometric test barely discriminates in Andean terrain, so
   the depth band decides nearly everything the tool selects. The scan measures a
   median of ~2×10⁶ g/cm² (about 7.5 km of standard rock) across accepted candidates;
   the default band is a placeholder. A band grounded in the τ production-and-escape
   optimum would make §4.8 physical rather than merely ordinal.
2. **Score composition** (§4.8): `product`, `mean` or `min` as the default? Now
   answerable from real score distributions rather than in the abstract.
3. **TAMBO muon shielding** (§5.1): is ">4 km rock" from ref. [2] Fig. 1 a
   site-selection criterion, or a description of the geometry?
4. **A published aperture curve as data** (§4.9). Either ref. [1] Fig. 25 or ref. [2]
   Fig. 3 as a two-column CSV would let the tool's shape be checked against a
   measured one, and would fix the normalisation.

Answered and folded into the plan: the differential acceptance table is unavailable
(§4.10 records how phase 1 proceeded without it), and TAMBO spacing starts at 100 m
(§5.2).

---

## References

[1] GRAND Collaboration (J. Alvarez-Muñiz et al.), *The Giant Radio Array for
Neutrino Detection (GRAND): Science and Design*, arXiv:1810.09994 (2018);
Sci. China Phys. Mech. Astron. (2020). Effective area: Fig. 25.
20 sub-arrays × 10,000 antennas over 10,000 km² each, E > 10⁸ GeV.

[2] TAMBO Collaboration, *Measuring the high-energy neutrino sky using the deep-valley
neutrino observatory TAMBO*, Nature Astronomy (2026),
doi:10.1038/s41550-026-02916-4. Exposure: Fig. 3.
