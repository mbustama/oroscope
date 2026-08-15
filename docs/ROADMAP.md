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

**✅ Fixed (§6.15).** Detector positions are now laid out in metres and only then
looked up in the pixel grid, so there is no stride to truncate. Measured density is
within 2% of analytic from 1000 m down to 60 m, on both square and triangular
lattices. Reported capacity fell by 7.3% on the Arequipa golden crop and 8.5% on the
synthetic ridge, both at 1 km; at TAMBO's 100 m the old count was 1.58× the analytic
density, so anything previously quoted there was overstated by more than half.

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

**Both published curves are now digitized** into `data/`, at the request of the
project. They were traced programmatically rather than read by eye: the plot frame and
major ticks were located from the rendered vector figure, and the curves extracted by
colour.

Ref. [1] Fig. 25 carries its own calibration check. The caption states the GRAND200k
curve is exactly 20x the GRAND10k one, and tracing both independently gives ratios of
**19.9-20.1** across the well-resolved range, which validates the axis calibration and
the tracing to about 0.5%. The TAMBO curve reproduces the paper's stated flattening
above 1 EeV, at 6.7e4 m^2 sr.

They remain transcriptions of figures, not tabulated values, and both are integral
over one array and one site, so they still cannot weight individual pixels.

**What they do enable is inferring the response.** Dividing a published curve by the
two factors this tool computes -- geometric aperture and the analytic decay
probability -- leaves everything else:

    response(E) = published(E) / (area * Omega * P_decay(E))

which is the neutrino-interaction, tau-exit and trigger physics that section 4.10 had
to leave out. Applied to the TAMBO curve it is well-conditioned above about 60 PeV and
rises roughly as **E^1.2**, consistent with a rising cross-section and lengthening tau
range. Below that the decay probability is small enough that the division is
meaningless, and `infer_response()` excludes that region rather than letting it
dominate. This is a better weight than a flat response, not a substitute for a
differential table.

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
| shower maturity | **1.07 × X_max** | showers are mature on arrival, but only just |
| footprint radius | 189 m | at 4335 m and 10.1 km |
| antennas across the footprint | **0.38** | a 1 km grid under-samples it ~2.6× |

Three of these are statements the tool could not previously make at all.

**(b) Geomagnetic weighting.** `--geomag_declination_deg` / `--geomag_inclination_deg`
weight each accepted cell by `sin(α)` to the field, giving a second solid angle
alongside the raw one, and `--no_geomagnetic` turns the effect off.

The defaults are for the Peruvian Andes, and their provenance differs:

| parameter | default | source |
| --- | --- | --- |
| declination | **-6.9°** | IGRF 2026 at Arequipa (16.4°S, 71.5°W) |
| inclination | **-14.0°** | centered-dipole estimate at the same point |

The inclination is an estimate because the IGRF value was not retrievable; the dipole
model reproduces it as -13.99° from Arequipa's -7.1° magnetic latitude.
`centered_dipole_inclination()` will do the same for any site. That approximation is
worth using for inclination only — the dipole declination at Arequipa is about -0.2°
against an IGRF -6.9°, so non-dipole terms dominate there and a dipole declination
would mislead. **Both defaults should be replaced with IGRF values per site.**

With the default field the azimuthal asymmetry is large:

| target azimuth | 0° (N) | 45° | 90° (E) | 135° | 180° (S) |
| --- | --- | --- | --- | --- | --- |
| `sin(α)` | 0.269 | 0.801 | **0.993** | 0.646 | 0.269 |

An east-facing target is worth 3.7× a north-facing one, for identical terrain.
Enabling the weighting drops the crop's median score from 0.167 to 0.109. Verified against the closed
form: near the magnetic equator a north-facing target retains under 5% of its raw
acceptance while an east-facing one keeps essentially all of it — two sites with
identical terrain statistics that no geometric measure can tell apart.

**(a) Atmospheric grammage.** The slant integral has a closed form, so no numerical
integration is needed, and it is checked against a 200,000-step numerical integral.
A given path at 4000 m carries 0.622 of its sea-level grammage — and the tool is
choosing between sites that differ by exactly that much altitude.

*Correction, after review:* this was first scored as a **band** around X_max, which is
wrong for radio. Radio emission comes from the region around shower maximum and then
simply propagates, and air is effectively transparent at 50–200 MHz, so a detector
well beyond maximum loses nothing on that account. The criterion is a **threshold**:
the shower must have matured, and past that there is no penalty here. The remaining
trade at greater distance is amplitude against footprint area, which belongs to the
footprint term (e) and not to grammage.

A particle array is the opposite case — charged-particle content peaks at maximum and
dies away after — so `--grammage_mode particle` keeps the band for TAMBO. One more
instance of criteria having to be per-channel. The measured 1.07 × X_max therefore
says these showers are mature on arrival with little margin, not that the sites are
optimally placed.

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

### 4.13 The production-and-escape optimum ✅ estimated, then corrected

The one substantive physics input still missing was a column-depth band with a basis.
It can be derived rather than assumed.

**A first version of this section was wrong and is corrected here.** It combined the
tau's decay and energy-loss lengths harmonically, which saturates at `1/β`, and used a
simple exponential for the tau's survival. Both are wrong, and the errors moved the
answer.

#### The corrected treatment

Decay and energy loss couple, because losing energy shortens the boosted decay length.
With `E(X) = E₀exp(-βX)` the decay probability per unit column depth is
`exp(βX)/X_decay(E₀)`, so

```
S(d)  = exp[ -(X_loss/X_decay)(exp(d/X_loss) - 1) ],   X_loss = 1/β
R_τ   = X_loss · ln(1 + X_decay/X_loss)
```

The range therefore **grows logarithmically** rather than saturating — the harmonic
form understated it by 2× at an EeV and 4× at 10 EeV — and survival falls as a *double*
exponential, far more sharply than `exp(-d/R)` beyond `1/β`.

The exit probability then has no closed form and is integrated numerically:

```
P(X) = ∫₀ˣ (dx/λ) exp(-x/λ) · S(X - x)
```

with the tau carrying `(1 - y)` of the neutrino energy, `⟨y⟩ ≈ 0.2`.

#### β, estimated

β is the radiative energy-loss coefficient in `-dE/dX = a + βE`, units cm²/g, and
`1/β` is the depth over which the tau's energy falls by `1/e`. Estimating from muon
coefficients in standard rock (brems 1.6, pair 2.0, photonuclear 0.4, ×10⁻⁶):
bremsstrahlung and pair production scale as `1/m²` and are suppressed by
`(m_μ/m_τ)² = 1/283`, contributing ~3%. Photonuclear depends on the lepton mass only
logarithmically and dominates entirely.

That gives **β_τ ≈ (0.4–1.0)×10⁻⁶ cm²/g**, cross-checked by `1/β` = 3.8–9.4 km of
rock, bracketing the ~10 km at which the tau range is usually quoted to saturate.
Since photonuclear grows with energy, β is modelled as `0.6×10⁻⁶ (E/1 EeV)^0.20`,
giving 0.38 at 100 PeV and 0.95 at 10 EeV. An estimate, not a fit — set the index to
zero for a constant.

#### The result

| E | β [10⁻⁶] | R_τ [g/cm²] | X_peak [g/cm²] | X_peak [km rock] |
| --- | --- | --- | --- | --- |
| 100 PeV | 0.38 | 1.1×10⁶ | 3.3×10⁶ | 12.5 |
| 1 EeV | 0.60 | 3.6×10⁶ | 5.7×10⁶ | 21.6 |
| 10 EeV | 0.95 | 5.1×10⁶ | 6.2×10⁶ | 23.5 |

The optimum **rises then flattens**, 12 to 23 km of standard rock — not the flat
18–30 km the erroneous version reported. It rises because the tau range grows
logarithmically, and flattens because β rises and tempers that growth.

Three consequences, and they are unchanged by the correction.

**(i) Topography sits at or below the optimum.** Arequipa's measured p90 column depth
is 4.9×10⁶ g/cm² (18.5 km of rock), giving **0.97–1.00 of peak yield across the whole
energy range**. Within what mountains can supply, more rock is never penalised, so the
depth criterion should be a rising ramp over the accessible range rather than a band.
The correction strengthened this: the optimum came *down* toward what terrain provides.

**(ii) The band is two and a half decades wide** at half maximum, 5×10⁵ to 6×10⁷ g/cm².
Column depth is an intrinsically weak discriminant and the criterion should not pretend
otherwise.

**(iii) The upper limit lives in the Earth chord, and narrows the window with energy.**
Setting `2R sin θ` equal to the upper band edge:

| E | chord [km] | elevation |
| --- | --- | --- |
| 100 PeV | 990 | −4.4° |
| 1 EeV | 452 | −2.0° |
| 10 EeV | 211 | −0.9° |

At 100 PeV the cut lies outside a ±3° window; by an EeV it has taken the bottom third,
by 10 EeV two thirds. **The effective arrival window is energy-dependent, its lower
edge climbing toward the horizon.** Worth checking against the collaboration's own
acceptance, since it predicts a specific narrowing.

#### Verification

Beyond the unit tests, the derivations are checked against limits and constructions
rather than against the code's own output: the numerical quadrature reproduces the
closed form when survival is made exponential (2×10⁻⁷ relative); survival reduces to
`exp(-d/X_decay)` as β→0; the decay length reproduces ref. [2]'s quoted 50 m – 5 km for
1–100 PeV; and the Earth chord is verified by direct construction, checking the
endpoint lands on the sphere.

One of those checks initially failed — and the check was wrong, not the code. It
compared the entry angle against `atan(sagitta/half-chord)`, which is `θ/2` rather than
`θ`. The chord formula `2R sin θ` is confirmed both by the tangent-chord angle theorem
and by explicit construction.

#### Remaining caveats

β is an estimate from scaling arguments, not a fit to published tables, and its energy
index is a guess consistent with photonuclear growth. Only charged-current attenuation
is counted; neutral-current regeneration would soften it. Neither affects (i).

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

### 5.3 First working TAMBO search, and GRAND+TAMBO combined ✅ delivered

The structural claim in §5 held up: nearly every TAMBO requirement turned out to be
expressible in configuration already — slope band, spacing, triangular grid, distance
window, arrival window, no Fresnel, no geomagnetic weighting, particle grammage mode.
One scan engine really does answer both. Only three code changes were needed.

**`min_width_km = 0` is now legal.** It was validated as strictly positive, but the
morphological opening it drives is precisely what deletes a strip-shaped array. A
"block-like array" is a GRAND assumption, not a general one.

**The shower band had to come from the energy** (§4.14 below). The default particle
band is (X_max, 4·X_max) = 700–2800 g/cm², but a canyon crossing supplies only the air
its own width contains. Every TAMBO candidate scored exactly zero until this was fixed.

**Two more capacity bugs, both found by cross-checking area against capacity.** See
§6.16; the second one inflated TAMBO's total by 38%.

Colca crop (1981 × 3061 at 1 arc-second, cut from the existing Arequipa DEM, which
already covers the canyon — verified 1673 m of incision against the published ~1.5 km):

| | area | sites | capacity | of its own area in the joint |
| --- | --- | --- | --- | --- |
| GRAND | 4580.2 km² | 1 | 5317 | 2.5% |
| TAMBO | 176.2 km² | 30 | 20385 | 65.3% |
| joint | 115.1 km² | | | Jaccard 0.025 |
| union | 4641.3 km² | | | |

**The co-location result is the interesting one.** Two thirds of TAMBO-viable ground is
also GRAND-viable, but that is only 2.5% of GRAND's, and the joint area is small in
absolute terms because the two deployable slope bands barely overlap — GRAND's 3–25°
against Colca's ~40° walls leaves only a 20–25° sliver. Co-location is possible but
marginal, and it is the *slope* criterion that decides, not the arrival geometry.

Visual check: the TAMBO selection traces the branching canyon network across the map
while GRAND takes the open plateau, which is what the physics should produce and is the
main evidence that the criteria are doing something real.

### 5.4 Per-role slope: the far wall ✅ delivered

The per-role criterion §5.2 called for. Slope was tested only at the candidate pixel —
the *near* wall the array stands on. The far wall, which is where the tau actually
exits, was never required to be a wall at all: the scan asked whether rock lay at the
right range and bearing, which on real Andean terrain is nearly always true somewhere.
That is why 92% of candidates passed before this existed.

The walk now records, per elevation bin, how fast the terrain was climbing along the
ray where it was first met — `dz/dd` between the previous sample and the intersection.
Measured along the arrival azimuth, so an obliquely-viewed wall counts as the tau would
actually cross it. It is reported as `target_slope_deg` and optionally bounded by
`min_target_slope_deg` / `max_target_slope_deg`, unset by default so GRAND is unchanged.

Verified against `synthetic.canyon`, whose wall slope is a fixture parameter. Filtered
to wall hits, the measured slope recovers the built value **exactly**:

| wall built | measured |
| --- | --- |
| 15° | 15.0° |
| 25° | 25.0° |
| 35° | 35.0° |
| 45° | 45.0° |

Two things the fixture made obvious and worth stating, because both broke a first draft
of the tests. The *unfiltered* mean is not the wall slope: rays aimed lower strike the
flat canyon floor, whose slope really is zero, so the mean over all accepted directions
is a mixture — filtering is what isolates the wall. And an *upper* bound does not empty
the result: a flat floor passes any ceiling, so a ceiling removes walls, not everything.

Effect at Colca, with a deliberately permissive 25° floor against ~40° walls:

| | before | after |
| --- | --- | --- |
| candidates accepted | 35.9% | **18.9%** |
| sites | 30 | 17 |
| capacity | 20385 | 10878 |
| usable area | 176.2 km² | 93.1 km² |
| joint with GRAND | 115.1 km² | 54.9 km² |

**The physical check that matters:** the surviving sites look at walls of 34.7–44.3°,
median **38.6°** — which is Colca's published wall steepness, recovered rather than
assumed. The criterion is selecting canyon-wall geometry, not terrain in general.

**Still open.** TAMBO's acceptance is a score cut (`min_score`), which is a knob rather
than a derivation; and the slope bands, the ±20° arrival window, the shower-band
fraction and the 25° far-wall floor are stated assumptions rather than collaboration
inputs. Every one of them is now a config knob, so re-running under different values
costs about six seconds.

---

## 6. Phase 3 — Performance (in progress)

Started ahead of phase 2, because at ~70 s per scan iterating on TAMBO criteria would
have been painful.

### 6.1 Delivered

**Separable morphology.** A rectangle of ones factorises into a column and a row, so
dilation or erosion by (h, w) is (h, 1) followed by (1, w) — O(N(h+w)) instead of
O(Nhw), and **bit-identical**, not an approximation. Verified against the direct
operation for several element shapes.

**The scan inner loop works in slope, not angle.** Every comparison the walk makes is
monotonic in elevation angle, so comparing `apparent/d` against pre-computed tangents
of the bin edges gives the same decisions without an arctangent per sample. At roughly
15 ns per sample the arctangent was about half the cost. The per-axis pixel steps are
hoisted out of the loop too.

Measured on the 2500² Arequipa crop, same configuration, 2 threads throughout:

| stage | before | after | |
| --- | --- | --- | --- |
| morphology | 10.50 s | **1.23 s** | 8.5× |
| arrival scan | 69.83 s | **52.46 s** | 1.33× |
| whole run | 82.0 s | **55.6 s** | 1.47× |

On all 12 cores — the CLI default, which the earlier figures did *not* use — the same
run is now **32.4 s**.

### 6.2 A bug the rewrite exposed

The original binned by `k = int((theta - elev_min)/bin)` and kept samples with `k >= 0`.
C-style `int()` truncates toward zero, so a sample up to one bin *below* `elev_min`
gave `k = 0` and was added to the lowest bin, inflating its column depth by as much as
a third. Working in slope against explicit tangent edges removes the boundary case
entirely. Selection barely moved at default settings — four candidates out of 69,000,
since nothing gates on depth unless a band or shielding floor is set — but the
*reported* depth was wrong, and would have propagated into any depth-gated run.

### 6.3 Balanced candidate ordering ✅ delivered

Numba's `prange` schedules statically, giving each thread one contiguous slice of the
index range. Candidates leave the topographic screen in *spatial* order, so that slice
is a contiguous patch of map — and walk cost varies enormously across a map, since rays
near an edge terminate early while interior ones run the full range. Some threads
therefore finished long before others: measured scaling was 2.4× on 12 cores against
4–5× for randomly scattered candidates.

Dealing *blocks* of neighbouring candidates round-robin keeps locality inside a block
while spreading each thread's slice across the whole map. Result-neutral, verified
against a single-threaded run.

| ordering | 8 threads |
| --- | --- |
| tile order | 31.3 s |
| blocks of 1024 | **23.3 s** |

Block size barely matters between 64 and 4096. Plain shuffling also balances, but
destroys locality and measured no better overall.

### 6.4 Two optimisations tried and rejected

Recorded because the negative results are as useful as the positive ones.

**Hoisting the per-sample division.** Most samples see terrain below the whole
acceptance window and need no slope at all, so the comparison can be made against
`tan(edge₀)·d` stepped incrementally. Tracking the horizon then needs an unreduced
fraction and two multiplies per sample. Net effect: **24.5 s → 26.1 s, slower.** The
division was already being pipelined; the extra multiplies were not free.

**Solving for the exit distance up front,** to remove the four bounds comparisons from
the inner loop. Worth about 5%, but it has to reproduce `int()`'s truncation exactly at
every edge, and the attempt altered 405 of 40,000 candidates — the bound is
`(cols − c0)/dc_per_m`, not `(cols − 1 − c0)/dc_per_m`. Silently truncating rays is not
worth 5%, and the branch predictor handles the test nearly as cheaply.

### 6.5 Where it stands

End-to-end on the 2500² Arequipa crop, 8 cores:

| stage | before phase 3 | now |
| --- | --- | --- |
| arrival scan | 69.8 s | **26.7 s** |
| morphology | 10.5 s | **1.7 s** |
| whole run | 82.0 s | **30.3 s** |

Extrapolating to the full Arequipa DEM, about 20× the pixels: roughly 10 minutes,
against something over an hour before.

### 6.6 The whole-raster sweep: premise measured, and wrong

This was billed as the main remaining structural win, on the reasoning that for
azimuths near north or south consecutive samples sit one DEM row apart — 10 KB — so
essentially every sample is a cache miss. **Measured, it is not true.** Holding
everything else fixed and varying only the azimuth:

| azimuth | row stride | seconds |
| --- | --- | --- |
| 90° (due east) | 0.00 | 2.27 |
| 45° | 0.71 | 2.29 |
| 0° (due north) | 1.00 | 2.26 |
| 270° | 0.00 | 2.01 |

Worst over best is **1.14×**. The hardware prefetcher handles a constant row stride
perfectly well, and with eight threads the memory-level parallelism hides what latency
remains. An earlier run appeared to show a large effect; that measurement had JIT
compilation folded into its first case.

So the sweep would buy at most ~14%, for a large rewrite which — as §4.2.2 already
noted — does not extend cleanly to a (θ, φ) scan carrying sub-surface chords, because
the running-max construction it relies on gives the horizon but not the per-bin first
crossings or the depth histogram. **Not done, and not recommended** unless profiling
on much larger DEMs contradicts this.

### 6.7 Bilinear profile sampling ✅ delivered

Nearest-neighbour sampling quantises the profile to pixel centres and treats terrain
as piecewise constant, so a ray is blocked by a whole pixel's worth of the nearby
maximum. It also samples through `int()`, which truncates toward zero and therefore
biases the sample point back toward the candidate by up to half a pixel —
asymmetrically, since the sign of the offset depends on the azimuth.

Interpolating removes both. It is exact on a plane, since bilinear reproduces linear
functions, which is how the tests check it.

The effect is not small. On the Arequipa crop, acceptance rises 13.4%, and 9.8% of
candidates change acceptance outright; site area rises 6.9% and capacity 5.5%. That is
the correction of a systematic pessimism, not noise. Cost is 1.44× on the scan.
`--nearest_sampling` restores the old behaviour.

### 6.8 Streaming DEM cache ✅ delivered

The cache was built with `tiff.imread(path).astype(np.float32)`, which materialises
the whole DEM and then a second full copy of it — defeating the out-of-core design the
rest of the pipeline rests on, and failing outright on the multi-gigabyte DEMs the
README targets. The page is now decoded straight into a native-dtype file and
converted a block of rows at a time.

Non-evictable (anonymous) memory while building the full Arequipa cache:

| | anonymous RSS |
| --- | --- |
| in-RAM conversion | 623 MiB |
| streaming | **132 MiB** |

and the streaming figure is bounded by the block size rather than the DEM, so a 20 GB
DEM stays at ~130 MiB instead of needing some 30 GB. Verified byte-identical to the
old conversion, and independent of block size.

Storing float32 with NaN rather than the DEM's own int16 is deliberate: NaN propagates
through the gradient and comparison chain in the screening stage, so nodata is excluded
without a sentinel test in every kernel. That costs twice the disk of an int16 cache
and buys correctness that would otherwise have to be re-established in half a dozen
places.

### 6.9 Topographic screen in gradient space ✅ delivered

The screen computed `arctan`, `sqrt` and `arctan2` over **every** pixel of every tile,
then used the slope only for a band comparison and read the aspect only at the pixels
that survived.

Both are avoidable. Slope rises monotonically with gradient magnitude, so

    min ≤ atan(√g) ≤ max   ⇔   tan²(min) ≤ g ≤ tan²(max)

which needs neither `sqrt` nor `arctan` (`slope_band_gradient_sq`, with bounds at or
beyond the vertical returned as "unbounded" since `tan` is singular there). And once
slope and altitude have reduced the tile to a subset, every remaining filter — aspect,
road distance, RFI circles and polygons — can work on that subset instead of on the
full tile. `arctan2` is then evaluated only where the aspect is actually read.

Measured on a 2500² Arequipa crop, against six filter configurations, all producing
**byte-identical candidate arrays**:

| configuration | before | after | |
| --- | --- | --- | --- |
| slope band only | 0.261 s | 0.209 s | 1.25× |
| aspect bounds | 0.218 s | 0.184 s | 1.18× |
| aspect bounds, wrapped | 0.218 s | 0.184 s | 1.19× |
| RFI circle + polygon | 0.334 s | 0.248 s | 1.35× |
| altitude band | 0.231 s | 0.119 s | 1.94× |
| narrow slope band | 0.201 s | 0.045 s | **4.47×** |

The gain scales with how much the filters reject, because the work now follows the
survivors rather than the tile. This also removed a redundant full-array pass in the
RFI branch, which computed `np.where(mask.ravel())` only to test whether it was empty.

### 6.10 Per-site passes made O(N) ✅ delivered, and a crash fixed

`analyze_sites_and_capacity` scanned the whole downsampled map once per site
(`labeled == site_id` for the mean aspect, then `labeled == original_id` per selected
site for the colouring, then `np.isin` for the mask), and
`summarize_observables_by_site` masked the whole accepted array once per site. All are
O(sites × pixels).

Replaced by one labelled pass (`scipy.ndimage.mean` over all candidate regions), one
lookup-table recolouring, and one stable `argsort` grouping. The sort is stable, so
each site's values keep the order a boolean mask would have produced and the
statistics are unchanged to the last bit.

On a 2000² map, `analyze_sites_and_capacity` / `summarize_observables_by_site`:

| sites | before | after | |
| --- | --- | --- | --- |
| 16 | 0.280 s | 0.226 s | 1.24× |
| 100 | 0.776 s | 0.223 s | 3.48× |
| 256 | *crashed* | 0.245 s | — |
| 900 | *crashed* | 0.232 s | — |

The new cost is flat in the site count, as it should be.

**The crash was real and pre-existing.** `labeled_viz` was `uint8`, so assigning the
256th site's colour raised `OverflowError` and took the run down *after* the physics
had been paid for. It is now sized from the label count via `np.min_scalar_type`.
This matters for phase 2: TAMBO needs a long strip of many small sub-arrays rather
than one blob (§5.2), which reaches 255 sites routinely. Pinned by
`tests/test_capacity.py::TestManySites`, which fails on the old code.

### 6.11 The scan is not compute-bound — two more arithmetic optimisations, both neutral

Recorded at length because the conclusion redirects any future work on the scan.

**Per-bin transcendentals hoisted out of the inner loop.** `cos(θ)`, `sin(θ)`,
`tan(θ)`, the solid-angle element and the Earth-chord term depend only on the
elevation bin, yet were recomputed for every (candidate × azimuth × bin); `cos(θ)` was
computed *before* the acceptance guards, so even rejected bins paid for it.
`exp(−z₀/H)` depends only on the candidate. `sin`/`cos` of the azimuth were recomputed
once per bin inside the geomagnetic term. Tabulating all of them is bit-identical and
measured **0.975×, i.e. 2.6% slower** — nothing, within noise.

**The histogram's bin search short-circuited.** A sample above the top bin edge fell
through the full linear scan over all 13 edges without ever breaking — the most common
case was the most expensive. Adding one comparison to catch it is bit-identical and
measured **0.992×**. Also nothing.

**Why.** The estimate that motivated both was wrong by more than two orders of
magnitude. Per (candidate, azimuth) the bin loop runs `n_bins` = 12 times, but the
profile *walk* runs `max_range/step` ≈ 2700 samples — a ratio of 225:1. Everything
outside the walk is under 1% of the kernel, so no amount of arithmetic saved there is
visible. Both changes were reverted: they are bit-identical and slightly slower, and
neutral complexity in the most delicate loop in the project is a net loss.

**What this means.** Arithmetic micro-optimisation of the scan is exhausted. Four
attempts have now failed for the same underlying reason: the two in §6.4 and the two
here. Anything further has to reduce the *number of samples* or the *memory traffic
per sample*, not the flops per sample.

### 6.12 What the machine contributes, and why timings here are noisy

Worth writing down, because two days of measurement were nearly wasted on it.

The workstation is a 13th-gen i5-1334U: a **hybrid** CPU with 2 P-cores (HT, 4.6 GHz
nominal, CPUs 0–3) and 8 E-cores (3.4 GHz, CPUs 4–11), 12 MiB of L3, and in practice
running at about 1.5 GHz under sustained load. Consequences:

- **Scaling saturates on the hardware, not the code.** Measured on 151k candidates:
  1.85× at 2 threads (92% efficient, both P-cores), 2.51× at 4 (the P-cores' HT
  ceiling), **3.70× at 8** (46%). The four E-cores added at the end contribute about
  0.3 of a P-core each.
- **Scheduling is not the problem.** `numba.set_parallel_chunksize` from the default
  static split down to 16 moved the 8-thread time only from 8.44 s to 8.17 s (~3%),
  so the existing block-dealt ordering (§6.3) is already doing its job.
- **Wall-clock A/B across processes is untrustworthy here.** The same unchanged code
  measured 43.6 s and 39.8 s in consecutive runs. An early cross-process comparison
  showed the §6.11 hoist at 1.15× — it is actually 0.975×. Measure A and B
  *alternating inside one process*, and prefer **single-threaded** runs on a
  subsample: one thread on a 12-thread box is essentially never descheduled.
- The 2500² crop is 25 MB of float32 against 12 MiB of L3, so it does not fit; the
  full DEM is ~20× larger again.

`bench/baseline.json` records `host.load_average_1min` for exactly this reason.

### 6.13 The Fresnel pass costs less than assumed

Suspected to be "a plausible chunk of the 1.44× bilinear cost". Measured on the 2500²
crop with the real Arequipa configuration: **43.55 s without it, 46.57 s with it** —
6.5% of the scan, not a chunk of it. Acceptance averages 14.7 accepted (azimuth,
elevation) cells per accepted candidate, and the second walk starts at the near-field
cut-off and stops at the shower offset, so it is much shorter than the main walk.
Fusing it into the main walk was therefore not attempted: the ceiling is 6.5% and the
change would entangle two passes with different termination conditions.

### 6.14 Where phase 3 stands

End-to-end on the 2500² Arequipa crop, 8 cores, default benchmark configuration:

| stage | before phase 3 | now |
| --- | --- | --- |
| arrival scan | 69.8 s | 15.3 s |
| topographic screen | 1.4 s | **0.25 s** |
| morphology | 10.5 s | **0.62 s** |
| capacity analysis | 0.6 s | **0.22 s** |

With the full physics configuration the scan is the whole story: 43.6 s of a 47.9 s
run, and §6.11 says that number is now bounded by memory traffic and by this machine's
cores rather than by arithmetic.

**Remaining ideas, in the order they are worth trying:**

1. **Fewer samples per walk.** The only lever with real headroom left. An early exit
   is not free — the depth histogram accumulates over the whole path and `horizon_deg`
   is a reported observable — so it needs an explicit decision about what those two
   mean when the walk stops early. See lead (k) in the handover.
2. **Memory traffic per sample.** Bilinear sampling touches two rows per sample rather
   than one; the DEM does not fit in L3. A tiled traversal that walks all azimuths of a
   block of candidates together is the shape of the idea, but note §6.6 measured the
   azimuth locality penalty at only 1.14×, which caps this.
3. `build_elevation_cache` still writes a raw file and then converts it (lead (j)). It
   is one-off per DEM and `load_dem` is 0.05 s on the crop, so this is about first-run
   latency on the full DEM, not throughput.

Not worth doing: multiprocessing (lead (i)) — screening and morphology are now 3% of
runtime combined, and the scan already uses every core the owner allows.

### 6.15 Capacity from metric grid placement ✅ delivered

The phase 2 blocker of §2.1. `count_grid_capacity` used to convert the ground spacing
into an integer pixel stride and step the array by it, truncating three times —
`spacing_r`, `spacing_c`, and the hex row pitch `int(spacing_r · sin60)`. Every
truncation shortens the spacing, so the count came out high, and worse as the spacing
approached the pixel size.

It now takes the pixel sizes and the spacing **in metres**, places positions in
continuous ground coordinates, and looks each one up in the pixel grid. There is no
stride left to truncate.

| requested spacing | old / analytic | new / analytic |
| --- | --- | --- |
| 1000 m (GRAND) | 1.074 | 1.00 |
| 150 m (TAMBO, published) | 1.423 | 1.00 |
| 100 m (TAMBO, starting) | 1.581 | 1.00 |
| 60 m | — | 1.00 |

Effect on the golden files, both at 1 km spacing: Arequipa crop **806 → 747 DUs
(−7.3%)**, synthetic ridge **611 → 559 (−8.5%)**. Site counts and areas are unchanged;
only the packing density moved. Regenerated deliberately.

The layout is still anchored at each site's bounding-box corner rather than fitted to
it, so this remains a capacity estimate for an arbitrarily-placed array rather than
the best achievable packing. A spacing finer than the DEM's pixels is now permitted
and yields several detectors per pixel — the honest continuum limit, though the
terrain mask cannot resolve whether those sub-pixel positions are usable.

`tests/test_capacity.py::TestDensityMatchesTheRequestedSpacing` replaces the
characterization tests that pinned the defect; the old ratios are now what a
regression looks like.

### 6.16 A second capacity bug: the bounding box is not the region ✅ fixed

Found by checking an invariant rather than by reading code. Capacity and area are
computed by different routes, so they must agree: `capacity × area-per-detector` should
reconcile with the summed site area. For GRAND it did (0.99). For TAMBO it did not
(0.725), and that gap was the bug.

`count_grid_capacity` was handed `final_map_disk[bbox]` — the region's bounding **box**,
not the region. A box also contains whatever else falls inside it: other sites, and
regions that failed the area threshold. Their pixels were counted as this site's, and
because the totals are summed over sites the same ground was sold more than once. One
compact site barely notices, which is why GRAND looked fine and the goldens moved by
only 2 DUs. A canyon network of thirty interleaved strips inflated its total by **38%**
(28054 → 20385). A synthetic L whose box spans the map over-counted by **2.07×**.

Restricting the count to `labeled[loc] == site_id` fixes it; both experiments now
reconcile to within 0.5%. Pinned by
`test_capacity.py::TestManySites::test_a_site_holds_no_more_detectors_than_its_own_area_allows`,
which states the physical invariant directly — no site holds more detectors than its
own area divides into — and fails on the previous code.

Worth noting the general lesson: **the area/capacity cross-check is a cheap invariant
that caught a bug neither the golden files nor 264 tests had noticed**, because both
quantities were wrong in the same direction only when many sites overlapped.

A related bias remains documented rather than fixed: per-site `area_km2` is measured on
the *downsampled* map while capacity is measured at full resolution, so at
`downsample_factor > 1` a feature only a few pixels wide loses area it keeps detectors
on. Both Colca configs therefore run at `downsample_factor: 1`.

### 6.17 What closing contributes, measured with a stride-1 run ✅ measured

§2.1 has said since the physics review that "closing, not the physics, determines most
of the reported site extent", but the size of the effect was never separated from the
`candidate_stride` sampling that closing partly compensates for. A stride-1 GRAND run
over the Colca crop settles it.

| | stride 5 | stride 1 |
| --- | --- | --- |
| candidates screened | 770,652 | 3,853,258 |
| directions accepted | 463,326 | 2,315,998 |
| **acceptance** | **60.1%** | **60.1%** |
| after closing | 5,080,873 | 5,296,905 |
| reported area | 4580 km² | 4809 km² |
| capacity | 5317 | 5610 |

Two results.

**Striding is unbiased.** The acceptance fraction is identical to three figures, and
the stride-corrected area estimate (accepted × 5) comes to 2120 km² against the
stride-1 truth of 2119 km² — agreement to 0.05%. `candidate_stride: 5` costs nothing
in accuracy and saves five sixths of the scan. The final reported figures differ by
only 5%, so stride-5 runs can be quoted directly.

**Closing inflates the area by 2.29×**, measured at stride 1 where no reconstruction is
needed. The 11× that a naive comparison suggests is an artefact of comparing against
un-reconstructed set pixels; the honest figure is 2.29×, and the honest physics-accepted
area at Colca is **~2120 km², not the reported 4580 km²**.

The closing element was tied to `antenna_spacing_km`, coupling two unrelated things and
hiding the effect. It is now `gap_close_km`, defaulting to the old behaviour, with 0 to
disable. Reproduced by `config/grand_colca_stride1.json`.

### 6.18 Tau decay probability in the score ✅ delivered

`decay_probability` existed but only `aperture.py` used it, so the per-candidate score
never asked whether the tau decays in the gap at all. GRAND gets this implicitly, since
its distance window is derived from the decay length. A canyon search does not: TAMBO's
window comes from the terrain.

The term is `1 − exp(−(d − shower_development)/L)` on the mean exit distance, present
only when `decay_energy_pev` is supplied — the probability is strongly energy-dependent
and one number cannot stand in for a spectrum:

| energy | decay length | P(decay within 3 km) |
| --- | --- | --- |
| 3 PeV | 147 m | 1.000 |
| 50 PeV | 2449 m | 0.706 |
| 100 PeV | 4898 m | 0.458 |
| 1 EeV | 48980 m | 0.059 |

At Colca with 55 PeV, the geometric midpoint of TAMBO's reach, this is the strongest
cut of any single criterion: acceptance 18.9% → **7.9%**, sites 17 → 5, capacity
10878 → **2056**, area 93.1 → 18.0 km².

**That result is worth stating plainly: under the current assumptions Colca supplies
about 2056 detector positions against TAMBO's 5000.** It is also the number most
sensitive to an assumption — a spectrum-folded treatment would replace the single
energy, and `min_score`, the far-wall floor and the shower band all move it. The right
next step is a sensitivity table over those four knobs rather than another criterion.

Note `shower_development_m` is set to 0 for TAMBO deliberately: the grammage band
already carries the requirement that enough air was traversed, and subtracting it here
too would count the same constraint twice.

### 6.19 Command line no longer loses to the config file ✅ fixed

The merge was `config > fallback > CLI`, and since `--generate_config` writes all 67
keys, a generated config made **every** command-line flag a silent no-op — no warning,
no error. The cause was that argparse cannot distinguish `--candidate_stride 5` from
its own default of 5, so honouring the command line would have let every default
overwrite the config.

`explicitly_passed()` answers that directly by re-parsing with `argparse.SUPPRESS` as
every default, so only options that actually appeared get set. An explicitly typed
option now wins and says so:

    Command line overrides config for 'min_slope_deg': 3.0 -> 44.0

Pinned by `test_screening.py::TestExplicitCommandLineDetection`, including the case
that made the old merge unfixable — a typed value that equals the default.

## Phase 4 — Usability *(sketch — to be scoped)*

Auto-detect `origin_lat`/`origin_lon` from the GeoTIFF tiepoint (verified present,
matching current configs to ~1e-4°); rename `src/setup.py`, which is not a packaging
file and whose name hijacks `pip install`; real packaging; rasterio/pyproj for CRS and
outputs; `--explain` funnel report; parameter sweeps.

## 7. Open questions

1. **~~A column-depth band for GRAND above 100 PeV~~** — estimated, §4.13.
2. **~~Score composition~~** — implemented; `product` is the default.
3. **~~TAMBO muon shielding~~** — implemented as `--muon_shielding_km`, default off.
4. **~~A published aperture curve as data~~** — both digitized into `data/`.
5. **IGRF field values per site.** Inclination now follows the DEM's own coordinates;
   declination still falls back to the Arequipa value and should be supplied per site.
6. **~~The tau survival model~~** — corrected in §4.13. Energy loss now shortens the
   decay length as the tau propagates, giving a double-exponential survival and a range
   that grows logarithmically rather than saturating.
7. **Beta, the tau energy-loss constant.** Now estimated at (0.4–1.0)×10⁻⁶ cm²/g from
   mass scaling, with an assumed energy dependence. Still worth pinning to whatever
   value or tabulation the collaboration uses; it moves the optimum in proportion,
   though not the siting conclusion.
8. **Neutral-current regeneration** is not modelled — only charged-current attenuation
   is counted, so the Earth-chord suppression of §4.13(iii) is somewhat overstated.

---

## References

[1] GRAND Collaboration (J. Alvarez-Muñiz et al.), *The Giant Radio Array for
Neutrino Detection (GRAND): Science and Design*, arXiv:1810.09994 (2018);
Sci. China Phys. Mech. Astron. (2020). Effective area: Fig. 25.
20 sub-arrays × 10,000 antennas over 10,000 km² each, E > 10⁸ GeV.

[2] TAMBO Collaboration, *Measuring the high-energy neutrino sky using the deep-valley
neutrino observatory TAMBO*, Nature Astronomy (2026),
doi:10.1038/s41550-026-02916-4. Exposure: Fig. 3.
