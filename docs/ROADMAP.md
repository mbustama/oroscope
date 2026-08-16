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

`src/oroscope/arrival_scan.py`, with 25 tests against terrain whose answer is known in closed
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

`src/oroscope/scoring.py`. Every component returns [0, 1] with a documented shape, and the
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

`src/oroscope/aperture.py` separates the estimate into three parts: the **geometric aperture**
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
`src/oroscope/physics.py` and the scan kernel, with 44 tests. All the analytic pieces are
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

### 6.20 Sensitivity of the TAMBO result ⚠️ the result is not robust

`src/oroscope/sensitivity.py` varies one parameter at a time about the Colca baseline. The
answer is that **2056 detector positions is not a number to quote.** Every criterion
sits near a cliff:

| parameter | value | capacity | vs baseline |
| --- | --- | --- | --- |
| `decay_energy_pev` | 3 | 10878 | 5.29× |
| | 10 | 10857 | 5.28× |
| | **55** | **2056** | 1.00× |
| | 100 | **0** | 0.00× |
| | 1000 | **0** | 0.00× |
| `min_score` | 0.0 | 45928 | 22.3× |
| | 0.2 | 15481 | 7.53× |
| | **0.35** | **2056** | 1.00× |
| | 0.5 | **0** | 0.00× |
| `min_target_slope_deg` | 0 | 7442 | 3.62× |
| | 15 | 5309 | 2.58× |
| | **25** | **2056** | 1.00× |
| | 35 | **0** | 0.00× |
| `grammage_band_fraction` | 0.05 / 0.1 / 0.2 | 2056 | 1.00× (inert) |

**The decay energy is the worst.** Across TAMBO's own 3 PeV – 1 EeV reach the answer
runs from 10878 to zero. A single energy standing in for a spectrum is not an
approximation here, it is the whole answer, and §6.18's number is an artefact of
picking 55 PeV. **This has to be folded over the spectrum before any TAMBO capacity is
quoted.**

**`min_score` does most of the remaining work,** and for a structural reason: the score
is a *product* of six components each in [0, 1], so it concentrates near zero and a
threshold anywhere in the middle sits on a cliff. Ranking sites and taking the best N
would be better behaved than thresholding a product; a weighted geometric mean would
also spread the distribution.

**`grammage_band_fraction` is inert** because the config sets `grammage_band_gcm2`
explicitly, and an explicit band correctly wins over a derived one. Correct, but a trap
worth a warning: setting the band silently disables the fraction.

The honest summary is that the pipeline now expresses TAMBO's geometry correctly — the
far-wall slopes it recovers are Colca's real ones — but the *capacity* it reports is
dominated by two modelling choices rather than by terrain. Fixing the decay treatment
is worth more than any further criterion.

### 6.21 Tau decay folded over a spectrum ✅ delivered

§6.20 found that the single-energy decay term *was* the answer rather than an
approximation to it: across TAMBO's own 3 PeV – 1 EeV reach the reported capacity ran
from 10878 to zero. `physics.spectrum_weighted_decay_probability` folds it instead,

    P(u) = ∫ E^-γ (1 − exp(−u/L(E))) dE  /  ∫ E^-γ dE

on a log-spaced grid, chunked over candidates so a 500k-candidate search costs about
1.1 s and tens of megabytes rather than the hundreds the full outer product would take.

**It worked.** The result is now robust where it was not:

| | single energy | folded (γ swept) |
| --- | --- | --- |
| range of the assumption | 3 – 1000 PeV | γ = 1.5 – 2.7 |
| capacity | 10878 → **0** | 7205 → **10495** |
| ratio across the range | ∞ | **1.46×** |

At the Colca baseline, γ = 2.0 over 3 PeV – 1 EeV: **15 sites, 9717 detector
positions, 17.5% acceptance.** The energy-range endpoints matter little (1.07× and
0.99× for the low and high edges), which is the point — the flux weighting is doing the
work rather than the choice of endpoint.

`min_score` is now the dominant remaining assumption, at 2.38× down to 0.20× across
0.2–0.5, and it is a threshold on a product. That remains the thing to fix next, and
§6.20's recommendation stands: rank sites and take the best N.

### 6.22 A memory leak, and safeguards against the next one ⚠️ found the hard way

Running the sweep killed the machine. The kernel's OOM killer took the process at
**6.9 GB anon RSS**, on a box with 15 GB where 9 were already in use.

**The leak was ours.** `generate_visualizations_and_outputs` created a 14×12-inch
figure and saved it without closing it. pyplot holds a global reference to every figure
it creates, so an unclosed one can never be collected — a single search never notices,
but a process that runs several accumulates the entire figure, canvas and artists, each
time. Reproduced in isolation at ~14 MB per iteration on a small image, far more on the
real one. Fixed with `plt.close('all')` in a `finally`, since the failure path leaked
as readily as the happy one.

**Three safeguards, because the leak will not be the last one.**

- `estimate_peak_memory_gb` predicts the anonymous allocations from the DEM size,
  `downsample_factor` and `candidate_stride`, and the run prints it against available
  memory and warns past 80%. The memory-mapped DEM is deliberately excluded: it is
  file-backed and evictable, and counting it would make every large search look
  impossible when the streaming design exists so that it is not.
- `apply_memory_cap` bounds `RLIMIT_AS`, defaulting to 80% of available. A search that
  outgrows the machine now fails with `MemoryError` naming itself, rather than letting
  the kernel pick a victim that may be the user's editor.
- **`sensitivity.py` runs each point in a subprocess.** Memory is reclaimed completely
  between points, and one failed point reports a failed row instead of ending the
  sweep. It costs a few seconds of JIT per point, which is the right trade for a sweep
  that otherwise cannot finish.

Measured after: the same sweep peaks at **1.2 GB**, against 6.9 GB before.

For the full DEM the estimator says 7.2 GiB at `downsample_factor: 1` and 5.1 GiB at 4,
against ~6–7 GiB typically available on this machine — so 4 is the right setting, which
is what §6.12 already recommended for a different reason.

**These two numbers were 4.5 and 2.3 until the full DEM was actually run**, and the
under-estimate cost that run 23 minutes before it hit its own cap. The estimator counted
only the arrays the scan returns, not the roughly three times as many live inside
`scoring.compose`, which is where the high-water mark really is. See §6.26a — including
why `downsample_factor` is the weaker of the two memory levers, contrary to what the
rest of this section implies.

### 6.23 The run explains itself ✅ delivered

Everything needed was already in the results JSON; what was missing was the reading of
it. `explain.explain_results(results)` takes the results dictionary and returns a
string — no files, no DEM, nothing re-run — so the pipeline, the library and a test all
get the same words, and a run from months ago can still be explained from its JSON.

**On by default**, per the owner; `--no_explain` suppresses it. Printed last, so it is
what a reader is left with, and saved as `explanation.txt` beside the results, because
these runs are meant to be handed to other people and a terminal scrollback is not.

It says four things:

- **Which constraint bound.** The funnel stage that removed the largest share of what
  reached it, plus the parameter behind it. Two stages are excluded by construction:
  `kept by stride N` is a deliberate subsample and would otherwise be named on nearly
  every run, and `after gap closing` adds pixels. A stage that leaves *nothing* wins
  outright over any ratio downstream of it.
- **What held each site back.** The score components were already named but were being
  dropped by `summarize_observables_by_site`, which kept a fixed field list and so
  stored the total and nothing else. Now stored, and the lowest median component is
  reported per site — under a product composition it bounds the total from above.
- **How much closing moved the area**, measured from this run rather than quoted.
- **Which numbers are assumptions**, with the sensitivity measured in §6.20 for each.

**Two things came out of running it on the existing results.**

*The attribution is unambiguous, and the same on every site.* `solid_angle` is the
weakest component at 15 of 15 TAMBO sites (median 0.57, everything else at 1.0 except
decay at 0.96), and at GRAND's single site the ranking is `solid_angle` 0.52,
`footprint` 0.61, `geomagnetic` 0.74, everything else 1.0. So the TAMBO result is set
almost entirely by `solid_angle_half_sr`, whose default 0.05 sr is flagged in
`assumptions.rst` as a GRAND-scale value. It does not saturate here, but it is what the
score is measuring.

*The closing factor is measurable per run, and disagrees with itself between the two
configs.* Closed pixels over stride-corrected accepted pixels gives **2.19× for GRAND**,
against the **2.29×** the stride-1 control measured in §6.11 — an independent
cross-check of that number, from a different quantity, agreeing to 4%. But **TAMBO gives
0.53×**: its closing element is 100 m, about three pixels, too small to bridge the gaps
`candidate_stride: 5` leaves, so its reported area *understates* the accepted set rather
than inflating it. ⚠️ **The §6.11 conclusion that striding is unbiased was measured at
GRAND's 1 km element and does not transfer to TAMBO's 100 m.** A stride-1 control at
TAMBO settings would settle it; until then TAMBO areas are a lower bound, not an upper
one, and the summary says so.

### 6.23b `oroscope-combine` was reading a stale mask ⚠️ found by accident, numbers changed

Found while checking that the run summary's area agreed with the combined report's: it
did not. The summary said TAMBO covered 83.6 km², the combined report said 44.5.

`combine_experiments.load_run` took `sorted(glob("*.tif"))[0]`. A directory re-run
since the rename holds both `oroscope_results_*.tif` and a stale
`grand_search_results_*.tif`, and **the legacy prefix sorts first** — so the overlay
used a superseded 48,663-pixel TAMBO mask against a current GRAND one. Nothing failed
and nothing warned; the report simply described a run that no longer existed. Now it
selects by prefix, current first, falling back to the legacy name and then to any
`.tif`, so a pre-rename directory still loads. `tests/test_combine.py` pins it, and
fails against the old behaviour.

**The Colca combination numbers change. The §5.1 table in the handover brief is wrong
for TAMBO:**

| | before (stale mask) | corrected |
| --- | --- | --- |
| TAMBO area | 44.5 km² | **83.6 km²** |
| joint (GRAND & TAMBO) | 26.4 km² | **50.1 km²** |
| union | 4598.3 km² | **4613.7 km²** |
| TAMBO's own area in the joint | 59.3% | **59.9%** |
| GRAND's own area in the joint | 0.6% | **1.1%** |

GRAND's own figures are unchanged (4580.2 km², 1 site, 5317) — its stale file happened
to hold an identical mask, which is why the fault stayed hidden. The conclusion is
unchanged too: co-location is still decided by slope, and the shared sliver is still a
percent of GRAND's area. But the shared area is nearly double what was reported.

### 6.23c The site list was longer than the result ⚠️ same check, second fault

The area check of §6.23b, repeated once `oroscope-combine` was fixed, found a second
disagreement — this one in the results file itself.

`analyze_sites_and_capacity` returns **every** site that cleared the area and capacity
thresholds, and `stop_at_target` then selects a prefix of that capacity-sorted list.
Only the selection reaches `total_sites`, `total_capacity` and the exported raster; the
full list is what goes into the JSON, with nothing marking the difference. So anything
totalling `sites` over-reported. Measured on a synthetic run at `target_antennas: 50`:

| | listed | selected |
| --- | --- | --- |
| sites | 2 | **1** |
| area | 243.9 km² | **215.7 km²** |
| capacity | 288 | **252** |

The first summary written in §6.23 said "2 sites, 243.9 km²" against a raster holding
one site of 215.7 — which is exactly the class of error the summary exists to prevent,
so it is worth being blunt about having made it.

Each record now carries **`selected`**, and `explain.selected_sites()` splits the list
(falling back, for files written before the flag, to the first `total_sites` entries —
exact, since selection walks the sorted list in order). The summary counts and sums the
selection, and reports the rest as what they are: the next best ground, not ground that
failed.

**The invariant, now tested:** the area the summary adds up equals the area of the
raster the run wrote, to 0.001% — checked in plain, `stop_at_target`, downsampled and
single modes. Both halves are pinned: removing the flag fails four subtests, and making
the summary total the raw list fails a fifth.

### 6.24 CLI/library parity ✅ delivered

Measured, then closed. The pipeline function is now what the command line calls, not a
subset of it.

- **`max_memory_gb` was a real gap** — applied only in `main()`, so the caller most
  likely to need it, a sweep, did not get it. The estimate, the warning and the cap are
  now `preflight_memory()`, called by the pipeline itself. `sensitivity.py` passes the
  ceiling as an ordinary parameter instead of the child re-applying it by hand.
- **`load_config`, `generate_config` and `default_config`** are ordinary functions. The
  67-key template used to be a literal inside `main()`, reachable only by running the
  CLI with `--generate_config`. A test now asserts the template names every parameter
  the pipeline accepts, since a template with holes falls back silently.
- **The pipeline returns its results dictionary.** It returned `None`, so every caller
  re-read the JSON it had just written; `tests/_support.py` no longer does. A run that
  finds no candidate at all now returns its funnel too — the case where the funnel is
  the entire answer. (`sensitivity.py` still reads the file: its caller is in another
  process, which is the point of running each point in one.)
- The origin-resolution print was already in the library, contrary to the handover note.

### 6.25 What the summary cost, measured

Only one stage does more work than before: `summarize_observables_by_site` now folds
the eight named score components as well as the twelve geometric observables, so a site
record carries 34 fields plus 3 per component.

Measured the way §6.12 says to — A/B alternating in one process, 15 pairs, 200k
candidates over 15 sites, medians: **76.6 ms without the components, 124.2 ms with,
1.62×**, or **+48 ms per run**. Against a real search's 19 s that is a quarter of a
percent, and it buys the attribution that made §6.23's two findings visible. Composing
and writing the summary itself is not measurable — it is string formatting over a
dictionary that is already in memory.

**`bench/baseline.json` is deliberately NOT refreshed.** The machine could not resolve
the difference: two consecutive benchmark passes over *identical* code reported
`synthetic_1800/ray_tracing` at +3.4% and then +72.1% against the same baseline, and
`arequipa_2500/ray_tracing` at +16.3% and +34.6%. Nothing in the scan path changed —
`arrival_scan.py`, `scoring.py` and `physics.py` are untouched by this work — so those
are noise, exactly as §6.12 warns. Updating the baseline from a run like that would
bake the noise in. Refresh it on a quiet machine, and expect `capacity_analysis` to be
the only stage that legitimately moved.

### 6.26 The full Arequipa DEM ✅ run

Run on 2026-08-16 at head `1fa8810`. Every number published before this came from
crops; this is the first search over the whole 10204 × 12603 DEM — 128.6 Mpx, a
117,430 km² footprint, 21.2× the Colca crop. Both configurations are the crop ones with
the three intended changes only (the full DEM, a null origin read from the tiepoint,
`downsample_factor: 4`); no criterion was touched, so the comparison below is like for
like.

| | area km² | sites | capacity | bound by | kept | closing | weakest |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GRAND, Colca crop | 4,580.2 | 1 | 5,317 | directions accepted | 60.1% | 2.19× | solid_angle 1/1 |
| **GRAND, full DEM** | **88,527.5** | **1** | **101,948** | directions accepted | **61.6%** | **2.10×** | solid_angle 1/1 |
| TAMBO, Colca crop | 83.6 | 15 | 9,717 | directions accepted | 17.5% | 0.53× | solid_angle 15/15 |
| **TAMBO, full DEM** | **111.9** | **26** | **9,024** | directions accepted | **9.7%** | **0.52×** | solid_angle 26/26 |
| joint, crop | 50.1 | | | | | | Jaccard 0.0109 |
| **joint, full DEM** | **50.2** | | | | | | **Jaccard 0.0006** |

**The four questions, answered.**

1. **The binding constraint is the same one the crops found** — `directions accepted`,
   for both experiments, as it was at Colca. So the crops were not misleading about
   *what* limits the answer, and the numbers derived from them do not need re-reading.
   But the *rate* differs sharply by experiment, and that is the real finding: GRAND
   keeps 61.6% at full scale against the crop's 60.1% — indistinguishable — while TAMBO
   keeps **9.7% against 17.5%**, barely half. Colca is ordinary ground for GRAND and
   exceptional ground for TAMBO, which is exactly what one would expect of a crop chosen
   for containing a canyon, and is now measured rather than suspected.
2. **The area.** GRAND: 88,527.5 km², against 96,946 km² for the crop scaled linearly by
   footprint — **0.91× the naive scale-up**, so the crop was very slightly better than
   typical and the linear extrapolation was close to right. It accepts 75.4% of the DEM
   against the crop's 82.5%. TAMBO: 111.9 km² against a naive 1,769 km² — **0.06×**, an
   order of magnitude short. Each run's own closing factor confirms both crop
   measurements: GRAND 2.10× (crop 2.19×, stride-1 control 2.29×), giving ~42,190 km² of
   physics-accepted ground; TAMBO **0.52×** (crop 0.53×). That TAMBO's closing factor
   reproduces to within 2% on twenty-one times the ground establishes it as structural —
   a 100 m element cannot bridge the gaps `candidate_stride: 5` leaves — and not a
   peculiarity of Colca. **TAMBO's area remains a lower bound.** §9.9 still stands.
3. **The spread, which is where the crop could say nothing.** GRAND is one contiguous
   region: 48 were labelled, one was large enough, and it holds 101,948 detectors against
   a target of 10,000 — at full scale GRAND is not capacity-limited at all, where on the
   crop it fell 47% short. TAMBO is the opposite: 15,105 regions labelled, 45 large
   enough, 26 selected, spread over **310 km** (151 N-S × 271 E-W). The good TAMBO ground
   is scattered across the whole region rather than concentrated in one canyon. That is a
   deployment finding, and it is the one a crop is structurally unable to produce.
4. **`solid_angle` is still the weakest component everywhere** — 26 of 26 TAMBO sites and
   GRAND's single site, as on the crops. Holding at 21× the area, it is a statement about
   the criterion rather than about Peru: the Arequipa result is set by
   `solid_angle_half_sr`. It belongs in §10 with the other assumptions, and it is the
   first thing to check before quoting any of these numbers.

**The programme-level number barely moved: joint 50.2 km² at full scale against 50.1 km²
in the crop.** Searching twenty-one times more ground found essentially no additional
co-locatable ground — the shared ground is Colca, and the DEM's other canyons do not
offer it. The Jaccard falls from 0.0109 to 0.0006 only because GRAND's union grew; the
shared area itself is unchanged. Co-location remains decided by slope, and 44.9% of
TAMBO's ground sits inside GRAND's while 0.1% of GRAND's sits inside TAMBO's.

**One comparison that must not be made naively.** TAMBO's full-DEM area (111.9 km²,
`downsample_factor: 4`) is not commensurable with the crop's (83.6 km², factor 1): area
is measured on the downsampled mask while capacity is measured at full resolution, and
the run puts the cost at ~30% for a canyon strip. TAMBO looking *larger* than the crop is
therefore not evidence that it found more ground at Colca — it found other ground
elsewhere. The conclusion in (1) survives this because it rests on the acceptance rate,
which is measured on the same grid in both runs, and because a 30% correction does not
close a 16× gap.

**What it cost, and what that corrected.** GRAND 24.2 minutes, TAMBO **1.2 minutes** —
not the "~25-30 minutes each" this section previously assumed. TAMBO is cheap because its
targets are 2-5 km away against GRAND's 10-40, so the profile walks are a fraction as
long. The whole store is ~26 minutes of compute, not ninety. The case for keeping
notebooks 7 and 8 out of CI still holds on GRAND's 25 minutes alone, but it should be
argued from the real number.

**The notebook reads the results; it does not produce them.** The runner writes the small
artefacts — results JSON, provenance, explanation, a few hundred KB — into
`results/arequipa_full/`, which is committed, and the notebook opens those. This is only
possible because `explain.explain_results()` is a pure function of the results
dictionary: no DEM, no re-run, no pipeline.

Notebooks 7 and 8 are therefore **excluded from the CI execution job**, which is a real loss of
coverage and is replaced deliberately: `tests/test_docs.py` asserts that every
`ss.<name>` and `explain.<name>` the generator writes still exists, which is the drift
an API rename produces and the one the execution would have caught.

`manifest.json` and the per-run provenance record when the store was built and from
what, so a stale store is detectable. The premise of storing rather than recomputing is
that nobody looks again, so the store has to say when it was last right.

### 6.26a The memory estimator was wrong by 2.5×, and the first attempt died of it

**The first run failed 23 minutes in**, at the scoring stage, with a clean `MemoryError`
against its own 5.5 GiB address-space cap. The cap did its job — a `MemoryError` naming
the run, not an OOM killer taking the user's session — but it fired because
`estimate_peak_memory_gb` had advertised **2.32 GiB** for a search that measured
**5.68 GiB peak RSS**.

The model's structure was right and its candidate count was accurate: it predicted
15.43M candidates against 15.14M actual. What it got wrong was `n_observables=12`, the
arrays `arrival_scan.scan` returns. **The peak is not at the end of the scan.** By the
time `scoring.compose` runs, the scan's 11 arrays are still live, a component has been
built for each criterion, `compose` has clipped a float64 *copy* of every component, and
the composition and scoring intermediates need several more — about **36** per-candidate
arrays at 115.5 MiB each, three times what was counted. Counting them by hand gives
5.75 GiB; the run measured 5.68. Fixed by splitting the term: `n_observables` keeps its
honest meaning and a new `n_scoring_arrays=24` names the rest.

Corrected, the estimator says **7.21 GiB at `downsample_factor: 1` and 5.08 GiB at 4**,
against 4.5 and 2.32 before.

**And the knob everything pointed at is the weaker one.** `downsample_factor` scales the
labelling and gradient arrays as its inverse square, but candidates are taken on the
*native* grid — `candidate_stride` subsamples the surviving-pixel list, not the map — so
downsampling never touches the term that dominates at this scale. Going from 1 to 4 cuts
the estimate by 1.4×, not the 16× the inverse-square reasoning suggests; going from
stride 5 to 10 at factor 4 cuts it from 5.08 to 2.83. Every document that explained the
choice of `downsample_factor: 4` by "the labelling arrays scale as its inverse square"
was giving a true reason for a conclusion it does not support, and has been corrected.

Two consequences worth carrying forward. The default cap is 80% of available and the
warning fires when the estimate exceeds 80% of available, so **the two are only coherent
if the estimate is accurate** — a wrong estimate silently disarms the warning while the
cap still bites, which is exactly the failure seen. And on a machine whose desktop
already holds half of RAM, this run needs its ceiling set explicitly:
`tools/run_full_dem.py` (then `run_arequipa_full.py`) gained `--max-memory-gb` for that, since a memory ceiling is
a property of the machine and not of the science, and the configs were left untouched.

### 6.27 Documentation drift, and a test for it

The README documented **34 of 83** command-line options, and one of the 34,
`--fresnel_buffer`, had not existed for months. It also described the *old* precedence
rule — config file over command line — which §4.x fixed long ago. A wrong number in a
docstring is caught by `test_doctests`; a missing option in the README was caught by
nothing.

`tests/test_docs.py` now pins the coverage rather than the prose: every CLI flag appears
in the README, no documented flag is imaginary, the precedence section states the rule
the code implements, every module has a docstring of more than a few words (three did
not, including `site_searcher` itself, whose `automodule` page therefore opened with a
bare list of functions), and every name in `__all__` exists.

Also fixed while there: the startup banner still described a single ray cast to a target
mountain and a "clearance buffer" over intervening terrain, both replaced by the arrival
scan some time ago; and `--output_directory_base_with_given_json` resolved *before* the
merge loop and so kept the old precedence — a config file silently beat an explicitly
typed command line, for that one flag only, while every other flag on the same command
line was honoured.

### 6.28 Coverage, and the three bugs raising it found

Coverage was **65%**, and the gaps were not evenly distributed: `combine_experiments`
was at 30% and `crop_dem` at 18% — the two modules that turn a search into a *place*,
and the two whose failures are silent by nature. `main()`'s configuration merge was
untested entirely, including the precedence rule fixed in §6.27.

Now **74%**: `combine_experiments` 88%, `crop_dem` 90%, `site_searcher` 79%,
`explain` 94%. 497 tests, up from 370 at the start of this work.

Writing them found three faults, which is the argument for having written them:

1. **`oroscope-combine` crashed on any `search_mode: single` run.** `capacity_of` did
   not catch the `ValueError` from `int('N/A')` — the string single mode writes — and
   the console summary asked for thousands-grouping on a string when a run reported no
   capacity. Either one ended the combination.
2. **`main()` leaked its log file and left its `TeeLogger` installed.** Running it
   twice in one process stacked interceptors and leaked a handle each time. Restored in
   a `finally` now, which also made the CLI testable at all.
3. A test of `crop_dem` asserting exact pixel boundaries was **the test's** fault, not
   the code's: `crop` floors the start and ceils the stop, so it returns the smallest
   pixel-aligned box *containing* the request, and asking exactly on a boundary leaves
   floating point to decide whether the neighbouring pixel comes too. The tests now
   assert that containment property, which is the actual contract.

**Left low deliberately:** `figures` (14%), `fetch_dem` (30%), `generate_env` (27%),
`sensitivity` (14%). Those draw pictures, download over the network, generate an
environment file and drive subprocesses; testing them buys little against what it
costs, and the first is exercised by the docs build and the notebooks on every push.

### 6.29 Why a site is *good*, not only why it is weak

§6.23 reported each site's weakest component, which answers "what is wrong with this
site". Once a site has been *selected*, the more useful question is the other one, and
the same named components answer it — but only with a table saying what a high score
means, because a component's name states what it measures and not what satisfying it
implies. `COMPONENT_MEANING` supplies that, and `site_strengths()` returns the criteria
a site satisfies with **the measurement that earned each**:

    Site 3555 — 16.14 km², 1,901 detectors, facing SE, centred -15.6366, -72.1703
        Measured: 1.08 sr of accepted sky, targets at 3,137 m, striking 39°
        terrain, 784,440 g/cm² of rock behind, at 2,640 m altitude.
        Satisfies 6 of 7 criteria: column depth, exit distance, footprint
        sampling, shower development, tau decay …
        Held back by accepted sky at 0.57.

**Site records now carry coordinates** — centre latitude/longitude and a bounding box —
because they carried area, capacity and facing but no position, so "which ground is
this?" meant opening the raster in a GIS. The bounding box was already computed and the
centroid costs one pass over a region that has just been scanned.

Also added, because they were in the results and unread: the energy the geometric
aperture favours, and a **what to try next** section of concrete commands chosen from
what the run did — sweep the constraint that bound it, replace an absolute score cut
with a rank, see the area without closing.

### 6.30 The combination explains itself, and a wrong explanation caught in the act

`oroscope-combine` reported a joint area and a Jaccard index; neither says *why*. The
answer is usually not about neutrinos: a pixel has one slope, one altitude and one
aspect, and every experiment deployed on it must accept those same values. Colca's
entire co-location result follows from GRAND's 3–25° deployable band against a canyon's
~40° walls — a 20–25° sliver, 23% of the narrower band.

`explain_combination()` computes that from the two runs' own recorded parameters and
names the band that limits the sharing.

⚠️ **The first version of it was confidently wrong, and the output caught it.** It
treated the *distance window* and the *arrival-elevation window* as shared constraints
too, and duly concluded that GRAND and TAMBO "cannot share ground at all: their target
distance bands are disjoint" — printed directly above the 50.1 km² they demonstrably
share. Those windows are asked of the **view** from a pixel, not of the pixel: two
experiments looking out from the same hillside at different ranges and different
elevations are in no conflict whatever. Only slope, altitude and aspect are properties
of the ground. Now separated, with the viewing windows reported explicitly as *not* an
obstacle.

Worth recording as a class of fault: a summary that reasons about the numbers can be
wrong in ways a summary that merely restates them cannot, and it will be wrong
persuasively. The defence is the same as everywhere else here — check it against a case
whose answer is already known.

### 6.31 A component that appeared when it had been switched off

Found while reading the new output: TAMBO runs listed `geomagnetic` at 1.00 among the
reasons their sites were good, with `use_geomagnetic: false` in the configuration.

Whether the weighting was applied is judged by comparing the weighted solid angle with
the plain one — but a candidate that accepted *no* directions has a ratio of zero by
construction, and the test looked at the whole array, so those zeros stood in as
evidence. The component was then created, identically 1 for every viable candidate.

Harmless under the default product composition, which is why it survived: TAMBO's
capacity is unchanged at 9717 and GRAND keeps its real component at 0.742. Wrong under
`mean` or `min`, and misleading in any summary. Judged on candidates with accepted sky
only, now.

### 6.32 The three sources of defaults now agree ✅ delivered

A parameter could state its default in three places — the pipeline's signature, the
argparse parser, and `default_config()` — and they disagreed on **ten** of them. So
omitting a parameter meant different things depending on which door you came in by:
`search_mode` was `single` from Python and `distributed` from a shell, and
`min_dist_km` was 30 km against 10.

§6.24 closed the parity gaps in *capability*. This is parity of **meaning**, and it is
the one a user actually trips over — it was found by writing notebook 7, which omitted
`search_mode`, quietly ran a single search at 30 km, and found nothing.

| parameter | was (signature) | now |
| --- | --- | --- |
| `search_mode` | `single` | `distributed` |
| `grid_type` | `square` | `hex` |
| `target_antennas` | 1000 | 10 000 |
| `min_dist_km` | 30.0 | 10.0 |
| `min_sub_array_size` | 100 | 500 |
| `max_road_dist_km` | `None` | 20.0 |

In every case the signature was the odd one out, so the CLI and the template — the
documented values, and the ones the bundled configurations use — were taken as correct
and the signature aligned to them.

Two of the ten were in the template rather than the signature, and one of those was a
real hazard:

- **`origin_lat`/`origin_lon` were `0.0`.** Zero is a *valid* coordinate — it is in the
  Gulf of Guinea — so a placeholder someone forgot to edit would georeference a run to
  the wrong continent rather than fail. They are `null` now, which means "read the DEM's
  own tiepoint" and is the recommended use anyway.
- **`generate_kml` was `true`** in the template against `False` everywhere else.

`dem_path` and `region_name` remain deliberate placeholders — a generated template is
meant to be edited, and those two say so by being obviously unreal.

Pinned by `tests/test_cli.py`, which compares the three sources pairwise over every
parameter, so this cannot drift back.

### 6.33 A package, not a path insert ✅ delivered

Every notebook used to open with ``sys.path.insert(0, '../src')`` and then
`import site_searcher`, `import physics`, `import explain`. That worked in a clone, and
it was the wrong thing to teach: an unconditional path insert *shadows* an installed
copy with whatever happens to be in the source tree, and the layout underneath it was
flat — ten top-level `py-modules`, so there was no `import oroscope` to write and a
dozen generic names like `physics` and `explain` landed on any installing user's path.

`src/oroscope/` is a real package now. `__init__.py` re-exports the whole public
surface — **131 names** — so `import oroscope` is the only setup step, while the
submodules stay importable when a narrower namespace reads better:

```python
import oroscope
results = oroscope.find_grand_regions_interactive(...)
print(oroscope.explain_results(results))

from oroscope import physics          # when that reads better
```

The move touched the intra-package imports (11 of them), `pyproject.toml`'s
`packages`/entry points/coverage source, every test, the notebooks' preamble, the
documentation examples, and **67 docstring examples**, which are executed by
`test_doctests` and so failed loudly until each said where its module came from.

**One real defect fell out of it, and it is the interesting part.** `import oroscope`
imports `combine_experiments`, which called `matplotlib.use("Agg")` at module level.
Harmless while it was a standalone module someone imported deliberately; not harmless
as a package front door, where it reached into every caller's session and overrode the
inline backend — so notebooks captured no figures at all. Trap 3, one level up from
where it bit twice before. The backend is chosen in `main()` now, where the command
line actually needs it. **A library must not decide how its user's figures are
rendered.** The same applied to `sensitivity`, which set `MPLBACKEND` at import for the
benefit of its subprocesses.

Settled for the release: the logo is a **PNG everywhere** — 1024×1024 RGBA, used by
Sphinx's `html_logo`, by the README and by PyPI alike. PyPI does not render SVG, and
carrying a vector for two of those and a raster for the third is how a project ends up
shipping two different logos.

### 6.34 The stride-1 control at TAMBO settings: acceptance is unbiased, area is not

§9.9 asked whether TAMBO's area is a lower bound. It is, **by 4.75×**, and the reason is
not the one the funnel's own closing factor suggested.

`config/tambo_colca_stride1.json` is `tambo_colca_config.json` with `candidate_stride: 1`
and nothing else moved. 26 seconds. Against the stride-5 run on the same crop:

| | stride 5 | stride 1 | |
| --- | --- | --- | --- |
| directions accepted | **17.491%** | **17.494%** | unbiased, 0.017% apart |
| accepted pixels | 83,343 (×5 = 416,715) | 416,776 | agree to 0.01% |
| after gap closing | 222,658 | 486,322 | |
| **area** | **83.6 km²** | **396.9 km²** | **4.75× under-report** |
| sites | 15 | 29 | |
| capacity | 9,717 | 45,856 | 4.72× |

**Acceptance is unbiased and area is not, and the two facts are not in tension.**
Striding decides which pixels are *tested*, and it tests a fair sample — 17.491 against
17.494 is as close as this measurement gets. But it also decides which pixels are
*marked*, and the mask is closed morphologically before any area is measured. Marking
one pixel in five leaves gaps of 5 px; at Colca's 30.7 m that is **154 m**. TAMBO's
closing element is `antenna_spacing_km` = 100 m, or 3.3 px. **A 100 m element cannot
bridge a 154 m gap**, so the mask never reconnects: it stays a scatter of isolated
pixels, most regions fall below `min_sub_array_size`, and the area collapses.

GRAND's element is 1 km — 32 px against the same 154 m gap — which is why the earlier
control found striding clean and why this went unnoticed for so long. The rule is one
line: **the closing element must outrun the stride gap.**

| | element | stride-5 gap | |
| --- | --- | --- | --- |
| GRAND | 1000 m (32.5 px) | 154 m | bridges |
| TAMBO | 100 m (3.3 px) | 154 m | **cannot bridge** |

`warn_stride_outruns_closing()` now checks this at the top of every run, printing the
comparison and naming the three ways out (raise `gap_close_km`, lower
`candidate_stride`, or read the area as a lower bound and say so). It is checked against
both experiments in the doctests, since it is exactly the sort of guard that is wrong in
the direction of silence.

**What this means for the published numbers.** The mechanism was already named correctly
— §5.1 said in as many words that TAMBO's 100 m element "cannot bridge the gaps
`candidate_stride: 5` leaves". What was missing was its size, and the 0.53× figure
standing in for it understated the problem badly. 0.53 is the *ratio of the closed
stride-5 mask to the stride-corrected accepted count*, which folds the closing and the
fragmentation into one number and so measures neither. Separated: closing alone inflates
by **1.17×**, mildly, as it does everywhere; fragmentation costs **4.75×**. A reader
seeing 0.53 would reasonably infer the area was low by about half. It is low by nearly a
factor of five.

So: **TAMBO's Colca area of 83.6 km² should be read as ~397 km², and its capacity of
9,717 as ~45,856** — nine times its 5,000 target rather than twice. The full-DEM figure
of 111.9 km² is under-reported for the same reason *and* by downsampling on top, and
cannot be corrected by simply applying 4.75× (the full run also uses
`downsample_factor: 4`). Measuring it directly is not possible on this machine: TAMBO at
stride 1 over the full DEM is 26.8M candidates, which the corrected estimator puts near
10 GiB. GRAND is unaffected throughout, at either scale.

The joint area is affected as well, and in the direction that matters: it is limited by
TAMBO's mask, so 50.1 km² at Colca is also a floor rather than an estimate.

### 6.35 One config→pipeline translation, and the silent bug the three copies hid

§9.10 called the triplicated config→pipeline mapping "where a new parameter gets
forgotten". It had already happened, and worse than forgetting.

`main()` translated a configuration into pipeline keywords across sixty explicit lines.
`sensitivity`'s child process splatted the payload straight in. `tools/run_full_dem.py`
re-derived a third version. Only the first was complete.

**The sweep child never resolved `rfi_zones`.** A preset name reached
`find_grand_regions_interactive`, which does `for item in rfi_zones: if item[0] ==
'circle'`. Given the string `"arequipa"` that iterates *characters*; `'a'[0]` is `'a'`,
never `'circle'`, so every zone was skipped. No exception and no warning — and because
the count is `len(rfi_zones)`, the run cheerfully printed **`RFI Zones: 8 active`**, one
per letter, while excluding nothing. A sweep on a GRAND config would have searched
straight through Arequipa city and reported that it had not.

The recorded sweeps (§6.20–6.21) are **not** affected: they were run on
`tambo_colca_config.json`, whose `rfi_zones` is `"none"`, and iterating `"none"` also
yields nothing — which is the right answer there, by luck rather than by design. The
child also never inverted `require_sky`, and never made the bands tuples.

`config_to_pipeline_kwargs()` is now the single translation and `run_from_config()` the
single entry point; `main()` is one call, the sweep child is one call, and the runner is
one call. Verified behaviour-preserving by re-running GRAND over the Colca crop and
diffing: funnel, results, aperture and every parameter identical. The sweep child now
records five zones and a funnel identical to the direct run.

Two things learned in the doing:

- **Unknown keys are now dropped *and named*.** A misspelled `min_slop_deg` used to be
  ignored in silence, because an explicit keyword-by-keyword mapping simply never reads
  it. The translation warns, which is how a typo becomes visible.
- **Bind the signature at import, not at call time.** The filter first read
  `inspect.signature(find_grand_regions_interactive)` on every call, which meant it
  followed whatever that name pointed at — so the existing CLI tests, which substitute a
  recorder, presented a bare `(*args, **kwargs)` and the translation dropped every
  parameter it was meant to pass. Eighteen tests failed at once and were right to.
  `_PIPELINE_PARAMS` is bound once, to the real function.

### 6.36 Configuration paths are relative to the configuration ✅ delivered

§9.11: the search resolved `dem_path` against the working directory, so the bundled
configurations ran only from `src/` and produced a `FileNotFoundError` anywhere else.

A configuration that says `"dem_path": "../input/dem/colca.tif"` is describing where the
DEM sits relative to *itself* — the only fixed point it can reason about. `load_config()`
now resolves `dem_path`, `road_map_path` and `resume_dir` that way, and `main()` applies
the same rule to the output base, since fixing the inputs alone would have left the
*outputs* landing wherever the caller happened to stand (from the repository root, the
default `../output/` writes a sibling of the repository).

**No shipped configuration had to change**, which is what made this safe: `config/` and
`src/` are both one level below the root, so `../input/dem/colca.tif` names the same file
read from either. Verified by running the TAMBO stride-1 control from the repository
root and getting the identical answer — 29 sites, 396.9 km², 45,856 detectors.

A path that resolves only against the working directory is left alone with a warning.
Silently breaking a setup that relied on the old behaviour would be a poor way to fix a
convenience wart. A base typed on the command line is likewise left relative to the
caller, because that is their own instruction rather than the configuration's.

**Not covered: `oroscope-fetch-dem`.** It writes `../input/dem/` and `../config/`
relative to the working directory and has no configuration file to be relative to, so it
still wants running from `src/`. Documented rather than fixed.

### 6.37 The benchmark baseline, and measuring the machine before trusting it

§9.12 asked for `bench/baseline.json` to be refreshed "on a quiet machine". There is no
quiet machine here, so the first job was to find out what this one can actually resolve.

**Two consecutive passes over identical code**, nothing changed between them:

| case/stage | pass 1 | pass 2 | spread |
| --- | --- | --- | --- |
| `arequipa_900/ray_tracing` | 1.228 s | 1.820 s | **48.2%** |
| `synthetic_1800/ray_tracing` | 5.833 s | 8.036 s | **37.8%** |
| `arequipa_2500/ray_tracing` | 17.553 s | 16.489 s | 6.5% |

Both of the first two exceed the 30% regression gate. A single-pass baseline on this
host does not record the cost of the code; it records one sample of the noise and then
gates on it. Note the pattern, which is the useful part: the 17-second case is stable
and the short ones are not. Scheduler placement — a thread landing on a 4.6 GHz P-core
or a 3.4 GHz E-core — is a roughly fixed cost, so it dominates a one-second stage and
averages out over a seventeen-second one.

**Three changes, and then the refresh.**

`--repeat N` runs each case N times and keeps the **minimum** per stage. The minimum
rather than the mean, for a reason worth stating: timing noise is one-sided. Nothing
makes a stage run faster than its true cost, while a great many things make it slower.
The minimum therefore converges on the real cost as N rises, where the mean wanders with
whatever else the machine was doing.

`spread_pct` is recorded per stage, so the baseline carries **what the machine could
resolve** alongside what it measured. This is the field that makes the rest honest.

The gate is now **spread-aware**: a stage whose recorded spread is at least half the
gate is reported with a `~` and never failed on. Verified on the refreshed baseline —
`synthetic_900/ray_tracing` came in at 1.20 s → 3.07 s, a 156% "regression", and was
correctly not gated, because its own baseline spread is 150%. Gating on it would fail
builds at random while telling nobody anything.

Refreshed with `--repeat 5` at a load average of 2.39 — *not* a quiet machine, and
recorded as such. What it is worth is now legible from the file itself:

| | resolvable (spread < 15%) | not resolvable |
| --- | --- | --- |
| `arequipa_2500` | ray_tracing 8.4%, capacity 5.1%, morphology 3.7%, outputs 1.7% | topographic_screen 18.3% |
| everything else | a few stages | most stages, up to 149.6% |

**So `arequipa_2500` is the case this host can gate on, and largely the only one.** A
re-run against the fresh baseline reproduces it to 0.0% on total and −0.3% on ray
tracing, which is the confirmation that the methodology works rather than that the
machine got quieter.

The expectation from §6.25 — that `capacity_analysis` would be the one stage
legitimately slower, by 1.62× — is now absorbed into the baseline rather than
outstanding. It sits at 0.340 s on `arequipa_2500` with a 5.1% spread, so it is
measurable, and future changes to it can be gated properly.

### 6.38 β is configurable, and does not affect a search — which the run had been claiming it did

§9.5 asked for β, the tau energy-loss constant, to stop being a source literal.
`physics.set_tau_energy_loss(reference=..., index=...)` adopts a value in one place for
every function that uses it, `tau_energy_loss_settings()` reports what is in force, and
`restore_tau_energy_loss()` puts the shipped estimate back. The `BETA_*` constants stay
as the documented default, and an explicit argument still overrides the module setting
for a single call.

**The more useful finding is where β does not reach.** Tracing it before plumbing it:
nothing in the search path uses β at all. It enters `tau_range_gcm2`, `tau_survival` and
`tau_exit_probability` — tau production and escape *through rock* — and the search does
not model those. What the search weights by is the decay length `L = (E/m_τ)·cτ`, which
is kinematics and carries no β. `aperture.py` likewise uses only that length.

So the run's own summary was wrong. `explain.py` listed β under **"WHICH OF THESE ARE
ASSUMPTIONS — Choices, not measurements. Check them before quoting a result"**, with a
hardcoded `0.6e-6 cm²/g`, for a quantity no reported number depends on. That is the
worse direction for an error of this kind: it invites a reader to discount a result over
a parameter that never touched it, and it pads a list whose whole value is that
everything on it matters. β has moved to the "not modelled at all" sentence, alongside
neutral-current regeneration and the trigger, and the text now says explicitly that the
decay length carries no β.

A test pins that, by asserting `tau_decay_length_m` is unchanged by a tenfold β. If β
ever does enter the search, the claim in the explanation becomes false and that test
fails — which is the only way a statement like this stays true.

### 6.39 Column depth is truncated by the walk, measured at 6.4× — and the fix has a cliff

§9.3: "column depth is bounded by the walk unless `max_range_km` is set". The parameter
already existed and defaults to `max_dist_km`, so the profile walk stops at the target
rather than continuing through the rock behind it. What was missing was the size of the
effect. TAMBO over the Colca crop, three walk lengths against an unchanged 2–5 km
distance window:

| `max_range_km` | sites | area km² | capacity | mean column depth g/cm² | directions accepted |
| --- | --- | --- | --- | --- | --- |
| 5 (the default) | 15 | 83.6 | 9,717 | 747,016 | 17.49% |
| 20 | 15 | 83.6 | 9,717 | **4,765,107** | 17.49% |
| 60 | **2** | **5.3** | **605** | 5,522,952 | **5.98%** |

**The reported depth is a 6.4× under-report at the default**, and correcting it costs
nothing: at 20 km the selection is byte-for-byte the same — 83,343 accepted pixels in
both, the same 15 sites, the same 83.6 km². The walk was simply stopping before it had
measured the thing it reports.

**But the knob is not monotone and must not be maximised.** At 60 km the same run keeps
5.98% of directions against 17.49%, and the result collapses to two sites. The mechanism
is not yet identified — the mean horizon barely moves (+13.6° → +14.0°), so it is not
simple distant blocking — and it is recorded here as measured behaviour rather than
explained. Anyone raising `max_range_km` should check the funnel, not assume.

The configurations are unchanged, deliberately: this is a reporting defect rather than a
selection one, and changing them would restate the published numbers for a reason
unrelated to the physics they describe. The run's own assumptions block now carries the
measured factor and the warning about the cliff, so the depth is read as the lower bound
it is.

### 6.40 Declination can follow the site — the model is a socket, not a shipped table

§9.6 and §10.1: inclination follows the DEM's coordinates through a centred dipole,
while declination falls back to Arequipa's −6.9° wherever the search happens to be,
because the dipole is unusable for it (−0.2° against a measured −6.9°).

`physics.set_declination_model(fn)` takes any callable `fn(lat, lon) -> degrees` and
`default_field_for_site` consults it before falling back, so declination now follows the
site exactly as inclination does. `declination_from_grid(lats, lons, values)` builds such
a callable by bilinear interpolation, which is the practical route: export a grid from
NOAA's geomagnetic calculator covering the DEM and hand it over.

**No IGRF implementation is shipped, and that is a decision rather than an omission.**
IGRF is a spherical-harmonic expansion with a couple of hundred coefficients per epoch.
Writing them from memory would produce declinations that look entirely plausible and are
wrong — which is the failure mode this project has now found half a dozen times, and the
one that is hardest to notice. Either install `ppigrf`/`pyIGRF` and pass its function, or
supply a grid. The socket is the part that was missing; the coefficients are somebody
else's published work and should arrive as data.

Nothing about a published number changes: with no model set the fallback is exactly as
before, which the tests pin. The run's assumptions block now says the declination is
constant *unless a model was supplied*, and names the two ways to supply one.

### 6.41 Neutral-current regeneration, to leading order and clearly labelled as such

§9.4: only charged-current attenuation was counted, so the Earth-chord suppression is
overstated. A CC interaction removes a neutrino from the beam; an NC one only degrades
its energy, and on a falling spectrum those degraded neutrinos scatter *down into* the
band and partially refill it.

`nc_regeneration_factor()` is the leading term, and it is derivable rather than
asserted. For `Φ(E) ~ E^-γ`, a neutrino seen at `E` after one NC scatter of inelasticity
`y` started at `E' = E/(1-y)`, where the flux is larger by `(1-y)^γ`; the Jacobian
`dE'/dE = 1/(1-y)` supplies one more power, giving `(1-y)^(γ-1)` per scatter. With
`τ_NC = (X_chord/X_CC)·(σ_NC/σ_CC)` scatters expected, the factor is
`1 + τ_NC·(1-y)^(γ-1)`.

Measured at 1 EeV, against the ±3° window that matters:

| elevation | absorption only | with regeneration | lift |
| --- | --- | --- | --- |
| −0.5° | 0.836 | 0.883 | 1.06× |
| −1.0° | 0.700 | 0.778 | 1.11× |
| −3.0° | 0.342 | 0.458 | **1.34×** |
| −5.0° | 0.168 | 0.262 | 1.56× |

So the correction is small at the top of the window and grows with the chord, which is
the expected shape: the deeper the crossing, the more NC scatters and the more refilling.

**Off by default, and it is an approximation.** It is the first term of a series, not a
solution of the cascade equations, and it omits the `ν_τ → τ → ν_τ` chain that makes the
Earth genuinely translucent to tau neutrinos at high energy. It is clipped at 1 so it can
offset absorption but never manufacture flux. A result that leans on it is a result that
needs a real transport code.

One thing to watch: `earth_absorption_cutoff_deg` — the prediction in §6 offered as the
cheapest external check (§9.7) — does *not* apply this correction, so that prediction
still describes absorption alone. Turning regeneration on would move its lower edge
upward, and if it is ever compared against a collaboration simulation the two treatments
must be matched.

**A correction to a wrong statement made while writing this.** The first draft of the
docstring said a steeper spectrum regenerates *more*. It is the opposite: the neutrinos
scattering into the band come from `E/(1-y)`, above it, where a steeper spectrum has
*less* flux. The formula was right and the prose was backwards — γ = 2.0 gives 1.315 and
γ = 2.7 gives 1.258 — and two of the doctest values were wrong on top of that. Caught by
running them, which is exactly what §8.2 exists for. A test now pins the direction.

### 6.42 The decay weighting is selectable: flux, acceptance, or both

§9.1: an event rate is `∫Φ(E)·A(E)·P(E)dE` and the weight used was the flux alone.
`spectrum_weighted_decay_probability` now takes `weight_by`, with the config parameter
`decay_weight_by` and the flag `--decay_weight_by`:

| | weight | asks |
| --- | --- | --- |
| `flux` *(default)* | `E^-γ` | of the neutrinos that **arrive**, what fraction decays usefully? |
| `acceptance` | `A(E)` | over the energies the **detector responds to**, what fraction decays usefully? No assumed spectrum — which is the point, since γ is an assumption. |
| `flux_times_acceptance` | `E^-γ·A(E)` | the event-rate integrand itself. |

`A(E)` arrives as `--decay_response_csv`, a two-column table; `aperture.infer_response()`
recovers one from a published integral curve by dividing out the geometric model.

The mechanism is pinned by a check worth stating: **`flux_times_acceptance` with a flat
response reproduces `flux` to twelve decimal places.** A constant must cancel in the
normalisation, and it does, so the acceptance factor is entering exactly where it should.
`acceptance` is likewise independent of γ, also to twelve places.

**Using it on TAMBO produced a warning, not a result.** Recovering `A(E)` from
`data/tambo_aperture_fig3.csv` against our own Colca configuration and re-running:

| `decay_weight_by` | sites | area km² | capacity | mean decay term |
| --- | --- | --- | --- | --- |
| `flux` | 15 | 83.6 | 9,717 | 0.9566 |
| `flux_times_acceptance` | 15 | 78.1 | 9,093 | 0.9401 |
| `acceptance` | **0** | **0.0** | **0** | — |

The inferred `A(E)` rises monotonically to its maximum at the *top* of the published
range, 8.8 EeV. Weighting by it alone puts all the weight where the tau decay length is
hundreds of kilometres against a 2–5 km canyon, so the decay term collapses, every
candidate falls below `min_score` and the search returns nothing.

**That is a statement about the inferred response, not about TAMBO.** `A(E)` here is
`published / (our geometric model)`, and our model's aperture peaks near 66 PeV while the
published curve keeps climbing. The ratio therefore absorbs every high-energy effect our
model does not reproduce — a different baseline distribution, a geometry the ±20° window
does not capture — and attributes all of it to "response". `infer_response`'s own
docstring says it is "a better weight than a flat response, not a substitute for a
differential table"; this is what that caveat looks like when it bites.

So: the selector is the deliverable, and it works. `flux` remains the default and every
published number is unmoved. `flux_times_acceptance` is usable now and costs ~6% of
capacity. **`acceptance` alone should not be used with an inferred response** — it wants
a real differential acceptance table, which remains the outstanding ask of §10.

### 6.43 `min_score` against `score_percentile`, measured — and left to the owner

§9.2 called `min_score` the dominant assumption and noted `--score_percentile` as the
scale-free alternative the configs do not use. Measured on TAMBO over the Colca crop:

| cut | sites | area km² | capacity | kept, of strided |
| --- | --- | --- | --- | --- |
| `min_score` 0.35 *(current)* | 15 | 83.6 | 9,717 | 17.49% |
| `score_percentile` 5 | 1 | 3.1 | 363 | 3.84% |
| `score_percentile` 10 | 5 | 18.3 | 2,121 | 7.68% |
| `score_percentile` 17.5 | 15 | 60.9 | 7,004 | 13.43% |
| `score_percentile` 25 | 17 | 94.2 | 10,976 | 19.19% |
| `score_percentile` 40 | 31 | 186.8 | 21,539 | 30.71% |

**`min_score` 0.35 is `score_percentile` 22.8 on this terrain**, by interpolation. Note
the percentile is taken over *viable* candidates rather than all strided ones, which is
why 17.5 keeps 13.4% of the strided set rather than 17.5%.

**There is no knee.** Area runs 3.1 → 18.3 → 60.9 → 94.2 → 186.8 km² across the range —
smooth, monotone, and close to linear above 10%. Nothing in the data marks 0.35 or any
other value as the natural cut. That is the strongest form of §5.1's claim that
`min_score` is an assumption rather than a measurement: if the terrain had a natural
threshold, a scan across the cut would show it, and it does not.

The case for switching is that a percentile means the same thing when the composition
changes, and an absolute cut on a product does not — add a component and every score
falls, so 0.35 silently becomes a harsher cut. The case against is that every published
number used 0.35, and switching restates all of them.

**Not switched.** The configs are unchanged. This changes published science and is the
owner's call, not a refactor to be slipped in; §10 already lists `min_score` among the
TAMBO assumptions to check with the collaboration, and this is the table to check it
against. If it is switched, `score_percentile: 22.8` reproduces the current selection
most closely.

### 6.44 `tau_exit_probability` under-resolves its own integral at large depth

Found while building `tau_in_rock` (§6.45), not looked for. `physics.tau_exit_probability`
integrates over the interaction depth with `x = np.linspace(0, X, samples)` — a fixed
number of points spread over the **whole** depth. But only interactions within roughly
one tau range of the far surface contribute anything, so as `X` grows the sample spacing
outruns the only region that matters and the trapezoid rule reports the area of a spike
it never resolved.

Measured at 3 PeV, where the tau range is 3.1×10⁴ g/cm² so the effect starts early:

| `samples` | X = 10⁷ | 10⁸ | 10⁹ |
| --- | --- | --- | --- |
| 2000 *(default)* | 2.335e−05 | 2.628e−05 | **8.884e−05** |
| 20,000 | 2.330e−05 | 2.181e−05 | 1.328e−05 |
| 200,000 | 2.330e−05 | 2.176e−05 | 1.103e−05 |
| 2,000,000 | 2.330e−05 | 2.176e−05 | 1.100e−05 |

At X = 10⁹ the default is **8× the converged value**, and worse than the magnitude: it
inverts the sign of the trend. Converged, the exit probability *falls* with depth, as it
must; at the default it *rises*, so the curve has a spurious maximum at the grid edge.

**The knock-on is `depth_band_from_energy`,** which locates the band by the half-maximum
of that curve. It returns (5.6×10⁷, 2.9×10⁸) for TAMBO's configured 3 PeV – 1 EeV range —
a band that excludes the true 1 EeV optimum of 5.7×10⁶ by more than an order of
magnitude, and whose low edge *rises* when the minimum energy is *lowered*
(100 PeV – 10 EeV gives 5.2×10⁵; 3 PeV – 10 EeV gives 5.6×10⁷). The docstring's stated
"roughly 5e5 to 2.6e8 across 100 PeV to 10 EeV" is the well-behaved case.

**No published number is affected.** Every config leaves `depth_band_gcm2` null, so runs
score against the default (10⁵, 10⁷) and never call this. `production_escape_optimum_gcm2`
searches to 10⁹ but its answers are in the resolved regime at ≥100 PeV: 3.30×10⁶ at
100 PeV, 5.71×10⁶ at 1 EeV, 6.23×10⁶ at 10 EeV, all reproducing published values.

**Not fixed.** The fix is a substitution — integrate in `u = X − x` on a log grid, or
truncate to a few tau ranges below the surface — and it changes a physics function's
outputs, so it wants its own change with a test that pins the converged values. Left to
the owner. `tools/make_animations.py` scores against the default band and says so.

### 6.45 Four more animations, and the two things building them measured

Added to `tools/make_animations.py`, taking it to eight: `the_azimuth_fan`,
`product_collapse`, `slope_criterion`, `tau_in_rock`. Notebook 9 builds all of them.
The filter was whether the *intermediate states carry the argument*; six candidates
were rejected because a static figure does them better, and those reasons are in the
handover rather than here.

Three of the eight now read `input/dem/colca.tif` when it is present, because they are
about what a criterion does to real ground, and fall back to synthetic terrain — saying
so on the figure — when it is not. The fallback is honest but degenerate for
`slope_criterion`: `synthetic.colca_like` has every wall at exactly 40.6°, so the
criterion flips between 39° and 42° instead of eroding.

**(a) The product collapse is real but is not evenly shared.** The premise was that each
component multiplied in drags the whole population toward zero. Measured over the central
40% of the Colca DEM, 119,788 candidates screened and 98,343 viable, against TAMBO's own
`min_score` of 0.35:

| after multiplying in | median score | above the cut |
| --- | --- | --- |
| *(nothing)* | 1.000 | 100.0% |
| `depth` | 1.000 | 100.0% |
| `solid_angle` | 0.233 | 35.9% |
| `distance` | 0.233 | 35.9% |
| `shower` | 0.228 | 34.3% |
| `decay` | 0.218 | 32.3% |
| `footprint` | 0.214 | 32.2% |

The collapse happens — 100% to 32.2%, median 1.000 to 0.214 — but **one component does
nine tenths of it** and two do nothing measurable. `depth` is ~1 everywhere because
Colca's canyon walls sit inside the default band. `distance` is **provably** inert: the
scan already applied the same 2–5 km window as a hard criterion, so every surviving
candidate scores 1 on it by construction. Scoring a criterion the scan has already
enforced is free, but it is also empty, and the funnel does not show it because the
candidates were gone before scoring saw them.

This does not weaken §6.43's case against thresholding a product — it sharpens it. The
cut moves from harmless to decisive on the addition of a single term, and which term
that is depends on the terrain, not on the configuration.

**(b) Where the wall-slope criterion bites.** `min_target_slope_deg` swept over the same
crop, against a wall-slope distribution whose quartiles are 22.2° / 29.7° / 34.6° and
whose 95th percentile is 41.2°:

| `min_target_slope_deg` | candidates accepting a direction | of no criterion at all |
| --- | --- | --- |
| unset | 106,926 | 100% |
| 15 | 103,725 | 97% |
| 25 *(TAMBO's value)* | 98,343 | 92% |
| 30 | 93,465 | 87% |
| 40 | 77,139 | 72% |
| 45 | 66,586 | 62% |
| 50 | 54,693 | 51% |
| 60 | 27,773 | 26% |
| 69 | 7,423 | 7% |

**The mask outlives its own median by 20°.** Half the candidates see a mean wall slope
below 29.7°, yet a 30° floor still keeps 87% of them, and the half-way point is not
reached until 50°. The criterion is applied *per direction* while the observable is a
mean over each candidate's accepted directions, so a candidate keeps its steepest
directions long after its average has fallen below the cut. Read the reported
`target_slope_deg` as a description of a candidate, never as a predictor of what a cut
will do to it. TAMBO's configured 25° costs 8%, which is consistent with §10 calling it
a deliberately permissive floor.

### 6.46 Peru, all of it, in four minutes — and what `--max_memory_gb` actually caps

The first search over a whole country. `config/grand_peru_survey.json`, GRAND's criteria
copied unchanged from `grand_arequipa_full.json` so the two are comparable, over
`input/dem/peru_SRTMGL3.tif` — 22,080 × 15,360 = **339.1 Mpx** at 3 arc-seconds, a
92.17 × 91.57 m pixel.

**3 arc-seconds is forced twice over, not chosen.** By memory: the same box is 3,052 Mpx
at 1 arc-second, and `estimate_peak_memory_gb` puts even 3 arc-seconds at 12.4 GiB at
`downsample_factor` 1 / `candidate_stride` 5. And by the API, which caps a request at
4,050,000 km² for 90 m datasets and 450,000 km² for every 30 m one; this box is about
2.86 million km², six times over the 30 m limit. Run at `downsample_factor` 4 and
`candidate_stride` 15, estimated 4.77 GiB.

| stage | time |
| --- | --- |
| load DEM | 4.6 s |
| topographic screen | 14.7 s |
| ray tracing | 200.1 s |
| morphology | 20.6 s |
| capacity | 8.3 s |

**Four minutes for a country**, against 26.8 minutes for the Arequipa DEM at 30 m —
2.6× the pixels but a ninth the candidates, and the candidates are the cost.

| | DEM pixels | slope band | after stride | directions accepted | rate |
| --- | --- | --- | --- | --- | --- |
| Peru, 90 m, stride 15 | 339,148,800 | 112,156,858 | 7,477,157 | 3,221,209 | **43.1%** |
| Arequipa, 30 m, stride 5 | 128,600,000 | — | — | — | 61.6% |

**17 sites, 563,411 km², 633,655 antenna positions** at 1 km hexagonal spacing. That is
43.8% of Peru's 1,285,216 km². One site holds 94.8% of it — 533,861 km² centred at
−10.23, −75.03, on the eastern Andean flank — and the next largest is 4,834 km².

**Read the area as a survey number, and here is the width of the bracket.** The 3,221,209
accepted strided pixels each stand for 15, which is 48.3 Mpx or **407,805 km²**; the
run reports 66.75 Mpx or 563,411 km² after closing, pruning and selection, a factor
**1.38** higher. Some of that gap is closing doing its job — filling holes inside
accepted ground, which the stride correction also does — and some is a 1.5 km element
bridging ground that was genuinely rejected. Both estimates are approximations and the
honest statement is **4–6 × 10⁵ km²**.

**"17 sites" is the number to distrust, not the area.** The largest site's bounding box
is the *entire DEM* — north −0.006, south −18.396, west −81.134, east −68.604 — so a
1.5 km closing element applied to a strided scatter across 339 Mpx has merged the whole
cordillera into one connected component, which `min_width_km: 2.0` then had no reason to
break up. Its accepted candidates are Andean (mean altitude 2,446 m, p50 2,256, p90
4,546) while the polygon enclosing them reaches the coast and the basin. Site *count*
and site *extent* from this run are artefacts of the element; the accepted-candidate
statistics inside them are not.

**One thing that was checked and turned out fine.** The obvious worry is that 3
arc-seconds changes the slope distribution itself — smoothing steep ground into the
3–25° band and roughening flat ground into it — which would make the whole screen a
resolution artefact. Measured on 20 Mpx of Arequipa, comparing the native grid against a
3×3 block mean of the same ground:

| grid | in 3–25° | slope quartiles |
| --- | --- | --- |
| 30 m | 67.6% | 5.2 / 11.7 / 21.4 |
| 90 m | 67.4% | 4.4 / 10.4 / 19.2 |

The quartiles do shift down, as smoothing must make them, but the *band fraction* moves
by 0.2 points: what is lost at the 25° ceiling is regained at the 3° floor. So the
screen's 74%-of-Peru is a fact about Peru and not about the grid. The prediction here
was wrong and the measurement is the useful part.

Peru accepts a *lower* fraction of its screened candidates than Arequipa does (43.1% vs
61.6%), which is the expected direction: the national box adds coastal desert below the
3° floor, high Andes above the 25° ceiling, and Amazon lowlands that are flat.

**No TAMBO counterpart, deliberately.** A 90 m pixel cannot resolve the geometry TAMBO
depends on — Colca's floor is ~1 km wide, so a canyon is ~11 pixels across and the wall
the array stands on is a handful, with the 20–60° band measured on a slope the grid has
already averaged. And TAMBO's 100 m closing element against this run's 1,382 m stride
gap is a ratio of 13.8; raising `gap_close_km` to 1.5 km to survive that is fifteen times
the array's own scale, so the mask would smear across the canyon instead of tracing the
wall — wrong in the opposite direction from §6.34. A national TAMBO answer needs 1
arc-second and tiling, and that is a different job.

**`--max_memory_gb` caps virtual address space, not resident memory.** The first attempt
ran the whole search successfully at `max_memory_gb 7.0` and then failed on the *map*
with `Unable to allocate 40.8 MiB`. Nothing was near exhausting 7 GiB of RAM. The cap is
`RLIMIT_AS`, which counts every mapping — including the 1.36 GB memory-mapped `.npy` and
the ping-pong buffers, all of it file-backed and evictable — while
`estimate_peak_memory_gb` deliberately estimates only *anonymous* memory and says so.
On a 339 Mpx DEM those two quantities differ by more than 2 GiB, so a cap set from the
estimate is tight by exactly the amount the estimate excludes. The run completed at 11.0.
The failure mode is benign but late: the search finishes, the JSON and the GeoTIFF are
written, and the picture is the thing you lose.

**Then that advice, followed once more, took the machine down** — and the pair is the
real lesson, so both halves are recorded here. A follow-up run at `candidate_stride` 10
(estimate 6.74 GiB) with the cap raised again to 13.0 died in the labelling stage and
killed the session with it. 13.0 was above the ~8.7 GiB the machine had available, and
**a cap above available memory is not a cap**: RLIMIT_AS aborts the process only if it
is reached before the kernel runs out of memory to give, and above that line the OOM
killer always gets there first. The limit was set, reported, and could never fire.

The two constraints pull opposite ways — clear the estimate by the mapped size of the
DEM, *and* stay under what the machine has — and nothing in the pre-flight said when
they could not both be met. `preflight_memory` now checks it, warns, and returns
`cap_exceeds_available` so a library caller can act on it:

```
⚙️  Estimated peak memory: 6.7 GiB, available 8.7 GiB
⚙️  Address space capped at 13.0 GiB (max_memory_gb=0 disables)
⚠️  That cap is above the 8.7 GiB currently available, so it cannot protect this
    machine: a runaway reaches the OOM killer before the limit fires.
```

Note that the *estimate* warning would not have fired: 6.74 against 8.71 available is
77%, just under the 80% threshold. The cap check is the one that catches this.

**When the two constraints cannot both be met, the configuration does not fit**, and the
answer is a coarser search — `candidate_stride` is the memory lever (§6.26a) — not a
bigger number in the cap. The committed config is now 8.0: above the measured ~7.0
virtual high-water mark, below available. The stride-10 variant was not retried and its
numbers are unknown; the committed survey is stride 15.

### 6.47 A joint site is not one polygon, and an optimiser told to work inside it finds nothing

The question asked was whether the code can lay out a *specific realization* at a
combined site — roughly 100 TAMBO units with a few GRAND antennas among them. The
answer is that it cannot yet, and the interesting part is not the missing placement
routine. It is that **the region a placement routine would obviously be pointed at is
the wrong region.**

Measured on the stored Colca combination, per pixel, against each experiment's own
deployable slope band:

| ground inside the 50.1 km² joint mask | px | share | km² |
| --- | --- | --- | --- |
| GRAND's band 3–25° only | 596 | 1.1% | 0.55 |
| both bands, 20–25° | 3,372 | 6.2% | 3.09 |
| TAMBO's band 20–60° only | 50,565 | 92.4% | 46.27 |
| neither | 195 | 0.4% | 0.18 |

The joint mask's slope quartiles are 31.1° / 36.3° / 41.8°. It is canyon wall. **Only
3.63 km² of it is ground a GRAND antenna could stand on**, and that 3.63 km² is
**1,702 disconnected fragments whose largest is 0.038 km²** — against the 0.866 km² a
single cell of a 1 km hexagonal lattice occupies. Not one fragment holds one antenna
position. `count_grid_capacity` reports 1 for the whole joint mask, and for once that is
not the anchored-not-fitted caveat understating things: the continuum limit is 4.2 and
even a perfectly fitted lattice could not reach it.

So an optimiser handed the intersection would report that a few GRAND antennas cannot be
placed at all, and it would be answering the question it was asked rather than the
question that was meant.

**The realization is buildable; it is just not one polygon.** The two halves are
unconstrained and adjacent rather than competing:

- **TAMBO.** 100 units at 100 m hexagonal spacing need 0.87 km². There are 49.4 km² of
  TAMBO-band ground inside the joint mask — 57× the room required. Siting them is not a
  constraint, it is a choice.
- **GRAND.** 10 antennas at 1 km need 8.66 km², which is not in the joint mask. It is
  next door: of GRAND's 4,580 km² Colca-crop mask, 2,995 km² is in the 3–25° band and
  2,744 km² of that lies in 48 patches large enough to hold a lattice cell. From joint
  ground the nearest such patch is a **median 0.92 km away; 53% of the joint mask is
  within 1 km of one and 84% within 3 km.**

Under one kilometre is *inside a single GRAND cell*. From GRAND's point of view the two
arrays are co-located; from the mask's point of view they never overlap. That is the
whole content of the "joint" idea, and it is why the Jaccard index of 0.0006 reported in
§6.26 reads as a failure and is not one.

**What this says about building the thing.** The three gaps named in the handover are
real and unchanged — the per-pixel score is aggregated to mean/median/p90 by
`summarize_observables_by_site` and never rastered; there is no placement routine; and
the score is a ranking proxy rather than an event rate, which still waits on `A(E)`
(§9.1, and see §6.42 for why an inferred table is unsafe). But the *domain* matters
before any of them:

1. Optimise over the **union** with a per-role band constraint, never over the
   intersection. The intersection is TAMBO ground that GRAND's region happens to
   enclose.
2. The coupling term is **shared line of sight**, not shared footprint. What makes a
   pairing joint is that both arrays watch the same wall, and that is a property the
   arrival scan already computes and the combination step currently discards.
3. A GRAND antenna's constraint is **patch size**, not area. 3.63 km² that never
   assembles 0.87 km² in one piece holds nothing, and a routine that reasons in total
   area will not notice.

Still analysis. Nothing here is implemented, and the placement routine remains unwritten
by choice rather than by oversight.

### 6.48 Ancash: the same question over steeper ground, and the terrain predicted the answer

The Arequipa pair repeated over Ancash — the Cordillera Blanca and the Callejón de
Huaylas — at the same 1 arc-second SRTMGL1 resolution, from the same source. 9,855 ×
6,958 = **68.6 Mpx**, 64,684 km², bounds from OpenStreetMap's administrative boundary
rather than eyeballed. Zero SRTM voids (OpenTopography serves the void-filled product,
which matters here: the glaciated ground above 5,000 m is 1.16% of the DEM and would
have been the first thing to drop out). Maximum elevation 6,744 m against Huascarán's
6,768.

**Every transferable criterion was held fixed**, and that is checkable rather than
assertable — the notebook diffs the configurations. TAMBO differs in **2 of 60**
settings, both bookkeeping: the file it reads and the name it prints. GRAND differs in
**3 of 52**, the only real one being `rfi_zones`. Arequipa's run excludes five
hand-curated circles and there is no Ancash preset; inventing one would have injected a
new assumption into a run whose entire purpose is comparison, so Ancash excludes nothing
and declares it. Arequipa's zones cover ~3,500 km² of a ~120,000 km² box, so read
Ancash's GRAND area as **at most ~3% flattered** on that account.

**The terrain made a prediction before any ray was traced.** Over land only:

| | Arequipa | Ancash |
| --- | --- | --- |
| median slope | 11.1° | **23.0°** |
| in GRAND's 3–25° band | 70.3% | **52.0%** |
| in TAMBO's 20–60° band | 24.1% | **58.0%** |

Ancash is twice as steep, so GRAND should do worse per unit area and TAMBO much better.
Naively that is 0.74× for GRAND and 2.41× for TAMBO.

**What the searches found.** Ancash is 0.533× Arequipa's pixel count, so that is the
ratio everything is read against: near 0.53× means "the same ground, less of it".

| | Arequipa | Ancash | ratio | per pixel |
| --- | --- | --- | --- | --- |
| GRAND, sites | 1 | 1 | | |
| GRAND, area km² | 88,527.5 | 43,091.2 | 0.49× | **0.91×** |
| GRAND, capacity | 101,948 | 49,447 | 0.49× | **0.91×** |
| TAMBO, sites | 26 | 35 | 1.35× | |
| TAMBO, area km² | 111.9 | 174.9 | 1.56× | **2.93×** |
| TAMBO, capacity | 9,024 | 14,290 | 1.58× | **2.97×** |
| joint, area km² | 50.20 | 75.25 | 1.50× | **2.81×** |
| Jaccard | 0.000567 | 0.001742 | 3.07× | |

**The prediction holds on both counts, and the direction is the whole result: Ancash is
worse for GRAND and about three times better for TAMBO.** GRAND's loss is milder than
the naive 0.74× because its 1 km closing element fills in around a fragmented mask;
TAMBO's gain exceeds the naive 2.41× because *both* of its stages improve — the slope
screen keeps 33.9M pixels against Arequipa's 26.8M **from a DEM half the size**, and
acceptance among those rises from 9.7% to 15.1%. GRAND's acceptance falls, 61.6% to
54.9%. (Compare those on the funnel rows `directions accepted / kept by stride 5`, not
positionally: Arequipa's GRAND funnel carries an extra `outside RFI zones` stage that
Ancash's does not.)

**GRAND's binding constraint moved, and TAMBO's did not.** At Arequipa GRAND binds at
`directions accepted` (61.6% kept) — plenty of deployable ground, and the arrival
geometry decides. At Ancash it binds at `slope 3.0-25.0 deg` (44.4% kept): the mountains
do not offer enough ground gentle enough to stand an array on, and the search never gets
as far as asking what that ground can see. TAMBO binds at `directions accepted` in both,
from opposite sides — 9.7% against 15.1%. **A criterion that binds is a statement about
the ground rather than about the configuration**, and here it moved when only the ground
moved, which is about as clean a demonstration as the funnel can give.

**One invariant worth keeping.** The joint region is **44.9% of TAMBO's mask at Arequipa
and 43.0% at Ancash** — essentially unchanged across two regions with very different
terrain. That is §6.47's finding arriving independently: the joint is TAMBO-limited, and
co-location costs GRAND almost nothing. The Jaccard index tripling is TAMBO's mask
growing, not the two experiments agreeing more.

**The best joint ground in the Ancash run is not in Ancash.** Asked where the joint
patches in the south-east of the map are, the answer is worth recording because it
generalises. The 91 joint patches were labelled and the largest reverse-geocoded:

| km² | lat, lon | nearest village | reverse-geocodes to |
| --- | --- | --- | --- |
| **5.37** | −10.5866, −77.0729 | Gorgor, 5.1 km | **Gorgor, Cajatambo, Lima** |
| 4.23 | −10.5345, −77.2167 | Aco, 5.0 km | Carhuapampa, Ocros, **Ancash** |
| 3.24 | −10.5945, −77.1550 | Manás, 1.3 km | **Manás, Cajatambo, Lima** |

The largest joint patch in the entire run — 3,442 m elevation, 3,539 m of relief within
8 km, median slope 26° — sits in **Lima region, not Ancash**. Regions are downloaded as
*bounding boxes* and departments are not rectangles, so Ancash's box reaches south past
the border. **37% of the joint ground (27.5 km² of 75.2) lies south of −10.45°**, in
that corner.

Two consequences. **File results by box, read them by geography** — a site named for the
run that found it may be administratively somewhere else, which matters for anything
involving permits, access or a collaboration's national footprint. And the Ancash and
Lima boxes overlap by **9,198 km²** (−10.79…−10.23 lat), so a Lima run will search this
same ground again — with **AW3D30 rather than SRTMGL1**. That is an accidental but
genuinely useful cross-check: the same terrain, two independent datasets, and a chance
to see whether the joint patches survive a change of DEM. Nothing else in this project
has had that test.

`results/region_comparison.md` holds the cross-region table and is regenerated by
`tools/compare_regions.py` from the stores, so it does not go stale when a region is
added.

Timing: GRAND 9.1 minutes, TAMBO under one, against Arequipa's 26.8 and ~1. Memory
estimate 2.92 GiB at `downsample_factor` 4 / `candidate_stride` 5, run with
`--max-memory-gb 5.0` against ~7 GiB available — comfortable, unlike Arequipa.

All the TAMBO caveats carry over unchanged: ~4.75× low from striding against a 100 m
closing element (§6.34) and ~30% again from downsampling. **Both regions are biased the
same way, which is why the ratio is the trustworthy number and the absolute areas are
not.**

### 6.49 The striding penalty for TAMBO is not 4.75×. On the Callejón de Huaylas it is 291×

A zoom-in over the Callejón de Huaylas and the Cañón del Pato — the Río Santa valley
between the Cordillera Blanca and the Cordillera Negra, cropped out of the Ancash DEM at
`−8.80…−9.90` lat, `−78.00…−77.20` lon, 3,961 × 2,881 = **11.4 Mpx**. Small enough to run
at `downsample_factor` **1** and `candidate_stride` **1**: 2.64 GiB estimated against
~7.8 GiB available. **The first unbiased run this project has done at scale.**

| | GRAND | TAMBO |
| --- | --- | --- |
| sites | 1 | **109** |
| area km² | 8,294.9 | **855.1** |
| capacity | 9,609 | **98,696** |
| joint | 637.1 km² | Jaccard 0.075, 74.5% of TAMBO |

Then the control: **the same crop, the same criteria, changing only the sampling** to the
`downsample_factor` 4 / `candidate_stride` 5 the department runs use.

| | ds 1 / stride 1 | ds 4 / stride 5 | ratio |
| --- | --- | --- | --- |
| GRAND sites / area / capacity | 1 / 8,294.9 / 9,609 | 1 / 7,537.9 / 8,658 | **1.1×** |
| TAMBO sites | 109 | 1 | **109×** |
| TAMBO area km² | 855.1 | 2.9 | **291×** |
| TAMBO capacity | 98,696 | 256 | **386×** |

**§6.34's 4.75× is not the size of this effect.** That measurement was taken on Colca,
varying the stride alone at `downsample_factor` 1. Here both levers move together, on
terrain whose accepted strips are numerous and individually small, and TAMBO loses
**two and a half orders of magnitude**.

**The mechanism is not the area measurement.** The funnels say exactly where it goes:

| stage | ds 1 / stride 1 | ds 4 / stride 5 |
| --- | --- | --- |
| slope 20–60° | 7,081,749 | 7,081,749 |
| directions accepted | 991,099 | 198,353 |
| after gap closing | 1,248,669 | 477,816 |
| **pixels in selected sites** | **912,320** | **3,136** |

Acceptance is **identical**: 991,099/7,081,749 = 14.0% at stride 1, and 198,353/1,416,351
= 14.0% at stride 5. Striding really is unbiased in acceptance, exactly as §6.34 says.
Closing differs by only 2.6×. **All 291× of the loss happens between closing and
selection**, in the region thresholds: at stride 5 the mask fragments into 7,954 labelled
regions of which 5 clear the area threshold and **1** clears `min_sub_array_size` (250
detectors). At stride 1 the mask is contiguous and 109 regions survive. The under-report
is fragmentation meeting a minimum-array-size cut, not pixels being miscounted.

GRAND is untouched for the reason it always was: a 1 km closing element bridges a 154 m
stride gap without noticing, so its mask never fragments.

**What this changes.** Every TAMBO number this project has published from a strided,
downsampled run is a lower bound by a factor that is **terrain-dependent and unbounded in
practice** — 4.75× at Colca, ~291× here. The honest reading is that strided TAMBO area
and capacity are *not* estimates of the true values at all; they are a different quantity
that happens to correlate. **Quote TAMBO numbers only from unbiased runs, or quote them
as "at least".** The Ancash and Arequipa department TAMBO figures (174.9 km² / 14,290 and
111.9 km² / 9,024) should be read in that light, and the ratio between them survives
better than either absolute.

It also means the Callejón de Huaylas was effectively invisible to the department run:
the ancash_full TAMBO mask contributes **1.2 km² inside this crop window** against the
crop's own 855.1 km².

Stored in `results/huaylas_full/`, controls included, with
`config/{grand,tambo}_huaylas.json` and their `_control` counterparts.

### 6.50 The memory pre-flight modelled the search and not the map, and the map is what kept dying

Three runs in one session finished their searches and then died drawing the picture,
having already written the JSON and the GeoTIFF — Peru at a 7.0 GiB cap, the Huaylas
combination at 5.5, the Huaylas TAMBO run at 5.5 and again at 6.0 and 7.0. A fourth,
given `max_memory_gb 0`, took the machine down.

`estimate_peak_memory_gb` models the candidate and scoring arrays, which is what §6.26a
calibrated it on and what a *search* allocates. The map is a separate peak landing on top
of them at the very end. It renders at `viz_ds = downsample_factor * 2`, so its raster is
`rows/(2d) × cols/(2d)` — at `downsample_factor` 1 on the Huaylas crop that is
1,981 × 1,441, exactly the array in the failure message.

Measured, shading through `LightSource.shade` and saving at 150 dpi:

| viz raster | peak RSS above idle |
| --- | --- |
| 700 × 500 (0.35 Mpx) | 126 MB |
| 1400 × 1000 (1.4 Mpx) | 263 MB |
| 1981 × 1441 (2.85 Mpx) | 539 MB |
| 2800 × 2000 (5.6 Mpx) | 959 MB |

**~190 bytes per viz pixel**, plus ~130 MB of matplotlib import and canvas.
`estimate_visualisation_memory_gb` carries both, `preflight_memory` adds it to the search
term and reports the split, and `preflight_memory(refuse=True)` now **raises rather than
warns** above `REFUSE_FRACTION` (0.8) of available. `tools/run_full_dem.py` passes
`refuse=True` on a real run and **rejects `--max-memory-gb 0` outright**.

**And the estimator is about 2× optimistic at `candidate_stride` 1.** Measured on the
Huaylas TAMBO run by polling `/proc/<pid>/status`:

| | estimate | measured |
| --- | --- | --- |
| search + map | 3.27 GiB | — |
| peak RSS | — | **5.31 GiB** |
| peak virtual (what `RLIMIT_AS` caps) | — | **6.51 GiB** |

So a cap set from the estimate is roughly half what the run needs at stride 1, which is
why 5.5, 6.0 and 7.0 all failed on the map while the search itself completed every time.
§6.26a calibrated the estimator on a *strided* run; nothing has calibrated it at stride 1,
and this is the one data point. **Treat the estimate as a lower bound, and size a cap
from measurement when the sampling is unusual.** The proper fix is a second calibration
point at stride 1 and a corrected `n_scoring_arrays`; not done.

The map for that run was never produced. It is not needed: `--reveal` on the combination
already renders TAMBO alone for the same crop, which is what notebook 10 shows.

### 6.51 Three departments, one question: the answer tracks median slope, and the joint share does not move

Lima completes the set. **Re-downloaded as SRTMGL1**, replacing the AW3D30 file that had
been there: the three department runs are compared against one another, and a dataset
difference would have sat inside every comparison as a confound indistinguishable from
a difference in the ground. 10,886 × 9,638 = **104.9 Mpx**, same 1 arc-second grid,
same `downsample_factor` 4 / `candidate_stride` 5, every transferable criterion copied.

| | Arequipa | Ancash | Lima |
| --- | --- | --- | --- |
| Mpx | 128.6 | 68.6 | 104.9 |
| **median slope** | **11.1°** | **23.0°** | **20.4°** |
| in GRAND's 3–25° | 70.3% | 52.0% | 54.6% |
| in TAMBO's 20–60° | 24.1% | 58.0% | 50.8% |
| GRAND area km² | 88,527.5 | 43,091.2 | 51,677.6 |
| GRAND capacity | 101,948 | 49,447 | 59,270 |
| GRAND acceptance | 61.6% | 54.9% | 49.6% |
| **GRAND per pixel** | 1.00× | **0.91×** | **0.72×** |
| TAMBO sites | 26 | 35 | 40 |
| TAMBO area km² | 111.9 | 174.9 | 190.9 |
| TAMBO capacity | 9,024 | 14,290 | 15,775 |
| **TAMBO per pixel** | 1.00× | **2.93×** | **2.09×** |
| joint km² | 50.2 | 75.2 | 88.3 |
| **joint as share of TAMBO's mask** | **44.9%** | **43.0%** | **46.2%** |

**Two results, and the second is the more interesting.**

**(1) The answer tracks median slope, in opposite directions for the two experiments.**
Order the regions by steepness — Arequipa 11.1°, Lima 20.4°, Ancash 23.0° — and GRAND
falls monotonically per pixel (1.00, 0.72, 0.91… not quite monotone, and Lima's dip below
Ancash is worth a look) while TAMBO rises (1.00, 2.09, 2.93). Two regions could be a
coincidence; three make it a property of the terrain rather than of the run. GRAND wants
ground gentle enough to stand an array on and steep country keeps giving it cliffs;
TAMBO wants exactly those cliffs.

**(2) The joint region is a near-constant share of TAMBO's mask — but the constant
depends on the sampling, and the strided value is not the real one.** Across the three
department runs it is 44.9%, 43.0%, 46.2%, on terrain that could hardly differ more. That
looks like a property, and §6.47 predicted it: the joint region is TAMBO-limited, so
co-location costs GRAND almost nothing.

**Then the two unbiased crops gave 74.5% (Huaylas) and 71.9% (Cajatambo).** Both at
`downsample_factor` 1 / `candidate_stride` 1, both on ground inside one of the department
boxes. So there are two constants, not one:

| sampling | joint as share of TAMBO's mask |
| --- | --- |
| 4 / 5 — Arequipa, Ancash, Lima | 44.9%, 43.0%, **46.2%** |
| 1 / 1 — Huaylas, Cajatambo | 74.5%, **71.9%** |

**The unbiased number is the true one.** Striding fragments TAMBO's mask and leaves
GRAND's untouched (§6.49), so the strided runs shrink the denominator's *quality* — what
survives is the scattered remainder, which overlaps GRAND's blob less. Roughly **three
quarters of TAMBO-viable ground is also GRAND-viable**, not four ninths.

The invariance itself survives, and is the finding worth carrying: at fixed sampling the
share barely moves across radically different terrain. But **quote ~72–75%, from the
unbiased runs, and never mix the two rows.** A Jaccard index that moves while the share
does not — 0.00057, 0.00174, 0.00170 for the departments against 0.075 and 0.138 for the
crops — is TAMBO's mask growing, not the two experiments agreeing more.

The full table is `results/region_comparison.md`, regenerated from the stores by
`tools/compare_regions.py`, and notebook 11 computes the whole comparison live.

**Caveats carried, unchanged.** All three department runs are strided and downsampled, so
every TAMBO area and capacity above is a lower bound by a terrain-dependent factor —
4.75× at Colca, **291× on the Huaylas crop** (§6.49). The ratios survive because all three
carry the bias equally; the absolute TAMBO numbers do not. Arequipa alone applies RFI
zones, worth ~2.9% of its box, which is why Ancash and Lima are both held at `none`.

### 6.52 A specific joint realization: what the code can answer today, and what it cannot

The question, asked again and worth answering properly: **can the code lay out ~100 TAMBO
units and ~10 GRAND antennas, co-located, at a chosen site?**

**It can already answer the feasibility half, unchanged.** Three things exist today:

1. :func:`~oroscope.site_searcher.count_grid_capacity` counts how many detectors of a
   given spacing and lattice fit on any boolean mask. So *"does this ground hold 100 TAMBO
   units?"* is a call, not a project.
2. The per-role slope bands are already computed, so *"where could each experiment
   stand?"* is a mask operation.
3. `oroscope-combine` gives the joint, union and membership rasters the masks come from.

So a **proposed** realization can be tested now. What cannot be done is **proposing** one:
there is no routine that places N detectors well.

**§6.47 measured this on the strided Colca joint mask. Repeated on the unbiased Cajatambo
crop — 805.1 km² of joint ground, sixteen times Colca's — the conclusion survives and
sharpens:**

| ground inside the 805.1 km² joint mask | share | km² |
| --- | --- | --- |
| GRAND's band 3–25° only | 3.9% | 31.5 |
| both bands, 20–25° | 10.9% | 87.8 |
| **TAMBO's band 20–60° only** | **84.7%** | **682.1** |

So **119.3 km² is GRAND-standable** — against Colca's 3.63 km², a very different number.
The continuum limit is 138 antennas at 1 km hexagonal spacing. And yet:

- The GRAND-standable ground inside the joint mask is **22,577 disconnected patches**.
- The largest is **1.252 km²**, against the 0.866 km² one lattice cell occupies.
- **Exactly one patch is large enough to hold a single antenna.**

**The binding constraint is patch size, not area**, and it survives removing the sampling
bias. A routine that reasons in total area would report 138 antennas where the ground
holds one. This is the single most important thing to know before writing an optimiser.

**What it would take, in order of difficulty.**

1. **Retain the score raster.** :func:`~oroscope.scoring.score_candidates` computes a score
   per candidate and only per-site aggregates reach the results file. Writing it out as a
   raster aligned with the mask is a small change and unlocks everything below.
2. **A placement routine.** With a score raster, greedy or blue-noise placement of N
   detectors maximising summed score subject to a minimum spacing is straightforward —
   **provided it is given a patch-aware feasibility test rather than an area budget.**
3. **A real objective.** The score is a *ranking proxy*, not an event rate. Optimising it
   is defensible but it is not maximising detected neutrinos, which needs the differential
   acceptance `A(E)` — still the outstanding physics ask (§9.1), and still unsafe to infer
   (§6.42).

**Two design points worth settling before any of that.**

**Optimise over the union, never the intersection.** The joint mask is TAMBO ground that
GRAND's region happens to enclose. An optimiser pointed at it will place the TAMBO units
easily and then report that the GRAND antennas cannot be placed at all — which is true of
the intersection and false of the site. At Colca the GRAND-deployable ground sat a median
0.92 km away, well inside one GRAND cell.

**The coupling term is shared line of sight, not shared footprint.** What makes a pairing
joint is that both arrays watch the same wall. That is a property the arrival scan already
computes per candidate and the combination step currently discards. Recovering it is
probably a better first move than the optimiser itself, because it is what an optimiser
would need as its objective.

Not implemented, and deliberately so: the measurement above says a naive formulation
would give a confidently wrong answer.

### 6.53 The score cut was invisible in its own funnel, and `--explain` blamed the geometry ✅ delivered

Found by audit, not by a failure. `run_arrival_scan` applies the score cut *before*
counting, and the pipeline then wrote that one post-cut number under **both** funnel
names:

```python
funnel.add("directions accepted", n_hits)
if min_score > 0:
    funnel.add(f"score >= {min_score:g}", n_hits)   # the same number
```

All twelve stored runs confirm it — the two rows are byte-identical wherever both
appear. TAMBO Arequipa: `directions accepted` 517,312, `score >= 0.35` 517,312.

**Three consequences, of rising severity.**

*The score row carries no information.* It is 100.000% of the stage above it by
construction, in every run that has one.

*The geometric acceptance is unrecorded whenever a cut is in force.* `directions
accepted` named the geometry and held the geometry-and-score product. Only runs with a
cut are affected: `min_score` 0 accepts every viable candidate, so `viable & (total >=
0)` is `viable` and every GRAND run in the store already held the right number.

*`--explain` gives advice that cannot work.* `binding_constraint` compares each stage
against the one before it, so a stage keeping exactly 100% can never be named however
much it removed. The `STAGE_KNOBS` entry pointing at `min_score` was unreachable code.
Driven with a funnel that `min_score` had emptied, the summary named `directions
accepted` and told the reader to change *the arrival window, the distance window,
`min_column_depth_gcm2` and `min_target_slope_deg`* — four knobs, none of them the one
that emptied the search. Under `score_percentile` it was worse: no row was written at
all, so a percentile keeping the top 22.8% made the arrival geometry look four times
less accepting than it is.

**Fixed.** `run_arrival_scan` takes the funnel and records the two counts it actually
decides: `directions accepted` from `viable` alone, and — only when a cut applies — a
row named for the cut that made it, `score >= 0.35` or `score in top 25%`. `STAGE_KNOBS`
gains the percentile label, so the binding constraint now names `score_percentile` when
that is what bound. Four tests pin it, including the two misdiagnoses above.

**A published figure was asserting the defect.** `figures.pipeline_stages` — the
schematic in `howitworks.rst` — hardcodes that Ancash run, so it drew *Arrival scan
1,022,530* above *Scoring 1,022,530*: two bars of identical width, the lower one
captioned "cut at `min_score`". The picture said the cut removed nothing; it removed a
great deal. The stored run cannot separate them and re-running it was not in scope, so
the two stages are drawn as one honestly-labelled bar rather than given an invented
number. A run made after this fix can be drawn as seven stages again.

**Stored results are not regenerated.** They carry the old meaning under `directions
accepted`, and `results/*/manifest.json` records the commit that produced them.
`tools/compare_regions.py` reads that key for its acceptance column, so a regeneration
against new runs would mix the two meanings — that column compares geometry-and-score
for TAMBO and geometry alone for GRAND, both before and after this change.

### 6.54 A mistyped score weight was accepted, dropped, and never mentioned ✅ delivered

Found by audit. `compose` merged the caller's weights with `if n in w` and nothing
behind it, so a key naming no component was discarded in silence. `parse_score_weights`
checked the *syntax* of `name=value` and never the name. Between them, a misspelling was
a request the tool accepted and ignored:

```
--score_weights geomag=0        # one character short of `geomagnetic`
```

Measured on a two-component product with `geomagnetic` 0.2 and `depth` 0.9: **0.18 with
the typo, 0.9 with the correct spelling.** The component the user had switched off ran
at full weight through the entire search. Nothing in the results, the funnel, the
explanation or the console recorded that the request had been dropped.

This is the same shape as §6.31 — a component that appeared when it had been switched
off — arriving by a different route, and it is the worst-behaved kind of input error in
this tool: the run completes, every number moves, and the output looks exactly like a
correct one. There is no downstream check that could catch it, because the score is not
independently predictable.

**Fixed at the choke point.** `parse_score_weights` handles both the CLI string and the
config mapping, so the name check goes there and covers both; it raises `SystemExit`
naming the offender and offering the nearest real component via `difflib`.
`scoring.SCORE_COMPONENTS` is the canonical set of ten, and a test asserts it agrees
with `explain.COMPONENT_MEANING` so the gate cannot drift into refusing legitimate
weights.

**And the quieter half.** A weight naming a *real* component that this run does not
have — `muon_shielding` on a run with shielding off — is equally inert and was equally
silent. `compose` now warns, because "that component is not in this composition" is a
different message from "that is not a component" and the user needs to hear it.

**A contradiction in the same docstring, resolved.** It claimed both "a weight of 0
excludes a component" and "ignored by `min`", and `min` did ignore weights entirely: a
component switched off by weight could still be the smallest, and so still decide the
score — the one outcome switching it off was meant to prevent. Zero now excludes in
every mode. `min` still ignores relative weights, which is correct rather than lazy: the
smallest component is the smallest however it is scaled. Excluding every component is
now an error rather than an empty `np.stack`.

613 → 622 tests.

### 6.55 The pre-flight was sized against the cheaper of the two searches ✅ delivered

Found by audit; latent, not live. `run_full_dem.py` runs one pre-flight and then two
searches, so the estimate has to cover whichever costs more. It read:

```python
downsample = max(int(c.get("downsample_factor") or 1) for c in sampling)
stride     = max(int(c.get("candidate_stride")  or 1) for c in sampling)
```

Both knobs scale memory **inversely** — a larger stride means fewer candidates, a
larger `downsample_factor` means smaller labelling arrays — so the costliest
configuration is the one with the *smallest* values and the answer is a `min`. As
written, a GRAND config at 4/5 beside a TAMBO config at 1/1 would have been pre-flighted
at 4/5 and the 1/1 run waved through unchecked.

It sits directly beneath a comment describing the same failure — the sampling hard-coded
at 4 and 5, "silently wrong for the huaylas crop, which runs at 1 and 1" — so the fix for
that bug reintroduced it in a new form. Every config pair in `config/` matches today, so
`min` and `max` agree on all of them and nothing has actually been mis-sized.

**Fixed**, and made testable: the choice is now `costliest_sampling`, a pure function
with its own examples. `tools/` had **no tests at all** before this, which is why a
regression in the one number standing between a run that does not fit and the OOM killer
went unnoticed. `tests/test_tools.py` adds nine, including one that asserts the property
rather than the implementation — no shipped config may run at a sampling finer than the
one estimated for it — and three that pin the monotonicity of
`estimate_peak_memory_gb` and `estimate_visualisation_memory_gb` in both knobs, since
`min` is only the right answer while those hold.

622 → 631 tests.

## Phase 4 — Usability *(sketch — to be scoped)*

Auto-detect `origin_lat`/`origin_lon` from the GeoTIFF tiepoint (verified present,
matching current configs to ~1e-4°); ~~rename `src/setup.py`~~ (done: `src/oroscope/fetch_dem.py`), which is not a packaging
file and whose name hijacks `pip install`; real packaging; rasterio/pyproj for CRS and
outputs; ~~`--explain` funnel report~~ (done, §6.23); parameter sweeps.

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
