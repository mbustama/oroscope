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

## 3. Phase 0 — Foundations

**Goal:** make every later change measurable. Nothing here changes results.

Phase 1 rewrites the scientific core; without a baseline we cannot tell a fix from a
regression. Phase 0 must land first.

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

Per-stage wall time and peak RSS on fixed inputs, committed as a baseline table.
Known starting point on the 2500² crop: **morphology is 6.8 s of a 9.4 s run** — the
cleanup dominates, not the physics. Phase 3 targets that number.

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

### 4.1 Terrain layer

Slope, aspect and roughness computed at an explicit `slope_baseline_m`, on a
correspondingly smoothed DEM. The scale dependence of §2.1 gets documented rather
than being an accident of `np.gradient`.

### 4.2 Azimuthal sweep engine

Replace per-candidate ray casting with a per-azimuth sweep across the raster. For a
fixed azimuth, the horizon for *every* pixel comes from one running-max walk —
sequential memory access, O(N) per azimuth, versus random access per candidate now.

Per pixel and azimuth it yields: horizon elevation angle, distance and height of the
horizon-defining terrain, and the qualifying targets within the distance window.

This is simultaneously the fix for §2.1 (fan instead of one ray) and §2.2 A
(occlusion), and it is faster than what it replaces — the rare case where the
physics fix pays for itself. Implemented with numba `prange`, bilinear sampling
instead of `int()` truncation (which currently biases every ray), and marching at DEM
resolution instead of 1 km jumps that can step over a narrow ridge.

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
column depth, energy). Getting that properly means either a table from the
collaborations' simulation chains or an explicit parameterisation with stated
assumptions. **This is the main open question for phase 1 — see §7.**

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

Two consequences worth flagging early:

- **Site shape differs fundamentally between experiments.** 5,000 units at 150 m
  triangular spacing need ~97 km² of usable canyon wall — a long strip along the
  canyon, not a compact blob. The current `min_width_km` opening with a square
  structuring element would destroy such a strip. Region/layout models must be
  per-experiment, not shared.
- **The 100 m vs 150 m spacing discrepancy** (§1.1 vs ref. [2]) should be confirmed;
  the roadmap assumes it is a deliberate starting point and keeps spacing configurable.

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
