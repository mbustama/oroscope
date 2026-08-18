# Outline — the co-location talk

**Phase 1 deliverable.** Structure and beats only; slide content comes next.
Built for **~20 minutes, 18 slides**. Slides marked ⏱ are the cuts for a 12-minute
version (drops to 13). Every number is sourced; nothing here is invented.

**Proposed titles**

1. *Slope competes, sight does not — siting a radio–particle hybrid in the Andes*
2. *Where can GRAND and TAMBO share a mountain?*
3. *One question, asked 128 million times: co-locating radio and particle detection in Peru*

---

## The spine

The talk carries **one claim** and **one map**.

> **Claim.** Two detectors that want opposite hillsides do not want opposite *views*.
> Co-location is limited by slope, not by sight — and even so, roughly three quarters
> of TAMBO-viable ground is also GRAND-viable.

> **Map.** Colca: GRAND alone → TAMBO alone → both, same frame, nothing moving.

Everything else exists to make those two land.

**The audience's prior**, which the talk must engage rather than ignore: *a radio array
wants flat ground and a particle array wants a canyon wall, so co-location is a
non-starter.* That prior is half right, and the talk's job is to say which half.

---

## Act I — Why bother (slides 1–6, ~6 min)

The brief is explicit that this comes first and is the part the audience will not
already believe. Do not rush it to get to the maps.

### 1. Title
Title, name, affiliation, date. Oroscope named as the tool, `v0.5.0`, public.

### 2. The messenger
**Lands:** the Earth is the target, the atmosphere is the calorimeter, and a mountain is
what puts them next to each other.

ν<sub>τ</sub> through rock → τ → escapes the surface → decays *in air* → extensive air
shower. The whole technique is that geometry. Say the word **Earth-skimming** once and
move.

*Visual:* `decay_and_shower` (tau decay length vs shower development).

### 3. One shower, two ways to see it
**Lands:** GRAND reads the radio emission, TAMBO reads the particles. Same shower.

The §2 comparison table, but shown as **two columns, not six rows** — the audience only
needs to retain that the numbers differ, and that they differ in kind:

| | GRAND | TAMBO |
| --- | --- | --- |
| channel | radio | particles |
| wants to see | 10–40 km away, ±3° about the horizon | 2–5 km across, ±20° |
| stands on | 3–25° slope, 1 km spacing | 20–60° slope, 150 m spacing |
| reaches | ≳100 PeV | ~3 PeV – 1 EeV |

Flag the last row deliberately — the energy ranges **overlap but are offset**. That sets
up slide 4.

### 4. Why co-locate — the argument the audience needs
**Lands:** four reasons, in decreasing order of how much a physicist will care.

1. **Cross-calibration.** Particle detection is established and its systematics are
   understood; radio detection of Earth-skimming taus is new. One site, the same sky,
   the same energy decade, two independent channels.
2. **A joint spectrum wider than either alone.** ~3 PeV to beyond 10 EeV from one site.
3. **Shared infrastructure.** Access, power, data, permits, community relations — at
   4,000 m in the Andes these are not a footnote in the budget.
4. **Hybrid events**, if the arrival directions and energies genuinely coincide.

**State the honesty line here, early:** point 4 is a motivation, not a result — Oroscope
answers a terrain question and computes no event rates. Saying this on slide 4 buys
credit that slide 14 will spend.

Then the pivot that sets up the whole rest of the talk:

> Whether the terrain *allows* co-location is a geometry question — and a geometry
> question is something a computer can answer over an entire country.

### 5. Why Peru
**Lands:** the Andes supply both landforms within tens of kilometres.

High plateau that GRAND wants; canyons among the deepest in the world that TAMBO needs.
Colca is TAMBO's chosen site. High, dry, sparsely populated.

**The counterpoint goes on this slide, not in a footnote.** Peru sits near the magnetic
equator, and the radio emission modelled here is geomagnetic only — so north–south
geometries are weakened, and a site with sin α → 0 scores exactly zero. This is the most
important caveat for a Peru-specific talk and it is stronger stated up front than
extracted in questions.

*Visual:* `canyon_geometry` — near wall, far wall, 41° walls 1.5 km deep, rim to rim
4.5 km, drawn to scale. (Strip the baked-in title; the slide carries it.)

### 6. ⏱ Why a canyon at all
**Lands:** what TAMBO is actually asking of the ground, geometrically.

Cut this if short — slide 5's figure already does most of the work.

---

## Act II — How the search works (slides 7–9, ~4 min)

### 7. One question, asked once per pixel
**Lands:** the structural claim about the method. One slide, no formulae.

> From this patch of ground, is there a **target surface** at the **right range**, in the
> **right direction**, at the **right relative orientation**, with the **right matter
> behind it**?

That is why one engine serves both experiments: **they differ in their numbers, not
their structure.** Adding an experiment is a configuration file, not a code path.

### 8. Two halves running in opposite directions
**Lands:** stages 1–4 *remove* candidates; stages 5–6 *rebuild a map*.

The brief calls this the single most misread thing in the pipeline, so it gets its own
slide and the figure already draws the split.

Say out loud only: *screen → stride → scan → score → close → prune*, and that **striding
is cost control, not a criterion**, and that **the count rises at closing**.

*Visual:* `pipeline_stages` — TAMBO over the full Ancash DEM, 68,571,090 pixels in,
789,552 out, log-width bars, the two halves bracketed.

### 9. ⏱ The arrival scan, live
**Lands:** what "walk outward and see what you hit" actually means.

`the_walk.mp4`, then `the_azimuth_fan.mp4` if time. The brief calls these the best
explanatory assets in the project for a live audience. `.gif` as fallback; test the
embed on the presenting machine before the room.

---

## Act III — The result (slides 10–15, ~7 min)

### 10. The trap
**Lands:** the confident wrong answer, and why it was wrong. This is the slide that makes
the talk memorable.

Treating both experiments' *viewing* windows as shared constraints gives the conclusion
that the two **cannot share ground at all** — a conclusion this project printed directly
above a map of the ground they demonstrably share.

Why it is wrong: GRAND scans ±3° at 10–40 km, TAMBO scans ±20° at 2–5 km. **Two
experiments looking out from the same hillside at different ranges and different
elevations are not in conflict.** Nothing is shared, so nothing competes.

Showing your own error is worth more than showing your own result. Keep it to 45 seconds
and do not be arch about it.

### 11. What actually competes
**Lands:** a pixel has one slope. Both experiments must accept that same value.

GRAND 3–25°, TAMBO 20–60° → they overlap only over **20–25°, which is 23% of the
narrower band**. Co-location is decided at the screening step, before any arrival
geometry is considered.

> **Slope competes. Sight does not.**

That line is the talk. Put it on the slide alone if the design allows.

### 12–14. The reveal — three slides, one frame
**Lands:** the map. Nothing shifts between the three; only the categories appear.

**Use Colca**, not Arequipa. It is TAMBO's chosen site, it is legible from the back of a
room, and the magenta traces the canyon rims where the roads also run — which is not a
coincidence and is worth saying aloud. Arequipa is busier and one place label is clipped;
it appears once, on slide 15, as the scale number.

1. `colca_GRAND_only.png` — a boundary enclosing most of the map. 4,446 km².
2. `colca_TAMBO_only.png` — a scatter of patches strung along the canyons. 80 km².
3. `colca_joint.png` — **123.3 km² magenta on the canyon rims. 60.7% of TAMBO's mask.**

Attribution on each: roads and place names © OpenStreetMap contributors (ODbL); DEM
SRTMGL1 via OpenTopography.

### 15. The numbers
**Lands:** the result is a range, and the range is set by array design, not by nature.

The §5 table, with **the two unbiased crops kept visually apart from the three strided
departments** — a rule, not a preference:

- Unbiased crops: **76.4%** (Cajatambo), **78.5%** (Huaylas)
- Departments, strided, therefore an under-report: **55.6 – 59.7%**

> Roughly **three quarters to four fifths** of TAMBO-viable ground is also GRAND-viable.

Arequipa as the one scale number: **88,208 km² of GRAND plateau, 1,037 km² of TAMBO
canyon, 619 km² shared.** And the asymmetry that makes the case: the joint is **0.7% of
GRAND's ground and 59.7% of TAMBO's** — co-location costs GRAND almost nothing and
consumes most of what TAMBO has.

### 16. ⏱ The partner does not have to stand on the joint
**Lands:** the practical result, and it is a more useful one than the headline.

What couples two arrays is a **shared line of sight to the same massif**, not a shared
footprint. Measured from each region's best TAMBO site: **100 GRAND antennas within
~10 km, 1,000 within 20–30 km, 5,000 within 60 km.** At Colca the nearest
GRAND-deployable ground sits a median **0.92 km** away — well inside one GRAND cell.

None of the three regions is limited by finding partner ground near a site.

---

## Act IV — What it does not tell you (slides 17–18, ~3 min)

Do not compress this. An audience of physicists will ask, the project's own
documentation is blunt about it, and volunteering it is worth more than surviving it.

### 17. Four caveats, on the slide, not in a footnote
1. **Reported area is an upper bound.** The mask is closed morphologically before area is
   measured — **2.35×** inflation at Colca against a stride-1 control.
2. **A strided TAMBO area is a lower bound**, by a terrain-dependent factor: **1.51× at
   Colca, 23.0× on the Callejón de Huaylas.** Striding is unbiased in *acceptance*
   (75.750% vs 75.736%); the whole loss is in reconstruction, where a fragmented mask
   fails the minimum-array-size cut.
3. **`min_score` is the dominant assumption.** A cut on a *product* of components, whose
   distribution piles up near zero. The shipped 0.35 is the **17.8th percentile** against
   a **median candidate score of 0.13**, and a sweep across the cut shows **no knee
   anywhere**.
4. **`solid_angle` is the weakest component at every selected site in every region** —
   16/16 Colca, 62/62 Ancash, 84/84 Lima, 85/85 Arequipa. The TAMBO result is set almost
   entirely by that one term.

*Visual:* `score_composition` for point 3; `striding_and_closing` for point 2 if there
is room. One of the two, not both.

### 18. What is missing, and what is next
- **Nothing has been checked against an external simulation.** Internally consistent,
  externally unvalidated. The cheapest available test is the Earth-absorption
  prediction — an arrival-window edge of **−4.4° at 100 PeV rising to −0.9° at 10 EeV** —
  and it has not been done. *Naming the specific test you have not run is stronger than
  admitting the general gap.*
- **Askaryan emission is not modelled** — geomagnetic only, which is why Peru's proximity
  to the magnetic equator hits north–south geometries hardest. (Callback to slide 5.)
- **Detector acceptance A(E) is not modelled**; two published *integral* curves are
  digitised, and the array-size correction cannot correct for the site.
- **Neutral-current regeneration is not modelled**, so Earth-chord suppression is
  somewhat overstated.
- **Declination** falls back to the Arequipa value rather than being per-site.
- A **national TAMBO answer** needs 1 arc-second over the whole country, and tiling.

### Closing line
Return to the claim, in one sentence, and stop:

> The terrain does not forbid it. Slope is the only thing that competes, and it still
> leaves three quarters of TAMBO's ground. The question worth asking next is not
> *whether* they can share, but *which* massif.

---

## Backup slides (after the last slide, for questions)

Cheap to make, and each answers a question that will be asked:

- The full §5 six-region table.
- `walk_mechanism` — one ray walked outward, how a direction is accepted.
- `striding_and_closing` — the cliff: a 3-px element recovers 0.04×, a 5-px one 0.68×.
- Arequipa and Huaylas joint maps.
- References: GRAND arXiv:1810.09994 (effective area, Fig. 25); TAMBO Nature Astronomy
  (2026), doi:10.1038/s41550-026-02916-4 (exposure, Fig. 3).

---

## Asset checklist

| Slide | Asset | Status |
| --- | --- | --- |
| 2 | `decay_and_shower.png` | rendered, `output/talk/` |
| 5 | `canyon_geometry.png` | rendered — strip figure title |
| 8 | `pipeline_stages.png` | rendered |
| 9 | `the_walk.mp4`, `the_azimuth_fan.mp4` | `output/animations/` — test embed |
| 12–14 | `colca_{GRAND_only,TAMBO_only,joint}.png` | `output/export_300dpi/`, 300 DPI |
| 17 | `score_composition.png`, `striding_and_closing.png` | rendered |
| backup | `walk_mechanism.png`, Arequipa + Huaylas maps | rendered / exported |

Nothing needs re-running. Everything on this list already exists.
