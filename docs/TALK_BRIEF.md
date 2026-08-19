# Handover brief — the co-location talk

**For a fresh session with no memory of this project.** Its job is a talk, not
development. Two phases: **first an outline**, then the contents and the slides in
**LibreOffice Impress**. The talk is tomorrow.

This is a **side project**. Do not audit, refactor or "improve" Oroscope. If you notice
something wrong in passing, say so in one line and carry on with the talk.

**Repository:** `~/Research/GRAND/oroscope`, public at `mbustama/oroscope`, version 0.5.0
on PyPI. **Python:** call `/home/mbustamante/anaconda3/envs/sssearch/bin/python` directly
— `conda activate sssearch` fails.

---

## 1. What the talk has to do

The owner's framing, in their order:

1. **Motivate why co-locating radio and particle detection is worth doing at all.**
   This comes first and it is the part the audience will not already believe.
2. **Motivate why Peru** is a good place to look for the opportunity.
3. **Convey the methods** — how the search actually works.
4. **Convey the results**, with plots.

The abstract the owner is writing is about **co-locating a particle detector and a radio
detector**. The talk should land one structural claim (§4) and one map (§7).

Plots are wanted throughout. There are 18 publication-resolution maps and 6 schematics
already rendered; see §8. **Do not invent numbers — every figure in §5–§7 is measured and
sourced here.**

---

## 2. The physics, in the order a slide needs it

**The messenger.** Ultra-high-energy tau neutrinos, Earth-skimming. A ντ travelling
through the Earth interacts in rock, produces a τ lepton, and if the τ escapes the surface
before decaying it decays *in air* and starts an extensive air shower. The whole
technique depends on that geometry: the Earth is the target, the atmosphere is the
calorimeter, and a mountain is what puts them next to each other.

**Two ways to see the same shower.**

| | GRAND | TAMBO |
| --- | --- | --- |
| detects | **radio** emission from the air shower | the **particles** themselves |
| wants | a wide view a few degrees below the horizon, targets 10–40 km away | a deep canyon: a near wall facing a far wall, 2–5 km across |
| deployable slope | **3–25°** | **20–60°** |
| arrival window | **±3°** around the horizon | **±20°** |
| detector spacing | 1 km antennas | **150 m** detection units |
| energy reach | ≳100 PeV | ~3 PeV – 1 EeV |

**Why co-locate.** They are sensitive to the same events by different channels, over
overlapping-but-offset energy ranges. One site serving both gives cross-calibration
between an established technique and a new one, shared access and power, and a joint
spectrum wider than either alone. The question is whether the terrain *allows* it — and
that is a geometry question a computer can answer over a whole country.

**Why Peru.** The Andes supply both landforms within a few tens of kilometres: high
plateau that GRAND wants, and canyons among the deepest in the world that TAMBO needs.
Colca Canyon is TAMBO's chosen site. The ground is high, dry and sparsely populated.
**State the honest counterpoint too:** Peru sits near the **magnetic equator**, and the
radio emission modelled here is geomagnetic, so north–south geometries are weakened —
see §9.

---

## 3. How Oroscope works — the methods slide

**One question, asked once per pixel of an elevation model:**

> From this patch of ground, is there a **target surface** at the **right range**, in the
> **right direction**, at the **right relative orientation**, with the **right matter
> behind it**?

That is what lets one engine serve both experiments: they differ in their *numbers*, not
their *structure*. Adding an experiment is a configuration, not a code path.

**Six stages.** Two halves that run in opposite directions — this is the single most
misread thing in the whole pipeline and is worth a slide of its own:

1. **Screening** — a cheap per-pixel test: slope band, altitude, aspect, distance to a
   road, radio-quiet zones. *Slope is where co-location is decided.*
2. **Striding** — keep one surviving pixel in N. **Cost control, not a criterion.**
3. **Arrival scan** — the expensive part. From each candidate, walk outward along a fan
   of bearings and find where each ray first meets terrain. A direction is *accepted* if
   the intersection lands in the distance window, the elevation window, and (optionally)
   strikes ground steep enough. One walk fills every elevation bin at once, so the
   **azimuth count sets the cost**.
4. **Scoring** — named components in [0,1], multiplied, then cut at `min_score`.
5. **Closing** — morphology fills the holes striding left. **The count rises here.**
6. **Pruning and selection** — regions too small or too poor in detectors are dropped;
   capacity is counted by actually packing a detector lattice, not by dividing area by
   spacing.

Stages 1–4 **remove** candidates. Stages 5–6 **rebuild a map** from what survived.

---

## 4. The result the talk is built around

**A pixel has one slope.** GRAND needs 3–25°, TAMBO needs 20–60°. Both must accept the
same ground, so **co-location is decided at the screening step** — and those bands barely
overlap.

**But what each asks of the *view* is in no conflict whatever.** GRAND scans ±3° around
the horizon at 10–40 km; TAMBO scans ±20° at 2–5 km. Two experiments looking out from the
same hillside at different ranges and different elevations do not compete. Treating the
viewing windows as shared constraints produced the confident and wrong conclusion that
the two "cannot share ground at all" — printed directly above the ground they
demonstrably share.

**So the answer is: they can share, and the limit is slope, not sight.**

Measured share of TAMBO's mask that is also GRAND-viable:

- **76.4% and 78.5%** on the two unbiased crops (Cajatambo, Huaylas)
- **55.6 – 59.7%** on the three departments, which are strided and therefore an
  under-report

**Roughly three quarters to four fifths of TAMBO-viable ground is also GRAND-viable, and
it is a range set by array design rather than a constant.** Co-location costs GRAND
almost nothing (the joint is under 1% of its mask) and consumes most of what TAMBO has.

**A second, practical result.** A partner array does not even have to stand on the joint
mask — what couples the two is a shared line of sight to the same massif, not a shared
footprint. Measured from each region's best TAMBO site, a GRAND array of 100 antennas
fits within ~10 km, 1,000 within 20–30 km, and 5,000 within 60 km.

---

## 5. The numbers — all measured, 2026-08-17

Six regions, both experiments and their combination, all on the same code.
`sampling` is `downsample_factor / candidate_stride`.

| region | sampling | GRAND sites / detectors / km² | TAMBO sites / detectors / km² | joint km² | share of TAMBO |
| --- | --- | --- | --- | --- | --- |
| colca | 1 / 5 | 1 / 5,315 / 4,569.4 | 16 / 10,437 / 203.0 | 123.3 | 60.7% |
| huaylas | **1 / 1** | 1 / 9,559 / 8,249.5 | 32 / 14,925 / 291.3 | 228.7 | **78.5%** |
| cajatambo | **1 / 1** | 1 / 6,457 / 5,573.8 | 44 / 39,658 / 774.5 | 591.7 | **76.4%** |
| ancash | 4 / 5 | 1 / 49,059 / 42,791.9 | 62 / 34,275 / 740.0 | 411.1 | 55.6% |
| lima | 4 / 5 | 1 / 58,669 / 51,209.0 | 84 / 42,549 / 915.4 | 509.8 | 55.7% |
| arequipa | 4 / 5 | 1 / 101,584 / 88,208.2 | 85 / 49,271 / 1,036.9 | 619.1 | 59.7% |

**Huaylas and Cajatambo are crops run unbiased** (every pixel a candidate). The three
departments are strided to fit in memory. **Never mix the two rows**, and for TAMBO
**quote the crops, never the departments** — see §6.

Arequipa is the flagship number: **88,208 km² of GRAND-viable plateau, 1,037 km² of TAMBO
canyon, 619 km² shared.**

---

## 6. Caveats that must appear on a slide, not in a footnote

The project's own documentation is deliberately blunt about these. A talk that quotes the
areas without them is over-claiming.

- **Reported area is not physics-accepted area.** The mask is closed morphologically
  before area is measured, which inflates it **2.35×** at Colca against a stride-1
  control. Read a reported area as an **upper bound**.
- **A strided TAMBO area is a lower bound**, by a terrain-dependent factor: **1.51× at
  Colca, 23.0× on the Callejón de Huaylas.** Striding is unbiased in *acceptance*
  (58.414% vs 58.415% for GRAND; 75.750% vs 75.736% for TAMBO) — the loss is entirely in
  reconstruction, where a fragmented mask fails the minimum-array-size cut.
- **`min_score` is the dominant assumption.** It is a cut on a *product* of components,
  whose distribution piles up near zero. The shipped 0.35 is the **17.8th percentile** on
  Colca terrain against a **median candidate score of 0.13**, and a sweep across the cut
  shows **no knee anywhere** — nothing in the data marks 0.35 as natural.
- **`solid_angle` is the weakest score component at every selected site in every region**
  — 16/16 at Colca, 62/62 Ancash, 84/84 Lima, 85/85 Arequipa. The TAMBO result is set
  almost entirely by that one term.

---

## 7. Plots — all already rendered, nothing to re-run

### The maps: `output/export_300dpi/` — 18 PNGs at 300 DPI (~3030×3280)

`<region>_GRAND_only.png`, `<region>_TAMBO_only.png`, `<region>_joint.png` for **colca,
huaylas, cajatambo, ancash, lima, arequipa**.

These are the three `--reveal` frames of the same overlay: everything except the
categories is identical across the three, so they **uncover a result progressively
without anything shifting between slides**. That is what they were made for. Each carries
terrain, roads, towns, a scale bar and a north arrow.

**Suggested use:** Arequipa or Colca as the three-slide reveal — GRAND alone (a boundary
enclosing most of the map) → TAMBO alone (a scatter of patches strung along the canyons)
→ both (the magenta that traces the canyon rims, where the roads also run, which is not a
coincidence).

To re-render at another size, no search is needed:

```bash
python -m oroscope.combine_experiments output/colca_full_grand output/colca_full_tambo \
    --labels GRAND TAMBO --out output/colca_render --reveal --dpi 300 \
    --roads input/roads/colca.geojson --settlements input/roads/colca_places.geojson
```

### The schematics: `oroscope.figures`, six builders

```python
from oroscope import figures
fig = figures.pipeline_stages()      # then fig.savefig("x.png", dpi=300, bbox_inches="tight")
```

| builder | what it shows | use it for |
| --- | --- | --- |
| `walk_mechanism` | one ray walked outward, how a direction is accepted | the methods slide |
| `canyon_geometry` | near wall, far wall, the angles TAMBO needs | why a canyon |
| `decay_and_shower` | tau decay length against shower development | the physics slide |
| `pipeline_stages` | the seven-stage funnel, log-width bars, two halves | **the methods slide** |
| `striding_and_closing` | the cliff: a 3-px element recovers 0.04×, a 5-px one 0.68× | the caveats slide |
| `score_composition` | why a product of components has no safe threshold | the `min_score` caveat |

### The animations: `output/animations/` — 8, each as `.mp4` and `.gif`

`the_walk`, `the_azimuth_fan`, `the_funnel`, `tau_in_rock`, `energy_window`,
`slope_criterion`, `stride_and_closing`, `product_collapse`.

**`the_walk` and `the_azimuth_fan` are the best explanatory assets in the project** for a
live audience — they show the arrival scan doing its work. Impress embeds `.mp4`; the
`.gif` is the fallback.

---

## 8. Suggested shape (the owner will revise — this is a starting point)

1. **The messenger** — Earth-skimming tau neutrinos; Earth as target, air as calorimeter.
2. **Two detectors, one shower** — radio vs particles; the comparison table of §2.
3. **Why co-locate** — cross-calibration, shared infrastructure, a wider joint spectrum.
4. **Why Peru** — plateau and deep canyon within tens of km; Colca; and the magnetic-equator caveat stated up front.
5. **The question** — the single structural question of §3, one slide, no formulae.
6. **How the search works** — `pipeline_stages`, and the two halves.
7. *(optional live asset)* — `the_walk` animation.
8. **The structural result** — a pixel has one slope; slope competes, sight does not.
9. **The reveal** — three map slides, GRAND → TAMBO → joint.
10. **The numbers** — the §5 table, crops and departments kept apart.
11. **What this does not tell you** — §6 and §9, plainly.
12. **What is next** — §9's open items.

---

## 9. What to say when asked what is missing

Be straight about these; the project is, and an audience of physicists will ask.

- **Nothing has been checked against an external simulation.** Everything is internally
  consistent and externally unvalidated. The cheapest available test is the
  Earth-absorption prediction — an arrival-window edge of **−4.4° at 100 PeV rising to
  −0.9° at 10 EeV** — and it has not been done.
- **Askaryan (charge-excess) emission is not modelled.** Only geomagnetic emission is, so
  a geometry with sin α → 0 scores exactly zero and a product composition rejects the
  site outright. **Peru is near the magnetic equator, so north–south geometries are hit
  hardest.** This is the most important caveat for a Peru-specific talk.
- **The detector acceptance *A(E)* is not modelled.** Two published *integral* curves are
  digitised in `data/`, and the array-size correction cannot correct for the site.
- **Neutral-current regeneration is not modelled**, so the Earth-chord suppression is
  somewhat overstated.
- **Declination** falls back to the Arequipa value rather than being per-site.
- A **national TAMBO answer** needs 1 arc-second over the whole country and tiling.

---

## 10. Practical notes for building the deck

- **LibreOffice Impress** is the target. Check availability with `soffice --version`.
  Building `.odp` programmatically is possible but fiddly; producing the content as
  structured markdown first and letting the owner paste, or generating via a template,
  are both reasonable — **ask the owner which they prefer before building slides.**
- **Figures:** no titles on the figures themselves (the slide carries the title); legend
  outside the axes at the top; scale bar and north arrow on every map; roads green;
  capitalise the first word of every label; attribution in the caption. The maps in
  `export_300dpi/` already follow all of this.
- **Attribution:** the roads and place names are **© OpenStreetMap contributors, ODbL**.
  The DEMs are **SRTMGL1** via OpenTopography. Both belong on the map slides.
- **References:** GRAND — arXiv:1810.09994, effective area Fig. 25. TAMBO — Nature
  Astronomy (2026), doi:10.1038/s41550-026-02916-4, exposure Fig. 3.
- The owner prefers **measured numbers over assertions**, and will ask where a figure
  came from. Every number in this brief is traceable to `results/<region>_full/` or to
  `docs/ROADMAP.md` §6.69–§6.74.
