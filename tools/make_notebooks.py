#!/usr/bin/env python
"""
Generates the tutorial notebooks in ``notebooks/``.

The notebooks are written from here rather than by hand for one reason: they must not
drift from the code they teach. Their prose lives in this file next to the calls it
describes, so changing an API and forgetting a notebook shows up as a failing execution
in CI rather than as a page that still reads plausibly and no longer works.

The generated notebooks are committed with their outputs, so they render on GitHub
without being run. Regenerate and re-execute with::

    python tools/make_notebooks.py
    jupyter nbconvert --execute --inplace notebooks/*.ipynb

Every notebook builds its own terrain. None of them needs the real DEM, which is
gitignored and a quarter of a gigabyte -- a tutorial that only its author can run is
not a tutorial.
"""

from __future__ import annotations

import pathlib

import nbformat as nbf

HERE = pathlib.Path(__file__).resolve().parents[1]
OUT = HERE / "notebooks"

# Put src/ on the path, so the notebooks work in a clone whether or not the package has
# been installed. Repeated verbatim at the top of every notebook so each one stands
# alone, which is how people actually open them.
PREAMBLE = """import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join('..', 'src')))

import numpy as np
import matplotlib.pyplot as plt
"""

FOOTER = """---

*Part of the [Oroscope](https://github.com/mbustama/oroscope) tutorials. \
{prev}{sep}{nxt}Full API reference: \
[oroscope docs](https://mbustama.github.io/oroscope/functions.html).*"""


def notebook(cells):
    nb = nbf.v4.new_notebook()
    nb.cells = [nbf.v4.new_markdown_cell(s) if kind == "md" else nbf.v4.new_code_cell(s)
                for kind, s in cells]
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
    }
    return nb


def footer(prev=None, nxt=None):
    p = f"Previous: [{prev[1]}]({prev[0]}). " if prev else ""
    n = f"Next: [{nxt[1]}]({nxt[0]}). " if nxt else ""
    return FOOTER.format(prev=p, sep="", nxt=n)


# --------------------------------------------------------------------------- 01
NB01 = [
("md", """# 1. Getting started

**Oroscope** searches digital elevation models for ground that could host a
particle-astrophysics observatory. It answers one question, and the fact that it is a
single question is what lets one engine serve experiments that look nothing alike:

> From this patch of ground, is there a target surface at the **right range**, in the
> **right direction**, at the **right relative orientation**, with the **right matter
> behind it**?

GRAND wants terrain a few degrees below the horizon and tens of kilometres away, to
catch radio from air showers started by Earth-skimming tau neutrinos. TAMBO wants a
canyon wall two to five kilometres across, to catch the particles themselves. They
differ in their *numbers*, not in their *structure*.

This notebook builds a piece of terrain, scans it, and reads the result. Nothing here
needs a real DEM."""),
("code", PREAMBLE + """
import arrival_scan
import site_searcher as ss

print("modules loaded")"""),
("md", """## A piece of terrain

A DEM is just a 2-D array of elevations plus a statement of how big a pixel is on the
ground. Here is a plain west-facing slope with a ridge to its east — the geometry a
detector on the slope would use to look *at* the ridge."""),
("code", """n = 400                      # pixels on a side
cell_deg = 1 / 3600.0        # 1 arc-second, as SRTM and AW3D30 supply

# Metric pixel sizes differ on each axis away from the equator, because a degree of
# longitude shrinks with the cosine of the latitude. resolve_grid_geometry does that.
grid = ss.resolve_grid_geometry("no-such-file.tif", -15.6, cell_size_deg=cell_deg)
print(f"pixel: {grid.cell_size_y:.1f} m north-south, {grid.cell_size_x:.1f} m east-west")
print(f"resolution source: {grid.source}")

cols = np.arange(n)[None, :].repeat(n, 0)
x_m = cols * grid.cell_size_x

# A ridge 1.2 km high centred 6.5 km east, on a plain at 2 km
z = (2000.0 + 1200.0 * np.exp(-((x_m - 6500.0) / 700.0) ** 2)).astype(np.float32)

fig, ax = plt.subplots(figsize=(7, 2.6))
ax.plot(x_m / 1000, z[0], color="#7A6A4F")
ax.fill_between(x_m[0] / 1000, 1900, z[0], color="#D9CDB8")
ax.set_xlabel("east (km)")
ax.set_ylabel("elevation (m)")
ax.set_title("a ridge to look at")
ax.set_ylim(1900, 3400)
ax.spines[["top", "right"]].set_visible(False)
plt.show()"""),
("md", """## Scanning one candidate

A *candidate* is a pixel the search is considering, given as `[row, col, aspect_deg]`.
The aspect is the downhill direction, and the scan fans its azimuths around it.

Put a detector on the plain west of the ridge, facing east."""),
("code", """candidate = np.array([[200.0, 60.0, 90.0]])     # row, col, aspect: due east

out = arrival_scan.scan(
    candidate, z, grid,
    n_azimuths=1, half_width_deg=0.0, use_aspect=True,   # look along the aspect only
    elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
    min_dist_km=1.0, max_dist_km=10.0, max_range_m=10000.0,
)

for key in ("cells", "solid_angle_sr", "mean_distance_m", "max_depth_gcm2",
            "horizon_deg", "target_slope_deg"):
    print(f"{key:>18}: {out[key][0]:,.3f}")"""),
("md", """What those mean:

- **`cells`** — how many (azimuth, elevation) directions were accepted.
- **`solid_angle_sr`** — the accepted solid angle. This is the quantity an aperture is
  proportional to, so it is the closest single number to "how good is this site".
- **`mean_distance_m`** — mean distance to the first terrain the ray strikes, which is
  where a tau would exit the rock.
- **`max_depth_gcm2`** — the most rock any accepted direction has behind it.
- **`horizon_deg`** — the highest terrain angle seen. Note it is **positive**: the
  ridge stands above the horizontal, so the whole ±3° window sits *below* the local
  horizon. That is the normal case in mountains, and it is why a simple hit/no-hit test
  discriminates so poorly.
- **`target_slope_deg`** — how steep the struck terrain is, along the arrival azimuth.

## Distance decides what is visible

Move the detector and the accepted geometry changes. The scan only accepts a direction
whose first intersection falls inside the distance window."""),
("code", """for col in (20, 60, 120, 180):
    cand = np.array([[200.0, float(col), 90.0]])
    r = arrival_scan.scan(cand, z, grid, n_azimuths=1, half_width_deg=0.0,
                          use_aspect=True, elev_min_deg=-3.0, elev_max_deg=3.0,
                          n_elev_bins=12, min_dist_km=1.0, max_dist_km=10.0,
                          max_range_m=10000.0)
    gap_km = (6500.0 - col * grid.cell_size_x) / 1000.0
    print(f"col {col:>4}  {gap_km:5.2f} km from the ridge   "
          f"cells {int(r['cells'][0]):>3}   Omega {r['solid_angle_sr'][0]:.3f} sr   "
          f"d {r['mean_distance_m'][0]:>7,.0f} m")"""),
("md", """The closest candidate accepts nothing: the ridge is inside the 1 km minimum, so every
direction that strikes it is rejected as too near. That minimum exists because the tau
needs room to decay and the shower room to develop.

## Where to go next

- **[2. The arrival scan](02_the_arrival_scan.ipynb)** — how one walk fills every
  elevation bin, and what the observables actually measure.
- **[3. The physics toolkit](03_physics_toolkit.ipynb)** — the closed-form physics,
  usable on its own with no terrain at all."""),
("md", footer(nxt=("02_the_arrival_scan.ipynb", "The arrival scan"))),
]

# --------------------------------------------------------------------------- 02
NB02 = [
("md", """# 2. The arrival scan

The engine at the centre of Oroscope. This notebook shows *how* it works, because the
mechanism explains both what it costs and what it can and cannot tell you.

Trace a ray **backwards** from a candidate along an arrival direction — azimuth φ,
elevation θ. Rays above the local horizon escape to the sky. Rays below it strike
terrain, and that first intersection is where the tau left the rock."""),
("code", PREAMBLE + """
import arrival_scan
import figures
import site_searcher as ss

grid = ss.resolve_grid_geometry("no-such-file.tif", -15.6, cell_size_deg=1/3600)"""),
("md", """## One walk fills every elevation bin

Writing the terrain's elevation angle at ground distance $d$ as

$$\\theta_{\\rm terrain}(d) = \\arctan\\!\\left(\\frac{z(d) - d^2/2R - z_0}{d}\\right)$$

a ray at angle θ first meets terrain at the smallest $d$ where
$\\theta_{\\rm terrain}(d) \\ge \\theta$. Since the **running maximum** of
$\\theta_{\\rm terrain}$ only increases, each new maximum claims a contiguous band of
elevation bins — so one pass fills them all.

That is why elevation binning is nearly free and the **azimuth count** is what sets the
cost of a search."""),
("code", """figures.walk_mechanism()
plt.show()"""),
("md", """The $d^2/2R$ term is the Earth's curvature dropping distant ground below a straight
line. It uses the **true** Earth radius: the neutrino and the tau are not refracted.
Only the radio path gets the inflated 4/3 radius, and only for the Fresnel term."""),
("code", """for k in (1.0, 4/3):
    R = arrival_scan.earth_radius_for_k(k)
    drop = (80_000.0 ** 2) / (2 * R)
    print(f"k = {k:.2f}:  R = {R/1000:>6,.0f} km   drop over 80 km = {drop:5.0f} m")"""),
("md", """## The ground at your own feet

A consequence that is easy to forget, and which broke four separate tests in this
project's history: **a detector standing on the ground has every steeply downward
direction blocked by the ground it stands on.**

On a uniform slope, a ray angled down more steeply than the slope intersects terrain
within a pixel or two."""),
("code", """# A plane sloping down toward the east at 20 degrees
n = 300
cols = np.arange(n)[None, :].repeat(n, 0)
slope_deg = 20.0
z_plane = (3000.0 - cols * grid.cell_size_x * np.tan(np.radians(slope_deg))).astype(np.float32)

cand = np.array([[150.0, 30.0, 90.0]])
r = arrival_scan.scan(cand, z_plane, grid, n_azimuths=1, half_width_deg=0.0,
                      use_aspect=True, elev_min_deg=-40.0, elev_max_deg=5.0,
                      n_elev_bins=45, min_dist_km=0.0, max_dist_km=6.0,
                      max_range_m=6000.0)
print(f"terrain slopes down at {slope_deg:.0f} deg")
print(f"horizon seen from the candidate: {r['horizon_deg'][0]:.1f} deg")
print(f"mean distance to first terrain: {r['mean_distance_m'][0]:,.0f} m")
print()
print("Rays steeper than the slope hit immediately; shallower ones never hit at all.")"""),
("md", """## Column depth accumulates over the whole walk

The ray at angle θ is underground wherever $\\theta_{\\rm terrain}(d) > \\theta$, so
binning the terrain angle and taking a suffix sum gives the underground path length for
every bin at once. Rays crossing several ridges accumulate **all** the rock they
traverse, not only the first chord.

There is a caveat worth knowing, and it is documented as a limitation: the walk stops at
`max_dist_km`, so for a short-range search the depth reported is bounded by how far is
left to walk rather than by the target's real thickness."""),
("code", """n = 400
cols = np.arange(n)[None, :].repeat(n, 0)
x_m = cols * grid.cell_size_x
z_two = (2000.0
         + 900.0 * np.exp(-((x_m - 4000.0) / 500.0) ** 2)
         + 1400.0 * np.exp(-((x_m - 9000.0) / 800.0) ** 2)).astype(np.float32)

cand = np.array([[200.0, 30.0, 90.0]])
for max_km in (6.0, 12.0, 20.0):
    r = arrival_scan.scan(cand, z_two, grid, n_azimuths=1, half_width_deg=0.0,
                          use_aspect=True, elev_min_deg=-3.0, elev_max_deg=3.0,
                          n_elev_bins=12, min_dist_km=1.0, max_dist_km=max_km,
                          max_range_m=max_km * 1000)
    print(f"walk to {max_km:>4.0f} km:  max depth {r['max_depth_gcm2'][0]:>12,.0f} g/cm^2   "
          f"cells {int(r['cells'][0]):>3}")"""),
("md", """## Azimuths are the cost

One profile walk per (candidate, azimuth), regardless of how finely the elevation
window is sampled. So doubling the elevation bins is nearly free; doubling the azimuths
is not."""),
("code", """import time

cands = np.column_stack([np.repeat(np.arange(50, 250), 20),
                         np.tile(np.arange(20, 40), 200),
                         np.full(4000, 90.0)]).astype(np.float64)

for n_az, n_bins in ((1, 12), (1, 48), (9, 12)):
    t = time.perf_counter()
    arrival_scan.scan(cands, z_two, grid, n_azimuths=n_az, half_width_deg=60.0,
                      elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=n_bins,
                      min_dist_km=1.0, max_dist_km=10.0, max_range_m=10000.0)
    dt = time.perf_counter() - t
    print(f"{n_az} azimuth(s), {n_bins:>2} bins: {dt*1000:7.1f} ms")

print()
print("Four times the bins costs almost nothing; nine times the azimuths costs nine times.")"""),
("md", """*(The first timing includes JIT compilation, so read the relative numbers of the last
two rather than the absolute ones.)*

## Where to go next

- **[3. The physics toolkit](03_physics_toolkit.ipynb)** — what the geometry gets
  multiplied by.
- **[4. Criteria and scoring](04_criteria_and_scoring.ipynb)** — turning observables
  into a comparable number."""),
("md", footer(prev=("01_getting_started.ipynb", "Getting started"),
              nxt=("03_physics_toolkit.ipynb", "The physics toolkit"))),
]

# --------------------------------------------------------------------------- 03
NB03 = [
("md", """# 3. The physics toolkit

`physics` has **no terrain in it**. It is closed-form physics that the scan needs but
cannot measure from a DEM: how much atmosphere a shower develops in, how much Earth a
neutrino crossed, how far a tau travels before it decays.

It is usable entirely on its own, which is the point of this notebook."""),
("code", PREAMBLE + """
import physics
import figures"""),
("md", """## The tau decay length

$L = (E/m_\\tau)\\,c\\tau$. This sets the scale of every useful detector-to-target
distance, and it is exactly analytic — no simulation input at all."""),
("code", """for e in (1.0, 3.0, 10.0, 100.0, 1000.0):
    print(f"{e:>7.0f} PeV   L = {physics.tau_decay_length_m(e):>10,.0f} m")

print()
print("TAMBO's published 50 m - 5 km tau range corresponds to about 1-100 PeV,")
print("and GRAND's 10-80 km window to roughly 0.2-1.6 EeV.")"""),
("md", """This is why a **single energy cannot stand in for a spectrum** when a search depends on
the tau decaying inside a gap. Across a 3 km canyon crossing the probability runs from
essentially one to a few per cent over one experiment's energy reach."""),
("code", """figures.decay_and_shower()
plt.show()"""),
("md", """## Shower development is measured in grammage, not metres

Air at 4000 m is a third thinner than at sea level, so a site search comparing
candidates at different altitudes while measuring path length in metres is comparing
unlike things."""),
("code", """for alt in (0.0, 2000.0, 4000.0):
    x = physics.slant_grammage_gcm2(alt, 0.0, 20000.0)
    rho = physics.air_density_kgm3(alt)
    print(f"20 km horizontal at {alt:>6,.0f} m:  rho = {rho:.3f} kg/m^3   X = {x:>6,.0f} g/cm^2")"""),
("md", """What counts as "enough" depends on **what is being detected**, and the two cases differ:

- **Radio** — emission comes from around shower maximum and then propagates through air
  that is transparent at 50–200 MHz. Being far past maximum costs nothing, so the
  criterion is a *threshold*.
- **Particles** — the charged-particle content peaks at maximum and dies after, so the
  criterion is genuinely a *band*.

The band follows from the primary energy through the shower profile."""),
("code", """x = np.linspace(1, 2000, 500)
fig, ax = plt.subplots(figsize=(7, 3.2))
for e, c in ((3.0, "#C8901A"), (55.0, "#B02A25"), (1000.0, "#7B2D6B")):
    xmax = float(physics.shower_maximum_gcm2(e))
    ax.plot(x, physics.shower_size_fraction(x, xmax), color=c,
            label=f"{e:g} PeV  ($X_{{\\\\rm max}}$ = {xmax:.0f})")
ax.axvspan(170, 390, color="#2C6E8F", alpha=0.12, lw=0)
ax.text(280, 0.92, "what a Colca\\ncrossing supplies", ha="center", fontsize=8,
        color="#2C6E8F")
ax.set_xlabel("atmospheric depth traversed (g/cm$^2$)")
ax.set_ylabel("particle content / peak")
ax.legend(frameon=False, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
plt.show()

lo, hi = physics.grammage_band_from_energy(3.0, 1000.0, fraction=0.1)
print(f"band at a tenth of peak content, 3 PeV - 1 EeV: {lo:.0f} - {hi:.0f} g/cm^2")"""),
("md", """That shaded strip is the whole siting problem for a canyon experiment: a crossing
supplies only the air its own width contains. The default particle band of
$(X_{\\rm max},\\,4X_{\\rm max})$ = 700–2800 g/cm² rejects **every** canyon.

## The Earth chord

For downgoing directions the neutrino has crossed a chord $2R\\sin\\theta$, which dwarfs
local topography."""),
("code", """lam = physics.neutrino_interaction_length_gcm2(1000.0)
print(f"interaction length at 1 EeV: {lam:.2e} g/cm^2\\n")
for theta in (-0.5, -1.0, -2.0, -3.0, -5.0):
    chord_km = physics.earth_chord_m(theta) / 1000
    surv = physics.neutrino_survival(theta, lam)
    print(f"{theta:>5.1f} deg:  chord {chord_km:>6,.0f} km   survival {surv:.3f}")"""),
("md", """This is a **falsifiable prediction**: the effective arrival window should narrow with
energy, its lower edge climbing toward the horizon."""),
("code", """for e in (100.0, 1000.0, 10000.0):
    cut = physics.earth_absorption_cutoff_deg(e)
    print(f"{e:>7,.0f} PeV:  window lower edge at {cut:+.1f} deg")

print()
print("If a collaboration's simulated window does not narrow this way,")
print("one of the two treatments has the absorption wrong.")"""),
("md", """## Geomagnetic emission

Radio emission goes as $|\\vec v \\times \\vec B|$: a shower travelling *along* the field
radiates almost none of it. Peru sits near the magnetic equator, where the field is
nearly horizontal and northward — so north–south showers are strongly suppressed."""),
("code", """dec, inc = physics.default_field_for_site(-16.4, -71.5)
B = physics.geomagnetic_unit_vector(dec, inc)
print(f"Arequipa field: declination {dec:.1f} deg, inclination {inc:.1f} deg\\n")

for name, az in (("north", 0.0), ("north-east", 45.0), ("east", 90.0)):
    print(f"{name:>10}-facing target:  sin(alpha) = "
          f"{physics.geomagnetic_sin_alpha(az, 0.0, B):.3f}")

print("\\nEast-facing targets are worth several times north-facing ones,")
print("which no purely geometric measure can see.")"""),
("md", """## Where to go next

- **[4. Criteria and scoring](04_criteria_and_scoring.ipynb)** — combining all of this
  into one comparable number, and the trap in doing so."""),
("md", footer(prev=("02_the_arrival_scan.ipynb", "The arrival scan"),
              nxt=("04_criteria_and_scoring.ipynb", "Criteria and scoring"))),
]

# --------------------------------------------------------------------------- 04
NB04 = [
("md", """# 4. Criteria and scoring

Every criterion becomes a component in $[0, 1]$ — a band, a saturating function, or a
ramp — and the components are composed into one number.

This notebook covers the shapes, the composition, and **two traps that are easy to fall
into and hard to notice**."""),
("code", PREAMBLE + """
import scoring"""),
("md", """## The three shapes"""),
("code", """x = np.linspace(0, 20, 400)
fig, axes = plt.subplots(1, 3, figsize=(10, 2.7))

axes[0].plot(x, scoring.band_score(x, 6, 14), color="#B02A25")
axes[0].set_title("band_score(x, 6, 14)", fontsize=9)

axes[1].plot(x, scoring.saturating_score(x, 4), color="#146B54")
axes[1].set_title("saturating_score(x, 4)", fontsize=9)

axes[2].plot(x, scoring.ramp_score(x, 4, 12), color="#2C6E8F")
axes[2].set_title("ramp_score(x, 4, 12)", fontsize=9)

for ax in axes:
    ax.set_ylim(-0.05, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.show()"""),
("md", """A **band** is for a quantity with an optimum on both sides. Column depth is the
motivating case: the tau must be *produced*, which needs rock, and must *escape*, which
limits how much.

A **saturating** score is for "more is better with diminishing returns and no natural
maximum" — accepted solid angle.

## Trap 1: a saturating score with the wrong scale stops discriminating

`solid_angle_half_sr` defaults to 0.05 sr, which suits GRAND. An experiment looking
across a canyon sees 0.2–1.5 sr, and against that the term saturates and carries no
information at all."""),
("code", """omega = np.linspace(0.0, 1.6, 300)
fig, ax = plt.subplots(figsize=(7, 3))
for half, c, lbl in ((0.05, "#B02A25", "0.05 sr (GRAND default)"),
                     (0.80, "#146B54", "0.80 sr (canyon scale)")):
    ax.plot(omega, scoring.saturating_score(omega, half), color=c, label=lbl)
ax.axvspan(0.2, 1.5, color="#999", alpha=0.15, lw=0)
ax.text(0.85, 0.25, "range actually observed\\nacross a canyon", ha="center", fontsize=8)
ax.set_xlabel("accepted solid angle (sr)")
ax.set_ylabel("score")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)
plt.show()

print("Over 0.2-1.5 sr the default spans "
      f"{scoring.saturating_score(1.5, 0.05) - scoring.saturating_score(0.2, 0.05):.3f} of score,")
print("against "
      f"{scoring.saturating_score(1.5, 0.8) - scoring.saturating_score(0.2, 0.8):.3f} "
      "at the canyon scale.")"""),
("md", """## Trap 2: a product score has no safe threshold

The default composition is a **product**. Multiply six components each in $[0,1]$ and
the result piles up near zero, so a threshold anywhere in the middle sits on a cliff."""),
("code", """rng = np.random.default_rng(0)
parts = {f"c{i}": rng.uniform(0.3, 1.0, 20000) for i in range(6)}

fig, axes = plt.subplots(1, 2, figsize=(10, 3))
for ax, mode in zip(axes, ("product", "mean")):
    total = scoring.compose(parts, mode)
    ax.hist(total, bins=60, color="#2C6E8F", alpha=0.75)
    ax.set_title(f"composition = '{mode}'", fontsize=9)
    ax.set_xlabel("score")
    ax.set_xlim(0, 1)
    ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.show()

total = scoring.compose(parts, "product")
for cut in (0.0, 0.1, 0.2, 0.35, 0.5):
    print(f"min_score {cut:.2f}  keeps {100 * (total >= cut).mean():5.1f}% of candidates")"""),
("md", """Measured on a real search, that cliff took the reported capacity from 45 928 detector
positions at `min_score` 0.0, to 2056 at 0.35, to **zero** at 0.5. A result that swings
by that much across a plausible range of a threshold is a result about the threshold.

**Prefer ranking sites and taking the best $N$** over thresholding a product. A weighted
geometric mean also spreads the distribution.

## Scoring candidates

`score_candidates` takes the scan's observables and returns both the total and the named
components — so a weak site can be *attributed* rather than merely ranked."""),
("code", """obs = {
    "cells": np.array([12, 40, 3]),
    "solid_angle_sr": np.array([0.10, 0.90, 0.05]),
    "mean_distance_m": np.array([2.0e4, 3.0e4, 6.0e4]),
    "max_depth_gcm2": np.array([8.0e5, 4.0e6, 2.0e5]),
}

total, parts = scoring.score_candidates(obs, {"solid_angle_half_sr": 0.3})
print("component scores:")
for name, vals in sorted(parts.items()):
    print(f"  {name:>14}: " + "  ".join(f"{v:5.3f}" for v in vals))
print("\\n         total: " + "  ".join(f"{v:5.3f}" for v in total))"""),
("md", """## Where to go next

- **[5. GRAND and TAMBO](05_grand_and_tambo.ipynb)** — the same engine, two
  experiments, and what actually differs between them."""),
("md", footer(prev=("03_physics_toolkit.ipynb", "The physics toolkit"),
              nxt=("05_grand_and_tambo.ipynb", "GRAND and TAMBO"))),
]

# --------------------------------------------------------------------------- 05
NB05 = [
("md", """# 5. GRAND and TAMBO

The claim this project rests on is that GRAND and TAMBO ask the **same structural
question** and differ only in their numbers. This notebook shows what those numbers are
and why each one is what it is."""),
("code", PREAMBLE + """
import arrival_scan
import figures
import site_searcher as ss

grid = ss.resolve_grid_geometry("no-such-file.tif", -15.6, cell_size_deg=1/3600)"""),
("md", """## What differs

| | GRAND | TAMBO | why |
|---|---|---|---|
| slope band | 3–25° | 20–60° | deployable ground vs a canyon wall |
| far-wall slope | unset | ≥ 25° | TAMBO's tau exits a *wall* |
| distance | 10–40 km | 2–5 km | to the horizon, or across a canyon |
| arrival window | ±3° | ±20° | Earth-skimming, vs a wall subtending tens of degrees |
| spacing | 1 km | 100 m | radio antennas vs particle detectors |
| `min_width_km` | 2.0 | 0.0 | a compact array vs a strip along a wall |
| Fresnel | 50 MHz | none | radio propagation, or none |
| geomagnetic | on | off | emission goes as $|v\\times B|$; particles do not care |
| grammage | threshold | band | radio propagates; particle content dies |

Every one of those is a configuration key. None is a code path.

## The canyon geometry

Drawn to scale — unlike a GRAND search, a canyon needs no vertical exaggeration, and
the arrival directions really do span tens of degrees."""),
("code", """figures.canyon_geometry()
plt.show()"""),
("md", """## Building a canyon, and scanning across it

The wall slope is a parameter, so the geometry the scan recovers can be checked against
the number the terrain was built with."""),
("code", """def canyon(n, cell_x, floor_width_m=1000.0, depth_m=1500.0, wall_slope_deg=40.6):
    \"\"\"Two opposing walls of known slope, with a flat floor between them.\"\"\"
    cols = np.arange(n, dtype=np.float64)[None, :].repeat(n, 0)
    x = cols * cell_x
    centre = (n * cell_x) / 2.0
    from_edge = np.abs(x - centre) - floor_width_m / 2.0
    rise = np.clip(from_edge, 0.0, None) * np.tan(np.radians(wall_slope_deg))
    return (3500.0 - depth_m + np.clip(rise, 0.0, depth_m)).astype(np.float32)


n = 400
for wall in (25.0, 40.6):
    z = canyon(n, grid.cell_size_x, wall_slope_deg=wall)
    col = int((n * grid.cell_size_x / 2.0 - 1200.0) / grid.cell_size_x)
    cand = np.array([[n // 2, float(col), 90.0]])       # part-way down the west wall
    r = arrival_scan.scan(cand, z, grid, n_azimuths=1, half_width_deg=0.0,
                          use_aspect=True, elev_min_deg=-25.0, elev_max_deg=25.0,
                          n_elev_bins=25, min_dist_km=0.3, max_dist_km=8.0,
                          max_range_m=8000.0, min_target_slope_deg=wall - 5.0)
    print(f"built a {wall:>4.1f} deg wall  ->  measured "
          f"{r['target_slope_deg'][0]:>4.1f} deg   "
          f"({int(r['cells'][0])} directions accepted)")"""),
("md", """## Why the far wall needs its own criterion

Without it, the scan asks only whether *rock* lies at the right range and bearing — and
on real mountainous terrain something nearly always does. Before this criterion existed,
**92% of Andean candidates passed** a canyon-shaped test.

Watch what the floor does: a shallow wall stops qualifying."""),
("code", """z = canyon(n, grid.cell_size_x, wall_slope_deg=15.0)      # a gentle valley
col = int((n * grid.cell_size_x / 2.0 - 1200.0) / grid.cell_size_x)
cand = np.array([[n // 2, float(col), 90.0]])

for floor in (None, 10.0, 30.0):
    kw = {} if floor is None else {"min_target_slope_deg": floor}
    r = arrival_scan.scan(cand, z, grid, n_azimuths=1, half_width_deg=0.0,
                          use_aspect=True, elev_min_deg=-25.0, elev_max_deg=25.0,
                          n_elev_bins=25, min_dist_km=0.3, max_dist_km=8.0,
                          max_range_m=8000.0, **kw)
    label = "no far-wall criterion" if floor is None else f"require >= {floor:.0f} deg"
    print(f"{label:>24}:  {int(r['cells'][0]):>3} directions accepted")

print("\\nA 15 deg valley satisfies a 10 deg floor and fails a 30 deg one,")
print("which is exactly the discrimination a canyon search needs.")"""),
("md", """## An honest caveat

Note that the *unfiltered* mean target slope is **not** the wall slope: rays aimed lower
strike the flat canyon floor, whose slope really is zero, so the mean over all accepted
directions is a mixture. Filtering is what isolates the wall.

This broke a first draft of the tests, and it is worth stating plainly rather than
discovering."""),
("code", """z = canyon(n, grid.cell_size_x, wall_slope_deg=40.6)
cand = np.array([[n // 2, float(col), 90.0]])
base = arrival_scan.scan(cand, z, grid, n_azimuths=1, half_width_deg=0.0,
                         use_aspect=True, elev_min_deg=-25.0, elev_max_deg=25.0,
                         n_elev_bins=25, min_dist_km=0.3, max_dist_km=8.0,
                         max_range_m=8000.0)
cut = arrival_scan.scan(cand, z, grid, n_azimuths=1, half_width_deg=0.0,
                        use_aspect=True, elev_min_deg=-25.0, elev_max_deg=25.0,
                        n_elev_bins=25, min_dist_km=0.3, max_dist_km=8.0,
                        max_range_m=8000.0, min_target_slope_deg=35.0)
print("wall built at 40.6 deg")
print(f"  unfiltered mean target slope: {base['target_slope_deg'][0]:5.1f} deg  <- a mixture")
print(f"  filtered to wall hits:        {cut['target_slope_deg'][0]:5.1f} deg  <- the wall")"""),
("md", """## Where to go next

- **[6. Combining experiments](06_combining_and_sensitivity.ipynb)** — where both are
  viable at once, and how firm any of it is."""),
("md", footer(prev=("04_criteria_and_scoring.ipynb", "Criteria and scoring"),
              nxt=("06_combining_and_sensitivity.ipynb", "Combining and sensitivity"))),
]

# --------------------------------------------------------------------------- 06
NB06 = [
("md", """# 6. Combining experiments, and how firm a result is

Two things that matter once a search produces numbers someone might act on: **where two
experiments can share ground**, and **how much the answer depends on assumptions**.

The second is the more important, and it is the one most easily skipped."""),
("code", """import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join('..', 'src')))

import numpy as np

import combine_experiments as combine
import physics"""),
("md", """## Combining is an overlay, and alignment is not optional

Each experiment is one run of the searcher with its own configuration, so combining them
is an overlay of the masks those runs produce. Three questions get different answers:

- **joint** — terrain satisfying *every* experiment. One site, one road, one power feed,
  two experiments.
- **union** — terrain satisfying *any*. How much of the region is useful to the
  programme as a whole.
- **each alone** — and what each would lose by being confined to the joint area.

The inputs must be pixel-aligned: same shape, same pixel size, same corner. That is
checked and **refused** rather than resampled, because two runs on differently-cropped
DEMs would silently compare the wrong ground."""),
("code", """# What the check actually compares: the six affine terms of the world file
world_a = (1/3600, 0.0, 0.0, -1/3600, -72.4, -15.3)
world_b = (1/3600, 0.0, 0.0, -1/3600, -72.1, -15.3)      # shifted 0.3 deg east

runs = [{"dir": "run_a", "mask": np.zeros((10, 10), bool), "world": world_a},
        {"dir": "run_b", "mask": np.zeros((10, 10), bool), "world": world_b}]

try:
    combine.check_alignment(runs)
except SystemExit as exc:
    print("refused, correctly:\\n")
    print(exc)"""),
("md", """Same shape, same pixel size, different ground — and it says so rather than overlaying
them.

## Reading a co-location result

On the Colca crop, with both experiments run over identical terrain:

| | area | sites | capacity | of its own area in the joint |
|---|---|---|---|---|
| GRAND | 4580.2 km² | 1 | 5317 | 1.2% |
| TAMBO | 93.1 km² | 17 | 10 878 | 58.9% |
| **joint** | 54.9 km² | | | Jaccard 0.012 |

The interesting part is *why* the joint is small. Two thirds of TAMBO-viable ground is
also GRAND-viable, but the two deployable **slope bands barely overlap** — GRAND's 3–25°
against Colca's ~40° walls leaves only a 20–25° sliver. Co-location is decided by slope,
not by arrival geometry."""),
("code", """grand = (3.0, 25.0)
tambo = (20.0, 60.0)
lo, hi = max(grand[0], tambo[0]), min(grand[1], tambo[1])
print(f"GRAND deployable band: {grand[0]:.0f}-{grand[1]:.0f} deg")
print(f"TAMBO wall band:       {tambo[0]:.0f}-{tambo[1]:.0f} deg")
print(f"overlap:               {lo:.0f}-{hi:.0f} deg  ({hi-lo:.0f} deg wide)")"""),
("md", """## How firm is a result?

`oroscope-sensitivity` varies one parameter at a time about a baseline and tabulates how
much each moves the answer. Run against a real TAMBO baseline, the verdict was blunt:

| parameter | | | | |
|---|---|---|---|---|
| `decay_energy_pev` | 3 → **10 878** | 55 → **2056** | 100 → **0** | 1000 → **0** |
| `min_score` | 0.0 → **45 928** | 0.2 → **15 481** | 0.35 → **2056** | 0.5 → **0** |
| `min_target_slope_deg` | 0° → **7442** | 15° → **5309** | 25° → **2056** | 35° → **0** |

Every criterion sits near a cliff. **The decay energy is the worst**: across TAMBO's own
3 PeV – 1 EeV reach the answer runs from 10 878 to zero, because a single energy cannot
stand in for a spectrum."""),
("code", """crossing_m = 3000.0
print(f"P(tau decays within a {crossing_m/1000:.0f} km crossing):\\n")
for e in (3.0, 10.0, 55.0, 100.0, 1000.0):
    L = physics.tau_decay_length_m(e)
    p = 1 - np.exp(-crossing_m / L)
    bar = "#" * int(round(p * 40))
    print(f"{e:>7.0f} PeV  {p:5.3f}  {bar}")"""),
("md", """That is a factor of seventeen inside one experiment's energy reach — and it is invisible
to every other term in the score.

**So: fold over the real spectrum before quoting a capacity.** A number computed at one
representative energy is a property of the energy chosen, not of the terrain.

## What to distrust in your own results

Three things this project measured about itself, worth checking in any search:

1. **Reported area is not physics-accepted area.** Morphological closing more than
   doubles it — measured at 2.29× with a stride-1 control run. Closing is not wrong; a
   site has to be a deployable region rather than a scatter of pixels. But the two
   numbers are different and should not be conflated.
2. **Candidate striding is unbiased** — acceptance is identical at strides 1 and 5, and
   the stride-corrected area matches the stride-1 truth to 0.05%. So that one *is* safe.
3. **Area and capacity are measured on different grids** at `downsample_factor > 1`, so
   a feature a few pixels wide loses area it keeps detectors on.

`docs/assumptions.rst` is the full list, and it is deliberately blunt."""),
("md", """## Where to go next

- The **[assumptions and limitations](https://mbustama.github.io/oroscope/assumptions.html)**
  page — what the numbers rest on.
- The **[physics](https://mbustama.github.io/oroscope/physics.html)** page — the
  derivation behind every criterion."""),
("md", footer(prev=("05_grand_and_tambo.ipynb", "GRAND and TAMBO"))),
]

NOTEBOOKS = {
    "01_getting_started.ipynb": NB01,
    "02_the_arrival_scan.ipynb": NB02,
    "03_physics_toolkit.ipynb": NB03,
    "04_criteria_and_scoring.ipynb": NB04,
    "05_grand_and_tambo.ipynb": NB05,
    "06_combining_and_sensitivity.ipynb": NB06,
}


def main():
    OUT.mkdir(exist_ok=True)
    for name, cells in NOTEBOOKS.items():
        nbf.write(notebook(cells), str(OUT / name))
        print(f"wrote {name}  ({len(cells)} cells)")


if __name__ == "__main__":
    main()
