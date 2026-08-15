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

# Notebook 7 draws nothing -- it is tables and prose -- so it must not import pyplot.
# `ruff check .` lints notebooks, and an unused import there fails CI like any other.
PREAMBLE_NO_PLOT = """import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join('..', 'src')))

import numpy as np
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
| GRAND | 4580.2 km² | 1 | 5317 | 1.1% |
| TAMBO | 83.6 km² | 15 | 9717 | 59.9% |
| **joint** | 50.1 km² | | | Jaccard 0.011 |

The interesting part is *why* the joint is small. Three fifths of TAMBO-viable ground is
also GRAND-viable, but the two deployable **slope bands barely overlap** — GRAND's 3–25°
against Colca's ~40° walls leaves only a 20–25° sliver. Co-location is decided by slope,
not by arrival geometry.

> An earlier version of this table reported TAMBO at 44.5 km² and the joint at 26.4.
> Both were wrong. `load_run` took the alphabetically first `.tif` in a run directory,
> and a directory re-run since the project was renamed holds both
> `oroscope_results_*.tif` and a stale `grand_search_results_*.tif` — the legacy prefix
> sorts first, so the overlay quietly used a superseded mask. Nothing failed; the
> report simply described a run that no longer existed. It is worth knowing that this
> class of fault produces a plausible number rather than an error."""),
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

Every criterion sits near a cliff. **The decay energy was the worst**: across TAMBO's
own 3 PeV – 1 EeV reach the answer ran from 10 878 to zero, because a single energy
cannot stand in for a spectrum.

That row is now history, and it is the reason the code changed. The decay term is
folded over a power-law spectrum instead, with the index pinned or marginalised
(`--decay_spectral_index`), and the same result then varies by **1.46×** across a
plausible range of index rather than without bound. `min_score` is what remains
dominant, at 2.38× to 0.20× about its baseline."""),
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
   the stride-corrected area matches the stride-1 truth to 0.05%. So that one *is*
   safe, *with the caveat below*.
3. **The closing element and the stride interact.** That striding result was measured
   with GRAND's 1 km closing element. Each run's own funnel gives the factor directly,
   and the two Colca configs disagree: GRAND's mask is 2.19× its stride-corrected
   accepted set — an independent check on the 2.29× above — while TAMBO's is **0.53×**,
   because a 100 m element is about three pixels and cannot bridge the gaps stride 5
   leaves. TAMBO's area is therefore a *lower* bound, not an upper one.
4. **Area and capacity are measured on different grids** at `downsample_factor > 1`, so
   a feature a few pixels wide loses area it keeps detectors on.
5. **Not every site in the results file is in the result.** `sites` lists everything
   that cleared the thresholds; with `stop_at_target`, only the first *n* were
   selected. Filter on each record's `selected` flag before totalling anything.

Every run reports 3 for itself, in its own summary. `docs/assumptions.rst` is the full
list, and it is deliberately blunt."""),
("md", """## Where to go next

- The **[assumptions and limitations](https://mbustama.github.io/oroscope/assumptions.html)**
  page — what the numbers rest on.
- The **[physics](https://mbustama.github.io/oroscope/physics.html)** page — the
  derivation behind every criterion."""),
("md", """## Where to go next

- **[7. Running a whole search](07_running_a_search.ipynb)** — driving the pipeline from
  Python, reading what it hands back, and the plan for the full Arequipa DEM."""),
("md", footer(prev=("05_grand_and_tambo.ipynb", "GRAND and TAMBO"),
              nxt=("07_running_a_search.ipynb", "Running a whole search"))),
]

# --------------------------------------------------------------------------- 07
NB07 = [
("md", """# 7. Running a whole search, and reading what comes back

The earlier notebooks drive the pieces: the scan kernel, the physics, the score shapes.
This one drives **the whole pipeline** — screen, scan, score, clean, label, pack, write
— as an ordinary Python call, and then reads the result properly.

It is also the place the **full Arequipa DEM** run belongs. That run is at the end,
guarded so this notebook still executes without the DEM, which is gitignored and a
quarter of a gigabyte.

Three things are worth knowing before the first call:

- **Everything the command line can do, the library can do.** Configuration files,
  the memory pre-flight, the run summary. There is no CLI-only behaviour left.
- **The pipeline returns its results.** It used to return `None` and leave callers to
  find and re-read the JSON it had just written.
- **It explains itself.** A plain-language account of what was found and why, printed
  and saved as `explanation.txt`, on by default."""),
("code", PREAMBLE_NO_PLOT + """
import contextlib
import io
import json
import tempfile

import explain
import site_searcher as ss

WORK = tempfile.mkdtemp(prefix="oroscope_nb07_")
print("working in", WORK)"""),
("md", """## Configuration is data, not a command-line concern

`default_config()` returns every knob the tool understands, with its default.
`generate_config(path, preset)` writes that as a template, and `load_config(path)` reads
one back. All three used to exist only inside `main()`, reachable by running the CLI.

A template naming **every** key matters more than it sounds: a config with holes in it
falls back silently for whatever it omits, and the fallback file is the least visible
input the tool has."""),
("code", """cfg = ss.default_config("arequipa")
print(f"{len(cfg)} keys, e.g.:")
for key in ("dem_path", "min_slope_deg", "max_slope_deg", "candidate_stride",
            "downsample_factor", "min_score", "explain"):
    print(f"   {key:>20}: {cfg[key]!r}")

path = os.path.join(WORK, "arequipa.json")
ss.generate_config(path, "arequipa")
print(f"\\nwritten and read back identically: {ss.load_config(path) == cfg}")"""),
("md", """## Before a big run: what will it cost?

`estimate_peak_memory_gb` predicts the *anonymous* allocations from the DEM's size and
two parameters. The memory-mapped DEM is deliberately excluded — it is file-backed and
the kernel can evict it, and counting it would make every large search look impossible
when the streaming design exists precisely so that it is not.

`downsample_factor` is the knob that matters: the labelling arrays scale as its inverse
square. Here is the real Arequipa DEM, 10204 × 12603 pixels."""),
("code", """rows, cols = 10204, 12603
print(f"Arequipa DEM: {rows} x {cols} = {rows*cols/1e6:.0f} Mpx\\n")
for ds in (1, 2, 4, 8):
    need = ss.estimate_peak_memory_gb(rows, cols, downsample_factor=ds)
    print(f"   downsample_factor {ds}:  {need:5.2f} GiB")

have = ss.available_memory_gb()
print(f"\\navailable right now: {have:.1f} GiB" if have else "\\n(memory not reportable here)")
print("\\nThis is why the full run uses downsample_factor 4.")"""),
("md", """The estimate is rough and says so — `survival_fraction` is the share of pixels passing
the topographic screen, which is terrain-dependent and unknown until the screen has run.
It is meant to catch the order-of-magnitude mistake, not to predict a number.

`preflight_memory` does the whole job: estimate, warn if it is close, and cap the
process's address space so a search that outgrows the machine fails with `MemoryError`
naming itself rather than letting the kernel's OOM killer pick a victim — which may be
your editor. A ten-point sweep once did exactly that at 6.9 GB."""),
("md", """## A complete run

Synthetic terrain, so this executes anywhere. A ridge with a slope in front of it: the
slope sees the ridge, and the ridge's own flank sees the terrain rising beyond."""),
("code", """def ridge_and_slope(n, cell_x):
    \"\"\"A valley between a ridge and a rising slope. Closed-form, no DEM needed.\"\"\"
    cols = np.arange(n, dtype=np.float64)[None, :].repeat(n, 0)
    x = cols * cell_x
    ridge = 1400.0 * np.exp(-((x - 0.30 * n * cell_x) / (0.05 * n * cell_x)) ** 2)
    rise = np.clip((x - 0.55 * n * cell_x) / (0.45 * n * cell_x), 0, 1) ** 2 * 1500.0
    return (2200.0 + ridge + rise).astype(np.float32)"""),
("code", """import tifffile as tiff

grid = ss.resolve_grid_geometry("no-such-file.tif", -15.6, cell_size_deg=1/3600)
z = ridge_and_slope(700, grid.cell_size_x)

dem = os.path.join(WORK, "ridge.tif")
tiff.imwrite(dem, z, extratags=[
    (33550, "d", 3, (1/3600, 1/3600, 0.0)),                    # ModelPixelScale
    (33922, "d", 6, (0.0, 0.0, 0.0, -72.3, -15.6, 0.0)),       # ModelTiepoint
])
print("wrote a GeoTIFF carrying its own resolution and corner:", os.path.basename(dem))"""),
("md", """Now the search itself. Note what is *not* passed: no origin — the DEM carries its own
corner in the tiepoint tag, and reading it removes the most error-prone input the tool
has. A supplied origin that disagrees with the file by more than ~100 m is reported
rather than silently honoured, because a wrong origin does not fail, it
mis-georeferences every output.

The run prints a great deal. It is captured here and unpacked below."""),
("code", """log = io.StringIO()
with contextlib.redirect_stdout(log), contextlib.redirect_stderr(io.StringIO()):
    results = ss.find_grand_regions_interactive(
        dem_path=dem,
        run_output_dir=os.path.join(WORK, "run"),
        # Note search_mode and grid_type. The function's own defaults are 'single'
        # and 'square'; the config template's are 'distributed' and 'hex'. Omitting
        # them here is not the same as omitting them from a config file -- see below.
        search_mode="distributed", grid_type="hex",
        target_antennas=200, min_sub_array_size=20,
        min_width_km=1.0, antenna_spacing_km=1.0,
        min_dist_km=3.0, max_dist_km=20.0,
        downsample_factor=2, tile_size=256, candidate_stride=5, num_cores=2,
    )

print(f"{len(log.getvalue().splitlines())} lines of output captured\\n")
print("returned:", ", ".join(sorted(results)))"""),
("md", """### One trap worth knowing

**The function's defaults are not the config template's defaults.** Five parameters
differ, and omitting one means different things depending on which door you came in by:

| parameter | `find_grand_regions_interactive` | `default_config()` |
|---|---|---|
| `search_mode` | `single` | `distributed` |
| `grid_type` | `square` | `hex` |
| `target_antennas` | 1000 | 10000 |
| `min_dist_km` | 30.0 | 10.0 |
| `min_sub_array_size` | 100 | 500 |

An earlier draft of this notebook omitted `search_mode` and quietly ran a *single*
search with a 30 km minimum distance, which on this small ridge found nothing at all.
The funnel said so plainly, which is the system working — but the safe habit when
driving the library is to start from `ss.default_config()` and override, rather than to
rely on the signature's defaults."""),
("code", """template = ss.default_config()
signature_defaults = {
    "search_mode": "single", "grid_type": "square",
    "target_antennas": 1000, "min_dist_km": 30.0, "min_sub_array_size": 100,
}
print(f"{'parameter':22} {'function':>12} {'config template':>18}")
for key, value in signature_defaults.items():
    print(f"{key:22} {str(value):>12} {str(template[key]):>18}")"""),
("md", """That dictionary is the same content the results JSON holds, plus the explanation and
the paths written. No re-reading the file it just wrote."""),
("code", """print(f"sites:    {results['results']['total_sites']}")
print(f"capacity: {results['results']['total_capacity']}")
print(f"stages:   {', '.join(results['timings_sec'])}")
print(f"files:    {len(results['output_files'])} written")
for f in results["output_files"]:
    print("   ", os.path.basename(f))"""),
("md", """## The funnel is the diagnostic

Every filter records how many pixels survived it. When a search returns little or
nothing, **the stage where the count collapses is the constraint responsible** — and
that is the single most useful thing anyone can be told about a disappointing run."""),
("code", """for stage, count in results["funnel"].items():
    print(f"   {stage:<34} {count:>12,}")

binding = explain.binding_constraint(results["funnel"])
print(f"\\nbinding constraint: {binding['stage']!r}")
print(f"   kept {100*binding['kept_fraction']:.1f}% of the {binding['before']:,} that reached it")
print(f"   change: {binding['knob']}")"""),
("md", """Two stages are excluded from that search by construction, and it is worth knowing why:

- **`kept by stride N`** is a deliberate subsample, not a filter. It removes four
  candidates in five and the acceptance is unchanged, so calling it the constraint
  would name the same answer on nearly every run.
- **`after gap closing`** *adds* pixels. A stage that grows the set cannot be what
  shrank it.

## Which sites are actually in the result

`sites` lists everything that cleared the area and capacity thresholds. With
`stop_at_target`, selection walks that capacity-sorted list until the target is met and
stops — so the list can be longer than the result. Only the selection is in
`total_sites`, `total_capacity` and the exported raster.

Each record says which it is."""),
("code", """log2 = io.StringIO()
with contextlib.redirect_stdout(log2), contextlib.redirect_stderr(io.StringIO()):
    truncated = ss.find_grand_regions_interactive(
        dem_path=dem, run_output_dir=os.path.join(WORK, "run2"),
        target_antennas=50, min_sub_array_size=5, stop_at_target=True,
        min_width_km=1.0, antenna_spacing_km=1.0,
        min_dist_km=3.0, max_dist_km=20.0,
        downsample_factor=2, tile_size=256, candidate_stride=5, num_cores=2,
    )

chosen, shortlisted = explain.selected_sites(truncated)
print(f"listed in the file: {len(chosen) + len(shortlisted)}")
print(f"selected:           {truncated['results']['total_sites']}\\n")
for site in chosen + shortlisted:
    mark = "selected" if site["selected"] else "not selected"
    print(f"   site {site['site_id']:>3}  {site['area_km2']:>8.2f} km²  "
          f"{site['capacity_exact']:>5} detectors   {mark}")

print(f"\\nsumming everything listed:  {sum(s['area_km2'] for s in chosen + shortlisted):8.2f} km²")
print(f"summing the selection:      {sum(s['area_km2'] for s in chosen):8.2f} km²  <- the raster")"""),
("md", """Totalling the wrong one over-reports, which is exactly the mistake this flag exists to
prevent. The sites that were not selected are the *next best ground*, not ground that
failed — worth keeping in the file, worth excluding from the totals.

## Attribution: what held each site back

The score is a product of **named** components, each in [0, 1], and each site's record
carries the distribution of every one. Under a product the lowest component bounds the
total from above, so naming it turns "this site scored 0.34" into something actionable."""),
("code", """site = chosen[0]
scan = site["arrival_scan"]
parts = {k[len("score_"):-len("_p50")]: v for k, v in scan.items()
         if k.startswith("score_") and k.endswith("_p50") and k != "score_p50"}

print(f"site {site['site_id']}, median score {scan['score_p50']:.3f}\\n")
for name, value in sorted(parts.items(), key=lambda kv: kv[1]):
    bar = "#" * int(round(value * 40))
    print(f"   {name:>14}  {value:5.3f}  {bar}")

name, value = explain.weakest_component(scan)
print(f"\\nweakest: {name} at {value:.3f}")"""),
("md", """On the real Colca configurations this is unambiguous: `solid_angle` is the weakest
component at **15 of 15** TAMBO sites, with everything else at 1.0 except the decay term
at 0.96. So that result is set almost entirely by `solid_angle_half_sr`, whose 0.05 sr
default is a GRAND-scale value.

That is the kind of statement the components make available and a single total does not.

## How much did closing move the area?

The reported area is not the physics-accepted area: the mask is closed morphologically
before areas are measured. The published figure is 2.29× at Colca, measured against a
stride-1 control — but each run has the number in it, as closed pixels over
stride-corrected accepted pixels."""),
("code", """ratio = explain.closing_inflation(results["funnel"],
                                 results["parameters"]["candidate_stride"])
print(f"this run: closing moved the mask by {ratio:.2f}x")
print("\\nOn the real configurations:")
print("   GRAND Colca  2.19x   (against 2.29x from a stride-1 control -- an independent check)")
print("   TAMBO Colca  0.53x   (a 100 m element cannot bridge the gaps stride 5 leaves,")
print("                        so its area is a LOWER bound, not an upper one)")"""),
("md", """## The run, explained — the whole summary, here

Everything above is assembled for you. `explain.explain_results` takes the results
dictionary and returns a string — it opens no files, runs nothing and needs no DEM, so
a run from months ago can still be explained from its JSON.

It is on by default, printed at the end of every run and saved as `explanation.txt`
beside the results, because these runs are meant to be handed to other people and a
terminal scrollback is not. `--no_explain` suppresses it.

This is the text in full. Its sections, in order:

| section | answers |
|---|---|
| **The run** | what was searched, at what resolution, by which commit |
| **The headline** | how many sites, how much area, how many detectors |
| **Where the candidates went** | the funnel, and **which constraint bound this run** |
| **From pixels to sites** | labelled regions → area threshold → capacity threshold |
| **The sites** | each one's area, capacity, facing, score and weakest criterion |
| **Why these sites qualify** | what the ground actually offers, criterion by criterion, with coordinates |
| **What energy this geometry favours** | where the geometric aperture peaks |
| **How to read these numbers** | the closing factor *for this run*, and what area is not |
| **Which of these are assumptions** | choices rather than measurements, with measured sensitivities |
| **What to try next** | concrete commands, chosen from what this run did |"""),
("code", """print(results["explanation"])"""),
("md", """## Provenance

Separate from the science outputs, and the answer to "what produced this number?"."""),
("code", """prov = results["provenance"]
print(f"commit:   {prov['git']['commit'][:10]} on {prov['git']['branch']}"
      f"  ({'dirty' if prov['git']['dirty'] else 'clean'} tree)")
print(f"DEM:      {os.path.basename(prov['dem']['path'])}")
print(f"          sha256 {prov['dem']['sha256'][:24]}...")
print(f"          {prov['dem']['cell_size_y_m']:.2f} m N-S x {prov['dem']['cell_size_x_m']:.2f} m E-W")
print(f"python:   {prov['platform']['python']} on {prov['platform']['system']}")
print(f"packages: {', '.join(f'{k} {v}' for k, v in list(prov['packages'].items())[:4])}, ...")"""),
("md", """---

## The full Arequipa DEM

Every number this project has published comes from **crops** — Colca, and small Arequipa
windows. The full DEM is the run that has never been done, and this is where it lands.

**Read, not run.** The cells below open results that were produced locally and stored in
`results/arequipa_full/`. They do not start a search. Each of these searches takes about
half an hour, CI executes notebooks on every push, and a tutorial costing ninety minutes
of compute per commit is a bill rather than a tutorial. The expensive half runs once, on
a machine that has the DEM; the notebook opens a few hundred kilobytes of JSON.

To produce or refresh the store:

```bash
python tools/run_arequipa_full.py --dry-run   # what it will do, and what it will cost
python tools/run_arequipa_full.py             # GRAND, TAMBO, then the combination
```

**Regenerate it when a configuration changes, and not otherwise.** The store carries a
manifest naming the configs and the time, so a stale one is detectable rather than
merely suspected.

Three searches, all at the same `downsample_factor` so their masks are pixel-aligned:

| | config | what it asks |
|---|---|---|
| **GRAND alone** | `config/grand_arequipa_full.json` | 3–25° deployable ground seeing a target 10–40 km away, within ±3° of the horizon |
| **TAMBO alone** | `config/tambo_arequipa_full.json` | a 20–60° near wall facing a ≥25° far wall, 2–5 km across |
| **Combined** | `combine_experiments` over both | joint, union, and how much of each sits inside the other |

**What it costs.** 10204 × 12603 pixels, about 129 Mpx. At `downsample_factor: 4` the
estimator says 2.3 GiB against the ~6 GiB typically free; at 1 it says 4.5 GiB, which is
why 4 is the setting. That choice has a price worth stating: area is measured on the
downsampled mask while capacity is measured at full resolution, so a feature a few
pixels wide keeps its detectors and loses area. **Read these areas as lower bounds**,
and more so for TAMBO's canyon strips than for GRAND's blobs."""),
("code", """STORE = os.path.abspath(os.path.join("..", "results", "arequipa_full"))

def load_stored(label):
    \"\"\"Reads one stored run, or returns None when the store does not have it yet.\"\"\"
    path = os.path.join(STORE, f"{label}_results.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)

manifest_path = os.path.join(STORE, "manifest.json")
if os.path.exists(manifest_path):
    with open(manifest_path) as f:
        manifest = json.load(f)
    print(f"store generated {manifest['generated']} by {manifest['generated_by']}")
    print(f"from {manifest['dem']}\\n")
    for name in manifest["files"]:
        print("   ", name)
else:
    manifest = None
    print("The full-DEM store is empty: these searches have not been run yet.\\n")
    print("Produce it with:")
    print("    python tools/run_arequipa_full.py --dry-run")
    print("    python tools/run_arequipa_full.py")

grand_full = load_stored("grand")
tambo_full = load_stored("tambo")"""),
("md", """### GRAND over the whole DEM

The four things worth reading, in this order:

1. **The funnel**, and specifically whether the binding constraint is the same one the
   crops found. If a full DEM is bound by a different stage than its crops were, the
   crops were not representative and every number derived from them needs re-reading.
2. **The area**, against the crop scaled up — and against the closing factor this run
   reports for itself, rather than the 2.29× quoted from Colca.
3. **The site count and their spread.** A crop cannot say whether the good ground is one
   region or fifty scattered ones, and that is a deployment question, not a physics one.
4. **The weakest score component**, which on the crops is `solid_angle` everywhere. If
   that holds at full scale it is a statement about the criterion, not about Peru."""),
("code", """def summarise(results, label):
    if results is None:
        print(f"{label}: not in the store yet.")
        return
    chosen, shortlisted = explain.selected_sites(results)
    area = sum(s["area_km2"] for s in chosen)
    binding = explain.binding_constraint(results["funnel"])
    ratio = explain.closing_inflation(results["funnel"],
                                      results["parameters"]["candidate_stride"])
    weakest = [explain.weakest_component(s.get("arrival_scan") or {}) for s in chosen]
    named = [w[0] for w in weakest if w]

    print(f"{label}")
    print(f"   sites          {results['results']['total_sites']:>10,}"
          f"   ({len(shortlisted)} more cleared the thresholds, not selected)")
    print(f"   capacity       {results['results']['total_capacity']:>10,}")
    print(f"   area           {area:>10,.1f} km²")
    if binding:
        print(f"   bound by       {binding['stage']}"
              f"  (kept {100*binding['kept_fraction']:.1f}%)")
    if ratio is not None:
        print(f"   closing moved  {ratio:>10.2f}x")
    if named:
        commonest = max(set(named), key=named.count)
        print(f"   weakest        {commonest} at {named.count(commonest)}/{len(named)} sites")

summarise(grand_full, "GRAND, full Arequipa DEM")"""),
("code", """if grand_full is not None:
    print(grand_full.get("explanation") or explain.explain_results(grand_full))
else:
    print("Nothing stored for GRAND yet -- see the cell above for how to produce it.")"""),
("md", """### TAMBO over the whole DEM

The more interesting of the two, because TAMBO's criteria are canyon-shaped and the
crop was *chosen* for containing a canyon. Over the whole DEM the question becomes: how
much other canyon is there, and is any of it as good?

Note that the areas here are the ones most affected by `downsample_factor: 4`: a strip
along a wall is exactly the feature that loses area to downsampling while keeping its
detectors."""),
("code", """summarise(tambo_full, "TAMBO, full Arequipa DEM")"""),
("code", """if tambo_full is not None:
    print(tambo_full.get("explanation") or explain.explain_results(tambo_full))
else:
    print("Nothing stored for TAMBO yet -- see above for how to produce it.")"""),
("md", """### Where both are viable

The overlay. On the Colca crop the answer was decided by slope: GRAND's 3–25° deployable
band against Colca's ~40° walls leaves only a 20–25° sliver, so the joint was about a
percent of GRAND's area and three fifths of TAMBO's.

Whether that survives at full scale is a real question. The crop contains one canyon
system; the DEM contains many, of varying wall slope, and the joint area is the
programme-level number — one site, one road, one power feed, two experiments."""),
("code", """report_path = os.path.join(STORE, "combined_report.json")
if not os.path.exists(report_path):
    print("The combination has not been produced yet.")
else:
    with open(report_path) as f:
        report = json.load(f)

    width = max(len(r["label"]) for r in report["runs"])
    print(f"   {'experiment'.ljust(width)} {'area km²':>12} {'sites':>7} "
          f"{'capacity':>10} {'in joint':>9}")
    print("   " + "-" * (width + 42))
    for r in report["runs"]:
        print(f"   {r['label'].ljust(width)} {r['area_km2']:>12,.1f} "
              f"{r['reported_sites']:>7,} {r['reported_capacity']:>10,} "
              f"{100*r['fraction_of_own_area_in_joint']:>8.1f}%")
    print("   " + "-" * (width + 42))
    print(f"   joint  {report['joint']['area_km2']:>10,.1f} km²")
    print(f"   union  {report['union']['area_km2']:>10,.1f} km²")
    for pair, stats in report["pairwise_overlap"].items():
        print(f"   {pair}: Jaccard {stats['jaccard']:.4f}")"""),
("md", """And the overlay explains itself too, the same way a search does — including the part
that is easy to get wrong. Co-location is decided by whichever *ground* property the
two experiments share least of, because a pixel has one slope and both have to accept
it. What each asks of the **view** — the distance window, the arrival elevations — may
differ freely: two experiments can look out from the same hillside at different ranges
without conflict.

`oroscope-combine` prints this and saves it as `combination_explanation.txt`."""),
("code", """if not os.path.exists(report_path):
    print("The combination has not been produced yet.")
else:
    runs = {label: res for label, res in (("GRAND", grand_full), ("TAMBO", tambo_full))
            if res is not None}
    print(explain.explain_combination(report, runs))"""),
("md", """### Comparing against the crop

The crop's numbers, for reference — GRAND 4580.2 km² in 1 site with 5317 detectors,
TAMBO 83.6 km² in 15 sites with 9717, joint 50.1 km². If the full DEM's binding
constraint or weakest component differs from these, the crop was not representative and
the comparison is the finding.

A crop is chosen because it is interesting. **A search over ground chosen for being
interesting is not a survey**, and that is the gap this run exists to close."""),
("code", """if grand_full is not None and tambo_full is not None:
    print("Both full-DEM runs are stored; compare their funnels against the crops':\\n")
    for label, res in (("GRAND", grand_full), ("TAMBO", tambo_full)):
        b = explain.binding_constraint(res["funnel"])
        print(f"   {label:>6}: bound by {b['stage']!r} "
              f"(kept {100*b['kept_fraction']:.1f}%)")
    print("\\n   Colca crop, for comparison:")
    print("   GRAND: bound by 'directions accepted' (kept 60.1%)")
    print("   TAMBO: bound by 'directions accepted' (kept 17.5%)")
else:
    print("Run both searches to make this comparison.")"""),
("md", footer(prev=("06_combining_and_sensitivity.ipynb", "Combining and sensitivity"))),
]

NOTEBOOKS = {
    "01_getting_started.ipynb": NB01,
    "02_the_arrival_scan.ipynb": NB02,
    "03_physics_toolkit.ipynb": NB03,
    "04_criteria_and_scoring.ipynb": NB04,
    "05_grand_and_tambo.ipynb": NB05,
    "06_combining_and_sensitivity.ipynb": NB06,
    "07_running_a_search.ipynb": NB07,
}


def cell_sources(nb):
    """A notebook's content, ignoring outputs and execution counts."""
    return [(c.cell_type, c.source) for c in nb.cells]


def main():
    """
    Writes only the notebooks whose *content* changed.

    Rewriting all of them unconditionally stripped the stored outputs from every
    notebook whenever one word of one notebook changed, which forced a full
    re-execution to restore them -- and re-executing is the expensive part. Comparing
    cell sources makes this idempotent: edit one notebook, re-execute one notebook.
    """
    OUT.mkdir(exist_ok=True)
    written = skipped = 0
    for name, cells in NOTEBOOKS.items():
        fresh = notebook(cells)
        path = OUT / name
        if path.exists():
            try:
                if cell_sources(nbf.read(str(path), as_version=4)) == cell_sources(fresh):
                    print(f"unchanged {name}  (outputs preserved)")
                    skipped += 1
                    continue
            except Exception:
                pass                      # unreadable or older format: just rewrite it
        nbf.write(fresh, str(path))
        print(f"wrote {name}  ({len(cells)} cells)  -- re-execute this one")
        written += 1

    if written:
        print(f"\n{written} rewritten, {skipped} unchanged. Re-execute what was "
              f"rewritten, so its stored outputs match its code:\n"
              f"    cd notebooks && jupyter nbconvert --execute --inplace <name>.ipynb")
    else:
        print(f"\nnothing to do: all {skipped} notebooks are up to date.")


if __name__ == "__main__":
    main()
