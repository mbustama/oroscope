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

# One import, the way a reader will actually use the library. Repeated verbatim at the
# top of every notebook so each one stands alone, which is how people open them.
#
# This used to be a sys.path insert pointing at src/, because the modules were flat and
# there was no package to import. Worse than ugly: an unconditional insert shadows an
# installed copy with whatever happens to be in the source tree. Both went with the
# package (roadmap 6.33); `pip install -e .` from a clone is now the only setup step.
PREAMBLE = """import numpy as np
import matplotlib.pyplot as plt
"""

# Notebooks 7 and 8 draw nothing -- tables and prose -- so they must not import pyplot.
# `ruff check .` lints notebooks, and an unused import there fails CI like any other.
PREAMBLE_NO_PLOT = """import numpy as np
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
import oroscope

print(f"oroscope {oroscope.__version__}")"""),
("md", """> **One import.** `import oroscope` is the whole setup — every function these notebooks
> use is on it, and the submodules stay available when a narrower namespace reads better
> (`from oroscope import physics`). Install it first, with `pip install oroscope`, or
> `pip install -e .` from a clone.

## A piece of terrain

A DEM is just a 2-D array of elevations plus a statement of how big a pixel is on the
ground. Here is a plain west-facing slope with a ridge to its east — the geometry a
detector on the slope would use to look *at* the ridge."""),
("code", """n = 400                      # pixels on a side
cell_deg = 1 / 3600.0        # 1 arc-second, as SRTM and AW3D30 supply

# Metric pixel sizes differ on each axis away from the equator, because a degree of
# longitude shrinks with the cosine of the latitude. resolve_grid_geometry does that.
grid = oroscope.resolve_grid_geometry("no-such-file.tif", -15.6, cell_size_deg=cell_deg)
print(f"pixel: {grid.cell_size_y:.1f} m north-south, {grid.cell_size_x:.1f} m east-west")
print(f"resolution source: {grid.source}")

cols = np.arange(n)[None, :].repeat(n, 0)
x_m = cols * grid.cell_size_x

# A ridge 1.2 km high centred 6.5 km east, on a plain at 2 km
z = (2000.0 + 1200.0 * np.exp(-((x_m - 6500.0) / 700.0) ** 2)).astype(np.float32)

fig, ax = plt.subplots(figsize=(7, 2.6))
ax.plot(x_m / 1000, z[0], color="#7A6A4F")
ax.fill_between(x_m[0] / 1000, 1900, z[0], color="#D9CDB8")
ax.set_xlabel("East (km)")
ax.set_ylabel("Elevation (m)")
ax.set_title("A ridge to look at")
ax.set_ylim(1900, 3400)
ax.spines[["top", "right"]].set_visible(False)
plt.show()"""),
("md", """## Scanning one candidate

A *candidate* is a pixel the search is considering, given as `[row, col, aspect_deg]`.
The aspect is the downhill direction, and the scan fans its azimuths around it.

Put a detector on the plain west of the ridge, facing east."""),
("code", """candidate = np.array([[200.0, 60.0, 90.0]])     # row, col, aspect: due east

out = oroscope.scan(
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
    r = oroscope.scan(cand, z, grid, n_azimuths=1, half_width_deg=0.0,
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
from oroscope import arrival_scan
from oroscope import figures
from oroscope import site_searcher as ss

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
from oroscope import physics
from oroscope import figures"""),
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
ax.set_xlabel("Atmospheric depth traversed (g/cm$^2$)")
ax.set_ylabel("Particle content / peak")
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
from oroscope import scoring"""),
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
ax.set_xlabel("Accepted solid angle (sr)")
ax.set_ylabel("Score")
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
    ax.set_xlabel("Score")
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
from oroscope import arrival_scan
from oroscope import figures
from oroscope import site_searcher as ss

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
("code", """import numpy as np

from oroscope import combine_experiments as combine, physics"""),
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

- **[7. Animating the mechanism](07_animating_the_mechanism.ipynb)** — the parts of
  all this that a still picture explains badly, as eight short films.
- **[8. Explaining a run](08_explaining_a_run.ipynb)** — driving the pipeline from
  Python, and reading what it says about a run that succeeds and one that does not."""),
("md", footer(prev=("05_grand_and_tambo.ipynb", "GRAND and TAMBO"),
              nxt=("07_animating_the_mechanism.ipynb", "Animating the mechanism"))),
]

# --------------------------------------------------------------------------- 08  explaining a run
SHOW_HELPER = """from IPython.display import Image, display


def show_figure(path, width=1100, caption=None):
    \"\"\"Displays a figure the pipeline wrote, downscaled so the notebook stays small.

    The searcher saves its map to disk and says where; a notebook that only prints the
    path makes the reader go and find it. Downscaled because the full-resolution PNG is
    a megabyte or two and these outputs are committed.
    \"\"\"
    if not os.path.exists(path):
        print(f"not here yet: {os.path.relpath(path)}")
        return
    try:
        from PIL import Image as PILImage
        import io as _io
        img = PILImage.open(path)
        if img.width > width:
            img = img.resize((width, round(img.height * width / img.width)),
                             PILImage.LANCZOS)
        buf = _io.BytesIO()
        img.convert("RGB").save(buf, format="PNG", optimize=True)
        data = buf.getvalue()
        # A shaded terrain map is photographic and does not compress as PNG: the full
        # Arequipa map is 1.1 MB at this width and 190 KB as JPEG, for no visible
        # difference. Line art stays PNG, where JPEG would ring around every edge.
        if len(data) > 250 * 1024:
            jpg = _io.BytesIO()
            img.convert("RGB").save(jpg, format="JPEG", quality=82, optimize=True)
            data = jpg.getvalue()
    except Exception:                       # pillow absent: show it full size
        with open(path, "rb") as f:
            data = f.read()
    if caption:
        print(caption)
    display(Image(data=data))"""


NB_EXPLAIN = [
("md", """# 8. Explaining a run

The earlier notebooks drive the pieces: the scan kernel, the physics, the score shapes.
This one drives **the whole pipeline** — screen, scan, score, clean, label, pack, write
— as an ordinary Python call, and then reads the result properly.

Two searches are run below — one that finds ground and one that finds none — because
the summary of an empty result is the one that matters most and is the easiest to
neglect. The full Arequipa DEM has its own notebook,
**[9](09_arequipa_dem.ipynb)**.

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
import os
import tempfile

import oroscope
from oroscope import explain, site_searcher as ss

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

**Two knobs, and the obvious one is the weaker.** `downsample_factor` scales the
labelling arrays as its inverse square, but candidates are taken on the *native* grid —
the stride subsamples the surviving-pixel list, not the map — so downsampling never
touches them. At full-DEM scale the per-candidate arrays dominate, and `candidate_stride`
is the lever on those. Here is the real Arequipa DEM, 10204 × 12603 pixels."""),
("code", """rows, cols = 10204, 12603
print(f"Arequipa DEM: {rows} x {cols} = {rows*cols/1e6:.0f} Mpx\\n")
for ds in (1, 2, 4, 8):
    need = ss.estimate_peak_memory_gb(rows, cols, downsample_factor=ds)
    print(f"   downsample_factor {ds}:  {need:5.2f} GiB")

print()
for stride in (1, 5, 10, 20):
    need = ss.estimate_peak_memory_gb(rows, cols, downsample_factor=4,
                                      candidate_stride=stride)
    print(f"   downsample_factor 4, candidate_stride {stride:2d}:  {need:6.2f} GiB")

have = ss.available_memory_gb()
print(f"\\navailable right now: {have:.1f} GiB" if have else "\\n(memory not reportable here)")
print("\\nThe full run uses downsample_factor 4 and candidate_stride 5,")
print("and measured 5.68 GiB peak RSS.")"""),
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
("md", """### The map it drew

The run writes a GeoTIFF, a world file, a KML and a PNG map, and prints where it put
them. Printing the path is not the same as showing the picture, so here it is — the
same file, read back from disk."""),
("code", SHOW_HELPER),
("code", """import glob

png = sorted(glob.glob(os.path.join(WORK, "run", "*.png")))
show_figure(png[0] if png else os.path.join(WORK, "run", "missing.png"),
            caption=f"{os.path.basename(png[0])} — altitude, the site outline, "
                    f"a scale bar and north" if png else None)"""),
("md", """### Where a default comes from

A parameter can state its default in three places — the function's signature,
`oroscope --help`, and `default_config()` — and **all three agree**. They did not
always: ten parameters disagreed, so omitting one meant different things depending on
which door you came in by. An earlier draft of this notebook omitted `search_mode`,
quietly ran a *single* search with a 30 km minimum distance, and found nothing at all
on this small ridge. The funnel said so plainly, which is the system working, but the
trap should not have existed.

It cannot come back: a test compares the three sources pairwise, for every parameter.

Starting from `default_config()` and overriding is still the clearer habit, because it
puts every knob in front of you rather than leaving them implicit."""),
("code", """import inspect

template = oroscope.default_config()
signature = {k: v.default for k, v in
             inspect.signature(oroscope.find_grand_regions_interactive).parameters.items()
             if v.default is not inspect.Parameter.empty}

print(f"{'parameter':22} {'signature':>12} {'template':>12}")
for key in ("search_mode", "grid_type", "target_antennas", "min_dist_km",
            "min_sub_array_size", "max_road_dist_km"):
    mark = "ok" if signature[key] == template[key] else "DIFFERS"
    print(f"{key:22} {str(signature[key]):>12} {str(template[key]):>12}   {mark}")"""),
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

## The other outcome: a search that finds nothing

A summary of a successful run is the easy case. The one that matters is the run that
comes back empty, because that is when a reader has no idea what to change — and it is
the case a bare results file serves worst: every section is zero and nothing says why.

Here is the same terrain asked an impossible question. The distance window is moved out
past anything this ridge can offer, so no arrival direction can be accepted."""),
("code", """empty = ss.find_grand_regions_interactive(
    dem_path=dem,
    run_output_dir=os.path.join(WORK, "empty"),
    search_mode="distributed", grid_type="hex",
    target_antennas=200, min_sub_array_size=20,
    min_width_km=1.0, antenna_spacing_km=1.0,
    min_dist_km=60.0, max_dist_km=90.0,      # further than this terrain reaches
    downsample_factor=2, tile_size=256, candidate_stride=5, num_cores=2,
    explain=False,                            # composed below instead, to keep this tidy
)

print(f"sites: {empty['results']['total_sites']}")
print(f"capacity: {empty['results']['total_capacity']}\\n")
for stage, count in empty["funnel"].items():
    print(f"   {stage:<34} {count:>12,}")"""),
("md", """The funnel is the whole answer, and the summary reads it: the stage where the count
reaches zero *is* the constraint, and everything downstream of it is zero for a reason
that is not its own. A stage that empties the map wins outright over any ratio below
it, which is why the report names that one rather than the largest percentage drop."""),
("code", """binding = explain.binding_constraint(empty["funnel"])
print(f"stage:      {binding['stage']}")
print(f"reached it: {binding['before']:,} pixels")
print(f"survived:   {binding['survivors']:,}")
print(f"fatal:      {binding['fatal']}")
print(f"change:     {binding['knob']}")"""),
("md", """And the summary in full. Note what it does *not* do: it does not apologise, and it does
not pad. It states the outcome, names the stage, names the parameter, and then tells you
what to try — which for an empty result is the only useful thing a report can say."""),
("code", """print(explain.explain_results(empty))"""),
("md", """Compare the two summaries. The successful one describes ground; this one describes a
constraint. Both are the same function reading the same shape of dictionary — which is
the point of keeping it a pure function of the results rather than something the
pipeline prints as it goes."""),
("md", """## Where to go next

- **[9. Arequipa, the full DEM](09_arequipa_dem.ipynb)** — the run that has never been
  done, and what to look at when it is.
- **[6. Combining and sensitivity](06_combining_and_sensitivity.ipynb)** — how firm any
  of this is."""),
("md", footer(prev=("07_animating_the_mechanism.ipynb", "Animating the mechanism"),
              nxt=("09_arequipa_dem.ipynb", "Arequipa, the full DEM"))),
]

# --------------------------------------------------------------------------- 09  arequipa
NB_AREQUIPA = [
("md", """# 9. Arequipa, the full DEM

Every number this project has published comes from **crops** — Colca, and small Arequipa
windows. The full DEM is the run that has never been done, and this notebook is where it
lands.

[Notebook 8](08_explaining_a_run.ipynb) covers how to drive the pipeline and how to read
what it says. This one is about one specific run, at a scale the crops cannot speak
for."""),
("code", """import json
import os

from oroscope import explain"""),
("md", """**Read, not run.** The cells below open results that were produced locally and stored in
`results/arequipa_full/`. They do not start a search. Each of these searches takes about
half an hour, CI executes notebooks on every push, and a tutorial costing ninety minutes
of compute per commit is a bill rather than a tutorial. The expensive half runs once, on
a machine that has the DEM; the notebook opens a few hundred kilobytes of JSON.

To produce or refresh the store:

```bash
python tools/run_full_dem.py --dry-run   # report the cost, then stop
python tools/run_full_dem.py             # GRAND, TAMBO, then the combination
```

**Start with `--dry-run`.** It begins nothing — no search, no file, no change to the
store — and prints the five things worth knowing before committing an hour of a
machine:

```text
DEM:       input/dem/arequipa_SRTMGL1.tif
estimate:  5.08 GiB at downsample_factor 4
available: 6.4 GiB
would run: grand, tambo, then combine
expected:  ~25 min for grand, ~1 min for tambo
store:     results/arequipa_full
```

`DEM` says whether the file is even present, so a missing DEM is reported before the
first search starts rather than after. `estimate` against `available` is what decides
`downsample_factor` and `candidate_stride` — the same DEM needs 7.2 GiB at 1 and
5.1 GiB at 4, and downsampling helps less than it looks because the candidates are taken
on the native grid regardless. `would run` honours `--only`, so `--only grand` runs one
search and skips the combination. And no memory cap is applied during a dry run, because
nothing is allocated.

**This is a run that needs its cap set.** The default ceiling is 80% of what the system
reports available, which on a machine whose desktop already holds half of RAM is below
what the search needs: the first attempt died 23 minutes in, at the scoring stage,
against a 5.5 GiB cap. Pass `--max-memory-gb` explicitly — the run measured 5.68 GiB
peak RSS.

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
estimator says 5.1 GiB against the ~6–7 GiB typically free; at 1 it says 7.2 GiB, which
is why 4 is the setting. That choice has a price worth stating: area is measured on the
downsampled mask while capacity is measured at full resolution, so a feature a few
pixels wide keeps its detectors and loses area — the run puts it at around 30% for a
canyon strip. **Read these areas as lower bounds**, and more so for TAMBO's canyon
strips than for GRAND's blobs. It is the reason TAMBO's full-DEM area below cannot be
compared directly against the crop's, which was measured at `downsample_factor: 1`."""),
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
    print("    python tools/run_full_dem.py --dry-run")
    print("    python tools/run_full_dem.py")

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
("md", """### The maps

The store holds the numbers, not the rasters — a GeoTIFF of a 129 Mpx mask is far too
large for a repository. But the PNG maps are in `output/`, and if you have run the
searches locally they are worth looking at rather than reading about: the funnel says
*how much* ground survived, the map says *where* it is.

If `output/` is empty here, produce it with `python tools/run_full_dem.py`. The
images below are stored in this notebook, so they are visible either way."""),
("code", SHOW_HELPER),
("code", """import glob

OUT = os.path.abspath(os.path.join("..", "output"))
for label, folder in (("GRAND", "arequipa_full_grand"),
                      ("TAMBO", "arequipa_full_tambo")):
    found = sorted(glob.glob(os.path.join(OUT, folder, "oroscope_results*.png")))
    show_figure(found[0] if found else os.path.join(OUT, folder, "none.png"),
                caption=f"{label} over the full Arequipa DEM")"""),
("md", """### The overlay, one experiment at a time

The combination is easier to read built up than all at once, so `--reveal` writes it as
three frames. Everything that is not a category — the terrain, the colour bar, the
roads, the towns, the scale bar, the legend — is identical in all three, so what
changes between them is only the result.

**First GRAND alone.** It is *outlined* rather than filled: it covers most of the
frame, and a translucent wash over that much of a map hides everything underneath.
Notice where the boundary runs — GRAND's constraint is deployable slope, so it accepts
the high plateau and declines the canyon walls and the steep coast."""),
("code", """def show_stage(name, caption):
    found = os.path.join(OUT, "arequipa_full_combined", name)
    show_figure(found, caption=caption)

show_stage("combined_overview_1_grand.png",
           "GRAND alone — 88,527 km², essentially the whole plateau")"""),
("md", """**Now TAMBO alone.** A different question entirely, and it shows: instead of a
boundary enclosing most of the map, a scatter of small patches strung along the canyon
systems. TAMBO needs a wall to stand on facing a wall to watch, and that exists only
where the ground is cut."""),
("code", """show_stage("combined_overview_2_tambo.png",
           "TAMBO alone — 112 km², following the canyons")"""),
("md", """**And both.** The magenta is the ground that satisfies the two at once: 50.2 km²,
0.1% of GRAND's and 44.9% of TAMBO's. The asymmetry is the finding — co-location costs
GRAND nothing and is most of what TAMBO has.

Look at where the magenta *is*. It traces the canyon rims, and the roads run along them
too, which is not a coincidence: a canyon rim is where a road goes in this terrain."""),
("code", """show_stage("combined_overview_3_both.png",
           "Both — the 50.2 km² that satisfies GRAND and TAMBO together")"""),
("md", """## What moves the answer

The numbers above come from one setting of every parameter. The honest question is how
much they would move under another, and that cannot be read off a single run.

The full DEM is far too slow to sweep — GRAND alone is 25 minutes a point — so what
follows searches a **small synthetic canyon** instead: a plateau cut by a gorge, a few
hundred pixels across, seconds per run. The absolute numbers mean nothing. The
*direction and steepness* of each response is the point, and those carry over."""),
("code", """import contextlib
import io
import tempfile

import numpy as np
import tifffile as tiff

from oroscope import site_searcher as ss

WORK = tempfile.mkdtemp(prefix="oroscope_nb08_")


def plateau_with_canyon(n=400, depth_m=1300.0, floor_px=30, wall_px=35):
    \"\"\"
    A high plateau cut by a gorge: ground a canyon search has an opinion about.

    Sized against the criteria rather than drawn freehand. At ~30 m pixels the rim-to-
    rim distance is 3.0 km, so the far wall sits inside a 1-4 km window, and the walls
    are ~51 deg, inside the 20-60 deg band. A narrower or steeper canyon returns
    nothing at all, which is a lesson about the criteria but a poor demonstration.
    \"\"\"
    col = np.arange(n)
    wall = np.clip((np.abs(col - n // 2) - floor_px) / wall_px, 0.0, 1.0)
    z = 3600.0 + wall[None, :] * depth_m
    # Gentle relief, so the plateau is not perfectly flat and slope has a distribution
    rr, cc = np.mgrid[0:n, 0:n]
    z = z + 90.0 * np.sin(rr / 47.0) + 60.0 * np.cos(cc / 61.0)
    return z.astype(np.float32)


dem = os.path.join(WORK, "canyon.tif")
tiff.imwrite(dem, plateau_with_canyon(), extratags=[
    (33550, "d", 3, (1 / 3600, 1 / 3600, 0.0)),
    (33922, "d", 6, (0.0, 0.0, 0.0, -72.0, -15.5, 0.0)),
])
print("synthetic canyon written:", os.path.basename(dem))"""),
("code", """BASE = dict(
    dem_path=dem, search_mode="distributed", grid_type="hex",
    min_slope_deg=20.0, max_slope_deg=60.0, min_target_slope_deg=25.0,
    min_dist_km=1.0, max_dist_km=4.0,
    elev_min_deg=-20.0, elev_max_deg=20.0, n_elev_bins=16,
    antenna_spacing_km=0.1, min_sub_array_size=20, min_width_km=0.0,
    target_antennas=5000, grammage_mode="particle",
    grammage_band_gcm2=(236.0, 1287.0),
    decay_energy_min_pev=3.0, decay_energy_max_pev=1000.0,
    solid_angle_half_sr=0.8, min_score=0.05,
    downsample_factor=1, tile_size=256, candidate_stride=3, num_cores=2,
    generate_kml=False, explain=False,
)


def search(tag, **overrides):
    \"\"\"One run, returning the few numbers worth comparing.\"\"\"
    params = dict(BASE, **overrides)
    quiet = io.StringIO()
    with contextlib.redirect_stdout(quiet), contextlib.redirect_stderr(io.StringIO()):
        res = ss.find_grand_regions_interactive(
            run_output_dir=os.path.join(WORK, tag), **params)
    sites = res["results"]["sites"]
    funnel = res["funnel"]
    strided = next(v for k, v in funnel.items() if k.startswith("kept by stride"))
    return dict(area=sum(s["area_km2"] for s in sites),
                sites=res["results"]["total_sites"],
                capacity=res["results"]["total_capacity"],
                accepted=100.0 * funnel["directions accepted"] / max(strided, 1))


baseline = search("baseline")
print(f"baseline: {baseline['sites']} sites, {baseline['area']:.1f} km², "
      f"{baseline['capacity']:,} detectors, {baseline['accepted']:.1f}% accepted")"""),
("md", """### One parameter at a time

Four knobs, each swept while everything else is held. Watch which ones bend the answer
and which barely touch it — and in particular watch where a response is *not* smooth,
because a cliff is where a result stops being a measurement and starts being a choice."""),
("code", """def sweep(label, key, values, fmt="{}"):
    print()
    print(label)
    print(f"   {'value':>12} {'sites':>6} {'area km²':>10} {'capacity':>10} {'accepted':>10}")
    print("   " + "-" * 52)
    for i, v in enumerate(values):
        r = search(f"{key}_{i}", **{key: v})
        mark = "  <- baseline" if v == BASE.get(key) else ""
        print(f"   {fmt.format(v):>12} {r['sites']:>6} {r['area']:>10.1f} "
              f"{r['capacity']:>10,} {r['accepted']:>9.1f}%{mark}")


sweep("Arrival window half-width (elev_min/max_deg)", "elev_max_deg",
      [5.0, 10.0, 20.0, 30.0], "±{:.0f}°")"""),
("code", """sweep("Far-wall slope floor (min_target_slope_deg)", "min_target_slope_deg",
      [25.0, 45.0, 50.0, 55.0], "{:.0f}°")"""),
("code", """sweep("Distance window, far edge (max_dist_km)", "max_dist_km",
      [2.0, 3.0, 4.0, 6.0], "{:.0f} km")"""),
("code", """sweep("Score cut (min_score)", "min_score",
      [0.0, 0.05, 0.10, 0.15], "{:.2f}")"""),
("md", """### What that shows

Three of these are physics and one is not, and they do not behave alike.

**The arrival window** widens the accepted set smoothly: more sky examined, more
directions accepted. It is the binding constraint on the real search too — both GRAND
and TAMBO are bound by `directions accepted` over Arequipa — so this is the response
that most directly sets the size of the answer.

**The far-wall floor** is a threshold compared against a physical distribution, and it
behaves exactly as that should: nearly flat while it sits below the bulk, then a cliff
as it crosses. These walls are ~51°, and the floor takes 52.7% of directions at 25°,
45.2% at 45°, 34.4% at 50° — and **nothing at all at 55°**. Ask for a wall steeper than
the terrain has and the search correctly returns an empty answer.

That is worth dwelling on, because it is the shape every terrain threshold has. The
question is never "is this value reasonable" but "where does it sit relative to the
distribution the ground actually offers", and only a scan like this answers it.

**The distance window** has to reach the far wall at all. At 2 and 3 km this canyon
returns *nothing*: the opposite wall is 3.0 km rim to rim, so a window that stops short
of it finds no target, however good the ground underfoot. Past that it opens up. A
window is two constraints, and the far edge is doing the work here.

**The score cut is not physics at all.** It is a threshold on a *product* of components,
and a product of numbers in [0, 1] piles up near zero — the site scores here have a
median near 0.10, so a cut of 0.15 is already biting into the body of the distribution
while 0.0 keeps everything. Over the real Colca crop `min_score` 0.35 is equivalent to
keeping the top 22.8% by rank, and a scan across the cut shows **no knee anywhere**:
nothing in the data marks 0.35 as the natural place to stand. It is the single most
consequential choice in the TAMBO configuration and the one least constrained by the
terrain.

Note also what `min_target_slope_deg: 0` means, since it is a trap: zero is *falsy*, so
it switches the far-wall criterion **off** rather than setting a 0° floor. Off is not
the permissive end of this sweep — it is a different search.

That last one is why `--score_percentile` exists: a rank means the same thing when the
composition changes, and an absolute cut on a product does not."""),
("md", """## Every assumption behind these numbers

A search produces authoritative-looking areas and detector counts, and the only defence
against those being over-read is writing down what is behind them. This is the complete
list for *this* run — not a generic one. Where an effect has been measured, the measured
size is given; where it has not, that is said instead.

The run says most of this about itself: `explanation.txt` carries a
**WHICH OF THESE ARE ASSUMPTIONS** block, and the cell below prints it. What follows adds
the ones that live in the configuration rather than the pipeline, and the ones measured
after the run."""),
("code", """for label, res in (("GRAND", grand_full), ("TAMBO", tambo_full)):
    if res is None:
        continue
    text = res.get("explanation") or explain.explain_results(res)
    block = text.split("WHICH OF THESE ARE ASSUMPTIONS")
    if len(block) > 1:
        print(f"===== {label} =====")
        print("WHICH OF THESE ARE ASSUMPTIONS" + block[1].split("WHAT TO TRY NEXT")[0])"""),
("md", """### 1. Choices written into the configuration

Both configs mark these `ASSUMPTION` in their own comments. They are the numbers a
collaboration should be asked about before any of this is quoted.

| | value | what it rests on |
|---|---|---|
| **`solid_angle_half_sr`** (TAMBO) | 0.8 sr | The scale at which accepted sky stops being scarce. The default 0.05 is GRAND-scale and saturates the term to ~1 across a canyon, so it stops discriminating. **This is what the TAMBO result actually turns on** — `solid_angle` is the weakest component at 26 of 26 sites. |
| **`min_score`** (TAMBO) | 0.35 | A cut on a *product* of components, whose distribution piles up near zero. Equivalent to `score_percentile` **22.8** on this terrain, and a scan across the cut shows **no knee** — area runs 3.1 → 186.8 km² across percentiles 5 → 40, smooth and near linear. Nothing marks 0.35 as natural. |
| **near-wall band** (TAMBO) | 20–60° | The wall the array stands on. Colca's walls are ~40°, far outside GRAND's 3–25°. |
| **`min_target_slope_deg`** (TAMBO) | 25° | The far wall the tau exits, measured along the arrival azimuth. Deliberately permissive against ~40° walls. Without it the scan only asks that rock is present at the right range, which is true almost everywhere in the Andes. |
| **arrival window** (TAMBO) | ±20° | What the far wall subtends from a detector on the near wall. |
| **`grammage_band_gcm2`** (TAMBO) | 236–1287 | = `grammage_band_from_energy(3, 1000, fraction=0.1)`. The 0.1 is a choice about detector capability — how far down the profile still counts as a usable shower — not a property of the shower. |
| **`decay_spectral_index`** (TAMBO) | 2.0 | The flux slope the decay term folds against. The canonical value. Costs 1.46× across a plausible range; may be given as a `(low, high)` pair to marginalise instead of choosing. |
| **arrival window** (GRAND) | ±3° | The Earth-skimming geometry. Note the ±3° window sits *below* the horizon almost everywhere: the median horizon is +7.3°. |
| **distance window** (GRAND) | 10–40 km | Far enough for the shower to develop, close enough for the signal. |

### 2. Choices forced by running at full scale

These are not physics. They are the price of searching 128.6 Mpx on one machine, and
each one moves a reported number.

| | value | what it costs |
|---|---|---|
| **`downsample_factor`** | 4 | Area is measured on the downsampled mask while capacity is measured at full resolution, so a feature a few pixels wide keeps its detectors and loses area — **~30% for a canyon strip**. Affects TAMBO far more than GRAND. |
| **`candidate_stride`** | 5 | Unbiased in *acceptance* — 17.494% against 17.491% at stride 1 — but the mask is closed before area is taken, and a 100 m element cannot bridge the 154 m gap stride 5 leaves. **TAMBO's area is 4.75× low** because of this. GRAND's 1 km element is unaffected. |
| **`gap_close_km`** | = detector spacing | Closing inflates GRAND's area **2.10×** in this run (2.29× measured against a stride-1 control). Reported area is not physics-accepted area. |
| **`max_range_km`** | unset | The profile walk stops at `max_dist_km`, so the reported **column depth is a property of where the walk stopped**. Measured on TAMBO: walking 4× further raised it **6.4×** with an identical selection. Read the depths as lower bounds. |
| **`min_width_km`** | 2.0 GRAND / 0.0 TAMBO | The opening step prunes tendrils. At GRAND's 2 km it would delete exactly the strip a canyon array is, which is why TAMBO sets it to 0. |

### 3. Physics not modelled at all

| | |
|---|---|
| **Detector acceptance `A(E)`** | An event rate is ∫Φ·A·P dE. `decay_weight_by` selects flux, acceptance, or both — **these numbers used `flux`**. Two published integral curves are supplied in `data/`, and what can and cannot be done with them is spelled out below. |
| **Tau production and escape through rock** | Not modelled, so β does not enter. The search weights by the decay length E/m·cτ, which is kinematics. |
| **Neutral-current regeneration** | Available (`nc_regeneration=True`) but **off here**, so Earth-chord suppression is overstated — by 1.06× at −0.5° rising to 1.56× at −5°. |
| **Shower simulation, detector response, trigger** | None of it. The scores rank sites; they are not apertures. |
| **Geology** | One standard rock density throughout. |
| **Geomagnetic declination** | Constant at Arequipa's −6.9° across the whole DEM. Inclination does follow the site. Right for southern Peru, and this DEM is southern Peru. |
| **External validation** | **Nothing here has been checked against an external simulation.** The Earth-absorption prediction — the window's lower edge climbing from −4.4° at 100 PeV to −0.9° at 10 EeV — is the cheapest such test and is ready for someone to run. |

### 4. The published effective areas: what we correct, and what we cannot

`data/` holds two curves from the collaborations. They are the closest thing to a real
`A(E)` this project has, and it is worth being exact about what they can carry.

| file | quantity | array it was simulated for |
|---|---|---|
| `tambo_aperture_fig3` | aperture, m² sr | 5,000 units at 150 m, **in Colca Canyon** |
| `grand_effective_area_fig25` | effective area, cm² (direction-averaged) | 10,000 antennas at 1 km, **at "HotSpot1"** |

Each is the output of a simulation of **one array at one site**. Oroscope changes both,
and only one of those can be fixed by arithmetic.

**The array size can be corrected.** An aperture scales with instrumented *ground*, so
`aperture.array_scale_factor` uses (N·s²)ₜₐᵣ𝒷ₑₜ / (N·s²)ₚᵤ𝒷ₗᵢₛₕₑ𝒹 — detector count *and*
spacing. Both matter: adding detectors at fixed spacing adds ground and scales the
aperture; adding them at fixed ground only makes the array denser, and past the point
where it already samples the Cherenkov cone that buys almost nothing. **TAMBO now runs
at 150 m to match the published simulation**, so the factor is a plain ratio of counts.

The linearity is checked rather than assumed: the GRAND paper states its 200k curve is
exactly 20× the 10k one, and tracing both independently from the figure gave 19.9–20.1×.

**The site cannot be corrected.** Folded into each curve, and inseparable from it, is
that site's distribution of column depth, target distance and arrival elevation, and its
trigger geometry. For TAMBO that terrain is Colca's walls. So a scaled curve says *"what
this many detectors would have achieved on the ground the simulation assumed"* — **not**
*"what they will achieve on this ground"*. Applying the TAMBO curve to an Ancash canyon
imports Colca's rock.

**Where this does and does not matter.** In the *score* it does not enter at all:
`spectrum_weighted_decay_probability` normalises its weights, so any constant on `A(E)`
cancels exactly, and only the shape in energy survives. Site *ranking* therefore does not
depend on any of this being right. It is only an absolute aperture that the scaling
normalises — and that is the number to treat with care.

**This is a workaround, not a simulation.** The right calculation is a full detector
simulation at the candidate site with the candidate layout. Everything above stands in
for that, and is defensible only for array size, only at fixed spacing, and only while
the candidate terrain is not wildly unlike the simulated terrain."""),
("code", """from oroscope import aperture

pub = aperture.PUBLISHED_ARRAYS["tambo_aperture_fig3"]
print(f"published: {pub['units']:,} units at {pub['spacing_km']*1000:.0f} m, {pub['site']}")

for units in (1000, 5000, 20000):
    f = aperture.array_scale_factor(units, "tambo_aperture_fig3")
    print(f"  {units:>6,} units at 150 m -> scale the published curve by {f:>5.2f}x")

# The trap the spacing term exists to avoid: the same ground, counted two ways.
same = aperture.array_scale_factor(11250, "tambo_aperture_fig3", target_spacing_km=0.10)
print("\\n11,250 units at 100 m cover the same ground as 5,000 at 150 m.")
print(f"  correct factor          : {same:.2f}x")
print(f"  by detector count alone : {11250/5000:.2f}x   <- wrong by the density ratio")"""),
("md", """### 5. And the one that is not an assumption at all

The layout is **anchored, not fitted**: detectors are placed from each site's bounding-box
corner rather than optimised. Capacity is an estimate for an arbitrarily placed array, and
a real deployment would do better.

---

**If you quote one number from this notebook, quote it with its caveat.** GRAND's
88,527.5 km² is a closed mask, 2.10× the accepted set. TAMBO's 111.9 km² is low by ~4.75×
from striding and ~30% again from downsampling. The joint 50.2 km² inherits both, and is
a floor rather than an estimate."""),
("md", """## Where to go next

- **[10. Ancash](10_ancash_dem.ipynb)** — the same question over the Cordillera
  Blanca, at the same resolution and with every criterion unchanged, so the difference
  in the answer is a difference in the ground.
- **[12. Peru, all of it](12_peru_dem.ipynb)** — the whole country at 3 arc-seconds,
  and how to read a coarse answer honestly."""),
("md", footer(prev=("08_explaining_a_run.ipynb", "Explaining a run"),
              nxt=("10_ancash_dem.ipynb", "Ancash"))),
]

# --------------------------------------------------------------------------- 07  animations
STILLS_HELPER = """from IPython.display import Image, display


def show_stills(stills, width=760):
    \"\"\"Stacks a few frames of one animation into one image, small enough to commit.

    A notebook whose outputs live in a repository cannot store a playable video, so what
    it stores is three frames read top to bottom: the beginning, the middle, the end.
    The MP4 written beside it plays.
    \"\"\"
    import io
    from PIL import Image as PILImage

    frames = [PILImage.open(io.BytesIO(s)).convert("RGB") for s in stills]
    scale = width / max(f.width for f in frames)
    frames = [f.resize((round(f.width * scale), round(f.height * scale)),
                       PILImage.LANCZOS) for f in frames]
    gap = 10
    sheet = PILImage.new("RGB", (width, sum(f.height for f in frames)
                                 + gap * (len(frames) - 1)), "white")
    y = 0
    for f in frames:
        sheet.paste(f, (0, y))
        y += f.height + gap
    buf = io.BytesIO()
    sheet.save(buf, format="JPEG", quality=78, optimize=True)
    display(Image(data=buf.getvalue()))"""

BUILD_HELPER = """def build(name, at=(0.0, 0.5, 1.0)):
    \"\"\"Builds one animation, writes its MP4, and shows stills from the same pass.\"\"\"
    fig, anim = ma.BUILDERS[name]()
    path, stills = ma.write_mp4_with_stills(name, fig, anim, OUT, at=at)
    print(f"{os.path.relpath(path, '..')}  "
          f"({os.path.getsize(path) / 1024:,.0f} KiB, {anim._save_count} frames)")
    show_stills(stills)"""

NB_ANIMATIONS = [
("md", """# 7. Animating the mechanism

Some of what this project does is a **process**, and a process is badly served by a
still picture. A ray sweeping down through the elevation window, a map draining stage by
stage, a score distribution walking under a threshold that never moves — in each case the
intermediate states *are* the argument, and a figure of the end state throws them away.

`tools/make_animations.py` builds eight of these. This notebook builds all eight, says
what each one is for, and shows the frames. It is the only notebook here whose output is
a set of files rather than a set of numbers.

**The filter was strict, and worth stating.** An animation is warranted when the
intermediate states carry the argument — not when the subject happens to be able to
move. Six candidates were rejected on exactly that test: a continuous sweep of the
closing element (`stride_and_closing` already shows the transition, and it is abrupt),
sin α against azimuth (a polar plot does it statically and better), decay probability
against distance and energy (that is a contour plot), the Earth chord against elevation
(one static diagram), Cherenkov footprint growth (a simple monotone relation), and an
array filling a site (decorative). Being able to animate something is not a reason
to."""),
# No numpy and no pyplot: this notebook computes nothing and draws nothing itself. The
# arrays and the axes are all inside the builders, and `ruff check .` lints notebooks,
# so an import here for the sake of a house style would fail CI.
("code", """import os
import sys

# tools/ is a directory of scripts, not an installed package, so it has to be put on the
# path. This is *not* the sys.path-insert that these notebooks otherwise avoid: that one
# pointed at src/ and shadowed an installed oroscope with whatever was in the tree. This
# one reaches a repository tool that is deliberately not shipped.
sys.path.insert(0, os.path.abspath(os.path.join("..", "tools")))
import make_animations as ma

OUT = os.path.abspath(os.path.join("..", "output", "animations"))
print(", ".join(ma.BUILDERS))"""),
("md", """## What is needed to run this

**`ffmpeg`, for the MP4s.** Every cell below writes an MP4 and nothing else; GIFs are a
conversion, covered at the end. Without `ffmpeg` the writes fail and the stills still
appear, which is enough to read the notebook but not enough to have the files.

**A DEM, optionally.** Five of the eight are built from committed code and synthetic
terrain and reproduce on any clone. Three — `the_azimuth_fan`, `product_collapse` and
`slope_criterion` — are about what a criterion does to *real* ground, which synthetic
terrain cannot honestly show, so they read `input/dem/colca.tif` when it is present and
fall back to a synthetic canyon when it is not. Each says on the figure which it used,
so a frame is never ambiguous about what it is a picture of.

Outputs land in `output/animations/`, which is gitignored."""),
("code", STILLS_HELPER),
("code", BUILD_HELPER),
("md", """---

## The mechanism, in two halves

Everything downstream — every criterion, every score, every map — rests on one
operation: from a candidate pixel, walk outward along a bearing and find where the ray
first meets terrain. The two animations below are the two axes of that walk.

### `the_walk` — sweeping the elevation

One backward ray sweeping down through the elevation window over a terrain profile. The
first intersection slides along the profile as the angle steepens, and the column depth
behind it accumulates. This is the single hardest thing in the project to convey in
prose, and the reason the elevation binning is nearly free: one pass over the profile
fills every bin at once, because the running maximum of the apparent terrain angle only
ever increases.

Watch the lower panel. Depth is plotted against the **angle being swept**, not against
distance — sharing the x-axis with the profile would make it look like a property of the
terrain, which it is not."""),
("code", """build("the_walk")"""),
("md", """### `the_azimuth_fan` — sweeping the bearing

The other half, and the one with no static counterpart. `the_walk` sweeps elevation at a
fixed bearing; this fixes the elevation window and sweeps the bearing through a full
360°, reporting what each one finds.

Three outcomes, and the whole point is that they are not evenly distributed. Bearings
near the candidate's aspect find the far canyon wall at two to three kilometres — inside
the accepted range, steep enough to be a wall. Bearings at right angles find the
candidate's *own* hillside a few hundred metres off, which is rock at the wrong range and
scores nothing. Bearings behind it are blocked at zero. The shaded wedge is the fan the
search actually tests: nine bearings within 60° of aspect, from
`arrival_scan.azimuth_fan(9, 60.0)`.

Acceptance against bearing is a polar quantity, and a polar quantity that changes is
exactly what a still figure renders badly — which is why this one exists and the sin α
version does not."""),
("code", """build("the_azimuth_fan")"""),
("md", """---

## What the pipeline does to a map

### `the_funnel` — where the candidates went

The funnel table in `--explain` says *how many* candidates each stage removed. It cannot
say **where on the ground** they were. Slope, then stride, then directions accepted, then
gap closing, then pruning — each stage drawn on the map with the surviving count.

The stage that surprises people is closing, which *adds* pixels. It is not a filter; it
is a repair of the holes striding left."""),
("code", """build("the_funnel")"""),
("md", """### `stride_and_closing` — why TAMBO's area was low by 4.75×

The measurement is in `docs/ROADMAP.md` §6.34; this is what it looks like. A strided mask
closed with an element that bridges the gaps, and the same mask closed with one that does
not.

The transition is **at** the gap and it is abrupt, not gradual: on this mask a 3-pixel
element recovers 0.04× of the accepted set and a 5-pixel one recovers 0.61×, fifteen
times more for two pixels of element. TAMBO's closing element is 3 pixels against a
5-pixel gap, which is how a real published area came out 4.75× low while every
intermediate number looked reasonable. Every run now warns
(`ss.warn_stride_outruns_closing`)."""),
("code", """build("stride_and_closing")"""),
("md", """### `slope_criterion` — where a criterion bites

The sensitivity sweeps in [notebook 6](06_combining_and_sensitivity.ipynb) say *how
much* a criterion costs. They cannot say *where*. This runs a full arrival scan for each
value of `min_target_slope_deg` and morphs the accepted mask over Colca as the cut climbs
through the wall-slope distribution beside it.

The surroundings go first and the canyon rims hold on longest, which is the criterion
doing exactly what it was written to do — separate a canyon from a hillside.

There is a lesson in the two panels disagreeing. Half the candidates see a mean wall
slope below 29.7°, yet a 30° floor still keeps 87% of them, and the half-way point is not
reached until 50°. **The mask outlives its own median by 20°**, because the criterion is
applied to each *direction* while the histogram is a mean over each candidate's accepted
directions — so a candidate keeps its steepest directions long after its average has
fallen under the cut. Read `target_slope_deg` as a description of a candidate, never as a
prediction of what a cut will do to it."""),
("code", """build("slope_criterion")"""),
("md", """---

## What the score does

### `product_collapse` — why a threshold on a product is treacherous

`min_score` is the dominant assumption in this project, and the reason a product
threshold is dangerous is **dynamic**: each component multiplied in drags the whole
population toward zero while the cut stays exactly where it was put. Six real components
of a real search, folded in one at a time, against TAMBO's own cut of 0.35. The weight of
the newest component ramps from 0 to 1 — under a product composition weights *are*
exponents, so every intermediate frame is a state `scoring.compose` can genuinely
produce, not a dissolve between two pictures.

Then measurement contradicted the premise, which is the useful part. The collapse is
real — 100% of viable candidates above the cut before any component, 32.2% after six —
but it is **not evenly shared**. `solid_angle` alone takes it from 100% to 35.9%.
`depth` does nothing, because Colca's walls sit inside the default band. And `distance`
does nothing *provably*: the scan already applied the same 2–5 km window as a hard
criterion, so every surviving candidate scores 1 on it by construction. Scoring a
criterion the scan has already enforced is free, but it is also empty.

That sharpens the case against thresholding a product rather than weakening it. The cut
goes from harmless to decisive on the addition of a single term — and which term that is
depends on the terrain, not on the configuration."""),
("code", """build("product_collapse")"""),
("md", """---

## What the physics says

### `tau_in_rock` — more rock is not better

The commonest misconception in this problem, corrected. A tau's energy and its survival
probability both falling as it burrows, against the column depth that maximises
production and escape *together*: `physics.production_escape_optimum_gcm2` gives
5.7×10⁶ g/cm² at 1 EeV, about 22 km of standard rock.

Production grows with the rock available — more target, more neutrino interactions.
Escape collapses once the depth passes the tau's range. Their product therefore has a
maximum, which is why the criterion is a **band** and not a floor. The code has always
had this right; until now it lived only in prose."""),
("code", """build("tau_in_rock")"""),
("md", """### `energy_window` — a prediction, animated over what it predicts

The arrival window narrowing as the energy rises, its lower edge climbing from −4.4° at
100 PeV to −0.9° at 10 EeV as Earth absorption removes the steeply upgoing directions.

This one is here for a different reason from the rest. It is not an explanation of
something the code does; it is a **falsifiable claim about the world** that the code
makes, and nothing in this project has yet been checked against an external simulation.
It is the cheapest such test available and it is ready for someone to run."""),
("code", """build("energy_window")"""),
("md", """---

## Turning an MP4 into an animated GIF

Everything above wrote MP4 only. MP4 is the right format to produce: it is a tenth the
size at the same quality, and it seeks. But a GIF plays inline in places an MP4 does not
— a GitHub README, a wiki, a chat message, a slide that has to survive someone else's
projector — so the conversion is worth having.

**Route one: ask the tool.** `make_animations.py` writes both formats by default and
takes `--format`, so the GIF never needs a conversion step at all:

```bash
python tools/make_animations.py --format gif                 # GIFs only
python tools/make_animations.py --format mp4,gif             # both, the default
python tools/make_animations.py --only the_walk --format gif # one of them
```

This route falls back to pillow when `ffmpeg` is absent, so it always produces
*something*. It is the right choice when the animation is being built now.

**Route two: convert an MP4 you already have.** This is the one to use for a file that
came from somewhere else, or when the MP4 took twenty seconds of terrain scanning to
produce and rebuilding it to change the format would be silly. Plain `ffmpeg`, in one
pass:

```bash
ffmpeg -i the_walk.mp4 \\
  -vf "fps=12,scale=900:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" \\
  -y the_walk.gif
```

The `palettegen`/`paletteuse` pair is what makes this worth doing properly. A GIF holds
256 colours; without a palette built from the actual footage, `ffmpeg` falls back to a
fixed web palette and shaded terrain bands horribly. With it, the same animation came out
**342 KiB against pillow's 761 KiB** and looked better.

Three knobs, in the order you will want them:

| | |
|---|---|
| `fps=12` | Match `ma.FPS`. Raising it above what was rendered duplicates frames and inflates the file for nothing. |
| `scale=900:-1` | Width in pixels; `-1` keeps the aspect ratio. This is the strongest size lever by far. |
| `-ss`/`-t` | Trim, if only part of the animation is wanted: `-ss 2 -t 4` takes four seconds from the two-second mark. |

To convert everything in one go:

```bash
cd output/animations
for f in *.mp4; do
  ffmpeg -loglevel error -i "$f" \\
    -vf "fps=12,scale=900:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" \\
    -y "${f%.mp4}.gif"
done
```

**If `ffmpeg` is unavailable**, route one still works — pillow is a hard dependency of
matplotlib, so `--format gif` cannot fail for want of an encoder. The files are larger
and the colours are worse, and for a talk that is usually a fine trade."""),
("md", """## Stills, for a paper

`ma.write_mp4_with_stills` returns the frames it grabbed along with the path it wrote, so
a figure for a paper or a slide costs nothing beyond the animation that was being built
anyway. It is one function rather than two for a reason worth knowing if you extend the
tool: **the builders accumulate.** The ray drawn at frame 30 is still on the axes at
frame 60, which is what makes the fan fill in — so the frames can be walked exactly once.
A second pass would start with everything already drawn, and its "first" frame would be a
lie."""),
("code", """fig, anim = ma.BUILDERS["energy_window"]()
path, stills = ma.write_mp4_with_stills("energy_window", fig, anim, OUT,
                                        at=(0.0, 0.35, 0.7, 1.0))
print(f"{len(stills)} stills at the video's own frame size, "
      f"{sum(len(s) for s in stills) / 1024:,.0f} KiB of PNG")"""),
("md", """## Where to go next

- **[2. The arrival scan](02_the_arrival_scan.ipynb)** — the walk that `the_walk` and
  `the_azimuth_fan` animate, in code.
- **[4. Criteria and scoring](04_criteria_and_scoring.ipynb)** — the component shapes
  that `product_collapse` multiplies together.
- **[9. Arequipa, the full DEM](09_arequipa_dem.ipynb)** — what the whole machine
  produced when it was finally pointed at everything.
- **[12. Peru, all of it](12_peru_dem.ipynb)** — the same machine pointed at a
  country, and how to read a coarse answer honestly."""),
("md", footer(prev=("06_combining_and_sensitivity.ipynb", "Combining and sensitivity"),
              nxt=("08_explaining_a_run.ipynb", "Explaining a run"))),
]

# --------------------------------------------------------------------------- 10  ancash
NB_ANCASH = [
("md", """# 10. Ancash

[Notebook 9](09_arequipa_dem.ipynb) ran the whole machine over Arequipa. This one runs
**the same three searches over different ground** — the Cordillera Blanca and the
Callejón de Huaylas, 300 km of the steepest tropical mountains on Earth, with Huascarán
at 6,768 m in the middle of them.

The point is not that Ancash is new terrain. It is that **every transferable criterion
is held fixed**, so a difference in the answer is a difference in the ground rather than
a difference in the question. That is a claim this notebook checks rather than asserts,
in the second cell below."""),
("code", """import json
import os

import numpy as np

from oroscope import explain
from oroscope import site_searcher as ss"""),
("md", """## What is the same, and what is not

Both configurations were derived from Arequipa's by changing the DEM path and the label,
and nothing else. Rather than ask you to believe that, here is the diff."""),
("code", """def settings(path):
    \"\"\"A configuration without its commentary: `_comment*` keys are prose.\"\"\"
    with open(path) as f:
        return {k: v for k, v in json.load(f).items() if not k.startswith("_")}


for experiment in ("grand", "tambo"):
    a = settings(os.path.join("..", "config", f"{experiment}_arequipa_full.json"))
    n = settings(os.path.join("..", "config", f"{experiment}_ancash_full.json"))
    keys = sorted(set(a) | set(n))
    differ = [k for k in keys if a.get(k) != n.get(k)]
    print(f"{experiment.upper()}: {len(keys)} settings, {len(differ)} differ")
    for k in differ:
        print(f"    {k:<14} arequipa={a.get(k)!r}")
        print(f"    {'':<14} ancash  ={n.get(k)!r}")"""),
("md", """**TAMBO is an exact repeat**: the two differences are the file it reads and the name it
prints.

**GRAND differs in one real setting, `rfi_zones`.** Arequipa's run excludes five
hand-curated circles — Arequipa city, Majes, Cerro Verde, La Joya, Mollendo — and there
is no Ancash preset. Inventing a five-circle Ancash list would have injected a new
assumption into a run whose entire purpose is comparison, so this run excludes nothing
and says so. The size of it: Arequipa's zones cover about 3,500 km² of a ~120,000 km²
box, so **read Ancash's GRAND area as at most ~3% flattered** on that account.

The real settlements are not ignored — 774 of them, from OpenStreetMap, are on the maps
below. They are context here, not a filter."""),
("md", """## Read, not run

The cells below open results produced locally and stored in `results/ancash_full/`. To
produce or refresh them:

```bash
export OPENTOPOGRAPHY_API_KEY=...          # free, see the CLI page
cd src && oroscope-fetch-dem --region ancash
python -m oroscope.fetch_roads --dem ../input/dem/ancash_SRTMGL1.tif --places

python tools/run_full_dem.py --region ancash --dry-run
python tools/run_full_dem.py --region ancash
```

Ancash is 69 Mpx against Arequipa's 129, so it costs roughly half as much: about
thirteen minutes for GRAND and under a minute for TAMBO."""),
("code", """STORE = os.path.abspath(os.path.join("..", "results", "ancash_full"))
AREQUIPA_STORE = os.path.abspath(os.path.join("..", "results", "arequipa_full"))


def load(store, label):
    path = os.path.join(store, f"{label}_results.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


runs = {label: load(STORE, label) for label in ("grand", "tambo")}
reference = {label: load(AREQUIPA_STORE, label) for label in ("grand", "tambo")}
for label, r in runs.items():
    print(f"ancash {label:<6} {'loaded' if r else 'NOT IN THE STORE YET'}")
for label, r in reference.items():
    print(f"arequipa {label:<4} {'loaded' if r else 'not present'}")"""),
("md", """## The terrain, before the searches

The two DEMs can be compared without running anything, and doing so makes a prediction
worth holding the results against. Slope is computed in tiles so that a 129 Mpx DEM does
not need 1.5 GB of gradients at once."""),
("code", """def slope_profile(tif, origin_lat, bands=(("GRAND 3-25", 3.0, 25.0),
                                          ("TAMBO 20-60", 20.0, 60.0)), block=2048):
    \"\"\"Median slope and band shares over a DEM, accumulated tile by tile.\"\"\"
    import tifffile
    grid = ss.resolve_grid_geometry(tif, origin_lat)
    z = tifffile.imread(tif)
    rows = z.shape[0]
    hist = np.zeros(9001)                       # 0.01 degree bins, 0-90
    counts = {name: 0 for name, _, _ in bands}
    total = 0
    for r0 in range(0, rows, block):
        r1 = min(rows, r0 + block)
        lo, hi = max(0, r0 - 1), min(rows, r1 + 1)
        tile = z[lo:hi].astype(np.float32)
        dy, dx = np.gradient(tile, grid.cell_size_y, grid.cell_size_x)
        slope = np.degrees(np.arctan(np.hypot(dy, dx)))[r0 - lo:r1 - lo]
        land = z[r0:r1] > 0                     # the sea is exactly 0 and has no slope
        s = slope[land]
        if not s.size:
            continue
        total += s.size
        hist += np.histogram(s, bins=9001, range=(0.0, 90.0))[0]
        for name, a, b in bands:
            counts[name] += int(((s >= a) & (s <= b)).sum())
    median = float((np.searchsorted(np.cumsum(hist), total / 2.0) + 0.5) / 100.0)
    return {"Mpx": z.size / 1e6, "median slope": median,
            **{name: 100.0 * counts[name] / total for name, _, _ in bands}}


DEMS = {"Arequipa": ("arequipa_SRTMGL1.tif", -14.5553),
        "Ancash":   ("ancash_SRTMGL1.tif",   -8.04958)}
terrain = {}
for name, (stem, lat) in DEMS.items():
    tif = os.path.abspath(os.path.join("..", "input", "dem", stem))
    if os.path.exists(tif):
        terrain[name] = slope_profile(tif, lat)

if terrain:
    cols = list(next(iter(terrain.values())))
    print(f"{'':<10}" + "".join(f"{c:>16}" for c in cols))
    for name, row in terrain.items():
        print(f"{name:<10}" + "".join(f"{row[c]:>16.1f}" for c in cols))
else:
    print("neither DEM is here; the measured values are quoted in the text below")"""),
("md", """**Ancash is twice as steep.** Median slope 23.0° against Arequipa's 11.1°, and the
consequence for the two experiments runs in opposite directions:

| share of land in the band | Arequipa | Ancash |
| --- | --- | --- |
| GRAND, 3–25° deployable | 70.3% | **52.0%** |
| TAMBO, 20–60° near wall | 24.1% | **58.0%** |

So before a single ray is traced: **Ancash should be worse for GRAND per unit area and
much better for TAMBO.** GRAND wants ground gentle enough to stand an array on and
Ancash keeps giving it cliffs; TAMBO wants exactly those cliffs. Hold that against what
the searches actually found."""),
("md", """## What the searches found"""),
("code", """def summary(results):
    if not results:
        return None
    r = results["results"]
    sites = r.get("sites") or []
    return {"sites": r.get("total_sites"),
            "area km2": sum(s.get("area_km2", 0.0) for s in sites),
            "capacity": r.get("total_capacity")}


print(f"{'':<22}{'sites':>8}{'area km2':>14}{'capacity':>12}")
for label in ("grand", "tambo"):
    for region, table in (("Arequipa", reference), ("Ancash", runs)):
        s = summary(table.get(label))
        if s:
            print(f"{label.upper() + ', ' + region:<22}{s['sites']:>8,}"
                  f"{s['area km2']:>14,.1f}{s['capacity']:>12,}")"""),
("code", """for label in ("grand", "tambo"):
    a, n = summary(reference.get(label)), summary(runs.get(label))
    if not (a and n):
        continue
    print(f"{label.upper()}  Ancash / Arequipa")
    for k in ("area km2", "capacity", "sites"):
        if a[k]:
            print(f"    {k:<10} {n[k] / a[k]:6.2f}x   "
                  f"({n[k]:,.0f} against {a[k]:,.0f})")
    print()"""),
("md", """Ancash is 0.533× Arequipa's pixel count, so that is the ratio to read everything
against: **anything near 0.53× means "the same ground, less of it", and a departure from
it is the terrain talking.**

And the terrain talks loudly. Per pixel, **GRAND is 0.91× and TAMBO is 2.93×** — the
prediction confirmed in both directions, from a comparison made before either search
ran.

GRAND's loss is milder than the naive 0.74× the band shares suggested, because its 1 km
closing element fills in around a mask that steep ground fragments. TAMBO's gain
*exceeds* its naive 2.41%, because **both** of its stages improve at once: the slope
screen keeps 33.9 million pixels in Ancash against 26.8 million in Arequipa — more
candidates from a DEM half the size — and the share of those that then accept a
direction rises from 9.7% to 15.1%. GRAND's acceptance moves the other way, 61.6% down
to 54.9%.

**So Ancash is a worse GRAND site and a much better TAMBO site than Arequipa, per unit
of ground.** That is the result, and it is a statement about mountains rather than about
software: TAMBO needs canyon walls and the Cordillera Blanca is made of them."""),
("md", """## Where each run lost its candidates"""),
("code", """for label in ("grand", "tambo"):
    print(f"=== {label.upper()} ===")
    for region, table in (("Arequipa", reference), ("Ancash", runs)):
        results = table.get(label)
        if not results:
            continue
        funnel = results["funnel"]
        binding = explain.binding_constraint(funnel)
        stages = list(funnel.items())
        print(f"  {region:<9} binds at: {binding['stage']}  "
              f"({100 * binding['kept_fraction']:.1f}% kept there)")
        print(f"  {'':<9} knob:     {binding['knob'].split(',')[0]}")
        print(f"  {'':<9} overall:  {stages[0][1]:>13,} px -> {stages[-1][1]:>12,} "
              f"({100 * stages[-1][1] / stages[0][1]:.2f}%)")
    print()"""),
("md", """**GRAND's binding constraint changed between the two regions, and TAMBO's did
not.** At Arequipa, GRAND was limited by `directions accepted` — plenty of deployable
ground, and the arrival geometry decided. At Ancash it is limited by
`slope 3.0-25.0 deg`: the mountains simply do not offer enough ground gentle enough to
stand an array on, and the search never gets as far as asking what that ground can see.

TAMBO binds at `directions accepted` in both, but from opposite sides — 9.7% of its
strided candidates accept a direction at Arequipa against 15.1% at Ancash.

This is the single most useful line in the comparison: **a criterion that binds is a
statement about the ground, not about the configuration**, and it moved when only the
ground moved."""),
("md", """## The maps"""),
("code", SHOW_HELPER),
("code", """OUT = os.path.abspath(os.path.join("..", "output"))
for label, title in (("grand", "GRAND over Ancash"), ("tambo", "TAMBO over Ancash")):
    import glob
    found = sorted(glob.glob(os.path.join(OUT, f"ancash_full_{label}",
                                          "oroscope_results*.png")))
    show_figure(found[0] if found else os.path.join(OUT, "missing.png"),
                caption=title)"""),
("md", """### Both at once

The combination, revealed one experiment at a time so the overlay is readable. Every
element that is not a category — terrain, colour bar, roads, towns, scale bar, legend —
is identical between the three frames, so what changes is only the result."""),
("code", """def show_stage(name, caption):
    show_figure(os.path.join(OUT, "ancash_full_combined", name), caption=caption)


for frame, caption in (("combined_overview_1_grand.png", "GRAND alone"),
                       ("combined_overview_2_tambo.png", "TAMBO alone"),
                       ("combined_overview_3_both.png", "Both — the joint ground")):
    show_stage(frame, caption)"""),
("code", """def joint(store):
    path = os.path.join(store, "combined_report.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        report = json.load(f)
    overlap = report["pairwise_overlap"]["GRAND & TAMBO"]
    return {"joint km2": report["joint"]["area_km2"],
            "jaccard": overlap["jaccard"],
            "share of TAMBO": 100.0 * overlap["fraction_of_TAMBO"]}


rows = {"Arequipa": joint(AREQUIPA_STORE), "Ancash": joint(STORE)}
rows = {k: v for k, v in rows.items() if v}
if rows:
    fmt = {"joint km2": "{:>18.1f}", "jaccard": "{:>18.5f}",
           "share of TAMBO": "{:>18.1f}"}
    cols = list(next(iter(rows.values())))
    print(f"{'':<10}" + "".join(f"{c:>18}" for c in cols))
    for name, row in rows.items():
        print(f"{name:<10}" + "".join(fmt[c].format(row[c]) for c in cols))"""),
("md", """**The joint area nearly doubles**, and per pixel it is 2.81× — tracking TAMBO
rather than GRAND, which is what you would expect if the joint region is limited by the
scarcer of the two.

The last column is the one to keep. **The joint region is 44.9% of TAMBO's mask at
Arequipa and 43.0% at Ancash** — essentially unchanged across two regions whose terrain
could hardly be more different. Co-location costs GRAND almost nothing and consumes
roughly half of what TAMBO has, wherever you look. The Jaccard index tripling is not the
two experiments agreeing more; it is TAMBO's mask growing while GRAND's shrinks."""),
("md", """---

## Zooming in: the Callejón de Huaylas and the Cañón del Pato

The department run above uses `downsample_factor` 4 and `candidate_stride` 5, because
69 Mpx will not fit in a desktop otherwise. Both cost area. **A crop is small enough to
run without either**, so the Río Santa valley between the Cordillera Blanca and the
Cordillera Negra — 11.4 million pixels, `−8.80…−9.90` lat, `−78.00…−77.20` lon — was cut
out and searched at 1 and 1.

```bash
cd src && oroscope-crop ../input/dem/ancash_SRTMGL1.tif ../input/dem/huaylas.tif \
    --north -8.80 --south -9.90 --west -78.00 --east -77.20
python tools/run_full_dem.py --region huaylas
```"""),
("code", """HUAYLAS = os.path.abspath(os.path.join("..", "results", "huaylas_full"))


def crop_summary(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        d = json.load(f)
    r = d["results"]
    sites = r.get("sites") or []
    return {"sites": r["total_sites"],
            "area km2": sum(s.get("area_km2", 0.0) for s in sites),
            "capacity": r["total_capacity"]}


print(f"{'':<40}{'sites':>8}{'area km2':>12}{'capacity':>12}")
for label in ("grand", "tambo"):
    for tag, path in (
            ("unbiased, ds 1 / stride 1",
             os.path.join(HUAYLAS, f"{label}_results.json")),
            ("control, ds 4 / stride 5",
             os.path.join(HUAYLAS, f"{label}_control_ds4_stride5.json"))):
        s = crop_summary(path)
        if s:
            print(f"{label.upper() + ', ' + tag:<40}{s['sites']:>8,}"
                  f"{s['area km2']:>12,.1f}{s['capacity']:>12,}")"""),
("md", """**The same ground, the same criteria, and only the sampling changed.** GRAND
moves by 1.1×. TAMBO moves by **291× in area and 386× in capacity**, and from 109 sites
to one.

That is not the 4.75× recorded in `docs/ROADMAP.md` §6.34. That figure was measured on
Colca varying the stride alone; here both levers move, on terrain whose accepted strips
are numerous and individually small.

The funnels say exactly where it goes, and it is not where you would guess."""),
("code", """def funnel_of(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)["funnel"]


a = funnel_of(os.path.join(HUAYLAS, "tambo_results.json"))
b = funnel_of(os.path.join(HUAYLAS, "tambo_control_ds4_stride5.json"))
if a and b:
    print(f"{'stage':<34}{'ds 1 / stride 1':>18}{'ds 4 / stride 5':>18}")
    for k in a:
        if k in b:
            print(f"{k:<34}{a[k]:>18,}{b[k]:>18,}")
    acc_a = 100 * a["directions accepted"] / a["slope 20.0-60.0 deg"]
    strided_b = next(v for k, v in b.items() if k.startswith("kept by stride"))
    acc_b = 100 * b["directions accepted"] / strided_b
    print(f"\\nacceptance: {acc_a:.1f}% at stride 1, {acc_b:.1f}% at stride 5")"""),
("md", """**Acceptance is identical — 14.0% either way.** Striding really is unbiased
there, exactly as §6.34 says. Closing differs by only 2.6×. *All* of the 291× happens
between closing and selection, in the region thresholds: at stride 5 the mask fragments
into 7,954 labelled regions, of which 5 clear the area threshold and **one** clears
`min_sub_array_size` of 250 detectors. At stride 1 the mask is contiguous and 109 regions
survive.

So the under-report is **fragmentation meeting a minimum-array-size cut**, not pixels
being miscounted. GRAND never suffers it because a 1 km closing element bridges a 154 m
stride gap without noticing.

**What to take from this.** Every TAMBO area and capacity in this project that came from
a strided, downsampled run is a lower bound by a factor that is terrain-dependent and, in
practice, unbounded — 4.75× at Colca, 291× here. Quote TAMBO numbers from unbiased runs,
or quote them as "at least". The Callejón de Huaylas was effectively invisible to the
department run: its TAMBO mask contributes **1.2 km² inside this crop window** against
the crop's own 855.1 km²."""),
("md", """### The zoom-in, drawn

The crop's own maps. `--reveal` writes the combination one experiment at a time, so the
same frame can be read three ways: GRAND alone, TAMBO alone, and the ground that
satisfies both. Everything that is not a category — terrain, colour bar, roads, towns,
scale bar — is identical between them."""),
("code", """CROP_OUT = os.path.abspath(os.path.join("..", "output", "huaylas_full_combined"))
for frame, caption in (
        ("combined_overview_1_grand.png",
         "GRAND over the crop — one site, 8,295 km²"),
        ("combined_overview_2_tambo.png",
         "TAMBO over the crop at stride 1 — 109 sites along the Río Santa"),
        ("combined_overview_3_both.png",
         "Both — 637 km² of joint ground, 74.5% of TAMBO's mask")):
    show_figure(os.path.join(CROP_OUT, frame), caption=caption)"""),
("md", """**And the same TAMBO search at the department run's sampling**, for the
comparison the numbers above make. One picture is 109 sites; the other is one."""),
("code", """show_figure(os.path.abspath(os.path.join(
                "..", "output", "huaylas_ctl_tambo", "oroscope_results_huaylas.png")),
            caption="TAMBO on the same crop at downsample 4 / stride 5 — the control")"""),
("md", """---

## The full explanation of each run

Every search writes a plain-language account of itself — what it found, where the
candidates went, why each site qualifies, which numbers are assumptions and how to read
them. Those are reproduced here in full rather than summarised, because the caveats
matter as much as the totals."""),
("code", """def explanation(store, name):
    path = os.path.join(store, name)
    if not os.path.exists(path):
        print(f"not in the store: {os.path.basename(path)}")
        return
    with open(path) as f:
        print(f.read())


explanation(STORE, "grand_explanation.txt")"""),
("code", """explanation(STORE, "tambo_explanation.txt")"""),
("code", """explanation(STORE, "combination_explanation.txt")"""),
("code", """explanation(HUAYLAS, "tambo_explanation.txt")"""),
("md", """## Where to go next

- **[9. Arequipa, the full DEM](09_arequipa_dem.ipynb)** — the run this one is held
  against, and the place to read the caveats that apply to both.
- **[12. Peru, all of it](12_peru_dem.ipynb)** — the whole country at a third of the
  resolution, and what that costs.
- **[6. Combining and sensitivity](06_combining_and_sensitivity.ipynb)** — how much of
  any of this survives a change of assumption."""),
("md", footer(prev=("09_arequipa_dem.ipynb", "Arequipa, the full DEM"),
              nxt=("11_lima_dem.ipynb", "Lima"))),
]

# --------------------------------------------------------------------------- 11  lima
NB_LIMA = [
("md", """# 11. Lima

The third department, and the one that closes the set. [Arequipa](09_arequipa_dem.ipynb)
is high plateau, [Ancash](10_ancash_dem.ipynb) is the Cordillera Blanca, and Lima is the
coastal contrast — desert shelf rising to the western Andean flank, with the Cordillera
Huayhuash in its north-east corner.

Three regions, one question, every transferable criterion held fixed. This notebook is
mostly about what the **three-way** comparison says, because a difference between two
regions can be a coincidence and a trend across three is harder to dismiss."""),
("code", """import glob
import json
import os

import numpy as np

from oroscope import site_searcher as ss"""),
("md", """## One dataset for all three

The Lima DEM used to be AW3D30 while Arequipa and Ancash were SRTMGL1. That would have
put a **dataset difference inside every comparison**, indistinguishable from a difference
in the ground, so Lima was re-downloaded as SRTMGL1 at the same 1 arc-second resolution.

```bash
export OPENTOPOGRAPHY_API_KEY=...
cd src && oroscope-fetch-dem --region lima
python -m oroscope.fetch_roads --dem ../input/dem/lima_SRTMGL1.tif --places

python tools/run_full_dem.py --region lima --dry-run
python tools/run_full_dem.py --region lima --max-memory-gb 6.0
```"""),
("code", """STORES = {name: os.path.abspath(os.path.join("..", "results", f"{name}_full"))
          for name in ("arequipa", "ancash", "lima")}


def load(region, label):
    path = os.path.join(STORES[region], f"{label}_results.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


runs = {r: {e: load(r, e) for e in ("grand", "tambo")} for r in STORES}
for r, table in runs.items():
    got = [e for e, v in table.items() if v]
    print(f"{r:<10} {', '.join(got) if got else 'NOT IN THE STORE YET'}")"""),
("md", """## The configurations, diffed

Same check as the Ancash notebook, extended to three. Anything that differs is either
bookkeeping or is called out."""),
("code", """def settings(path):
    with open(path) as f:
        return {k: v for k, v in json.load(f).items() if not k.startswith("_")}


for experiment in ("grand", "tambo"):
    base = settings(os.path.join("..", "config", f"{experiment}_arequipa_full.json"))
    print(f"{experiment.upper()}, against arequipa:")
    for region in ("ancash", "lima"):
        other = settings(os.path.join("..", "config",
                                      f"{experiment}_{region}_full.json"))
        differ = [k for k in sorted(set(base) | set(other))
                  if base.get(k) != other.get(k)]
        bookkeeping = {"dem_path", "region_name"}
        real = [k for k in differ if k not in bookkeeping]
        print(f"    {region:<9} {len(differ)} differ "
              f"({len(differ) - len(real)} bookkeeping)"
              + (f" -> {real}" if real else ""))"""),
("md", """`rfi_zones` is the only real difference, and it is Arequipa that is the odd one
out: it excludes five hand-curated circles around its own towns, and neither Ancash nor
Lima has an equivalent list. Holding both at `none` keeps *those two* exactly comparable
and leaves one stated caveat against Arequipa — its zones cover ~3,500 km² of a
~120,000 km² box, so about 2.9%."""),
("md", """## The terrain, before any search

The same tiled slope profile the Ancash notebook uses, over all three."""),
("code", """def slope_profile(tif, origin_lat, block=2048):
    import tifffile
    grid = ss.resolve_grid_geometry(tif, origin_lat)
    z = tifffile.imread(tif)
    hist = np.zeros(9001)
    bands = {"GRAND 3-25": (3.0, 25.0), "TAMBO 20-60": (20.0, 60.0)}
    counts = dict.fromkeys(bands, 0)
    total = 0
    for r0 in range(0, z.shape[0], block):
        r1 = min(z.shape[0], r0 + block)
        lo, hi = max(0, r0 - 1), min(z.shape[0], r1 + 1)
        dy, dx = np.gradient(z[lo:hi].astype(np.float32),
                             grid.cell_size_y, grid.cell_size_x)
        slope = np.degrees(np.arctan(np.hypot(dy, dx)))[r0 - lo:r1 - lo]
        s = slope[z[r0:r1] > 0]                    # the sea is exactly 0
        if not s.size:
            continue
        total += s.size
        hist += np.histogram(s, bins=9001, range=(0.0, 90.0))[0]
        for name, (a, b) in bands.items():
            counts[name] += int(((s >= a) & (s <= b)).sum())
    return {"Mpx": z.size / 1e6,
            "median slope": float((np.searchsorted(np.cumsum(hist), total / 2.0)
                                   + 0.5) / 100.0),
            **{k: 100.0 * v / total for k, v in counts.items()}}


DEMS = {"arequipa": ("arequipa_SRTMGL1.tif", -14.5553),
        "ancash": ("ancash_SRTMGL1.tif", -8.04958),
        "lima": ("lima_SRTMGL1.tif", -10.2283)}
terrain = {}
for name, (stem, lat) in DEMS.items():
    tif = os.path.abspath(os.path.join("..", "input", "dem", stem))
    if os.path.exists(tif):
        terrain[name] = slope_profile(tif, lat)

if terrain:
    cols = list(next(iter(terrain.values())))
    print(f"{'':<10}" + "".join(f"{c:>15}" for c in cols))
    for name, row in terrain.items():
        print(f"{name:<10}" + "".join(f"{row[c]:>15.1f}" for c in cols))"""),
("md", """## The three-way result"""),
("code", """def summary(results):
    if not results:
        return None
    r = results["results"]
    sites = r.get("sites") or []
    funnel = results.get("funnel") or {}
    strided = next((v for k, v in funnel.items()
                    if k.startswith("kept by stride")), None)
    accepted = funnel.get("directions accepted")
    return {"px": next(iter(funnel.values()), None),
            "sites": r.get("total_sites"),
            "area": sum(s.get("area_km2", 0.0) for s in sites),
            "capacity": r.get("total_capacity"),
            "acceptance": 100.0 * accepted / strided if accepted and strided else None}


base_px = (summary(runs["arequipa"]["grand"]) or {}).get("px")
for experiment in ("grand", "tambo"):
    print(f"=== {experiment.upper()} ===")
    print(f"{'region':<10}{'Mpx':>8}{'sites':>8}{'area km2':>12}"
          f"{'capacity':>11}{'accept':>9}{'area /px':>11}")
    ref = summary(runs["arequipa"][experiment])
    for region in ("arequipa", "ancash", "lima"):
        s = summary(runs[region][experiment])
        if not s:
            continue
        scale = s["px"] / base_px if base_px else None
        per = (f"{s['area'] / ref['area'] / scale:.2f}x"
               if ref and scale and ref["area"] else "-")
        print(f"{region:<10}{s['px'] / 1e6:>8.1f}{s['sites']:>8,}{s['area']:>12,.1f}"
              f"{s['capacity']:>11,}{s['acceptance']:>8.1f}%{per:>11}")
    print()"""),
("md", """## Both at once, three ways"""),
("code", """print(f"{'region':<10}{'joint km2':>12}{'jaccard':>11}{'share of TAMBO':>17}")
for region in ("arequipa", "ancash", "lima"):
    path = os.path.join(STORES[region], "combined_report.json")
    if not os.path.exists(path):
        continue
    with open(path) as f:
        report = json.load(f)
    overlap = report["pairwise_overlap"]["GRAND & TAMBO"]
    print(f"{region:<10}{report['joint']['area_km2']:>12,.1f}"
          f"{overlap['jaccard']:>11.5f}"
          f"{100 * overlap['fraction_of_TAMBO']:>16.1f}%")"""),
("md", """**The share of TAMBO's mask is the invariant to watch.** It held near 44% across
Arequipa and Ancash, whose terrain could hardly differ more. A third region is what turns
that from a coincidence into a property of the two experiments: the joint region is
TAMBO-limited, and co-location costs GRAND almost nothing.

The full cross-region table, regenerated from the stores whenever a region is added, is
in [`results/region_comparison.md`](../results/region_comparison.md) — produced by
`python tools/compare_regions.py`."""),
("md", """## The maps"""),
("code", SHOW_HELPER),
("code", """OUT = os.path.abspath(os.path.join("..", "output"))
for label, title in (("grand", "GRAND over Lima"), ("tambo", "TAMBO over Lima")):
    found = sorted(glob.glob(os.path.join(OUT, f"lima_full_{label}",
                                          "oroscope_results*.png")))
    show_figure(found[0] if found else os.path.join(OUT, "missing.png"), caption=title)"""),
("code", """for frame, caption in (("combined_overview_1_grand.png", "GRAND alone"),
                       ("combined_overview_2_tambo.png", "TAMBO alone"),
                       ("combined_overview_3_both.png", "Both — the joint ground")):
    show_figure(os.path.join(OUT, "lima_full_combined", frame), caption=caption)"""),
("md", """---

## Zooming in: Cajatambo and the upper Pativilca

Lima's TAMBO sites are not spread evenly. Seventeen of the forty, including the largest,
sit in the north-west corner around Cajatambo — **the same ground the Ancash box reaches
into from the other side**, where the largest joint patch of the whole Ancash run turned
out to be (notebook 10). So this crop does two jobs: it is Lima's densest TAMBO ground,
and running it from the Lima side checks that two independently downloaded, differently
aligned crops of one terrain agree.

8.3 Mpx, so it runs unbiased at `downsample_factor` 1 / `candidate_stride` 1:

```bash
cd src && oroscope-crop ../input/dem/lima_SRTMGL1.tif ../input/dem/cajatambo.tif \
    --north -10.30 --south -11.10 --west -77.60 --east -76.80
python tools/run_full_dem.py --region cajatambo --max-memory-gb 6.0
```"""),
("code", """CROP = os.path.abspath(os.path.join("..", "results", "cajatambo_full"))

print(f"{'':<34}{'sites':>8}{'area km2':>12}{'capacity':>12}")
for label in ("grand", "tambo"):
    for tag, store in ((f"{label.upper()}, Lima department (4/5)", STORES["lima"]),
                       (f"{label.upper()}, Cajatambo crop (1/1)", CROP)):
        path = os.path.join(store, f"{label}_results.json")
        if not os.path.exists(path):
            continue
        with open(path) as f:
            r = json.load(f)["results"]
        sites = r.get("sites") or []
        print(f"{tag:<34}{r['total_sites']:>8,}"
              f"{sum(x.get('area_km2', 0.0) for x in sites):>12,.1f}"
              f"{r['total_capacity']:>12,}")"""),
("md", """**The crop alone finds six times the TAMBO area of the entire Lima department
run, and eight times the detector positions — from 8% of its pixels.** That is not new
terrain; it is the same terrain measured without striding and downsampling, and it is the
291× of notebook 10 arriving independently on different ground.

It also revises something. The joint region's share of TAMBO's mask sat near 44% across
all three department runs, which looked like a constant of the two experiments."""),
("code", """print(f"{'':<28}{'sampling':>10}{'joint km2':>12}{'share of TAMBO':>17}")
for tag, store, samp in (("Arequipa", STORES["arequipa"], "4 / 5"),
                         ("Ancash", STORES["ancash"], "4 / 5"),
                         ("Lima", STORES["lima"], "4 / 5"),
                         ("Huaylas crop", os.path.abspath(os.path.join(
                             "..", "results", "huaylas_full")), "1 / 1"),
                         ("Cajatambo crop", CROP, "1 / 1")):
    path = os.path.join(store, "combined_report.json")
    if not os.path.exists(path):
        continue
    with open(path) as f:
        report = json.load(f)
    o = report["pairwise_overlap"]["GRAND & TAMBO"]
    print(f"{tag:<28}{samp:>10}{report['joint']['area_km2']:>12,.1f}"
          f"{100 * o['fraction_of_TAMBO']:>16.1f}%")"""),
("md", """**There are two constants, not one, and the strided value is not the real one.**
At 4 / 5 the share is ~44%; at 1 / 1 it is ~73%, on two separate crops. Striding
fragments TAMBO's mask and leaves GRAND's untouched, so what survives a strided run is
the scattered remainder — which overlaps GRAND's contiguous blob less than the intact
mask does.

So the honest statement is that **roughly three quarters of TAMBO-viable ground is also
GRAND-viable**, not four ninths. The *invariance* survives — at fixed sampling the share
barely moves across radically different terrain, which is the real content of the
TAMBO-limited finding — but the number to quote comes from the unbiased runs, and the two
rows must never be mixed."""),
("code", """for frame, caption in (
        ("combined_overview_1_grand.png", "GRAND over the Cajatambo crop"),
        ("combined_overview_2_tambo.png", "TAMBO at stride 1 — 97 sites"),
        ("combined_overview_3_both.png", "Both — 805 km², 71.9% of TAMBO's mask")):
    show_figure(os.path.abspath(os.path.join(
        "..", "output", "cajatambo_full_combined", frame)), caption=caption)"""),
("md", """---

## Co-locating the two arrays

The searches answer *where each experiment could go*. Deployment asks something else:
**given a site chosen for one, is there enough ground nearby for the other?**

The tempting way to ask is "how much of the joint mask can host GRAND", and it gives the
wrong answer. On the unbiased Cajatambo crop the GRAND-viable ground *inside* the joint
mask is 22,577 fragments of which exactly **one** is large enough for a single 1 km
lattice cell — while 976 km² of perfectly good GRAND ground lies within 20 km. An
optimiser pointed at the intersection would call that site impossible.

**A partner array does not have to stand on the joint mask.** What couples two arrays is
a shared line of sight to the same massif, not a shared footprint — and GRAND's own
targets are 10–40 km away, so an antenna 20 km from a TAMBO strip is watching the same
wall. `ce.colocation_capacity` measures that directly."""),
("code", """from oroscope import combine_experiments as ce

OUTDIRS = {r: os.path.abspath(os.path.join("..", "output", f"{r}_full_grand"))
           for r in ("arequipa", "ancash", "lima")}


def best_tambo_site(region):
    path = os.path.join(STORES[region], "tambo_results.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        sites = json.load(f)["results"]["sites"]
    return max(sites, key=lambda s: s["capacity_exact"]) if sites else None


print(f"{'region':<11}{'anchor TAMBO site':>19}{'GRAND within 10 km':>21}"
      f"{'within 20 km':>15}{'within 40 km':>15}")
for region in ("arequipa", "ancash", "lima"):
    site = best_tambo_site(region)
    if not site or not os.path.isdir(OUTDIRS[region]):
        continue
    run = ce.load_run(OUTDIRS[region])
    rows = ce.colocation_capacity(run["mask"], run["world"],
                                  site["center_lat"], site["center_lon"],
                                  radii_km=(10.0, 20.0, 40.0), spacing_km=1.0)
    caps = "".join(f"{r['capacity']:>15,}" for r in rows)
    print(f"{region:<11}{site['capacity_exact']:>13,} units{caps[:21]}{caps[21:]}")"""),
("md", """And the question a proposal actually asks — **how far would the partner array
have to spread?**"""),
("code", """WANTED = (100, 1000, 5000)
print(f"{'region':<11}" + "".join(f"{'r for ' + str(w) + ' GRAND':>20}" for w in WANTED))
for region in ("arequipa", "ancash", "lima"):
    site = best_tambo_site(region)
    if not site or not os.path.isdir(OUTDIRS[region]):
        continue
    run = ce.load_run(OUTDIRS[region])
    cells = []
    for want in WANTED:
        r = ce.smallest_radius_for(run["mask"], run["world"], site["center_lat"],
                                   site["center_lon"], want, spacing_km=1.0)
        cells.append(f"{r:.0f} km" if r else "> 80 km")
    print(f"{region:<11}" + "".join(f"{c:>20}" for c in cells))"""),
("md", """**The answer is consistent across all three regions**, which is what makes it
useful: about **10 km for 100 antennas, 20–30 km for 1,000, and 40–60 km for 5,000**.
None of the three is limited by finding GRAND ground near a TAMBO site — the limit is
simply that 1,000 antennas at 1 km spacing need 866 km², and steep country does not offer
that within a few kilometres.

Two caveats, and the first is the one to keep.

**The anchor is each region's largest TAMBO site, and at ``4 / 5`` those capacities are
badly understated** — 551, 726 and 990 units against 6,357 on the unbiased Cajatambo
crop. The *locations* survive striding even where the capacities do not, so the radii
above are sound while the "units" column is a floor.

**The capacity is an anchored lattice**, laid from the region's bounding-box corner
rather than fitted, so it is an estimate for an arbitrarily placed array. Fitting would
do better."""),
("md", """---

## The full explanation of each run

Reproduced in full rather than summarised: the caveats matter as much as the totals, and
they are otherwise sitting in files nobody opens."""),
("code", """def explanation(region, name):
    path = os.path.join(STORES[region], name)
    if not os.path.exists(path):
        print(f"not in the store: {name}")
        return
    with open(path) as f:
        print(f.read())


explanation("lima", "grand_explanation.txt")"""),
("code", """explanation("lima", "tambo_explanation.txt")"""),
("code", """explanation("lima", "combination_explanation.txt")"""),
("md", """## Where to go next

- **[9. Arequipa, the full DEM](09_arequipa_dem.ipynb)** and
  **[10. Ancash](10_ancash_dem.ipynb)** — the other two thirds of the comparison.
- **[12. Peru, all of it](12_peru_dem.ipynb)** — the whole country at 3 arc-seconds,
  and why a coarse answer has to be read differently.
- **[6. Combining and sensitivity](06_combining_and_sensitivity.ipynb)** — how much of
  any of this survives a change of assumption."""),
("md", footer(prev=("10_ancash_dem.ipynb", "Ancash"),
              nxt=("12_peru_dem.ipynb", "Peru, all of it"))),
]

# --------------------------------------------------------------------------- 12  peru
NB_PERU = [
("md", """# 12. Peru, all of it

Every other notebook here works on a crop, or on one department. This one is the search
run over **a whole country** — 22,080 × 15,360 pixels, 339 million of them, from the
Ecuadorian border to the Chilean one and from the Pacific to the Brazilian lowlands.

It is a **survey**, and the word is doing real work. The criteria are copied unchanged
from the full-Arequipa run so the two are comparable, but the grid is 3 arc-seconds
rather than 1, and almost everything worth knowing about the answer is a consequence of
that. This notebook is as much about reading a coarse result honestly as it is about
Peru."""),
("code", PREAMBLE_NO_PLOT + """import json
import os

from oroscope import site_searcher as ss"""),
("md", """## Read, not run

The cells below open results produced locally and stored in
`output/grand_peru_survey/`. They do not start a search — it needs a 302 MB DEM that is
not in the repository, and four minutes of eight cores.

To produce them yourself:

```bash
# 1. A free OpenTopography key: register at portal.opentopography.org/myopentopo,
#    then open "myOpenTopo Authorizations and API Key" and copy it.
export OPENTOPOGRAPHY_API_KEY=...

# 2. The DEM (302 MB) and a template config
cd src && oroscope-fetch-dem --region peru

# 3. The survey, about four minutes
oroscope --config_path ../config/grand_peru_survey.json
```

Everything below degrades gracefully if the store is absent: the figures are saved in
this notebook, so it reads either way."""),
("code", """STORE = os.path.abspath(os.path.join("..", "output", "grand_peru_survey"))
RESULTS = os.path.join(STORE, "oroscope_results_peru_SRTMGL3.json")

run = json.load(open(RESULTS)) if os.path.exists(RESULTS) else None
print("stored run found" if run else f"not here yet: {os.path.relpath(RESULTS)}")"""),
("md", """## Why 3 arc-seconds, and why that is not a preference

Two independent limits force it, and it is worth seeing both because they bracket what
is possible on a desktop.

**The API.** OpenTopography caps a single request by area, per dataset: 450,000 km² for
every 30 m dataset and 4,050,000 km² for the 90 m ones. Peru's bounding box is about
2.86 million km² — comfortably inside the 90 m limit and six times over the 30 m one.
A 1 arc-second national DEM is not one download.

**Memory.** `estimate_peak_memory_gb` answers this directly, and the shape of its answer
is the useful part."""),
("code", """rows, cols = 22080, 15360            # Peru at 3 arc-seconds
print(f"{rows:,} x {cols:,} = {rows*cols/1e6:,.0f} Mpx\\n")
print("estimated peak anonymous memory, GiB")
print(f"{'':>8}" + "".join(f"{'stride ' + str(s):>12}" for s in (10, 15, 30)))
for ds in (1, 2, 4):
    row = "".join(f"{ss.estimate_peak_memory_gb(rows, cols, ds, s):>12.2f}"
                  for s in (10, 15, 30))
    print(f"  ds {ds:<4}" + row)"""),
("md", """Read down a column rather than across a row. **`candidate_stride` is the memory lever,
not `downsample_factor`** — candidates are taken on the *native* grid, so downsampling
scales the labelling arrays as its inverse square but barely touches the dominant term.
That is the single most useful fact for planning a large run, and it is the opposite of
what most people reach for first.

At 1 arc-second the same country is 3,052 Mpx, nine times this table. The chosen cell is
`downsample_factor` 4 with `candidate_stride` 15: **4.77 GiB**, against roughly 8 GiB
available on the machine this ran on."""),
("md", """## Raising the stride costs area unless the closing element keeps up

This is the trap that made a published TAMBO area 4.75× too low
([notebook 7](07_animating_the_mechanism.ipynb) animates it as `stride_and_closing`).
Striding marks one surviving pixel in N; the mask is then closed morphologically before
area is measured. If the closing element cannot bridge the gap the stride leaves, the
mask never reconnects — it stays a scatter of isolated pixels, small regions fall under
the size thresholds, and the area collapses.

At 3 arc-seconds the pixel is ~92 m, so the gap grows with it."""),
("code", """grid = ss.resolve_grid_geometry("no-such-file.tif", -9.2, cell_size_deg=3/3600)
print(f"pixel: {grid.cell_size_y:.1f} m N-S, {grid.cell_size_x:.1f} m E-W\\n")

for stride in (10, 15, 30):
    gap = ss.stride_gap_m(stride, grid.cell_size_y)
    print(f"stride {stride:>2}  gap {gap:>6,.0f} m")
    for element_km in (1.0, 1.5, 3.0):
        verdict = ss.warn_stride_outruns_closing(stride, grid.cell_size_y,
                                                 element_km, 1.0, quiet=True)
        state = "ok" if verdict is None else f"WARNS, {verdict['ratio']:.2f}x short"
        print(f"    closing {element_km:>4} km : {state}")"""),
("md", """So stride 15 needs an element of at least ~1.4 km. The run uses **1.5 km**, which is
1.5× GRAND's own 1 km antenna spacing.

That is a declared bias, and its direction is the opposite of the TAMBO failure: an
element larger than the array's own scale will bridge some gaps that are *real*, so the
area comes out slightly high rather than catastrophically low. For a survey that is the
right way to be wrong, and saying so is the price of using it."""),
("md", """## What the run found"""),
("code", """if run:
    t = run["timings_sec"]
    for stage, seconds in t.items():
        print(f"  {stage:<22} {seconds:>7.1f} s")
    print(f"  {'TOTAL':<22} {sum(t.values()):>7.1f} s")"""),
("md", """**Four minutes for a country**, against 26.8 minutes for the Arequipa DEM at 1
arc-second. Peru is 2.6× the pixels but a ninth of the candidates, and the candidates
are what the ray tracing costs — which is the same lesson the memory table taught,
arriving from the other direction."""),
("code", """if run:
    funnel = run["funnel"]
    first = next(iter(funnel.values()))
    prev = None
    for stage, n in funnel.items():
        step = f"{100*n/prev:6.1f}% of previous" if prev else ""
        print(f"  {stage:<34} {n:>13,}  {step}")
        prev = n"""),
("md", """Two rows deserve attention.

**`directions accepted` over `kept by stride 15` is the acceptance rate: 43.1%.** The
full Arequipa DEM gave 61.6% for the same criteria. Lower is the expected direction —
the national box adds coastal desert below the 3° slope floor, high Andes above the 25°
ceiling, and Amazon lowlands that are simply flat.

**`after gap closing` is larger than `directions accepted` by a factor of 24.** Closing
is not a filter; it is the step that undoes striding, and most of that factor is the 15
pixels each kept candidate stands for."""),
("code", """if run:
    r = run["results"]
    sites = sorted(r["sites"], key=lambda s: -s["area_km2"])
    print(f"{r['total_sites']} sites, "
          f"{sum(s['area_km2'] for s in sites):,.0f} km2, "
          f"{r['total_capacity']:,} antenna positions at 1 km\\n")
    print(f"{'area km2':>12} {'capacity':>10}   centre")
    for s in sites[:6]:
        print(f"{s['area_km2']:>12,.0f} {s['capacity_exact']:>10,}   "
              f"{s['center_lat']:7.2f}, {s['center_lon']:8.2f}")
    print(f"{'...':>12}   and {len(sites)-6} smaller")"""),
("md", """## The map"""),
("code", SHOW_HELPER),
("code", """show_figure(os.path.join(STORE, "oroscope_results_peru_SRTMGL3.png"),
            caption="GRAND over the whole of Peru, 3 arc-seconds")"""),
("md", """The Pacific is the flat ground bottom-left, the Amazon basin the low blue to the
north-east, and the accepted terrain follows the cordillera between them — which is what
it should do, and is the first thing to check on any result of this size."""),
("md", """## Now read it honestly

Three caveats, in descending order of how much they should change what you say.

### 1. The area is a bracket, not a number

The stride-corrected accepted set and the reported area disagree by 38%, and both are
approximations."""),
("code", """if run:
    px_km2 = grid.cell_size_y * grid.cell_size_x / 1e6
    accepted = run["funnel"]["directions accepted"]
    selected = run["funnel"]["pixels in selected sites (est.)"]
    stride_corrected = accepted * 15 * px_km2
    reported = selected * px_km2
    print(f"  accepted, stride-corrected : {stride_corrected:>10,.0f} km2")
    print(f"  reported after closing     : {reported:>10,.0f} km2")
    print(f"  ratio                      : {reported/stride_corrected:>10.2f}")
    print(f"\\n  quote this as 4-6 x 10^5 km2, or {100*reported/1_285_216:.0f}% of Peru "
          f"with the bracket attached")"""),
("md", """Some of that factor is closing doing its job — filling holes *inside* accepted ground,
which the stride correction also does — and some is the 1.5 km element bridging ground
that was genuinely rejected. Nothing in this run separates them.

### 2. "17 sites" is the number to distrust — more than the area

A site is a connected component of the closed mask. At this scale, with a 1.5 km
element applied to a strided scatter across 339 Mpx, components merge."""),
("code", """if run:
    big = max(run["results"]["sites"], key=lambda s: s["area_km2"])
    b = big["bounds"]
    print(f"largest site: {big['area_km2']:,.0f} km2 "
          f"({100*big['area_km2']/sum(s['area_km2'] for s in run['results']['sites']):.1f}%"
          f" of the total)")
    print(f"  its bounding box: {b['north']:.2f} to {b['south']:.2f} lat, "
          f"{b['west']:.2f} to {b['east']:.2f} lon")
    print("  the DEM's own:    0.00 to -18.40 lat, -81.40 to -68.60 lon")
    print(f"\\n  mean altitude of its accepted candidates: "
          f"{big['arrival_scan']['altitude_m_mean']:,.0f} m")"""),
("md", """**The largest site's bounding box is the entire DEM.** One connected component holds
94.8% of the area, and `min_width_km: 2.0` had no reason to break it up. Its accepted
candidates are Andean — mean altitude 2,446 m — while the polygon enclosing them reaches
the coast and the basin.

So the *candidate* statistics inside a site are meaningful and its *extent* is not. Site
count and site geometry from a run at this stride are artefacts of the closing element.

### 3. The one that was checked and came back fine

The obvious worry about 3 arc-seconds is that it moves the slope screen itself —
smoothing steep ground down into the 3–25° band and roughening flat ground up into it —
which would make the whole result an artefact of resampling.

That is testable directly: take real terrain at its native resolution, block-mean it to
three times the pixel, and compare."""),
("code", """def band_fraction(z, cell_y, cell_x, lo=3.0, hi=25.0):
    \"\"\"Share of a DEM inside a slope band, and its slope quartiles.\"\"\"
    dy, dx = np.gradient(np.asarray(z, dtype=np.float64), cell_y, cell_x)
    slope = np.degrees(np.arctan(np.hypot(dy, dx)))
    return 100.0 * ((slope >= lo) & (slope <= hi)).mean(), np.percentile(slope, [25, 50, 75])


AREQUIPA = os.path.abspath(os.path.join("..", "input", "dem", "arequipa_SRTMGL1.npy"))
if os.path.exists(AREQUIPA):
    z30 = np.asarray(np.load(AREQUIPA, mmap_mode="r")[3000:7000, 4000:9000], dtype=np.float64)
    # The DEM's own geometry, not a nominal arc-second, so this agrees to the digit
    # with the same measurement quoted in ROADMAP 6.46.
    g = ss.resolve_grid_geometry(AREQUIPA.replace(".npy", ".tif"), -14.5553)
    source = "20 Mpx of the Arequipa DEM"
else:
    g = ss.resolve_grid_geometry("no-such-file.tif", -14.56, cell_size_deg=1/3600)
    # No DEM here. Rough terrain built in place rather than imported, so this cell runs
    # on a bare clone -- the numbers will not be Peru's, but the comparison is the point.
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[0:1200, 0:1200] * g.cell_size_x
    z30 = (900.0 * np.sin(xx / 7000.0) * np.cos(yy / 9000.0)
           + 250.0 * np.sin(xx / 1300.0) + rng.normal(0.0, 6.0, (1200, 1200)))
    source = "synthetic terrain (no DEM present)"

h, w = (z30.shape[0] // 3) * 3, (z30.shape[1] // 3) * 3
z90 = z30[:h, :w].reshape(h // 3, 3, w // 3, 3).mean(axis=(1, 3))

f30, q30 = band_fraction(z30, g.cell_size_y, g.cell_size_x)
f90, q90 = band_fraction(z90, g.cell_size_y * 3, g.cell_size_x * 3)
print(f"{source}\\n")
print(f"{'grid':<8}{'in 3-25 deg':>14}   slope quartiles")
print(f"{'30 m':<8}{f30:>13.1f}%   {np.round(q30, 1)}")
print(f"{'90 m':<8}{f90:>13.1f}%   {np.round(q90, 1)}")
print(f"\\nband moves by {f90-f30:+.1f} points")"""),
("md", """On the real DEM the quartiles do shift down — smoothing must make them — but the *band
fraction* moves by 0.2 points: 67.6% at 30 m against 67.4% at 90 m. **What is lost at the
25° ceiling is regained at the 3° floor.**

So the screen's reach is a fact about Peru rather than about the grid. The prediction
here was that resampling would dominate, and it was wrong; that is the useful part, and
it is why the check is in the notebook rather than the assumption."""),
("md", """## Why there is no TAMBO counterpart

`config/` holds `grand_peru_survey.json` and deliberately no TAMBO equivalent. Two
reasons, both about the grid rather than about Peru.

**The geometry is below the pixel.** Colca's floor is about 1 km wide and its rim-to-rim
separation about 4.5 km. At 92 m a canyon is roughly eleven pixels across the floor, and
the wall the array stands on — where the 20–60° slope band is measured — is a handful.
The grid has averaged away the thing being screened for.

**The closing element cannot be made to work.** TAMBO's 100 m antenna spacing sets a
100 m element against this run's 1,382 m stride gap."""),
("code", """gap = ss.stride_gap_m(15, grid.cell_size_y)
print(f"stride 15 leaves a {gap:,.0f} m gap\\n")

for label, element_km in (("TAMBO's own, 100 m", 0.1), ("raised to 1.5 km", 1.5)):
    v = ss.warn_stride_outruns_closing(15, grid.cell_size_y, element_km, 0.1, quiet=True)
    state = "ok" if v is None else f"WARNS, {v['ratio']:.1f}x short"
    print(f"  {label:<20} {state}")

print(f"\\n...but 1.5 km is {1500/100:.0f}x TAMBO's own array scale, so a mask closed")
print("with it smears across the canyon instead of tracing the wall.")"""),
("md", """Wrong in the opposite direction from the 4.75× under-report, and no less wrong. **A
national TAMBO answer needs 1 arc-second and tiling**, and that is a different job.

## One more thing this run taught, the hard way

`--max_memory_gb` is applied as `RLIMIT_AS`, which caps **virtual** address space — so it
counts every mapping, including the memory-mapped DEM and the scratch buffers. But
`estimate_peak_memory_gb` estimates **anonymous** memory and deliberately excludes them,
because the kernel can evict a file-backed page.

On a 339 Mpx DEM those two quantities differ by more than 2 GiB. The first attempt set
the cap from the estimate, ran the entire search successfully, and then died on the
*map* with `Unable to allocate 40.8 MiB` — with the JSON and the GeoTIFF already
written, so only the picture was lost.

The obvious fix is to raise the cap. Raising it too far is worse than not setting it:
**a cap above available memory is not a cap**, because the OOM killer arrives before the
limit can fire. A later run set it at 13 GiB on a machine with 8.7 GiB available and
took the machine down.

`preflight_memory` now reports both numbers and says when they cannot be reconciled."""),
("code", """report = ss.preflight_memory("no-such-file.tif", downsample_factor=4,
                             candidate_stride=15, max_memory_gb=0)
print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in report.items()})
print("\\n(max_memory_gb=0 disables capping, which is what makes this safe to run here)")"""),
("md", """When the cap cannot both clear the mapped DEM and stay under available memory, **the
configuration does not fit** — and the answer is a coarser search, raising
`candidate_stride`, rather than a bigger number in the cap.

## Where to go next

- **[9. Arequipa, the full DEM](09_arequipa_dem.ipynb)** — the same machine at 1
  arc-second, where the areas are quotable rather than bracketed.
- **[7. Animating the mechanism](07_animating_the_mechanism.ipynb)** —
  `stride_and_closing` is the trap in this notebook, animated.
- **[6. Combining and sensitivity](06_combining_and_sensitivity.ipynb)** — what a
  criterion costs, which is the sweep this survey never ran."""),
("md", footer(prev=("11_lima_dem.ipynb", "Lima"),
              nxt=("13_turning_the_knobs.ipynb", "Turning the knobs"))),
]

# Reading order, which is also the numbering. The regional runs are 09 onward and are
# numbered by region rather than by the order they happened to be written, so a new one
# slots in without renumbering its neighbours -- 11 is reserved for Lima.
# --------------------------------------------------------------------------- 13  knobs
NB_KNOBS = [
("md", """# 13. Turning the knobs

[Notebook 6](06_combining_and_sensitivity.ipynb) says a result is only as firm as its
assumptions. This one **turns each knob and measures what moves**, on terrain small
enough that every cell runs in seconds.

The magnitudes here are not the ones to quote — synthetic terrain is not Peru, and the
real numbers are in notebooks 9 to 12. What transfers is the **shape** of each response:
which knobs are gentle, which sit on cliffs, and which change one reported number while
leaving another alone."""),
("code", PREAMBLE + """import contextlib
import io
import time

from oroscope import arrival_scan, scoring
from oroscope import site_searcher as ss
from scipy.ndimage import binary_closing, label"""),
("md", """## A canyon to work on

Built in place rather than imported, so this notebook runs on a bare clone. Colca's
published geometry: ~1.5 km deep, ~4.5 km rim to rim, which implies ~40.6° walls."""),
("code", """CELL_DEG = 1 / 3600.0
grid = ss.resolve_grid_geometry("no-such-file.tif", -15.6, cell_size_deg=CELL_DEG)

n = 700
cols = np.arange(n, dtype=np.float64)[None, :].repeat(n, 0)
x = cols * grid.cell_size_x
centre = (n * grid.cell_size_x) / 2.0
from_edge = np.abs(x - centre) - 500.0                 # 1 km floor
rise = np.clip(from_edge, 0.0, None) * np.tan(np.radians(40.6))
z = (3500.0 - 1500.0 + np.clip(rise, 0.0, 1500.0)).astype(np.float32)
# A gentle along-axis tilt, so the canyon is not perfectly translation-invariant
z = z + (np.arange(n, dtype=np.float32)[:, None] * 0.4)

print(f"{n} x {n} px, pixel {grid.cell_size_y:.1f} x {grid.cell_size_x:.1f} m")
print(f"elevation {z.min():,.0f} - {z.max():,.0f} m")

plt.figure(figsize=(4.2, 3.6))
plt.imshow(z, cmap="Greys_r")
plt.colorbar(label="Elevation (m)")
plt.xticks([])
plt.yticks([])
plt.tight_layout()"""),
("md", """## The funnel, live

Every run records how many pixels survived each stage. This is the same table
`--explain` prints, on this canyon, for TAMBO's criteria."""),
("code", """TAMBO = dict(min_slope_deg=20.0, max_slope_deg=60.0,
             min_dist_km=2.0, max_dist_km=5.0,
             elev_min_deg=-20.0, elev_max_deg=20.0, n_elev_bins=20,
             n_azimuths=9, half_width_deg=60.0,
             min_target_slope_deg=25.0, max_range_m=5000.0)


def screen(stride):
    \"\"\"Stage 1 and 2: the slope band, then striding.

    The screen's progress bar goes to stderr, where a notebook stores it as output. It
    reports a step that takes a fraction of a second here, so it is carrying nothing.
    \"\"\"
    with contextlib.redirect_stderr(io.StringIO()):
        return ss.get_candidates_chunked(
            z, grid, None, -15.6, -72.4,
            min_slope_deg=TAMBO["min_slope_deg"],
            max_slope_deg=TAMBO["max_slope_deg"], candidate_stride=stride)


def scan(cand, **over):
    \"\"\"Stage 3: the arrival scan.\"\"\"
    kw = {k: v for k, v in TAMBO.items()
          if k not in ("min_slope_deg", "max_slope_deg")}
    kw.update(over)
    return arrival_scan.scan(cand, z, grid, **kw)


cand = screen(1)
obs = scan(cand)
accepted = int((obs["cells"] > 0).sum())
print(f"{'DEM pixels':<28}{z.size:>12,}")
print(f"{'slope 20-60 deg, stride 1':<28}{len(cand):>12,}"
      f"   {100*len(cand)/z.size:5.1f}% of the DEM")
print(f"{'directions accepted':<28}{accepted:>12,}"
      f"   {100*accepted/len(cand):5.1f}% of candidates")"""),
("md", """---

## Knob 1: `candidate_stride`

Keeps one surviving pixel in N. **Cost control, not a criterion** — so the thing to
check is whether it biases the answer. Acceptance says no."""),
("code", """print(f"{'stride':>7}{'candidates':>13}{'accepted':>11}{'acceptance':>13}{'scan time':>12}")
for stride in (1, 2, 5, 10):
    c = screen(stride)
    t = time.perf_counter()
    o = scan(c)
    dt = time.perf_counter() - t
    a = int((o["cells"] > 0).sum())
    print(f"{stride:>7}{len(c):>13,}{a:>11,}{100*a/len(c):>12.1f}%{dt:>11.2f}s")"""),
("md", """**Acceptance barely moves; the cost falls as 1/N.** That is the whole case for
striding, and it is why the funnel labels it "not a constraint".

The catch is not here. It is that striding leaves a **gap**, and what happens to that gap
decides the reported area."""),
("md", """## Knob 2: `gap_close_km`, and why the penalty depends on how thin the ground is

Striding leaves a gap. Closing repairs it — but how *well* depends on something the
parameter does not mention: **the width of the accepted feature**. A wide blob survives a
marginal element; a narrow strip does not.

That is the difference between a 4.75× under-report at Colca and a 291× one on the
Callejón de Huaylas, and it is worth seeing on its own, so here it is on synthetic
strips of controlled width rather than on the canyon above."""),
("code", """def recovery(width_px, stride=5, elements=(3, 5, 9)):
    \"\"\"How much of a strip of a given width survives striding and closing.\"\"\"
    m = 400
    rr, cc = np.mgrid[0:m, 0:m]
    truth = np.abs(cc - 200 - 40 * np.sin(rr / 50.0)) < width_px / 2
    strided = np.zeros_like(truth)
    strided[::stride, ::stride] = truth[::stride, ::stride]
    out = {}
    for k in elements:
        closed = binary_closing(strided, np.ones((k, k)))
        out[k] = (closed.sum() / truth.sum(), label(closed)[1])
    return out


print("Fraction of the accepted strip recovered, after stride 5 and closing:")
print(f"{'strip width':>12}{'element 3 px':>14}{'element 5 px':>14}{'element 9 px':>14}")
for width in (3, 6, 12, 30, 80):
    r = recovery(width)
    print(f"{width:>10} px" + "".join(f"{r[k][0]:>13.2f}x" for k in (3, 5, 9)))"""),
("md", """Two things, and the second is the one that is easy to miss.

**An element below the gap never recovers anything** — 0.04× at every width, because the
marks simply never touch. That is the failure the warning catches.

**An element at or above the gap recovers a fraction that depends on the width.** At 5 px
against a 5 px gap: 0.10× for a 3-pixel strip, 0.91× for an 80-pixel one. So "the element
outruns the gap" is necessary and **not sufficient**. A canyon wall a few pixels across
loses most of itself to a nominally adequate element, and that is exactly the terrain
TAMBO selects.

And the loss does not show up where you would look for it — in the area — but in the
**region count**, because pruning then discards what fragmented:"""),
("code", """print(f"{'strip width':>12}{'regions after closing at 5 px':>32}")
for width in (3, 6, 12, 30, 80):
    print(f"{width:>10} px{recovery(width)[5][1]:>32,}")"""),
("md", """A strip that should be one region becomes dozens or hundreds, and
`min_sub_array_size` then throws away every piece too small to hold an array. **That is
where the factor of 291 went** — not into the area measurement.

The warning that fires on the first condition:"""),
("code", """for stride in (1, 5, 15):
    verdict = ss.warn_stride_outruns_closing(stride, grid.cell_size_y,
                                             gap_close_km=0.1,      # TAMBO's 100 m
                                             antenna_spacing_km=0.1, quiet=True)
    gap = ss.stride_gap_m(stride, grid.cell_size_y)
    state = "ok" if verdict is None else f"WARNS, {verdict['ratio']:.2f}x short"
    print(f"stride {stride:>2}: gap {gap:>6.0f} m against a 100 m element -> {state}")"""),
("md", """---

## Knob 3: `min_score`, the dominant assumption"""),
("code", """cfg = dict(solid_angle_half_sr=0.8, grammage_mode="particle",
           grammage_band_gcm2=(236.0, 1287.0), shower_development_m=0.0,
           decay_energy_min_pev=3.0, decay_energy_max_pev=1000.0,
           decay_spectral_index=2.0, composition="product")
obs["altitude_m"] = z[cand[:, 0].astype(int), cand[:, 1].astype(int)].astype(float)
total, parts = scoring.score_candidates(obs, cfg, distance_window_m=(2000.0, 5000.0))
viable = obs["cells"] > 0

print(f"components: {', '.join(parts)}")
print(f"viable candidates: {viable.sum():,}")
print()
print(f"{'min_score':>10}{'kept':>10}{'of viable':>12}")
for cut in (0.0, 0.1, 0.2, 0.35, 0.5, 0.7):
    kept = int((total[viable] >= cut).sum())
    print(f"{cut:>10.2f}{kept:>10,}{100*kept/viable.sum():>11.1f}%")"""),
("md", """On this terrain the fall is gentle until 0.5 and then steep — but **where the
steep part sits is not a property of the terrain**. It moves when a component is added,
because a product of numbers in [0, 1] concentrates near zero and every extra factor
pushes the whole population down. On the real Colca run the same cut at 0.0, 0.35 and 0.5
gave 45,928, 2,056 and **zero** detector positions.

`score_percentile` asks for a rank instead, which is scale-free. It is not the default
only because every published number here used `min_score`."""),
("md", """## Knob 4: `min_target_slope_deg`, which describes the far wall"""),
("code", """print(f"{'min_target_slope':>18}{'accepted':>11}{'of candidates':>15}")
for floor in (None, 15.0, 25.0, 35.0, 38.0, 42.0, 45.0):
    o = scan(cand, min_target_slope_deg=floor)
    a = int((o["cells"] > 0).sum())
    label_ = "unset" if floor is None else f"{floor:.0f} deg"
    print(f"{label_:>18}{a:>11,}{100*a/len(cand):>14.1f}%")"""),
("md", """These walls were built at 40.6°, and the criterion finds them: everything below
that passes, everything above fails, with the step in between. On real terrain the wall
slopes are a *distribution* rather than one value, so the step becomes a slope — and it
sits further right than the median, because the criterion is applied per direction while
the reported `target_slope_deg` is a mean over each candidate's accepted directions. See
[notebook 10](10_ancash_dem.ipynb).

The important point is that this is the **far** wall, the one the tau exits.
`min_slope_deg`/`max_slope_deg` describe the **near** ground the array stands on. A single
slope band cannot express both, which is why there are two."""),
("md", """## Knob 5: `downsample_factor`, which moves one number and not the other"""),
("code", """def accepted_mask(stride=1, element_px=3):
    \"\"\"The accepted set at a stride, closed with an element of a given size.\"\"\"
    c = screen(stride)
    o = scan(c)
    ok = o["cells"] > 0
    mask = np.zeros(z.shape, dtype=bool)
    mask[c[ok, 0].astype(int), c[ok, 1].astype(int)] = True
    return binary_closing(mask, np.ones((element_px, element_px)))


m = accepted_mask()
px_km2 = grid.cell_size_y * grid.cell_size_x / 1e6
full_km2 = m.sum() * px_km2
cap = ss.count_grid_capacity(np.ascontiguousarray(m), grid.cell_size_y,
                             grid.cell_size_x, 100.0, 1)
print(f"{'downsample':>11}{'area km2':>11}{'vs full':>9}{'capacity':>11}")
for ds in (1, 2, 4, 8):
    area = m[::ds, ::ds].sum() * px_km2 * ds * ds
    print(f"{ds:>11}{area:>11.2f}{area/full_km2:>8.2f}x{cap:>11,}")"""),
("md", """**Area is measured on the downsampled mask; capacity is counted at full
resolution.** On a feature this wide the area survives downsampling almost exactly — which
is the honest result here, and worth stating, because the ~30% loss quoted for real canyon
strips is a statement about *thin* features and about the interaction with pruning, not
about subsampling arithmetic.

The trap is the one that survives regardless: **the ratio of capacity to area is not a
property of the terrain**, because the two were measured on different grids. And
`candidate_stride` is the memory lever, not this one — candidates are taken on the native
grid, so downsampling scales the labelling arrays as its inverse square and barely touches
the dominant term."""),
("md", """---

## What to take away

| Knob | Response | Watch for |
| --- | --- | --- |
| `candidate_stride` | acceptance flat, cost 1/N | the gap it leaves |
| `gap_close_km` | **cliff at the stride gap** | areas collapsing, region counts exploding |
| `min_score` | **cliff**, moving with component count | prefer `score_percentile` |
| `min_target_slope_deg` | sharp at the true wall slope | it is the *far* wall, not the near one |
| `downsample_factor` | area falls, capacity does not | comparing the two across grids |

Two of the five sit on cliffs, and a third moves one reported number without moving
another. **Read this notebook before quoting any single number from a run** — and read
[notebook 6](06_combining_and_sensitivity.ipynb) for what a sensitivity sweep does with
the whole pipeline rather than one stage."""),
("md", """## Where to go next

- **[7. Animating the mechanism](07_animating_the_mechanism.ipynb)** — several of the
  responses above, as films.
- **[6. Combining and sensitivity](06_combining_and_sensitivity.ipynb)** — sweeps over
  the whole pipeline.
- **[10. Ancash](10_ancash_dem.ipynb)** — where the striding cliff was measured at 291×
  on real ground."""),
("md", footer(prev=("12_peru_dem.ipynb", "Peru, all of it"))),
]

NOTEBOOKS = {
    "01_getting_started.ipynb": NB01,
    "02_the_arrival_scan.ipynb": NB02,
    "03_physics_toolkit.ipynb": NB03,
    "04_criteria_and_scoring.ipynb": NB04,
    "05_grand_and_tambo.ipynb": NB05,
    "06_combining_and_sensitivity.ipynb": NB06,
    "07_animating_the_mechanism.ipynb": NB_ANIMATIONS,
    "08_explaining_a_run.ipynb": NB_EXPLAIN,
    "09_arequipa_dem.ipynb": NB_AREQUIPA,
    "10_ancash_dem.ipynb": NB_ANCASH,
    "11_lima_dem.ipynb": NB_LIMA,
    "12_peru_dem.ipynb": NB_PERU,
    "13_turning_the_knobs.ipynb": NB_KNOBS,
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
