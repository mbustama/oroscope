#!/usr/bin/env python
"""
Animations of the things in this search that a still picture explains badly.

Eight, chosen because each one shows a *process* whose intermediate states are the
point. Anything that is a single state is already a figure in ``oroscope.figures``, and
turning it into a movie would be decoration rather than explanation.

``the_walk``
    One backward ray sweeping down through the elevation window over real terrain: the
    first intersection sliding along the profile, and the column depth accumulating
    behind it. This is the mechanism the whole search rests on and the hardest thing in
    the project to convey in prose.
``the_azimuth_fan``
    The other half of that mechanism. ``the_walk`` sweeps elevation at one bearing;
    this sweeps the bearing and reports what each one finds -- a wall at the right
    range, the candidate's own hillside a few hundred metres away, or sky. Acceptance
    against bearing is a polar quantity, which is what a still picture renders badly.
``the_funnel``
    The map draining stage by stage -- slope, stride, directions accepted, closing,
    pruning -- with the surviving count. The funnel table says where the candidates
    went; this shows *where on the ground* they went.
``stride_and_closing``
    Why TAMBO's area is 4.75x low. A strided mask closed with an element that bridges
    the gaps, and the same mask closed with one that does not. The measurement is in
    ROADMAP 6.34; this is what it looks like.
``product_collapse``
    Why a threshold on a product is treacherous, and the reason is dynamic: every
    component multiplied in drags the whole population toward zero while the cut stays
    where it was put. Six real components, one real cut, and the fraction surviving
    read off each time.
``slope_criterion``
    The sensitivity sweeps say *how much* a criterion costs. This says *where*. The
    accepted mask over Colca as ``min_target_slope_deg`` climbs through the wall-slope
    distribution, and the rims letting go.
``tau_in_rock``
    More rock is not better -- the commonest misconception in the problem. A tau's
    energy and survival falling as it burrows, against the depth that maximises
    production and escape together (5.7e6 g/cm^2 at 1 EeV).
``energy_window``
    The arrival window narrowing as energy rises, its lower edge climbing from -4.4
    degrees at 100 PeV to -0.9 at 10 EeV. A falsifiable prediction, animated over the
    quantity it is a prediction about.

Everything is built from committed code. Six of the eight use synthetic terrain and so
reproduce on any clone; ``the_azimuth_fan``, ``product_collapse`` and
``slope_criterion`` are about what happens on *real* ground, so they read the Colca DEM
when it is present and fall back to synthetic terrain -- saying which they used -- when
it is not. See :func:`_colca_ground`.

    python tools/make_animations.py                  # all eight, MP4 then GIF
    python tools/make_animations.py --only the_walk
    python tools/make_animations.py --format gif --out docs/source/_static

MP4 needs ffmpeg; GIF falls back to pillow, which is always available. Outputs land in
``output/animations/`` by default, which is gitignored -- pass ``--out`` to place the
small ones somewhere they can be committed.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")                                # a tool, not a library: see trap 3
import matplotlib.animation as animation             # noqa: E402
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
from scipy.ndimage import (                          # noqa: E402
    binary_closing, binary_dilation, label,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "tests"))

from oroscope import arrival_scan, physics, scoring   # noqa: E402
from oroscope import site_searcher as ss              # noqa: E402
from oroscope.figures import (                       # noqa: E402
    DETECTOR, INK, MUTED, ROCK_EDGE, ROCK_FILL, RULE, SEQUENCE, WINDOW,
    _styled, _tidy, _minus,
)
import synthetic                                     # noqa: E402

FPS = 12

# Colca as the two published configurations describe it. Quoted here rather than read
# from config/*.json so an animation cannot silently change when a configuration does:
# these are illustrations of a mechanism, not a reproduction of a run.
COLCA_ORIGIN_LAT, COLCA_ORIGIN_LON = -15.299861, -72.400139
TAMBO_SLOPE_BAND = (20.0, 60.0)                      # the near wall the array stands on
TAMBO_SCAN = dict(n_azimuths=9, half_width_deg=60.0,
                  elev_min_deg=-20.0, elev_max_deg=20.0, n_elev_bins=20,
                  min_dist_km=2.0, max_dist_km=5.0, max_range_m=5000.0)
TAMBO_MIN_TARGET_SLOPE = 25.0
TAMBO_MIN_SCORE = 0.35

ACCEPT, TOO_NEAR, SKY = "#0F6B54", "#B02A25", "#9AA0A6"


def _colca_ground(fraction=0.4):
    """
    Colca: the real DEM when it is present, synthetic terrain when it is not.

    Three of these animations are about where a criterion bites on *actual* ground,
    which synthetic terrain cannot honestly show. The DEM is gitignored, so they fall
    back to the canyon :func:`synthetic.colca_like` builds -- Colca's published depth
    and rim-to-rim separation, and nothing else about it -- and say which they used.
    The fallback's walls are all at exactly 40.6 degrees, so a criterion sweeping
    across that value flips rather than erodes; the real DEM has a distribution.

    Only the middle ``fraction`` of each axis is taken. The canyon runs through it, and
    the slope sweep runs a full arrival scan per frame: over all 6.1 Mpx that is about
    a minute a frame rather than half a second.

    Returns
    -------
    elevation : ndarray
    map_grid : MapGrid
    label : str
        What the frames should say the terrain is.
    """
    tif = os.path.join(REPO, "input", "dem", "colca.tif")
    npy = tif.replace(".tif", ".npy")
    if os.path.exists(tif) or os.path.exists(npy):
        if not os.path.exists(npy):
            ss.build_elevation_cache(tif, npy)
        z = np.load(npy, mmap_mode="r")
        grid = ss.resolve_grid_geometry(tif, COLCA_ORIGIN_LAT)
        rows, cols = z.shape
        lo, hi = (1.0 - fraction) / 2.0, (1.0 + fraction) / 2.0
        z = np.ascontiguousarray(z[int(lo * rows):int(hi * rows),
                                   int(lo * cols):int(hi * cols)])
        return z, grid, "Colca, SRTM 1-arcsec"

    grid = ss.resolve_grid_geometry("no-such-file.tif", COLCA_ORIGIN_LAT,
                                    cell_size_deg=1.0 / 3600.0)
    z = synthetic.colca_like(700, grid.cell_size_x)
    return z, grid, "Synthetic canyon — no DEM present"


def _colca_candidates(z, grid, stride=5):
    """
    TAMBO's topographic screen over that ground, as the pipeline runs it.

    The screen's tqdm bar goes to stderr, where a notebook stores it as output. It
    reports a step that takes a fifth of a second on a crop this size, so it is carrying
    nothing, and three of these animations are built inside a notebook whose outputs are
    committed. Silenced here rather than in the library, which is right to be noisy when
    it is chewing through a 6 Mpx DEM.
    """
    with contextlib.redirect_stderr(io.StringIO()):
        return ss.get_candidates_chunked(
            z, grid, None, COLCA_ORIGIN_LAT, COLCA_ORIGIN_LON,
            min_slope_deg=TAMBO_SLOPE_BAND[0], max_slope_deg=TAMBO_SLOPE_BAND[1],
            candidate_stride=stride)


def the_walk():
    """One ray sweeping the elevation window, with the column depth it accumulates."""
    cell_x = synthetic.cell_sizes(-15.6)[1]
    # Reversed, and the detector lifted onto the ridge shoulder. Taken the other way
    # round the ground rises immediately in front of the detector, so every ray exits
    # at 0.0 km and the sweep shows nothing -- which is the "a detector on the ground
    # has every steep downward direction blocked at its own feet" effect, true but
    # degenerate. The interesting geometry is a detector looking out across ground that
    # falls away toward a distant target.
    z = synthetic.ridge_and_slope(700, cell_x)[350, ::-1]
    d = np.arange(z.size) * cell_x
    z0 = z[0] + 250.0
    elevations = np.linspace(2.0, -7.0, 90)

    with _styled():
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.6, 5.6), height_ratios=[3, 1])
        ax.set_xlim(0, d[-1] / 1000.0)
        ax.fill_between(d / 1000.0, z, z.min() - 200, color=ROCK_FILL, lw=0)
        ax.plot(d / 1000.0, z, color=ROCK_EDGE, lw=1.1)
        ax.plot([0], [z0], "o", color=DETECTOR, ms=6, zorder=5)
        ax.annotate("Detector", (0, z0), textcoords="offset points", xytext=(6, 10),
                    color=DETECTOR, fontsize=9)
        ray, = ax.plot([], [], color=WINDOW, lw=1.6, zorder=4)
        hit, = ax.plot([], [], "o", color=INK, ms=5, zorder=6)
        under, = ax.plot([], [], color="#B02A25", lw=2.6, alpha=0.75, zorder=3)
        title = ax.set_title("")
        ax.set_ylabel("Elevation (m)")
        ax.set_ylim(z.min() - 200, z.max() + 400)
        _tidy(ax)

        # Depth against the angle being swept, not against distance: the shared x-axis
        # would make it look like a profile of the terrain, which it is not.
        depths = []
        line, = bx.plot([], [], color="#B02A25", lw=1.4)
        bx.set_xlabel("Arrival elevation (°)")
        bx.set_ylabel("Column depth\n(10³ g/cm²)")
        bx.set_xlim(elevations[0], elevations[-1])
        ax.set_xlabel("Ground distance (km)")
        _tidy(bx)

        # Pre-compute so the axis can be scaled once rather than jumping every frame.
        for e in elevations:
            ray_z = z0 + np.tan(np.radians(e)) * d - d ** 2 / (2 * 6.371e6)
            below = ray_z < z
            first = np.argmax(below) if below.any() else -1
            depth = (below.sum() * cell_x * 100.0 * physics.CRUST_DENSITY_GCM3
                     if first >= 0 else 0.0)
            depths.append(depth / 1000.0)
        bx.set_ylim(0, max(depths) * 1.1 + 1)

        def frame(i):
            e = elevations[i]
            ray_z = z0 + np.tan(np.radians(e)) * d - d ** 2 / (2 * 6.371e6)
            below = ray_z < z
            ray.set_data(d / 1000.0, ray_z)
            if below.any():
                first = int(np.argmax(below))
                hit.set_data([d[first] / 1000.0], [z[first]])
                under.set_data(d[below] / 1000.0, ray_z[below])
                note = (f"exit at {d[first]/1000.0:.1f} km, "
                        f"{depths[i]*1000:,.0f} g/cm² behind")
            else:
                hit.set_data([], [])
                under.set_data([], [])
                note = "escapes to the sky"
            title.set_text(f"Arrival elevation {_minus(e)}°  —  {note}")
            line.set_data(elevations[: i + 1], depths[: i + 1])
            return ray, hit, under, title, line

        return fig, animation.FuncAnimation(fig, frame, frames=len(elevations),
                                            interval=1000 // FPS, blit=False)


def the_funnel():
    """The map draining stage by stage, with what survives each one."""
    cell_y, cell_x = synthetic.cell_sizes(-15.6)
    z = synthetic.ridge_and_slope(400, cell_x)
    dy, dx = np.gradient(z, cell_y, cell_x)
    slope = np.degrees(np.arctan(np.hypot(dy, dx)))

    stages = []
    everything = np.ones_like(z, dtype=bool)
    stages.append(("DEM pixels", everything))
    band = (slope >= 3.0) & (slope <= 25.0)
    stages.append(("Slope 3–25°", band))
    strided = np.zeros_like(band)
    strided[::5, ::5] = band[::5, ::5]
    stages.append(("Kept by stride 5", strided))
    rng = np.random.default_rng(0)
    accepted = strided & (rng.random(z.shape) < 0.60)
    stages.append(("Directions accepted", accepted))
    closed = binary_closing(accepted, structure=np.ones((17, 17)))
    stages.append(("After gap closing", closed))
    lab, n = label(closed)
    if n:
        sizes = np.bincount(lab.ravel())
        sizes[0] = 0
        pruned = lab == int(sizes.argmax())
    else:                                            # pragma: no cover - defensive
        pruned = closed
    stages.append(("After pruning", pruned))

    total = z.size
    with _styled():
        fig, ax = plt.subplots(figsize=(6.4, 6.0))
        ax.imshow(z, cmap="Greys_r", interpolation="nearest")
        overlay = ax.imshow(np.zeros(z.shape + (4,)), interpolation="nearest")
        title = ax.set_title("")
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_color(RULE)

        hold = 8                                     # frames each stage stays up

        def frame(i):
            name, mask = stages[min(i // hold, len(stages) - 1)]
            rgba = np.zeros(z.shape + (4,))
            rgba[..., 0], rgba[..., 1], rgba[..., 2] = 0.06, 0.42, 0.33
            rgba[..., 3] = np.where(mask, 0.62, 0.0)
            overlay.set_data(rgba)
            kept = int(mask.sum())
            title.set_text(f"{name}\n{kept:,} px   ({100*kept/total:.2f}% of the DEM)")
            return overlay, title

        return fig, animation.FuncAnimation(fig, frame, frames=len(stages) * hold,
                                            interval=1000 // FPS, blit=False)


def stride_and_closing():
    """Why an element smaller than the stride gap loses the mask entirely."""
    rng = np.random.default_rng(1)
    n = 200
    truth = np.zeros((n, n), dtype=bool)
    rr, cc = np.mgrid[0:n, 0:n]
    truth |= np.abs(cc - 60 - 18 * np.sin(rr / 26.0)) < 9      # a canyon-wall strip
    truth |= np.abs(cc - 140 - 12 * np.cos(rr / 18.0)) < 7
    truth &= rng.random((n, n)) < 0.96

    strided = np.zeros_like(truth)
    strided[::5, ::5] = truth[::5, ::5]

    # The gap striding leaves, in pixels. The transition is at the gap and it is sharp:
    # measured on this mask, a 3 px element recovers 0.04x of the accepted set and a
    # 5 px one recovers 0.61x -- fifteen times more for two pixels of element.
    gap = 5
    base = int(truth.sum())
    panels = [("Accepted, every pixel", truth, None),
              (f"Marked one pixel in {gap}", strided, None)]
    for k in (3, 5, 9):
        panels.append((f"Closed, element {k} px", binary_closing(strided, np.ones((k, k))), k))

    with _styled():
        fig, ax = plt.subplots(figsize=(5.6, 5.8))
        img = ax.imshow(truth, cmap="Greens", vmin=0, vmax=1.4, interpolation="nearest")
        title = ax.set_title("")
        note = ax.text(0.5, -0.06, "", transform=ax.transAxes, ha="center",
                       va="top", fontsize=9, color=MUTED)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_color(RULE)

        hold = 14

        def frame(i):
            name, mask, element = panels[min(i // hold, len(panels) - 1)]
            img.set_data(mask)
            kept = int(mask.sum())
            title.set_text(f"{name}\n{kept:,} px   ({kept/base:.2f}× the accepted set)")
            if element is None:
                note.set_text(f"The gap left is {gap} px — 154 m at 30.7 m pixels.\n"
                              f"TAMBO's closing element is 3 px; GRAND's is 33.")
            elif element < gap:
                note.set_text(f"{element} px cannot bridge a {gap} px gap. The mask "
                              f"never reconnects.")
            else:
                note.set_text(f"{element} px outruns the gap, and the mask comes back. "
                              f"The transition is at the gap,\nnot gradual: 3 px gives "
                              f"0.04×, 5 px gives 0.61×.")
            return img, title, note

        return fig, animation.FuncAnimation(fig, frame, frames=len(panels) * hold,
                                            interval=1000 // FPS, blit=False)


def energy_window():
    """The arrival window's lower edge climbing as Earth absorption bites."""
    energies = np.logspace(np.log10(30.0), np.log10(30000.0), 80)
    cutoffs = [physics.earth_absorption_cutoff_deg(e) for e in energies]

    with _styled():
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
        ax.axhline(0.0, color=RULE, lw=0.9)
        ax.plot(energies, cutoffs, color=RULE, lw=1.0, zorder=1)
        trace, = ax.plot([], [], color="#B02A25", lw=1.8, zorder=3)
        head, = ax.plot([], [], "o", color="#B02A25", ms=6, zorder=4)
        band = ax.axhspan(-3.0, 3.0, color=WINDOW, alpha=0.12, zorder=0)
        title = ax.set_title("")
        ax.set_xscale("log")
        ax.set_xlabel("Tau energy (PeV)")
        ax.set_ylabel("Arrival elevation (°)")
        ax.set_ylim(min(cutoffs) - 1.0, 3.6)
        ax.annotate("Nominal ±3° window", (energies[2], 3.0), fontsize=9,
                    color=WINDOW, va="bottom")
        _tidy(ax)

        def frame(i):
            trace.set_data(energies[: i + 1], cutoffs[: i + 1])
            head.set_data([energies[i]], [cutoffs[i]])
            # axhspan gives a Rectangle in (axes-x, data-y), so move its edges rather
            # than replacing vertices -- set_xy on a Rectangle takes a corner, not a
            # polygon, and silently means something else.
            band.set_y(cutoffs[i])
            band.set_height(3.0 - cutoffs[i])
            title.set_text(f"{energies[i]:,.0f} PeV — Earth absorption cuts below "
                           f"{_minus(cutoffs[i])}°, leaving "
                           f"{3.0 - cutoffs[i]:.1f}° of window")
            return trace, head, band, title

        return fig, animation.FuncAnimation(fig, frame, frames=len(energies),
                                            interval=1000 // FPS, blit=False)


def the_azimuth_fan():
    """What each bearing finds: a wall at the right range, its own hillside, or sky."""
    z, grid, ground = _colca_ground()
    cand = _colca_candidates(z, grid, stride=25)

    # Stand on the candidate the search itself would rank highest, not one picked by
    # eye -- with room around it for the plan view. `probe` is the ordinary scan.
    probe = arrival_scan.scan(cand, z, grid,
                              min_target_slope_deg=TAMBO_MIN_TARGET_SLOPE, **TAMBO_SCAN)
    half = int(6000.0 / grid.cell_size_x)
    room = ((cand[:, 0] > half) & (cand[:, 0] < z.shape[0] - half)
            & (cand[:, 1] > half) & (cand[:, 1] < z.shape[1] - half))
    if not room.any():                               # pragma: no cover - terrain dependent
        raise RuntimeError("no screened candidate sits far enough from the DEM edge")
    idx = np.flatnonzero(room)
    best = int(idx[np.argmax(probe["solid_angle_sr"][idx])])
    r0, c0, aspect = float(cand[best, 0]), float(cand[best, 1]), float(cand[best, 2])

    # One bearing per frame. Setting the candidate's aspect column *is* how you aim a
    # single-azimuth scan: azimuth_fan(1, 0.0) is [0.0], an offset of nothing from the
    # aspect, so the walk goes exactly where the third column points.
    bearings = arrival_scan.azimuth_fan(72, None)
    one = dict(TAMBO_SCAN, n_azimuths=1, half_width_deg=0.0)
    n_bins = TAMBO_SCAN["n_elev_bins"]
    accepted, struck, reach = [], [], []
    for b in bearings:
        c = np.array([[r0, c0, float(b)]])
        a = arrival_scan.scan(c, z, grid,
                              min_target_slope_deg=TAMBO_MIN_TARGET_SLOPE, **one)
        # The same walk with the range window and the wall-slope floor lifted. Without
        # it "found nothing" and "found the wrong thing" are the same zero.
        s = arrival_scan.scan(c, z, grid, **dict(one, min_dist_km=0.0))
        accepted.append(int(a["cells"][0]))
        struck.append(int(s["cells"][0]))
        reach.append(float(s["mean_distance_m"][0]))
    accepted, struck = np.array(accepted), np.array(struck)
    reach = np.array(reach)

    outcomes = [(ACCEPT, "Wall accepted"), (SEQUENCE[0], "Wall too gentle"),
                (TOO_NEAR, "Own hillside"), (SKY, "Sky")]

    def verdict(i):
        """Which of ``outcomes`` this bearing found, and how to say it."""
        if accepted[i] > 0:
            return 0, f"wall at {reach[i]/1000:.1f} km — {accepted[i]} of {n_bins} bins"
        if struck[i] == 0:
            return 3, "sky — nothing within 5 km"
        if reach[i] < TAMBO_SCAN["min_dist_km"] * 1000.0:
            return 2, f"its own hillside, {reach[i]/1000:.1f} km away"
        return 1, f"rock at {reach[i]/1000:.1f} km, but no wall steep enough"

    # Only the outcomes this candidate actually produces. A legend naming a colour that
    # never appears asks the reader to hunt for something that is not there -- and on a
    # canyon wall with a 2 km floor, most bearings are its own hillside and nothing is
    # sky.
    present = sorted({verdict(i)[0] for i in range(len(bearings))})

    win = z[int(r0) - half:int(r0) + half, int(c0) - half:int(c0) + half]
    ex_x, ex_y = half * grid.cell_size_x / 1000.0, half * grid.cell_size_y / 1000.0
    theta = np.radians(bearings)

    with _styled():
        fig = plt.figure(figsize=(10.4, 5.4))
        ax = fig.add_subplot(1, 2, 1)
        bx = fig.add_subplot(1, 2, 2, projection="polar")

        ax.imshow(win, cmap="Greys_r", origin="upper",
                  extent=(-ex_x, ex_x, -ex_y, ex_y), interpolation="nearest")
        ax.set_xlabel("East (km)")
        ax.set_ylabel("North (km)")
        ax.set_aspect("equal")
        rays = []
        for i in range(len(bearings)):
            line, = ax.plot([], [], lw=1.1, alpha=0.75, zorder=3,
                            color=outcomes[verdict(i)[0]][0])
            rays.append(line)
        live, = ax.plot([], [], lw=2.4, color=INK, zorder=5)
        # The range window, drawn. Without it the rays that stop after 200 m look like
        # a rendering failure rather than the criterion doing its job.
        ring = np.linspace(0, 2 * np.pi, 200)
        for radius in (TAMBO_SCAN["min_dist_km"], TAMBO_SCAN["max_dist_km"]):
            ax.plot(radius * np.sin(ring), radius * np.cos(ring), ls=":", lw=0.9,
                    color=WINDOW, zorder=4)
        ax.text(0.03, 0.965, "Dotted: the accepted range, 2–5 km",
                transform=ax.transAxes, fontsize=8.5, color=WINDOW, va="top")
        ax.plot([0], [0], "o", color=DETECTOR, ms=7, zorder=6)
        ax.set_xlim(-ex_x, ex_x)
        ax.set_ylim(-ex_y, ex_y)
        ss.add_scale_bar(ax, 1.0)
        ss.add_north_arrow(ax)

        # The configured fan, so the sweep is read against what the search actually
        # tests: nine bearings within 60 degrees of the aspect, and nothing outside.
        fan = aspect + arrival_scan.azimuth_fan(TAMBO_SCAN["n_azimuths"],
                                                TAMBO_SCAN["half_width_deg"])
        wedge = np.radians(np.linspace(fan[0], fan[-1], 60))
        bx.fill_between(wedge, 0, n_bins, color=WINDOW, alpha=0.10, lw=0, zorder=0)
        bx.plot(np.radians(fan), np.full(fan.size, n_bins), "|", color=WINDOW,
                ms=9, zorder=2)
        bx.set_theta_zero_location("N")
        bx.set_theta_direction(-1)
        bx.set_rlim(0, n_bins)
        bx.set_rlabel_position(112.5)
        bx.set_yticks([0, n_bins // 2, n_bins])
        bx.grid(color=RULE, lw=0.6)
        rock, = bx.plot([], [], color=ROCK_EDGE, lw=1.2, zorder=3)
        good, = bx.plot([], [], color=ACCEPT, lw=2.0, zorder=4)
        spot, = bx.plot([], [], "o", color=INK, ms=5, zorder=5)
        # The two polar curves need naming where they are: the four-colour legend
        # describes the rays on the map, and grey there means sky, which would read as
        # a contradiction against a grey curve meaning the opposite.
        bx.text(0.5, -0.10, "Elevation bins against bearing. Outer curve: bins that "
                "found rock at all.\nInner: bins accepted. The gap is rock at the "
                "wrong range; outside is sky.",
                transform=bx.transAxes, ha="center", va="top", fontsize=8.5,
                color=MUTED, linespacing=1.4)

        fig.legend([plt.Line2D([], [], color=outcomes[k][0], lw=2.2) for k in present],
                   [outcomes[k][1] for k in present],
                   loc="upper center", ncol=min(4, len(present)),
                   bbox_to_anchor=(0.5, 1.0), frameon=False)
        title = fig.text(0.5, 0.925, "", ha="center", fontsize=10)
        fig.text(0.29, 0.045, f"{ground}. Candidate aspect {aspect:.0f}°.",
                 ha="center", fontsize=8.5, color=MUTED)
        fig.subplots_adjust(top=0.86, bottom=0.17)

        def frame(i):
            d = reach[i] if struck[i] else TAMBO_SCAN["max_range_m"]
            x, y = d * np.sin(theta[i]) / 1000.0, d * np.cos(theta[i]) / 1000.0
            rays[i].set_data([0, x], [0, y])
            live.set_data([0, x], [0, y])
            rock.set_data(theta[: i + 1], struck[: i + 1])
            good.set_data(theta[: i + 1], accepted[: i + 1])
            spot.set_data([theta[i]], [accepted[i]])
            off = (bearings[i] - aspect + 180.0) % 360.0 - 180.0
            title.set_text(f"Bearing {bearings[i]:.0f}°  "
                           f"({_minus(off, '{:+.0f}')}° from aspect)"
                           f"  —  {verdict(i)[1]}")
            return (*rays, live, rock, good, spot, title)

        return fig, animation.FuncAnimation(fig, frame, frames=len(bearings),
                                            interval=1000 // FPS, blit=False)


def product_collapse():
    """Six components multiplied in under a cut that never moves."""
    z, grid, ground = _colca_ground()
    cand = _colca_candidates(z, grid, stride=5)
    obs = arrival_scan.scan(cand, z, grid,
                            min_target_slope_deg=TAMBO_MIN_TARGET_SLOPE, **TAMBO_SCAN)
    # Altitude is the one observable the scan does not return: it is a property of the
    # candidate pixel, not of what the walk found, and the footprint term needs it.
    obs["altitude_m"] = z[cand[:, 0].astype(int), cand[:, 1].astype(int)].astype(float)

    # TAMBO's configuration, restricted to the parameters that switch a component on.
    cfg = dict(solid_angle_half_sr=0.8, grammage_mode="particle",
               grammage_band_gcm2=(236.0, 1287.0), shower_development_m=0.0,
               decay_energy_min_pev=3.0, decay_energy_max_pev=1000.0,
               decay_spectral_index=2.0, spacing_m=100.0, composition="product")
    window = (TAMBO_SCAN["min_dist_km"] * 1000.0, TAMBO_SCAN["max_dist_km"] * 1000.0)
    _, parts = scoring.score_candidates(obs, cfg, distance_window_m=window)

    viable = obs["cells"] > 0
    parts = {k: v[viable] for k, v in parts.items()}
    names = list(parts)
    n = int(viable.sum())

    # Weights are exponents under a product composition, so ramping the newest one from
    # 0 to 1 is a continuous walk between two states compose() itself can produce --
    # not a dissolve between two pictures.
    hold = 15
    edges = np.linspace(0.0, 1.0, 61)
    frames = len(names) * hold

    def population(i):
        k, t = divmod(i, hold)
        k = min(k, len(names) - 1)
        active = names[: k + 1]
        w = {nm: 1.0 for nm in active}
        w[active[-1]] = t / (hold - 1) if k or hold > 1 else 1.0
        return k, active, scoring.compose({nm: parts[nm] for nm in active}, "product", w)

    # Where each component leaves the population. Kept as a running tally because the
    # measurement is not the one the shape of the histogram suggests: two components do
    # essentially all of the work, and `distance` provably does none, since the scan
    # already applied the same window as a hard criterion before scoring saw it.
    survivors = [100.0]
    for k in range(len(names)):
        run = scoring.compose({nm: parts[nm] for nm in names[: k + 1]}, "product")
        survivors.append(100.0 * float((run >= TAMBO_MIN_SCORE).mean()))
    survivors = np.array(survivors)

    with _styled():
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.2, 6.6),
                                     height_ratios=[2.1, 1])
        counts, _ = np.histogram(np.ones(n), bins=edges)
        bars = ax.bar(edges[:-1], counts, width=np.diff(edges), align="edge",
                      color=WINDOW, alpha=0.85, lw=0)
        ax.axvline(TAMBO_MIN_SCORE, color="#B02A25", lw=1.6, zorder=5)
        ax.text(TAMBO_MIN_SCORE + 0.012, 0.95, f"min_score = {TAMBO_MIN_SCORE}",
                transform=ax.get_xaxis_transform(), color="#B02A25", fontsize=9,
                va="top")
        med = ax.axvline(1.0, color=INK, lw=1.0, ls="--", zorder=4)
        med_label = ax.text(1.0, 0.95, " median", transform=ax.get_xaxis_transform(),
                            color=INK, fontsize=9, va="top")
        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_ylim(0.5, n * 2.0)
        ax.set_xlabel("Composed score")
        ax.set_ylabel("Candidates")
        title = ax.set_title("")
        _tidy(ax)

        bx.plot(np.arange(len(survivors)), survivors, color=RULE, lw=1.0, zorder=1)
        tally, = bx.plot([], [], color="#B02A25", lw=1.9, marker="o", ms=4, zorder=3)
        bx.set_xticks(np.arange(len(survivors)))
        bx.set_xticklabels(["none"] + names, rotation=20, ha="right")
        bx.set_ylim(0, 105)
        bx.set_ylabel("Above the cut (%)")
        bx.set_xlabel("Components multiplied in")
        _tidy(bx)

        fig.text(0.5, 0.025,
                 f"{n:,} candidates on {ground}. The cut does not move; the population "
                 f"walks under it — but not evenly.\nsolid_angle does most of the work "
                 f"and distance does none, the scan having already applied that same "
                 f"window as a hard criterion.",
                 ha="center", fontsize=8.5, color=MUTED, linespacing=1.5)
        fig.subplots_adjust(top=0.90, bottom=0.20, hspace=0.55)

        def frame(i):
            k, active, total = population(i)
            counts, _ = np.histogram(total, bins=edges)
            for bar, h in zip(bars, counts):
                bar.set_height(max(h, 0.5))
            kept = 100.0 * float((total >= TAMBO_MIN_SCORE).mean())
            median = float(np.median(total))
            med.set_xdata([median, median])
            med_label.set_x(median)
            x = np.append(np.arange(k + 1), k + (i % hold) / (hold - 1.0))
            tally.set_data(x, np.append(survivors[: k + 1], kept))
            title.set_text(" × ".join(active) + f"\n{kept:.1f}% still above the cut"
                           f"   (median {median:.3f})")
            return (*bars, med, med_label, tally, title)

        return fig, animation.FuncAnimation(fig, frame, frames=frames,
                                            interval=1000 // FPS, blit=False)


def slope_criterion():
    """min_target_slope_deg climbing through the wall-slope distribution, on the ground."""
    z, grid, ground = _colca_ground()
    cand = _colca_candidates(z, grid, stride=5)
    rr, cc = cand[:, 0].astype(int), cand[:, 1].astype(int)

    # The distribution the cut crosses, measured with the cut off. This is a *mean*
    # over each candidate's accepted directions, whereas the criterion is applied to
    # each direction separately -- which is why the mask survives well past the median.
    base = arrival_scan.scan(cand, z, grid, **TAMBO_SCAN)
    wall = base["target_slope_deg"][base["cells"] > 0]

    thresholds = np.arange(0.0, 72.0, 3.0)
    masks, kept = [], []
    block = np.ones((5, 5), dtype=bool)              # one stride cell, so dots read as area
    for t in thresholds:
        r = arrival_scan.scan(cand, z, grid,
                              min_target_slope_deg=(None if t == 0 else float(t)),
                              **TAMBO_SCAN)
        ok = r["cells"] > 0
        m = np.zeros(z.shape, dtype=bool)
        m[rr[ok], cc[ok]] = True
        masks.append(binary_dilation(m, block))
        kept.append(int(ok.sum()))
        print(f"   {t:4.0f}° -> {kept[-1]:>8,} candidates", flush=True)
    kept = np.array(kept)

    # The map keeps square pixels, so its panel is as wide as the DEM is wide. Sizing
    # the figure from the crop rather than fixing it means the histogram beside it
    # comes out the same height whatever shape the DEM turns out to be.
    map_h = 3.7
    map_w = map_h * z.shape[1] / z.shape[0]
    hist_w = 4.2

    with _styled():
        fig, (ax, bx) = plt.subplots(1, 2, figsize=(map_w + hist_w + 1.4, map_h + 1.7),
                                     gridspec_kw={"width_ratios": [map_w, hist_w]})
        ax.imshow(z, cmap="Greys_r", interpolation="nearest")
        overlay = ax.imshow(np.zeros(z.shape + (4,)), interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_color(RULE)
        ss.add_scale_bar(ax, grid.cell_size_x / 1000.0)
        ss.add_north_arrow(ax)

        bx.hist(wall, bins=np.arange(0.0, 72.0, 1.5), color=ROCK_FILL,
                edgecolor=ROCK_EDGE, lw=0.5)
        cut = bx.axvline(0.0, color="#B02A25", lw=1.8, zorder=5)
        bx.set_xlabel("Mean wall slope a candidate sees (°)")
        bx.set_ylabel("Candidates")
        bx.set_xlim(0, 71)
        _tidy(bx)
        title = fig.text(0.5, 0.945, "", ha="center", fontsize=10.5)
        fig.text(0.5, 0.035,
                 f"{ground}. The criterion is applied per direction, the histogram is a "
                 f"mean over each candidate's accepted ones — so the mask outlives its "
                 f"own median.", ha="center", fontsize=8.5, color=MUTED)
        fig.subplots_adjust(top=0.89, bottom=0.14)

        hold = 3

        def frame(i):
            k = min(i // hold, len(thresholds) - 1)
            rgba = np.zeros(z.shape + (4,))
            rgba[..., 0], rgba[..., 1], rgba[..., 2] = 0.06, 0.42, 0.33
            rgba[..., 3] = np.where(masks[k], 0.70, 0.0)
            overlay.set_data(rgba)
            cut.set_xdata([thresholds[k], thresholds[k]])
            share = 100.0 * kept[k] / kept[0]
            title.set_text(f"min_target_slope_deg = {thresholds[k]:.0f}°   —   "
                           f"{kept[k]:,} candidates, {share:.0f}% of what no criterion "
                           f"at all accepts")
            return overlay, cut, title

        return fig, animation.FuncAnimation(fig, frame, frames=len(thresholds) * hold,
                                            interval=1000 // FPS, blit=False)


def tau_in_rock():
    """More rock is not better: production rises, escape collapses, the product peaks."""
    energy_pev = 1000.0
    depth = np.logspace(4.0, 8.0, 90)
    optimum = physics.production_escape_optimum_gcm2(energy_pev)
    exit_p = physics.tau_exit_probability(depth, energy_pev)
    survival = physics.tau_survival(depth, energy_pev)
    # The band a run actually scores against, not physics.depth_band_from_energy() --
    # which returns (5.6e7, 2.9e8) for TAMBO's 3 PeV - 1 EeV range and so excludes this
    # optimum by more than an order of magnitude. See ROADMAP 6.44.
    band = scoring.DEFAULT_SCORE_CONFIG["depth_band_gcm2"]

    # The tau's own energy, integrated rather than assumed: dE/dx = -beta(E) E, and
    # beta itself rises with energy, so the loss is not a single exponential.
    e_rel = np.empty_like(depth)
    e, x_prev = energy_pev, 0.0
    for i, x in enumerate(depth):
        e *= float(np.exp(-physics.tau_energy_loss_beta(max(e, 1e-3)) * (x - x_prev)))
        e_rel[i], x_prev = e / energy_pev, x

    km = depth / physics.CRUST_DENSITY_GCM3 / 1e5

    with _styled():
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.6, 6.0), sharex=True)
        ax.plot(depth, survival, color=RULE, lw=1.0)
        ax.plot(depth, e_rel, color=RULE, lw=1.0)
        surv, = ax.plot([], [], color="#B02A25", lw=1.9, zorder=3)
        loss, = ax.plot([], [], color=WINDOW, lw=1.9, zorder=3)
        ax.set_ylabel("Fraction remaining")
        ax.set_ylim(-0.03, 1.05)
        ax.set_xscale("log")
        ax.annotate("Survival", (depth[0], 1.0), xytext=(4, -12), fontsize=9,
                    textcoords="offset points", color="#B02A25")
        ax.annotate("Energy", (depth[0], 1.0), xytext=(4, -26), fontsize=9,
                    textcoords="offset points", color=WINDOW)
        _tidy(ax)

        bx.axvspan(*band, color=WINDOW, alpha=0.10, lw=0, zorder=0)
        bx.plot(depth, exit_p, color=RULE, lw=1.0)
        trace, = bx.plot([], [], color=DETECTOR, lw=1.9, zorder=3)
        head, = bx.plot([], [], "o", color=INK, ms=5, zorder=4)
        bx.axvline(optimum, color=INK, lw=1.0, ls="--", zorder=2)
        bx.annotate(f"Optimum {optimum/1e6:.1f}×10⁶ g/cm²\n"
                    f"= {optimum/physics.CRUST_DENSITY_GCM3/1e5:.0f} km of rock",
                    (optimum, exit_p.max()), xytext=(-12, -34), fontsize=9,
                    textcoords="offset points", color=INK, va="top", ha="right")
        bx.annotate("The scored depth band", (np.sqrt(band[0] * band[1]), 0.0),
                    xytext=(0, 6), fontsize=9, ha="center",
                    textcoords="offset points", color=WINDOW)
        bx.set_xlabel("Column depth traversed (g/cm²)")
        bx.set_ylabel("Exit probability")
        bx.set_ylim(0, exit_p.max() * 1.18)
        _tidy(bx)

        title = fig.text(0.5, 0.955, "", ha="center", fontsize=10.5)
        fig.text(0.5, 0.055,
                 "A 1 EeV tau. Production grows with the rock available; escape "
                 "collapses once the depth passes the tau's range.\nTheir product has a "
                 "maximum, so the criterion is a band and not a floor: past 22 km, more "
                 "rock is worse.",
                 ha="center", fontsize=8.5, color=MUTED, linespacing=1.5)
        fig.subplots_adjust(top=0.90, bottom=0.20, hspace=0.12)

        def frame(i):
            surv.set_data(depth[: i + 1], survival[: i + 1])
            loss.set_data(depth[: i + 1], e_rel[: i + 1])
            trace.set_data(depth[: i + 1], exit_p[: i + 1])
            head.set_data([depth[i]], [exit_p[i]])
            side = "short of" if depth[i] < optimum else "past"
            title.set_text(f"{km[i]:,.1f} km of rock — energy {100*e_rel[i]:.0f}%, "
                           f"survival {100*survival[i]:.0f}%, {side} the optimum")
            return surv, loss, trace, head, title

        return fig, animation.FuncAnimation(fig, frame, frames=len(depth),
                                            interval=1000 // FPS, blit=False)


BUILDERS = {"the_walk": the_walk, "the_azimuth_fan": the_azimuth_fan,
            "the_funnel": the_funnel, "stride_and_closing": stride_and_closing,
            "product_collapse": product_collapse, "slope_criterion": slope_criterion,
            "tau_in_rock": tau_in_rock, "energy_window": energy_window}


def write_mp4_with_stills(name, fig, anim, out_dir, at=(0.0, 0.5, 1.0)):
    """
    Writes one MP4 and returns still frames from the *same* pass over the animation.

    Two reasons this is one function and not two. It is half the work — several of these
    spend twenty seconds scanning terrain before a frame is drawn. And the builders
    **accumulate**: the ray drawn at frame 30 is still on the axes at frame 60, which is
    what makes the fan fill in. So the frames can be walked exactly once; a second pass
    would begin with everything already drawn, and its "first" frame would be a lie.

    For showing an animation somewhere that cannot play one — a notebook whose outputs
    are committed to a repository, a printed page, a slide that has to survive a
    projector.

    Parameters
    ----------
    name : str
        Basename for the file, without extension.
    fig, anim : Figure, FuncAnimation
        As a builder in :data:`BUILDERS` returns them.
    out_dir : str
        Directory to write into; created if absent.
    at : sequence of float, optional
        Where to take stills, as fractions of the way through.

    Returns
    -------
    path : str
        The MP4 written.
    stills : list of bytes
        PNG bytes, one per entry in ``at``, at the video's own frame size.
    """
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{name}.mp4")
    n = getattr(anim, "_save_count", None) or 1
    wanted = sorted({min(n - 1, max(0, round(f * (n - 1)))) for f in at})
    stills = []
    # Driving the writer by hand means anim.save() never runs, and matplotlib's
    # finaliser warns that an animation was discarded without being rendered. It was
    # rendered; it just was not rendered by the method the finaliser watches.
    anim._draw_was_started = True
    writer = animation.FFMpegWriter(fps=FPS, bitrate=2400)
    with writer.saving(fig, path, dpi=100):
        for i in range(n):
            # FuncAnimation offers no public way to ask for frame i, so this is the
            # callable it was built with. Walked in order, for the reason above.
            anim._func(i)
            writer.grab_frame()
            if i in wanted:
                buf = io.BytesIO()
                fig.savefig(buf, format="png")     # fig.dpi is the writer's: same size
                stills.append(buf.getvalue())
    plt.close(fig)
    return path, stills


def write(name, fig, anim, out_dir, formats):
    """Writes one animation in each requested format, reporting what it cost."""
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for fmt in formats:
        path = os.path.join(out_dir, f"{name}.{fmt}")
        try:
            if fmt == "mp4":
                anim.save(path, writer=animation.FFMpegWriter(fps=FPS, bitrate=2400))
            else:
                anim.save(path, writer=animation.PillowWriter(fps=FPS))
        except Exception as e:                       # pragma: no cover - env dependent
            print(f"   {fmt}: skipped ({type(e).__name__}: {e})")
            continue
        size = os.path.getsize(path) / 1024.0
        print(f"   {fmt}: {os.path.relpath(path, REPO)}  ({size:,.0f} KiB)")
        written.append(path)
    plt.close(fig)
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=sorted(BUILDERS), help="build just one")
    ap.add_argument("--out", default=os.path.join(REPO, "output", "animations"),
                    help="where to write (default: output/animations/, gitignored)")
    ap.add_argument("--format", default="mp4,gif",
                    help="comma-separated: mp4, gif (default: both)")
    args = ap.parse_args()

    formats = [f.strip() for f in args.format.split(",") if f.strip()]
    names = [args.only] if args.only else list(BUILDERS)
    for name in names:
        print(f"\n{name} ...", flush=True)
        fig, anim = BUILDERS[name]()
        write(name, fig, anim, args.out, formats)
    print(f"\nwrote to {os.path.relpath(args.out, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
