#!/usr/bin/env python
"""
Animations of the things in this search that a still picture explains badly.

Four, chosen because each one shows a *process* whose intermediate states are the
point. Anything that is a single state is already a figure in ``oroscope.figures``, and
turning it into a movie would be decoration rather than explanation.

``the_walk``
    One backward ray sweeping down through the elevation window over real terrain: the
    first intersection sliding along the profile, and the column depth accumulating
    behind it. This is the mechanism the whole search rests on and the hardest thing in
    the project to convey in prose.
``the_funnel``
    The map draining stage by stage -- slope, stride, directions accepted, closing,
    pruning -- with the surviving count. The funnel table says where the candidates
    went; this shows *where on the ground* they went.
``stride_and_closing``
    Why TAMBO's area is 4.75x low. A strided mask closed with an element that bridges
    the gaps, and the same mask closed with one that does not. The measurement is in
    ROADMAP 6.34; this is what it looks like.
``energy_window``
    The arrival window narrowing as energy rises, its lower edge climbing from -4.4
    degrees at 100 PeV to -0.9 at 10 EeV. A falsifiable prediction, animated over the
    quantity it is a prediction about.

Everything is built from committed code and synthetic terrain, so these reproduce on
any clone with no DEM present.

    python tools/make_animations.py                  # all four, MP4 then GIF
    python tools/make_animations.py --only the_walk
    python tools/make_animations.py --format gif --out docs/source/_static

MP4 needs ffmpeg; GIF falls back to pillow, which is always available. Outputs land in
``output/animations/`` by default, which is gitignored -- pass ``--out`` to place the
small ones somewhere they can be committed.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")                                # a tool, not a library: see trap 3
import matplotlib.animation as animation             # noqa: E402
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
from scipy.ndimage import binary_closing, label      # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "tests"))

from oroscope import physics                         # noqa: E402
from oroscope.figures import (                       # noqa: E402
    DETECTOR, INK, MUTED, ROCK_EDGE, ROCK_FILL, RULE, WINDOW, _styled, _tidy, _minus,
)
import synthetic                                     # noqa: E402

FPS = 12


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


BUILDERS = {"the_walk": the_walk, "the_funnel": the_funnel,
            "stride_and_closing": stride_and_closing, "energy_window": energy_window}


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
