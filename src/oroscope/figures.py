"""
Publication-quality schematics of the geometry the search computes.

Each routine returns a :class:`matplotlib.figure.Figure` built from the same physics
the code uses, so a diagram cannot drift away from the implementation it illustrates.
They are here rather than in the documentation tree so they can be imported, restyled
and reused --- in a talk, a proposal, or a paper --- without copying code out of an
``.rst`` file.

Styling is applied per-figure through a context manager rather than by mutating global
``rcParams``, so importing this module does not change the appearance of anybody
else's plots.

Examples
--------
>>> from oroscope import figures
>>> fig = figures.walk_mechanism()
>>> fig.savefig("walk.pdf", bbox_inches="tight")   # doctest: +SKIP
"""

from __future__ import annotations

import contextlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

__all__ = ["walk_mechanism", "canyon_geometry", "decay_and_shower"]

# House style for anything a reader sees on a figure: **the first word of every axis
# label, title, legend entry and annotation is capitalised.** Applies to the notebooks
# and to the maps the pipeline writes, not only to this module. Units stay as they are
# ("Elevation (km)", not "Elevation (KM)"), and a label that starts with a function
# name or a symbol is left alone -- `band_score(x, 6, 14)` and `$d$` are code and
# mathematics, not prose.

# A restrained palette, fixed here so every figure in a set matches.
INK, MUTED, RULE = "#1A1A1A", "#6B6B6B", "#C9CCC8"
ROCK_FILL, ROCK_EDGE = "#D9CDB8", "#7A6A4F"
DETECTOR = "#0F6B54"
WINDOW = "#2C6E8F"
# Ordered low elevation angle to high, and warm so they read against the terrain
SEQUENCE = ["#C8901A", "#D2621B", "#B02A25", "#7B2D6B"]

_STYLE = {
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.labelsize": 9.5,
    "axes.titlesize": 10,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.frameon": False,
    "mathtext.fontset": "dejavusans",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
}


@contextlib.contextmanager
def _styled():
    """Applies the figure style locally, leaving global rcParams untouched."""
    with mpl.rc_context(_STYLE):
        yield


def _tidy(ax):
    """Drops the top and right spines and greys what remains."""
    ax.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(RULE)


def _minus(value, fmt="{:+.1f}"):
    """Formats with a typographic minus rather than a hyphen."""
    return fmt.format(value).replace("-", "−")


def walk_mechanism(earth_radius_km=6371.0, detector_elevation_km=1.05,
                   bin_edges_deg=(-1.2, -0.4, 0.4, 1.2), figsize=(10.2, 4.15)):
    r"""
    One profile walk fills every elevation bin at once.

    The central algorithmic claim, drawn. Panel (a) is the terrain profile with the
    first terrain met for each elevation bin; panel (b) is the quantity that makes it
    work, the apparent elevation angle

    .. math::

        \theta_{\rm terrain}(d) = \arctan\!\left(
            \frac{z(d) - d^2/2R - z_0}{d}\right),

    together with its running maximum. Because that maximum only increases, each new
    value claims a contiguous band of bins, so one pass fills them all --- which is why
    the elevation binning is nearly free and the azimuth count sets the cost.

    Parameters
    ----------
    earth_radius_km : float, optional
        Earth radius used for the :math:`d^2/2R` curvature drop. The true radius, not
        the inflated radio one: the particle trajectory is not refracted.
    detector_elevation_km : float, optional
        Elevation of the candidate pixel the walk starts from.
    bin_edges_deg : sequence of float, optional
        Elevation-bin edges to trace, in degrees. Drawn in the order given, coloured
        low to high.
    figsize : tuple of float, optional
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
        The finished two-panel figure.

    Examples
    --------
    >>> from oroscope import figures
    >>> fig = figures.walk_mechanism()
    >>> len(fig.axes)
    2
    """
    edges = np.asarray(bin_edges_deg, dtype=float)
    d = np.linspace(0.05, 90, 3000)
    # A synthetic profile with a near rise, an intermediate ridge and a far massif, so
    # that the running maximum has somewhere to plateau and somewhere to jump
    z = (2.35 * np.exp(-((d - 52) / 13) ** 2)
         + 1.30 * np.exp(-((d - 30) / 7) ** 2)
         + 0.85 * np.exp(-((d - 76) / 10) ** 2)
         + 0.30 * np.exp(-((d - 12) / 4) ** 2))
    z0 = float(detector_elevation_km)
    z_app = z - d ** 2 / (2.0 * earth_radius_km)

    theta = np.degrees(np.arctan((z_app - z0) / d))
    running = np.maximum.accumulate(theta)
    first = [int(np.argmax(running >= e)) for e in edges]

    with _styled():
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 2, width_ratios=[1.22, 1.0], wspace=0.26,
                              left=0.065, right=0.985, top=0.845, bottom=0.135)
        ax0, ax1 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

        # -- (a) the profile and the rays
        ax0.fill_between(d, -0.35, z_app, color=ROCK_FILL, lw=0, zorder=1)
        ax0.plot(d, z_app, color=ROCK_EDGE, lw=1.1, zorder=3)
        for edge, i, colour in zip(edges, first, SEQUENCE):
            ax0.plot([0, d[i]], [z0, z_app[i]], color=colour, lw=1.5, zorder=4,
                     solid_capstyle="round")
            ax0.plot(d[i], z_app[i], "o", color=colour, ms=5.0, zorder=5,
                     mec="white", mew=0.8)
            ax0.annotate(_minus(edge) + "$\\degree$", (d[i], z_app[i]),
                         textcoords="offset points", xytext=(6, -9),
                         fontsize=8.0, color=colour, weight="bold")
        ax0.plot(0, z0, "^", color=DETECTOR, ms=10, zorder=6, mec="white", mew=0.9)
        ax0.text(1.6, z0 + 0.30, "detector", color=DETECTOR, fontsize=9, weight="bold")
        ax0.text(46, 2.72, "first terrain met, per elevation bin",
                 fontsize=8.4, color=INK, ha="center")
        ax0.set_xlabel("Ground distance  $d$  (km)")
        ax0.set_ylabel("Elevation (km)")
        ax0.set_title("(a)  one profile walk", loc="left", weight="bold")
        ax0.set_xlim(0, 90)
        ax0.set_ylim(-0.35, 3.0)
        # Measured from the axes box rather than asserted, so it cannot go stale
        box = ax0.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        exaggeration = (90.0 / box.width) / (3.35 / box.height)
        ax0.text(0.985, 0.045,
                 f"vertical exaggeration $\\approx$ {exaggeration:.0f}:1",
                 transform=ax0.transAxes, ha="right", fontsize=7.8, color=MUTED,
                 style="italic")
        _tidy(ax0)

        # -- (b) the mechanism
        ax1.axhspan(edges[0], edges[-1], color=WINDOW, alpha=0.07, lw=0, zorder=0)
        keep = d > 6
        ax1.plot(d[keep], theta[keep], color=MUTED, lw=0.85, zorder=2,
                 label=r"$\theta_{\rm terrain}(d)$")
        ax1.plot(d[keep], running[keep], color=INK, lw=1.9, zorder=3,
                 label="Running maximum")
        for edge, i, colour in zip(edges, first, SEQUENCE):
            ax1.axhline(edge, color=colour, lw=0.8, ls=(0, (4, 3)), zorder=1)
            ax1.plot(d[i], edge, "o", color=colour, ms=5.0, zorder=5,
                     mec="white", mew=0.8)
            ax1.annotate("", (d[i], edge), (d[i], -2.05),
                         arrowprops=dict(arrowstyle="-", color=colour, lw=0.6,
                                         alpha=0.5))
        ax1.text(88, edges[-1] + 0.10, "accepted elevation window", ha="right",
                 va="bottom", fontsize=8.2, color=WINDOW)
        ax1.annotate("Each new maximum claims\nthe bins it has risen past",
                     (d[first[2]], running[first[2]]), (41, -0.55), fontsize=8.4,
                     color=INK, ha="left",
                     arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8,
                                     connectionstyle="arc3,rad=0.25"))
        ax1.set_xlabel("Ground distance  $d$  (km)")
        ax1.set_ylabel("Apparent elevation angle  (deg)")
        ax1.set_title("(b)  why it fills every bin at once", loc="left", weight="bold")
        ax1.set_xlim(6, 90)
        ax1.set_ylim(-2.05, 1.95)
        ax1.legend(loc="lower left", fontsize=8.3, handlelength=1.6,
                   borderpad=0.2, labelspacing=0.35)
        _tidy(ax1)

    return fig


def canyon_geometry(depth_m=1500.0, floor_width_m=1000.0, wall_slope_deg=40.6,
                    figsize=(7.6, 4.0)):
    r"""
    The canyon-crossing geometry a particle array such as TAMBO selects on.

    Drawn to scale, unlike the long-range figure: a canyon is a few kilometres across
    and one and a half deep, so no exaggeration is needed and none is applied. That is
    itself worth showing --- the arrival directions really do span tens of degrees
    here, where GRAND's span three.

    Two criteria are marked because they are separate and are easy to conflate. The
    **near wall** is the ground the array stands on and must be deployable; the **far
    wall** is where the tau exits and must be steep. A single slope band cannot express
    both.

    Parameters
    ----------
    depth_m : float, optional
        Rim-to-floor depth. Colca is about 1500 m.
    floor_width_m : float, optional
        Width of the flat valley floor.
    wall_slope_deg : float, optional
        Slope of both walls. Colca's published depth and ~4.5 km rim separation imply
        about 40.6 degrees, which is far outside GRAND's 3-25 degree deployable band.
    figsize : tuple of float, optional
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from oroscope import figures
    >>> fig = figures.canyon_geometry()
    >>> round(fig.get_figwidth(), 1)
    7.6
    """
    run = depth_m / np.tan(np.radians(wall_slope_deg))
    rim = floor_width_m / 2.0 + run
    x = np.linspace(-rim * 1.35, rim * 1.35, 1200)
    from_edge = np.abs(x) - floor_width_m / 2.0
    z = np.clip(np.clip(from_edge, 0, None) * np.tan(np.radians(wall_slope_deg)),
                0, depth_m)

    # A detector part way up the near (left) wall, looking across
    x_det = -(floor_width_m / 2.0 + run * 0.55)
    z_det = float(np.interp(x_det, x, z))

    with _styled():
        fig, ax = plt.subplots(figsize=figsize)
        ax.fill_between(x, -120, z, color=ROCK_FILL, lw=0, zorder=1)
        ax.plot(x, z, color=ROCK_EDGE, lw=1.2, zorder=3)

        # Rays to the far wall, spanning what the detector can see of it
        targets = np.linspace(floor_width_m / 2.0 + run * 0.12,
                              floor_width_m / 2.0 + run * 0.98, 4)
        for xt, colour in zip(targets, SEQUENCE):
            zt = float(np.interp(xt, x, z))
            ax.plot([x_det, xt], [z_det, zt], color=colour, lw=1.4, zorder=4,
                    solid_capstyle="round")
            ax.plot(xt, zt, "o", color=colour, ms=4.6, zorder=5, mec="white", mew=0.8)
        lo = np.degrees(np.arctan((np.interp(targets[0], x, z) - z_det)
                                  / (targets[0] - x_det)))
        hi = np.degrees(np.arctan((np.interp(targets[-1], x, z) - z_det)
                                  / (targets[-1] - x_det)))

        ax.plot(x_det, z_det, "^", color=DETECTOR, ms=11, zorder=6, mec="white", mew=0.9)
        # Every label sits in the void or outside the rock, so nothing overlaps a ray
        ax.annotate("Array on the near wall\n(must be deployable)", (x_det, z_det),
                    (-rim * 1.30, depth_m * 1.13), fontsize=8.4,
                    color=DETECTOR, ha="left", va="center",
                    arrowprops=dict(arrowstyle="->", color=DETECTOR, lw=0.9,
                                    connectionstyle="arc3,rad=0.25"))
        ax.annotate("Far wall:\nwhere the tau exits\n(must be steep)",
                    (targets[-1], float(np.interp(targets[-1], x, z))),
                    (rim * 0.98, depth_m * 0.52), fontsize=8.4, color=INK, ha="right",
                    va="center",
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9,
                                    connectionstyle="arc3,rad=0.25"))
        ax.annotate("", (-rim, depth_m * 1.16), (rim, depth_m * 1.16),
                    arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.9))
        ax.text(0, depth_m * 1.20, f"rim to rim  $\\approx$ {2 * rim / 1000:.1f} km",
                ha="center", fontsize=8.6, color=MUTED)
        ax.text(0.0, depth_m * 0.12,
                f"arrival directions span $\\approx$ {_minus(lo, '{:.0f}')}$\\degree$ "
                f"to {hi:+.0f}$\\degree$",
                ha="center", fontsize=8.6, color=WINDOW, weight="bold")

        ax.set_xlabel("Horizontal distance (m)")
        ax.set_ylabel("Elevation above the valley floor (m)")
        ax.set_title(f"Across a canyon: {wall_slope_deg:.0f}$\\degree$ walls, "
                     f"{depth_m / 1000:.1f} km deep, drawn to scale",
                     loc="left", weight="bold")
        ax.set_ylim(-120, depth_m * 1.34)
        ax.set_aspect("equal", adjustable="box")
        _tidy(ax)
        fig.tight_layout()

    return fig


def decay_and_shower(energies_pev=(3.0, 10.0, 55.0, 100.0, 1000.0),
                     crossing_m=3000.0, figsize=(7.4, 3.9)):
    r"""
    Why a single energy cannot stand in for a spectrum.

    The tau must decay inside the gap for a shower to reach the detector, with
    probability :math:`1 - \exp(-d/L)` for a boosted decay length :math:`L`. Across a
    canyon that probability runs from essentially one to a few per cent over a single
    experiment's energy reach, which is why a capacity computed at one representative
    energy is an artefact of the energy chosen rather than a property of the terrain.

    Parameters
    ----------
    energies_pev : sequence of float, optional
        Energies to mark, in PeV.
    crossing_m : float, optional
        Gap the tau must decay within, in metres.
    figsize : tuple of float, optional
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from oroscope import figures
    >>> fig = figures.decay_and_shower()
    >>> len(fig.axes)
    1
    """
    from oroscope import physics

    grid = np.logspace(np.log10(1.0), np.log10(3000.0), 400)
    lengths = np.array([physics.tau_decay_length_m(e) for e in grid])
    prob = 1.0 - np.exp(-crossing_m / lengths)

    with _styled():
        fig, ax = plt.subplots(figsize=figsize)
        ax.semilogx(grid, prob, color=INK, lw=2.0, zorder=3)
        ax.fill_between(grid, 0, prob, color=WINDOW, alpha=0.10, lw=0, zorder=1)

        for energy, colour in zip(energies_pev, SEQUENCE + [WINDOW]):
            p = 1.0 - np.exp(-crossing_m / physics.tau_decay_length_m(energy))
            ax.plot(energy, p, "o", color=colour, ms=6, zorder=5, mec="white", mew=0.9)
            ax.annotate(f"{p:.2f}", (energy, p), textcoords="offset points",
                        xytext=(4, 7), fontsize=8.2, color=colour, weight="bold")

        ax.axvspan(3.0, 1000.0, color=MUTED, alpha=0.07, lw=0, zorder=0)
        ax.text(55, 1.04, "TAMBO's energy reach", ha="center", fontsize=8.4,
                color=MUTED)
        ax.set_xlabel("Tau energy (PeV)")
        ax.set_ylabel(f"P(decays within {crossing_m / 1000:.0f} km)")
        ax.set_title("A single energy cannot stand in for a spectrum",
                     loc="left", weight="bold")
        ax.set_ylim(0, 1.14)
        ax.set_xlim(1, 3000)
        _tidy(ax)
        fig.tight_layout()

    return fig


def pipeline_stages(figsize=(9.2, 5.4)):
    r"""
    How a DEM becomes a list of sites: the stages, and what each one removes.

    The vocabulary this project uses --- *screening*, *striding*, *the arrival scan*,
    *scoring*, *closing*, *pruning* --- is introduced nowhere in one place, and the
    terms are not guessable. This is that place, drawn.

    The widths are proportional to the survivors at each stage, on a logarithmic scale
    because the range is six orders of magnitude and a linear funnel would show one
    visible bar and six slivers. The numbers are a real run: TAMBO over the full Ancash
    DEM, 68.6 Mpx.

    Read it as two halves. Everything down to the arrival scan **removes** candidates;
    everything below **rebuilds a map from them**, which is why the count rises again at
    closing. Confusing those two halves is the single commonest way to misread a funnel
    table.

    The arrival scan and the scoring share one bar because that run cannot separate
    them. Its funnel recorded the post-cut count under both names --- the defect fixed
    in :func:`~oroscope.site_searcher.run_arrival_scan` --- so drawing them as two
    stages showed a scoring bar exactly as wide as the one above it, which asserts that
    the ``min_score`` cut removed nothing. It removed a great deal. A run made after the
    fix records the two counts separately and can be drawn as separate stages; this one
    is quoted as stored rather than re-run, so it is drawn as what it actually measured.

    Parameters
    ----------
    figsize : tuple of float, optional
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from oroscope import figures
    >>> fig = figures.pipeline_stages()
    >>> len(fig.axes)
    1
    """
    stages = [
        ("DEM pixels", 68_571_090, "Every pixel of the elevation model."),
        ("Slope screen", 33_943_993, "Ground steep enough to stand the array on,\n"
                                     "and not so steep it cannot be built."),
        ("Striding", 6_788_807, "One surviving pixel in N is kept as a candidate.\n"
                                "Cost control, not a criterion."),
        ("Arrival scan\n+ scoring", 1_022_530,
         "Walk outward along several bearings: is there a\n"
         "target at the right range? Score what there is on\n"
         "named components in [0, 1], and cut at min_score."),
        ("Closing", 2_543_406, "Morphology fills the holes striding left.\n"
                               "The count RISES here."),
        ("Pruning + selection", 186_704, "Regions too small or too poor in detectors\n"
                                         "are dropped."),
    ]
    colours = [MUTED, ROCK_EDGE, ROCK_EDGE, WINDOW, DETECTOR, DETECTOR]

    with _styled():
        fig, ax = plt.subplots(figsize=figsize)
        top = np.log10(stages[0][1])
        centre, span = 1.30, 0.62          # bars live here; labels sit either side
        for i, ((name, n, note), colour) in enumerate(zip(stages, colours)):
            half = np.log10(max(n, 10)) / top * span
            y = -i
            ax.add_patch(Rectangle((centre - half, y - 0.30), 2 * half, 0.60,
                                   facecolor=colour, alpha=0.9, lw=0))
            ax.text(centre - span - 0.06, y, f"{name}\n{n:,}", ha="right",
                    va="center", fontsize=9, color=INK, linespacing=1.35)
            ax.text(centre + span + 0.06, y, note, ha="left", va="center",
                    fontsize=8.5, color=MUTED, linespacing=1.35)
            if i:
                ax.annotate("", xy=(centre, y + 0.30), xytext=(centre, y + 0.70),
                            arrowprops=dict(arrowstyle="-|>", color=RULE, lw=1.2))

        # Which half of the pipeline each stage belongs to.
        bracket = centre - span - 0.78
        ax.plot([bracket, bracket], [-0.35, -3.35], color=ROCK_EDGE, lw=2.0, alpha=0.6)
        ax.text(bracket - 0.05, -1.85, "Removes\ncandidates", ha="right", va="center",
                fontsize=9, color=ROCK_EDGE, weight="bold", linespacing=1.35)
        ax.plot([bracket, bracket], [-3.65, -5.35], color=DETECTOR, lw=2.0, alpha=0.6)
        ax.text(bracket - 0.05, -4.5, "Rebuilds\na map", ha="right", va="center",
                fontsize=9, color=DETECTOR, weight="bold", linespacing=1.35)

        ax.set_xlim(bracket - 0.95, centre + span + 2.05)
        ax.set_ylim(-6.0, 0.55)
        ax.axis("off")
        ax.text(centre, -5.85, "TAMBO over the full Ancash DEM. Bar widths are "
                               "logarithmic in the survivor count.",
                ha="center", va="center", fontsize=8.5, color=MUTED)
        fig.tight_layout()
    return fig


def striding_and_closing(stride=5, element_px=(3, 5, 9), figsize=(8.6, 3.1)):
    r"""
    Why the closing element has to outrun the gap that striding leaves.

    Striding keeps one surviving pixel in ``stride``, so the accepted set becomes a
    lattice of isolated marks. Morphological *closing* is what turns that back into a
    region --- but only if its structuring element is larger than the gap. Below the
    gap the marks never touch and the mask stays a scatter; above it the region
    reappears almost intact.

    The transition is **at** the gap and it is abrupt, not gradual. That is the whole
    content of the figure, and it is the mechanism behind a real 4.75x under-report of
    TAMBO's area --- and a 291x one on steeper ground, where the accepted strips are
    narrower still.

    Parameters
    ----------
    stride : int, optional
        Keeps every Nth surviving pixel; also the gap it leaves, in pixels.
    element_px : tuple of int, optional
        Closing element sizes to draw, in pixels. One below the gap, one at it, one
        above.
    figsize : tuple of float, optional
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from oroscope import figures
    >>> fig = figures.striding_and_closing()
    >>> len(fig.axes)
    5
    """
    from scipy.ndimage import binary_closing

    n = 120
    rng = np.random.default_rng(1)
    rr, cc = np.mgrid[0:n, 0:n]
    truth = np.abs(cc - 55 - 16 * np.sin(rr / 20.0)) < 11
    truth &= rng.random((n, n)) < 0.97

    strided = np.zeros_like(truth)
    strided[::stride, ::stride] = truth[::stride, ::stride]

    panels = [("Accepted, every pixel", truth, None),
              (f"Marked one pixel in {stride}", strided, None)]
    for k in element_px:
        panels.append((f"Closed, element {k} px",
                       binary_closing(strided, np.ones((k, k))), k))

    base = int(truth.sum())
    with _styled():
        fig, axes = plt.subplots(1, len(panels), figsize=figsize)
        for ax, (name, mask, element) in zip(axes, panels):
            ax.imshow(mask, cmap="Greens", vmin=0, vmax=1.45,
                      interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            for side in ax.spines.values():
                side.set_color(RULE)
            verdict = ("" if element is None
                       else "  ✗" if element < stride else "  ✓")
            ax.set_xlabel(f"{name}{verdict}\n{mask.sum() / base:.2f}× the accepted set",
                          fontsize=8.5)
        fig.tight_layout()
    return fig


def score_composition(cut=0.35, figsize=(7.8, 3.4)):
    r"""
    Why a threshold on a product of components sits on a cliff.

    Each component scores a candidate in [0, 1] against one named criterion --- depth,
    accepted solid angle, exit distance, and so on. They are combined by
    **multiplication**, so a candidate has to be good at everything, and the composed
    score of several components piles up near zero however good the terrain is.

    A cut placed in the middle of that pile is therefore not a mild preference: it is a
    cliff, and where it lands depends on how many components happen to be enabled.
    Adding a component moves every score down and so silently tightens the cut.

    Parameters
    ----------
    cut : float, optional
        Where ``min_score`` is placed, for illustration.
    figsize : tuple of float, optional
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from oroscope import figures
    >>> fig = figures.score_composition()
    >>> len(fig.axes)
    2
    """
    rng = np.random.default_rng(0)
    n = 40_000
    parts = [rng.beta(5.0, 2.0, n) for _ in range(6)]

    with _styled():
        fig, (ax, bx) = plt.subplots(1, 2, figsize=figsize)

        running = np.ones(n)
        for i, part in enumerate(parts, start=1):
            running = running * part
            ax.hist(running, bins=60, range=(0, 1), histtype="step",
                    color=SEQUENCE[min(i - 1, len(SEQUENCE) - 1)], lw=1.3,
                    density=True)
        ax.axvline(cut, color="#B02A25", lw=1.6)
        ax.text(cut + 0.02, ax.get_ylim()[1] * 0.92, f"min_score = {cut}",
                color="#B02A25", fontsize=9, va="top")
        ax.set_xlabel("Composed score, as components are multiplied in")
        ax.set_ylabel("Density")
        ax.set_xlim(0, 1)
        _tidy(ax)

        kept = []
        running = np.ones(n)
        for part in parts:
            running = running * part
            kept.append(100.0 * (running >= cut).mean())
        bx.plot(range(1, len(parts) + 1), kept, marker="o", color="#B02A25", lw=1.8)
        bx.set_xlabel("Components multiplied in")
        bx.set_ylabel("Above the cut (%)")
        bx.set_ylim(0, 105)
        _tidy(bx)
        fig.tight_layout()
    return fig
