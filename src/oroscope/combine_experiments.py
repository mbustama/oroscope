#!/usr/bin/env python
"""
Combines the results of two or more experiment searches over the same ground.

GRAND and TAMBO ask the same structural question -- from this patch of ground, is
there a target surface at the right range, in the right direction, with the right
matter behind it? -- and differ in their numbers rather than their structure. So each
experiment is one run of the searcher with its own configuration, and combining them
is an overlay of the masks those runs produce.

Three questions get different answers, and all three are worth reporting:

  joint   terrain that satisfies *every* experiment at once. This is the co-location
          case: one site, one road, one power feed, two experiments.
  union   terrain that satisfies *any* of them. This is the coverage case: how much of
          the region is useful to the programme as a whole.
  each    what each experiment gets on its own, and what it would lose by being
          confined to the joint area.

The inputs must be pixel-aligned: same shape, same pixel size, same corner. That is
not a detail to paper over -- two runs on differently-cropped DEMs would silently
overlay the wrong ground -- so it is checked and refused rather than resampled.

    python combine_experiments.py ../output/grand_colca_config ../output/tambo_colca_config \\
        --labels GRAND TAMBO --out ../output/combined_colca
"""

from __future__ import annotations

import argparse
import glob
import json
import os

import numpy as np

__all__ = ["load_run", "check_alignment", "read_world_file",
           "pixel_area_km2", "capacity_of", "main",
           "dem_for_run", "relief_for_mask", "elevation_for_mask", "geographic_extent"]
import tifffile as tiff

from oroscope import explain as explain_mod
from oroscope import site_searcher as ss

# Matplotlib is only needed for the overview image. The backend is chosen in main(),
# NOT here: this module is imported by `import oroscope`, and forcing Agg at import
# reached into every caller's session -- it overrode the inline backend in a notebook,
# so figures were captured as nothing at all. A library must not decide how its user's
# figures are rendered. Trap 3, one level up from where it bit before.
import matplotlib
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.patches import Patch                 # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402


def geographic_extent(world: tuple[float, ...], shape: tuple[int, int]) -> tuple:
    """
    The ``(left, right, bottom, top)`` an overlay needs to sit in degrees.

    Without it ``imshow`` puts a map in *pixel* coordinates, which is why the overview
    used to carry no axes at all: pixel indices are not worth labelling, so the ticks
    were switched off rather than made meaningful.

    The world file's terms are the centre of the top-left pixel, so half a pixel is
    added back at each edge to give the outer bounds of the raster.

    Parameters
    ----------
    world : tuple of float
        The six affine terms, from :func:`read_world_file`.
    shape : tuple of int
        ``(rows, cols)`` of the raster.

    Returns
    -------
    tuple of float
        ``(left, right, bottom, top)`` in degrees, for ``imshow(extent=...)``.

    Examples
    --------
    >>> from oroscope import combine_experiments as ce
    >>> world = (0.01, 0.0, 0.0, -0.01, -72.0, -15.0)     # 0.01 deg pixels
    >>> left, right, bottom, top = ce.geographic_extent(world, (100, 200))
    >>> round(left, 3), round(right, 3)
    (-72.005, -70.005)
    >>> round(bottom, 3), round(top, 3)
    (-16.005, -14.995)
    """
    cell_x, _, _, cell_y, x0, y0 = world
    rows, cols = shape
    left = x0 - 0.5 * cell_x
    right = x0 + (cols - 0.5) * cell_x
    top = y0 - 0.5 * cell_y                          # cell_y is negative going south
    bottom = y0 + (rows - 0.5) * cell_y
    return (left, right, bottom, top)


def dem_for_run(run: dict) -> str | None:
    """
    The DEM a run searched, if it can still be found.

    The combiner works from the masks the runs wrote, so it never needed the DEM --
    which is why the overview lost its terrain. The path is recorded twice, in the
    results parameters and in ``provenance.json``, so it can be recovered when the file
    is still there and reported as absent when it is not.

    Parameters
    ----------
    run : dict
        A run, as :func:`load_run` returns.

    Returns
    -------
    str or None
        Absolute path to the DEM, or ``None`` when it is not recorded or not present.
    """
    results = run.get("results") or {}
    candidates = [(results.get("parameters") or {}).get("dem")]
    prov = os.path.join(run.get("dir", ""), "provenance.json")
    if os.path.exists(prov):
        try:
            with open(prov) as f:
                candidates.append((json.load(f).get("dem") or {}).get("path"))
        except (OSError, ValueError):                # pragma: no cover - defensive
            pass
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return None


def elevation_for_mask(dem_path: str, shape: tuple[int, int]) -> np.ndarray | None:
    """
    The run's DEM, read strided onto the mask grid.

    The mask is the DEM downsampled, so it is read strided to match rather than
    resampled. The float32 ``.npy`` cache the searcher leaves beside the DEM is
    memory-mapped when present, which makes a strided read of a 129 Mpx DEM cost almost
    nothing; otherwise the GeoTIFF is read whole.

    Parameters
    ----------
    dem_path : str
        Path to the DEM the run searched.
    shape : tuple of int
        ``(rows, cols)`` of the mask to match.

    Returns
    -------
    ndarray or None
        Elevation in metres with the given shape, or ``None`` if the DEM could not be
        read or is smaller than the mask.
    """
    cache = os.path.splitext(dem_path)[0] + ".npy"
    try:
        if os.path.exists(cache):
            elevation = np.load(cache, mmap_mode="r")
        else:
            elevation = tiff.imread(dem_path)
    except Exception:                                # pragma: no cover - unreadable DEM
        return None

    rows, cols = shape
    step = max(1, int(round(elevation.shape[0] / float(rows))))
    small = np.asarray(elevation[::step, ::step], dtype=np.float64)
    if small.shape[0] < rows or small.shape[1] < cols:
        return None
    return small[:rows, :cols]


def relief_for_mask(dem_path: str, shape: tuple[int, int]) -> np.ndarray | None:
    """
    Hillshaded relief for the mask grid, as a greyscale array in [0, 1].

    Shaded relief rather than a colour elevation ramp, deliberately: the overlay needs
    the colour for its categories, and a terrain colourmap underneath would compete
    with them for the reader's attention. Relief carries the shape of the ground
    without spending any hue on it.

    The mask is the DEM downsampled, so the DEM is read strided to match rather than
    resampled. The float32 ``.npy`` cache the searcher leaves beside the DEM is
    memory-mapped when present, which makes a strided read of a 129 Mpx DEM cost
    almost nothing; otherwise the GeoTIFF is read whole.

    Parameters
    ----------
    dem_path : str
        Path to the DEM the run searched.
    shape : tuple of int
        ``(rows, cols)`` of the mask to match.

    Returns
    -------
    ndarray or None
        Relief in [0, 1] with the given shape, or ``None`` if the DEM could not be
        read or does not divide onto the mask grid.
    """
    cache = os.path.splitext(dem_path)[0] + ".npy"
    try:
        if os.path.exists(cache):
            elevation = np.load(cache, mmap_mode="r")
        else:
            elevation = tiff.imread(dem_path)
    except Exception:                                # pragma: no cover - unreadable DEM
        return None

    rows, cols = shape
    step = max(1, int(round(elevation.shape[0] / float(rows))))
    small = np.asarray(elevation[::step, ::step], dtype=np.float64)
    if small.shape[0] < rows or small.shape[1] < cols:
        return None
    small = small[:rows, :cols]

    finite = np.isfinite(small)
    if not finite.any():                             # pragma: no cover - empty DEM
        return None
    small = np.where(finite, small, np.nanmedian(small[finite]))

    shaded = matplotlib.colors.LightSource(azdeg=315, altdeg=45).hillshade(
        small, vert_exag=1.5, dx=1.0, dy=1.0)
    return np.clip(shaded, 0.0, 1.0)


def read_world_file(tfw_path: str) -> tuple[float, ...]:
    """
    Reads the six affine terms of an ESRI world file.

    Parameters
    ----------
    tfw_path : str
        Path to the ``.tfw`` file.

    Returns
    -------
    tuple of float
        ``(pixel_size_x, rot_y, rot_x, pixel_size_y, upper_left_x, upper_left_y)``.

    Raises
    ------
    ValueError
        If the file does not hold exactly six terms.
    """
    with open(tfw_path) as f:
        terms = [float(line.strip()) for line in f if line.strip()]
    if len(terms) != 6:
        raise ValueError(f"{tfw_path}: expected 6 terms, found {len(terms)}")
    return tuple(terms)


def load_run(run_dir: str) -> dict:
    """
    Loads one search's mask, georeferencing and results JSON.

    The searcher writes the mask as a downsampled GeoTIFF beside a world file and a
    results JSON, all sharing a base name.

    Parameters
    ----------
    run_dir : str
        A run's output directory.

    Returns
    -------
    dict
        ``dir``, ``tif``, ``mask``, ``world`` and ``results``. ``results`` is ``None``
        when no results JSON is present.

    Raises
    ------
    SystemExit
        If no mask GeoTIFF is found, or if the world file beside it is missing --
        without which alignment cannot be confirmed.
    """
    # By prefix, not alphabetically. A directory re-run since the rename holds both
    # `oroscope_results_*.tif` and a stale `grand_search_results_*.tif`, and the
    # legacy name sorts first -- so taking tifs[0] silently overlaid an old mask
    # against a current one. Measured: it reported TAMBO at Colca as 44.5 km2 from a
    # superseded 48,663-pixel mask, against 83.6 km2 for the run it claimed to
    # describe. Nothing failed; the number was simply wrong.
    tif = None
    for prefix in (ss.RESULTS_PREFIX, ss.LEGACY_RESULTS_PREFIX):
        found = sorted(glob.glob(os.path.join(run_dir, prefix + "*.tif")))
        if found:
            tif = found[0]
            break
    if tif is None:
        # An unprefixed mask is still a mask; only fall back once the named ones fail.
        loose = sorted(glob.glob(os.path.join(run_dir, "*.tif")))
        if not loose:
            raise SystemExit(f"no mask GeoTIFF found in {run_dir}")
        tif = loose[0]
    base = os.path.splitext(tif)[0]
    tfw = base + ".tfw"
    if not os.path.exists(tfw):
        raise SystemExit(f"no world file beside {tif}; cannot confirm alignment")

    mask = tiff.imread(tif).astype(bool)
    world = read_world_file(tfw)

    results = None
    found = ss.find_results_json(run_dir)
    if found:
        with open(found) as f:
            results = json.load(f)

    return {"dir": run_dir, "tif": tif, "mask": mask, "world": world, "results": results}


def check_alignment(runs: list[dict]) -> None:
    """
    Refuses to overlay masks that do not describe the same ground.

    Comparing shapes is not enough: two crops of the same size taken from different
    corners would overlay cleanly and mean nothing. The world file's pixel size and
    upper-left corner are what actually pin the ground down.

    Parameters
    ----------
    runs : list of dict
        Loaded runs, from :func:`load_run`. The first is taken as the reference.

    Raises
    ------
    SystemExit
        If any run differs from the reference in shape, pixel size or corner. Refusing
        is deliberate: resampling would silently compare the wrong terrain.
    """
    ref = runs[0]
    for other in runs[1:]:
        if other["mask"].shape != ref["mask"].shape:
            raise SystemExit(
                f"masks are not the same shape: {ref['dir']} is {ref['mask'].shape}, "
                f"{other['dir']} is {other['mask'].shape}.\n"
                f"Both runs must use the same DEM crop and the same downsample_factor.")
        for i, name in enumerate(("pixel size x", "rot y", "rot x", "pixel size y",
                                  "upper-left x", "upper-left y")):
            a, b = ref["world"][i], other["world"][i]
            # Georeferencing is in degrees; 1e-9 deg is well under a millimetre
            if abs(a - b) > 1e-9:
                raise SystemExit(
                    f"masks do not cover the same ground: {name} is {a!r} in "
                    f"{ref['dir']} and {b!r} in {other['dir']}.")


def pixel_area_km2(world: tuple[float, ...], reference_latitude_deg: float) -> float:
    """
    Ground area of one mask pixel, in km^2.

    The world file is in degrees, so the east-west size shrinks with the cosine of the
    latitude. This uses the same convention as the searcher's own grid geometry.

    Parameters
    ----------
    world : tuple of float
        The six affine terms, from :func:`read_world_file`.
    reference_latitude_deg : float
        Latitude at which to evaluate the east-west pixel size, in degrees. Normally
        the centre of the map.

    Returns
    -------
    float
        Area of one pixel, in km^2.
    """
    cell_deg_x, _, _, cell_deg_y, _, _ = world
    km_y = abs(cell_deg_y) * 110.6
    km_x = abs(cell_deg_x) * 111.32 * np.cos(np.radians(reference_latitude_deg))
    return km_x * km_y


def capacity_of(results):
    """
    The capacity a run reported, or ``None`` when it did not report one.

    ``search_mode: single`` writes the string ``'N/A'`` rather than a number, so
    ``int()`` raises ``ValueError`` -- which was not caught, and took the whole
    combination down with it. A run that does not report a capacity is an ordinary
    case here, not an error: the overlay is about ground, not detectors.

    Parameters
    ----------
    results : dict or None
        A run's parsed results JSON, or ``None`` when the directory held none.

    Returns
    -------
    int or None
        The reported capacity, or ``None`` when the run did not report one.
    """
    try:
        return int(results["results"]["total_capacity"])
    except (TypeError, KeyError, ValueError):
        return None


def main():
    # Headless, because this is the command line and it writes a PNG to disk.
    # Set here rather than at import, so importing the module leaves a
    # notebook's own backend alone.
    matplotlib.use("Agg")
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+", help="output directories of the searches to combine")
    ap.add_argument("--labels", nargs="+", default=None,
                    help="name for each run, in the same order (default: directory names)")
    ap.add_argument("--out", default=None,
                    help="directory for the combined outputs (default: alongside the first run)")
    ap.add_argument("--mode", choices=["joint", "union", "all"], default="all",
                    help="which combination to write as the headline mask (default: all)")
    ap.add_argument("--require", nargs="+", default=None,
                    help="labels that must all be satisfied for the joint mask. Defaults to "
                         "every run, but allows a joint of a subset when combining three or more.")
    ap.add_argument("--no_image", action="store_true", help="skip the overview PNG")
    ap.add_argument("--reveal", action="store_true",
                    help="Also write the overlay as three frames -- first experiment, "
                         "second experiment, both -- for uncovering a result "
                         "progressively in a talk. Everything but the categories is "
                         "identical in all three, so nothing shifts between slides.")
    ap.add_argument("--roads", type=str, default=None,
                    help="GeoJSON of road geometry to draw as context, from "
                         "oroscope-fetch-roads. Access is the question a joint area "
                         "cannot answer on its own.")
    ap.add_argument("--settlements", type=str, default="auto",
                    help="Named places to mark on the overview: 'auto' (the default) "
                         "uses whichever curated list has points inside the map, or "
                         "give a preset name ('arequipa', 'lima') or 'none'.")
    ap.add_argument("--no_explain", action="store_false", dest="explain",
                    help="Skip the plain-language account of the overlay. It is printed "
                         "by default and saved as combination_explanation.txt: what each "
                         "experiment brings, how much ground they can share, and which "
                         "screening band decides that. Co-location is usually settled by "
                         "slope rather than by anything about the physics.")
    args = ap.parse_args()

    if len(args.run_dirs) < 2:
        raise SystemExit("combining needs at least two runs")

    runs = [load_run(d) for d in args.run_dirs]
    labels = args.labels or [os.path.basename(os.path.normpath(d)) for d in args.run_dirs]
    if len(labels) != len(runs):
        raise SystemExit(f"got {len(labels)} labels for {len(runs)} runs")

    check_alignment(runs)

    world = runs[0]["world"]
    top_lat = world[5]
    rows = runs[0]["mask"].shape[0]
    centre_lat = top_lat + 0.5 * rows * world[3]      # world[3] is negative going south
    px_km2 = pixel_area_km2(world, centre_lat)

    masks = {label: run["mask"] for label, run in zip(labels, runs)}

    required = args.require or labels
    unknown = [r for r in required if r not in masks]
    if unknown:
        raise SystemExit(f"--require names unknown labels: {unknown}. Known: {labels}")

    joint = np.logical_and.reduce([masks[r] for r in required])
    union = np.logical_or.reduce(list(masks.values()))

    report = {
        "runs": [],
        "pixel_area_km2": px_km2,
        "joint_requires": required,
        "joint": {"pixels": int(joint.sum()), "area_km2": float(joint.sum() * px_km2)},
        "union": {"pixels": int(union.sum()), "area_km2": float(union.sum() * px_km2)},
        "pairwise_overlap": {},
    }

    for label, run in zip(labels, runs):
        m = run["mask"]
        own = int(m.sum())
        shared = int(np.logical_and(m, joint).sum())
        report["runs"].append({
            "label": label,
            "dir": run["dir"],
            "pixels": own,
            "area_km2": float(own * px_km2),
            "reported_capacity": capacity_of(run["results"]),
            "reported_sites": (run["results"]["results"]["total_sites"]
                               if run["results"] else None),
            "area_in_joint_km2": float(shared * px_km2),
            "fraction_of_own_area_in_joint": (shared / own) if own else 0.0,
        })

    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            both = int(np.logical_and(masks[a], masks[b]).sum())
            either = int(np.logical_or(masks[a], masks[b]).sum())
            report["pairwise_overlap"][f"{a} & {b}"] = {
                "area_km2": float(both * px_km2),
                "jaccard": (both / either) if either else 0.0,
                "fraction_of_" + a: (both / masks[a].sum()) if masks[a].any() else 0.0,
                "fraction_of_" + b: (both / masks[b].sum()) if masks[b].any() else 0.0,
            }

    out_dir = args.out or os.path.join(os.path.dirname(os.path.normpath(args.run_dirs[0])),
                                       "combined")
    os.makedirs(out_dir, exist_ok=True)

    written = []
    to_write = {"joint": joint, "union": union}
    if args.mode != "all":
        to_write = {args.mode: to_write[args.mode]}
    for name, m in to_write.items():
        path = os.path.join(out_dir, f"combined_{name}.tif")
        tiff.imwrite(path, m.astype(np.uint8))
        tfw_path = os.path.join(out_dir, f"combined_{name}.tfw")
        with open(tfw_path, "w") as f:
            f.write("\n".join(f"{t:.12f}" for t in world) + "\n")
        written += [path, tfw_path]

    # A membership map: which combination of experiments each pixel satisfies
    code = np.zeros(runs[0]["mask"].shape, dtype=np.uint8)
    for bit, label in enumerate(labels):
        code |= (masks[label].astype(np.uint8) << bit)
    code_path = os.path.join(out_dir, "combined_membership.tif")
    tiff.imwrite(code_path, code)
    written.append(code_path)
    report["membership_encoding"] = {
        "description": "bit i is set where run i's mask is set",
        "bits": {str(i): label for i, label in enumerate(labels)},
    }

    if not args.no_image and len(labels) == 2:
        a, b = labels
        extent = geographic_extent(world, code.shape)
        # Size the figure to the ground it shows. set_aspect on its own is not enough:
        # attach_colorbar's divider pins the axes position, so matplotlib satisfies the
        # aspect by widening the *data* limits instead, padding the map with empty
        # degrees until it looks like the raster failed to fill its frame. Giving the
        # figure the right shape up front leaves nothing to adjust.
        _centre_lat = 0.5 * (extent[2] + extent[3])
        _aspect = 1.0 / np.cos(np.radians(_centre_lat))
        _ratio = abs(extent[3] - extent[2]) * _aspect / abs(extent[1] - extent[0])
        only_a = masks[a] & ~masks[b]
        only_b = masks[b] & ~masks[a]
        both = masks[a] & masks[b]
        joint_km2 = float(both.sum()) * px_km2

        dem_path = next((p for p in (dem_for_run(r) for r in runs) if p), None)
        elevation = elevation_for_mask(dem_path, code.shape) if dem_path else None
        base_alpha = 0.55 if elevation is not None else 1.0

        # A class covering most of the frame is outlined, not filled: a translucent
        # wash over 80% of a map turns the whole thing to mud and makes the other
        # classes look like terrain. Decided ONCE, from the full set, so a stage that
        # omits a class still draws and labels the rest exactly as the final figure
        # does -- otherwise the legend would change shape between slides.
        fill_limit = 0.40 * code.size
        layers = (("a", only_a, "#1B6CA8", base_alpha, 1),
                  ("b", only_b, "#C8621B", 0.95, 2),
                  ("both", both, "#E8189B", 1.0, 3))
        outlined = {colour for _, mask, colour, _, _ in layers
                    if mask.any() and mask.sum() > fill_limit and elevation is not None}

        def key(colour, label, alpha=1.0):
            """A filled patch, or a line when that class is outlined rather than filled."""
            if colour in outlined:
                return Line2D([0], [0], color=colour, lw=1.6, label=f"{label} (outline)")
            return Patch(facecolor=colour, alpha=alpha, label=label)

        places = ss.resolve_settlements(
            args.settlements, (extent[2], extent[3], extent[0], extent[1]))

        def to_deg(lat, lon):
            """This overlay's axes are already degrees."""
            return (lon, lat)

        def draw(visible, path):
            """
            Draws the overlay showing only ``visible`` categories, and saves it.

            Everything outside the categories -- the terrain, the colour bar, the
            roads, the towns, the scale bar, the legend -- is built identically every
            time. That is the point: these are frames of a progressive reveal for a
            talk, so nothing may shift as one is replaced by the next.
            """
            fig, ax = plt.subplots(
                figsize=(11.0, max(3.5, min(13.0, 11.0 * _ratio))))

            im = None
            if elevation is not None:
                vmin, vmax = ss.altitude_limits(elevation)
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
                # Light greys only, so the darkest ground still sits well behind the
                # overlay rather than competing with it for contrast.
                grey = matplotlib.colors.LinearSegmentedColormap.from_list(
                    "altitude_grey", ["#3C3C3C", "#F2F2F2"])
                shaded = matplotlib.colors.LightSource(azdeg=315, altdeg=45).shade(
                    np.nan_to_num(elevation, nan=vmin), cmap=grey, norm=norm,
                    blend_mode="soft", vert_exag=2.0, dx=1.0, dy=1.0)
                # Water is water, not the bottom of the ramp. And nodata is neither:
                # filling NaN with vmin before shading painted the empty corner of the
                # Arequipa DEM as dark low ground.
                water = elevation <= ss.SEA_LEVEL_M
                missing = ~np.isfinite(elevation)
                shaded[water & ~missing, :3] = matplotlib.colors.to_rgb(ss.WATER_COLOUR)
                shaded[missing, :3] = matplotlib.colors.to_rgb(ss.NODATA_COLOUR)
                ax.imshow(shaded, extent=extent, origin="upper",
                          interpolation="bilinear", zorder=0)
                im = matplotlib.cm.ScalarMappable(norm=norm, cmap=grey)

            for name, mask, colour, alpha, order in layers:
                if name not in visible or not mask.any():
                    continue
                if colour in outlined:
                    ax.contour(mask.astype(float), levels=[0.5], colors=[colour],
                               linewidths=1.4, extent=extent, origin="upper",
                               zorder=order)
                    continue
                rgba = np.zeros(code.shape + (4,), dtype=np.float64)
                rgba[..., :3] = matplotlib.colors.to_rgb(colour)
                rgba[..., 3] = np.where(mask, alpha, 0.0)
                ax.imshow(rgba, extent=extent, origin="upper",
                          interpolation="nearest", zorder=order)

            ax.set_xlabel("Longitude (°)")
            ax.set_ylabel("Latitude (°)")
            ax.tick_params(direction="out", length=4)
            # "auto" rather than set_aspect(1/cos(lat)): attach_colorbar's divider pins
            # the axes position, so a fixed aspect is satisfied by widening the data
            # limits rather than shrinking the box, padding the map with empty degrees.
            # The figure was shaped to the ground above, so letting the axes fill it
            # gives those proportions with nothing left over.
            ax.set_aspect("auto")
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
            ss.add_scale_bar(ax, 111.32 * np.cos(np.radians(_centre_lat)))
            ss.add_north_arrow(ax)

            roads_drawn = 0
            if args.roads:
                from oroscope import fetch_roads as fetch_roads_mod
                roads_drawn = ss.add_roads(
                    ax, fetch_roads_mod.load_roads(args.roads), to_deg)
                # No attribution on the figure: it belongs in the caption, where a
                # reader of the paper or the slide will see it. The loader still keeps
                # it inside the GeoJSON, so it is not lost with the picture.
            ss.add_settlements(ax, places, to_deg)

            if im is not None:
                ss.attach_colorbar(fig, ax, im, "Altitude (m)")

            handles = [
                key("#1B6CA8", f"{a} only — {masks[a].sum() * px_km2:,.0f} km²",
                    base_alpha),
                key("#C8621B", f"{b} only — {masks[b].sum() * px_km2:,.0f} km²", 0.95),
                key("#E8189B", f"Both — {joint_km2:,.1f} km²"),
            ]
            if roads_drawn:
                handles.append(Line2D([0], [0], color=ss.ROAD_COLOUR, lw=1.0,
                                      alpha=0.7, label="Roads"))
            ax.legend(handles=handles, framealpha=0.9, fontsize="small", ncol=4,
                      loc="lower left", bbox_to_anchor=(0.0, 1.01), borderaxespad=0.0,
                      columnspacing=1.6)
            fig.savefig(path, dpi=140, bbox_inches="tight")
            plt.close(fig)
            written.append(path)

        everything = {"a", "b", "both"}
        draw(everything, os.path.join(out_dir, "combined_overview.png"))
        if args.reveal:
            # Frames for uncovering a result progressively in a talk. The legend, the
            # axes, the colour bar and the roads are identical in all three, so only
            # the categories appear as one frame replaces the next.
            for suffix, visible in ((f"1_{a.lower()}", {"a"}),
                                    (f"2_{b.lower()}", {"b"}),
                                    ("3_both", everything)):
                draw(visible, os.path.join(out_dir, f"combined_overview_{suffix}.png"))

        if elevation is None:
            print(f"   {ss.C.WARN}{ss.Icon.WARN}No DEM found for either run, so the overview "
                  f"has no relief. The path is recorded in each run's results and "
                  f"provenance; it is only missing if the DEM has moved.{ss.C.RESET}")

    report_path = os.path.join(out_dir, "combined_report.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    written.append(report_path)

    # ---- console summary
    print(f"\nCombining {len(runs)} searches over {rows}x{runs[0]['mask'].shape[1]} "
          f"pixels of {px_km2:.4f} km² each\n")
    width = max(len(label) for label in labels + ["joint", "union"]) + 2
    print(f"   {'experiment'.ljust(width)} {'area km²':>12} {'sites':>7} {'capacity':>10}"
          f" {'in joint':>10}")
    print("   " + "-" * (width + 43))
    for entry in report["runs"]:
        # Formatted before the f-string, not inside it. Thousands grouping is invalid
        # for a string, so `{cap if cap is not None else '-':>10,}` raised ValueError
        # the moment a run had no capacity to report -- which is every search_mode
        # 'single' run, since those write 'N/A', and any directory missing its results
        # JSON. Combining crashed rather than printing a dash.
        cap = entry["reported_capacity"]
        cap_s = f"{cap:,}" if cap is not None else "-"
        sites = entry["reported_sites"]
        sites_s = f"{sites:,}" if sites is not None else "-"
        print(f"   {entry['label'].ljust(width)} {entry['area_km2']:>12,.1f}"
              f" {sites_s:>7}"
              f" {cap_s:>10}"
              f" {entry['fraction_of_own_area_in_joint']*100:>9.1f}%")
    print("   " + "-" * (width + 43))
    print(f"   {('joint (' + ' & '.join(required) + ')').ljust(width)} "
          f"{report['joint']['area_km2']:>12,.1f}")
    print(f"   {'union (any)'.ljust(width)} {report['union']['area_km2']:>12,.1f}")

    if report["pairwise_overlap"]:
        print("\n   co-location")
        for pair, stats in report["pairwise_overlap"].items():
            print(f"      {pair}: {stats['area_km2']:,.1f} km² shared, "
                  f"Jaccard {stats['jaccard']:.3f}")

    # The overlay, explained. Same reasoning as the search's own summary: the report
    # gives a joint area and a Jaccard index, and neither says why either is what it
    # is. Failures are reported and swallowed -- a summary that cannot be written is
    # not a reason to lose a combination that already succeeded.
    if args.explain:
        try:
            text = explain_mod.explain_combination(
                report, {label: run["results"] for label, run in zip(labels, runs)
                         if run.get("results")})
            print("\n" + text)
            explanation_path = os.path.join(out_dir, "combination_explanation.txt")
            with open(explanation_path, "w", encoding="utf-8") as f:
                f.write(text)
            written.append(explanation_path)
        except Exception as e:                       # pragma: no cover - defensive
            print(f"   could not compose the combination summary: {e}")

    print("\n   written:")
    for path in written:
        print(f"      {path}")


if __name__ == "__main__":
    main()
