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

import argparse
import glob
import json
import os

import numpy as np

__all__ = ["load_run", "check_alignment", "read_world_file",
           "pixel_area_km2", "capacity_of", "main"]
import tifffile as tiff

# Matplotlib is only needed for the overview image, and the searcher already forces a
# headless backend when it is imported in a pipeline context
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.patches import Patch                 # noqa: E402


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
    tifs = sorted(glob.glob(os.path.join(run_dir, "*.tif")))
    if not tifs:
        raise SystemExit(f"no mask GeoTIFF found in {run_dir}")
    tif = tifs[0]
    base = os.path.splitext(tif)[0]
    tfw = base + ".tfw"
    if not os.path.exists(tfw):
        raise SystemExit(f"no world file beside {tif}; cannot confirm alignment")

    mask = tiff.imread(tif).astype(bool)
    world = read_world_file(tfw)

    results = None
    jsons = sorted(glob.glob(os.path.join(run_dir, "grand_search_results_*.json")))
    if jsons:
        with open(jsons[0]) as f:
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
    try:
        return int(results["results"]["total_capacity"])
    except (TypeError, KeyError):
        return None


def main():
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
        fig, ax = plt.subplots(figsize=(10, 7))
        # 0 neither, 1 a only, 2 b only, 3 both
        img = np.zeros(code.shape, dtype=np.uint8)
        img[masks[a] & ~masks[b]] = 1
        img[masks[b] & ~masks[a]] = 2
        img[masks[a] & masks[b]] = 3
        cmap = matplotlib.colors.ListedColormap(
            ["#EDEFEC", "#2C6E8F", "#B0781E", "#7B2D8E"])
        ax.imshow(img, cmap=cmap, vmin=0, vmax=3, interpolation="nearest")
        ax.set_title(f"{a} and {b}: where each is viable, and where both are")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(handles=[
            Patch(facecolor="#2C6E8F", label=f"{a} only"),
            Patch(facecolor="#B0781E", label=f"{b} only"),
            Patch(facecolor="#7B2D8E", label="both (co-located)"),
        ], loc="lower right", framealpha=0.9)
        png = os.path.join(out_dir, "combined_overview.png")
        fig.savefig(png, dpi=140, bbox_inches="tight")
        plt.close(fig)
        written.append(png)

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
        cap = entry["reported_capacity"]
        print(f"   {entry['label'].ljust(width)} {entry['area_km2']:>12,.1f}"
              f" {entry['reported_sites'] if entry['reported_sites'] is not None else '-':>7}"
              f" {cap if cap is not None else '-':>10,}"
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

    print("\n   written:")
    for path in written:
        print(f"      {path}")


if __name__ == "__main__":
    main()
