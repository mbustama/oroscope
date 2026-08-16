#!/usr/bin/env python
"""
Tabulates every full-DEM region against every other, and writes it down.

The per-region results live in ``results/<region>_full/``, one store each, and a store
answers questions about *its own* region only. The interesting questions are
comparative — is this ground better or worse, and for which experiment — and nothing
was holding the answer except a roadmap entry that goes stale the moment a region is
added.

So this reads the stores and regenerates ``results/region_comparison.md``. Run it after
``run_full_dem.py`` finishes a new region.

    python tools/compare_regions.py                 # rewrite the comparison
    python tools/compare_regions.py --print         # to stdout, write nothing

**Everything is normalised per DEM pixel**, because the DEMs are different sizes and a
bigger box trivially finds more ground. A ratio near the pixel ratio means "the same
terrain, less of it"; a departure from it is the ground talking.

The terrain block needs the DEMs, which are gitignored. Without them that section is
omitted and the rest still writes.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))

import numpy as np                                          # noqa: E402

from oroscope import site_searcher as ss                    # noqa: E402
from oroscope.fetch_dem import REGIONS as DEM_REGIONS       # noqa: E402

EXPERIMENTS = ("grand", "tambo")
OUT = os.path.join(REPO, "results", "region_comparison.md")

# The reference every ratio is quoted against. Arequipa is the region every published
# number in this project came from, so it is the sensible baseline rather than the
# first one alphabetically.
BASELINE = "arequipa"


def stores():
    """Every ``results/<region>_full/`` that holds at least one experiment."""
    root = os.path.join(REPO, "results")
    found = {}
    for name in sorted(os.listdir(root)) if os.path.isdir(root) else []:
        if not name.endswith("_full"):
            continue
        region = name[: -len("_full")]
        path = os.path.join(root, name)
        if any(os.path.exists(os.path.join(path, f"{e}_results.json"))
               for e in EXPERIMENTS):
            found[region] = path
    return found


def read(store, label):
    path = os.path.join(store, f"{label}_results.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def summarise(results):
    """Sites, area and capacity, plus the two funnel rows worth comparing."""
    if not results:
        return None
    r = results["results"]
    sites = r.get("sites") or []
    funnel = results.get("funnel") or {}
    # By name, never by position: a run with RFI zones carries an extra stage, so the
    # same index means different things in two regions. That is a real trap -- it made
    # GRAND's acceptance look like 20% at Arequipa when it is 61.6%.
    strided = next((v for k, v in funnel.items() if k.startswith("kept by stride")), None)
    accepted = funnel.get("directions accepted")
    band = next((v for k, v in funnel.items() if k.startswith("slope ")), None)
    return {
        "sites": r.get("total_sites"),
        "area_km2": sum(s.get("area_km2", 0.0) for s in sites),
        "capacity": r.get("total_capacity"),
        "dem_px": next(iter(funnel.values()), None),
        "slope_band_px": band,
        "strided": strided,
        "accepted": accepted,
        "acceptance": (100.0 * accepted / strided) if (accepted and strided) else None,
    }


def combination(store):
    path = os.path.join(store, "combined_report.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        report = json.load(f)
    pair = report.get("pairwise_overlap") or {}
    overlap = next(iter(pair.values()), {})
    return {"joint_km2": report.get("joint", {}).get("area_km2"),
            "jaccard": overlap.get("jaccard"),
            "share_of_tambo": 100.0 * overlap["fraction_of_TAMBO"]
            if "fraction_of_TAMBO" in overlap else None}


def terrain(region, block=2048):
    """Median slope and per-experiment band shares, over land, tile by tile."""
    spec = DEM_REGIONS.get(region)
    if not spec:
        return None
    tif = os.path.join(REPO, "input", "dem", spec["filename"])
    if not os.path.exists(tif):
        return None
    import tifffile

    grid = ss.resolve_grid_geometry(tif, spec["north"])
    z = tifffile.imread(tif)
    hist = np.zeros(9001)
    bands = {"grand_band": (3.0, 25.0), "tambo_band": (20.0, 60.0)}
    counts = dict.fromkeys(bands, 0)
    total = 0
    for r0 in range(0, z.shape[0], block):
        r1 = min(z.shape[0], r0 + block)
        lo, hi = max(0, r0 - 1), min(z.shape[0], r1 + 1)
        tile = z[lo:hi].astype(np.float32)
        dy, dx = np.gradient(tile, grid.cell_size_y, grid.cell_size_x)
        slope = np.degrees(np.arctan(np.hypot(dy, dx)))[r0 - lo:r1 - lo]
        s = slope[z[r0:r1] > 0]                     # the sea is exactly 0
        if not s.size:
            continue
        total += s.size
        hist += np.histogram(s, bins=9001, range=(0.0, 90.0))[0]
        for name, (a, b) in bands.items():
            counts[name] += int(((s >= a) & (s <= b)).sum())
    if not total:
        return None
    median = float((np.searchsorted(np.cumsum(hist), total / 2.0) + 0.5) / 100.0)
    return {"median_slope": median,
            **{k: 100.0 * v / total for k, v in counts.items()}}


def box_overlaps():
    """Pairs of regions whose bounding boxes intersect, and by how much."""
    out = []
    names = sorted(DEM_REGIONS)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ra, rb = DEM_REGIONS[a], DEM_REGIONS[b]
            lat_lo, lat_hi = max(ra["south"], rb["south"]), min(ra["north"], rb["north"])
            lon_lo, lon_hi = max(ra["west"], rb["west"]), min(ra["east"], rb["east"])
            if lat_lo >= lat_hi or lon_lo >= lon_hi:
                continue
            mid = (lat_hi + lat_lo) / 2.0
            km_ns = (lat_hi - lat_lo) * 111.32
            km_ew = (lon_hi - lon_lo) * 111.32 * math.cos(math.radians(mid))
            out.append({"a": a, "b": b, "km2": km_ns * km_ew,
                        "lat": (lat_lo, lat_hi), "lon": (lon_lo, lon_hi),
                        "datasets": (ra["demtype"], rb["demtype"])})
    return out


def row(values, fmt="{}"):
    return "| " + " | ".join(fmt.format(v) if v is not None else "—"
                             for v in values) + " |"


def build():
    found = stores()
    if not found:
        return "No `results/<region>_full/` stores found. Run `tools/run_full_dem.py`.\n"

    data = {}
    for region, store in found.items():
        data[region] = {
            **{e: summarise(read(store, e)) for e in EXPERIMENTS},
            "combined": combination(store),
            "terrain": terrain(region),
        }

    regions = sorted(data, key=lambda r: (r != BASELINE, r))
    base = data.get(BASELINE)
    px = {r: (data[r]["grand"] or data[r]["tambo"] or {}).get("dem_px") for r in regions}

    L = []
    L.append("# Full-DEM regions, compared\n")
    L.append("*Generated by `tools/compare_regions.py` from the stores in "
             "`results/<region>_full/`. Regenerate it when a region is added; do not "
             "edit by hand.*\n")
    L.append("Every criterion is held fixed across regions except where a table says "
             "otherwise, so a difference here is a difference in the **ground**. "
             f"Ratios are against **{BASELINE}**, and the honest way to read them is "
             "**per DEM pixel** — a bigger box trivially finds more ground.\n")

    L.append("\n## The DEMs\n")
    L.append(row(["region", "grid", "Mpx", "vs " + BASELINE, "dataset"]))
    L.append(row(["---"] * 5))
    for r in regions:
        spec = DEM_REGIONS.get(r, {})
        mpx = (px[r] / 1e6) if px[r] else None
        ratio = (px[r] / px[BASELINE]) if (px[r] and px.get(BASELINE)) else None
        L.append(row([r, "1 arc-sec" if spec.get("demtype", "").endswith("GL1")
                      or spec.get("demtype") == "AW3D30" else "3 arc-sec",
                      f"{mpx:,.1f}" if mpx else None,
                      f"{ratio:.3f}×" if ratio else None,
                      spec.get("demtype")]))

    if any(data[r]["terrain"] for r in regions):
        L.append("\n## The terrain, before any search\n")
        L.append("Over land only. This predicts the result: GRAND wants ground gentle "
                 "enough to stand an array on, TAMBO wants the walls GRAND cannot use.\n")
        L.append(row(["region", "median slope", "in GRAND's 3–25°", "in TAMBO's 20–60°"]))
        L.append(row(["---"] * 4))
        for r in regions:
            t = data[r]["terrain"]
            L.append(row([r,
                          f"{t['median_slope']:.1f}°" if t else None,
                          f"{t['grand_band']:.1f}%" if t else None,
                          f"{t['tambo_band']:.1f}%" if t else None]))

    for exp in EXPERIMENTS:
        L.append(f"\n## {exp.upper()}\n")
        L.append(row(["region", "sites", "area km²", "capacity", "acceptance",
                      "area /px", "capacity /px"]))
        L.append(row(["---"] * 7))
        b = base.get(exp) if base else None
        for r in regions:
            s = data[r][exp]
            if not s:
                L.append(row([r, None, None, None, None, None, None]))
                continue
            scale = (px[r] / px[BASELINE]) if (px[r] and px.get(BASELINE)) else None
            def per(v, bv):
                if not (b and bv and scale):
                    return None
                return f"{v / bv / scale:.2f}×"
            L.append(row([r, f"{s['sites']:,}", f"{s['area_km2']:,.1f}",
                          f"{s['capacity']:,}",
                          f"{s['acceptance']:.1f}%" if s["acceptance"] else None,
                          per(s["area_km2"], b["area_km2"] if b else None),
                          per(s["capacity"], b["capacity"] if b else None)]))
        L.append("\n*Acceptance is `directions accepted / kept by stride`, read by "
                 "stage **name**: a run with RFI zones carries an extra funnel stage, "
                 "so the same index means different things in two regions.*\n")

    if any(data[r]["combined"] for r in regions):
        L.append("\n## Both at once\n")
        L.append(row(["region", "joint km²", "Jaccard", "share of TAMBO's mask",
                      "joint /px"]))
        L.append(row(["---"] * 5))
        for r in regions:
            c = data[r]["combined"]
            scale = (px[r] / px[BASELINE]) if (px[r] and px.get(BASELINE)) else None
            bc = base.get("combined") if base else None
            per = (f"{c['joint_km2'] / bc['joint_km2'] / scale:.2f}×"
                   if (c and bc and scale and bc.get("joint_km2")) else None)
            L.append(row([r,
                          f"{c['joint_km2']:,.1f}" if c else None,
                          f"{c['jaccard']:.5f}" if c and c.get("jaccard") else None,
                          f"{c['share_of_tambo']:.1f}%"
                          if c and c.get("share_of_tambo") else None,
                          per]))
        L.append("\n**The share of TAMBO's mask is the number to watch.** It has stayed "
                 "near 44% across regions whose terrain could hardly differ more, which "
                 "says the joint region is TAMBO-limited and that co-location costs "
                 "GRAND almost nothing (ROADMAP §6.47). A Jaccard index that moves "
                 "while that share does not is TAMBO's mask growing, not the two "
                 "experiments agreeing more.\n")

    overlaps = box_overlaps()
    if overlaps:
        L.append("\n## Where the boxes overlap\n")
        L.append("Regions are downloaded as **bounding boxes**, and departments are not "
                 "rectangles. Two boxes can cover the same ground — with different "
                 "datasets — so ground found in one region may be found again in "
                 "another, and some of it lies outside the department it is filed "
                 "under.\n")
        L.append(row(["regions", "overlap km²", "latitudes", "longitudes", "datasets"]))
        L.append(row(["---"] * 5))
        for o in overlaps:
            L.append(row([f"{o['a']} ∩ {o['b']}", f"{o['km2']:,.0f}",
                          f"{o['lat'][0]:.2f}…{o['lat'][1]:.2f}",
                          f"{o['lon'][0]:.2f}…{o['lon'][1]:.2f}",
                          " vs ".join(o["datasets"])]))
        L.append("\nMeasured consequence: **37% of Ancash's joint ground (27.5 km² of "
                 "75.2) lies south of −10.45°**, in the corner of its box that reaches "
                 "into Lima region. The largest joint patch in the whole Ancash run — "
                 "5.37 km² at −10.5866, −77.0729 — reverse-geocodes to **Gorgor, "
                 "Cajatambo, Lima**, not Ancash. File results by box, read them by "
                 "geography.\n")

    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--print", dest="to_stdout", action="store_true",
                    help="write to stdout and leave the file alone")
    args = ap.parse_args()

    text = build()
    if args.to_stdout:
        print(text)
        return 0
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write(text)
    print(f"wrote {os.path.relpath(OUT, REPO)} ({len(text):,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
