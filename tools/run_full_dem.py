#!/usr/bin/env python
"""
Runs a full DEM locally, for every experiment, and stores the artefacts a notebook reads.

Notebook 8 reads results rather than producing them, for one reason: GRAND over the
Arequipa DEM takes about 25 minutes and CI executes every notebook on every push. A
tutorial that costs half an hour of compute per commit is not a tutorial, it is a bill.
So the expensive part runs here, on a machine that already has the DEM, and what it
leaves behind is a few hundred kilobytes of JSON that the notebook opens instantly.

(TAMBO is the cheap one at about a minute: its targets are 2-5 km away against GRAND's
10-40, so the profile walks are a fraction as long. The cost is GRAND's alone.)

**Run this when the configuration changes, not otherwise.** The stored results record
which commit and which parameters produced them, so a stale store is detectable rather
than merely suspected.

    python tools/run_full_dem.py --dry-run              # what it would cost, then stop
    python tools/run_full_dem.py                        # Arequipa, all three
    python tools/run_full_dem.py --region ancash        # the same over Ancash
    python tools/run_full_dem.py --only grand           # one of them

This was ``run_arequipa_full.py``, which hard-coded one DEM, one store and one pair of
configs. It handles a region table now, because the second region asked for was not a
reason to copy two hundred lines.

``--dry-run`` starts nothing, writes nothing and touches no store. It reports the DEM
it would search, the pre-flight memory estimate against what the system reports free,
which searches would run, the wall time to expect for each, and where the artefacts
would land -- the five things worth knowing before committing an hour of a machine.
No memory cap is applied in a dry run, because nothing is allocated.

Three searches per region, all over the same DEM at the same ``downsample_factor`` so
that their masks are pixel-aligned and can be overlaid:

1. **GRAND alone** -- ``config/grand_<region>_full.json``
2. **TAMBO alone** -- ``config/tambo_<region>_full.json``
3. **The combination** -- ``combine_experiments`` over the two, giving joint, union and
   co-location.

What is stored, per run, in ``results/<region>_full/``: the results JSON, the provenance
record and the explanation. Not the rasters -- a GeoTIFF of a 129 Mpx mask is far too
large for a repository, and the notebook does not need it. The full outputs stay in
``output/``, which is gitignored.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "src")
sys.path.insert(0, SRC)

from oroscope import combine_experiments as ce  # noqa: E402
from oroscope import site_searcher as ss  # noqa: E402

# One entry per region with a full-DEM pair of configurations. `expect` is wall time
# measured on this machine, not predicted: Arequipa's GRAND run is the 25-minute one.
REGIONS = {
    "arequipa": {
        "dem": os.path.join(REPO, "input", "dem", "arequipa_SRTMGL1.tif"),
        "expect": "~25 min for grand, ~1 min for tambo",
    },
    "ancash": {
        "dem": os.path.join(REPO, "input", "dem", "ancash_SRTMGL1.tif"),
        "expect": "~13 min for grand, under a minute for tambo (69 Mpx against 129)",
    },
}

EXPERIMENTS = ("grand", "tambo")


def region_paths(region):
    """DEM, config pair, store and output prefix for one region."""
    spec = REGIONS[region]
    return {
        "dem": spec["dem"],
        "expect": spec["expect"],
        "store": os.path.join(REPO, "results", f"{region}_full"),
        "prefix": f"{region}_full",
        "configs": {e: os.path.join(REPO, "config", f"{e}_{region}_full.json")
                    for e in EXPERIMENTS},
    }


def default_map_context(dem):
    """
    Road and place files for this DEM, when oroscope-fetch-roads has produced them.

    Both are optional context. A search does not need them and none of its numbers
    move; a map without them cannot say whether the good ground is reachable or what
    the nearest town is called, which is the question a site count leaves open.
    """
    stem = os.path.splitext(os.path.basename(dem))[0]
    roads = os.path.join(REPO, "input", "roads", f"{stem}.geojson")
    places = os.path.join(REPO, "input", "roads", f"{stem}_places.geojson")
    return (roads if os.path.exists(roads) else None,
            places if os.path.exists(places) else None)


def run_one(paths, label, out_root, max_memory_gb=None,
            roads=None, places=None):
    """
    Runs one configuration over the full DEM and returns its results dictionary.

    The configuration-to-pipeline translation used to be written out here as well, a
    third copy alongside ``main()`` and the sensitivity child. It is
    ``ss.run_from_config`` now, so a parameter added to the pipeline reaches this
    runner without anyone having to remember that this file exists.

    ``dem_path`` is overridden because the configs spell it relative to ``src/``.
    """
    overrides = {"dem_path": paths["dem"]}
    if max_memory_gb is not None:
        overrides["max_memory_gb"] = max_memory_gb
    if roads:
        overrides["roads_geojson"] = roads
    if places:
        overrides["settlements"] = places

    config_path = paths["configs"][label]
    out_dir = os.path.join(out_root, f"{paths['prefix']}_{label}")
    print(f"\n=== {label.upper()} ===")
    print(f"config: {os.path.relpath(config_path, REPO)}")
    print(f"output: {os.path.relpath(out_dir, REPO)}")

    started = time.time()
    results = ss.run_from_config(config_path, run_output_dir=out_dir, **overrides)
    print(f"\n{label}: finished in {(time.time() - started) / 60:.1f} minutes")
    return results, out_dir


def store(paths, label, out_dir):
    """Copies the small, readable artefacts into the committed results store."""
    os.makedirs(paths["store"], exist_ok=True)
    kept = []

    found = ss.find_results_json(out_dir)
    if found:
        shutil.copy(found, os.path.join(paths["store"], f"{label}_results.json"))
        kept.append(f"{label}_results.json")
    for name, target in (("provenance.json", f"{label}_provenance.json"),
                         ("explanation.txt", f"{label}_explanation.txt")):
        path = os.path.join(out_dir, name)
        if os.path.exists(path):
            shutil.copy(path, os.path.join(paths["store"], target))
            kept.append(target)
    return kept


def combine(paths, out_root, roads=None, places=None, reveal=True):
    """Overlays the two runs and stores the combined report."""
    dirs = [os.path.join(out_root, f"{paths['prefix']}_{label}")
            for label in EXPERIMENTS]
    missing = [d for d in dirs if not os.path.exists(d)]
    if missing:
        print(f"\ncannot combine: {missing} not present. Run both searches first.")
        return None

    out_dir = os.path.join(out_root, f"{paths['prefix']}_combined")
    print("\n=== COMBINED ===")
    argv = sys.argv
    sys.argv = ["combine_experiments.py", *dirs,
                "--labels", "GRAND", "TAMBO", "--out", out_dir]
    if roads:
        sys.argv += ["--roads", roads]
    if places:
        sys.argv += ["--settlements", places]
    if reveal:
        sys.argv += ["--reveal"]
    try:
        ce.main()
    finally:
        sys.argv = argv

    os.makedirs(paths["store"], exist_ok=True)
    report = os.path.join(out_dir, "combined_report.json")
    if os.path.exists(report):
        shutil.copy(report, os.path.join(paths["store"], "combined_report.json"))
        return ["combined_report.json"]
    return []


def write_manifest(paths, region, kept):
    """
    Records what the store holds and what produced it.

    Without this a stale store is indistinguishable from a current one, and the whole
    point of storing rather than recomputing is that nobody looks again.
    """
    manifest = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "generated_by": "tools/run_full_dem.py",
        "region": region,
        "dem": os.path.relpath(paths["dem"], REPO),
        "configs": {k: os.path.relpath(v, REPO)
                    for k, v in paths["configs"].items()},
        "files": sorted(kept),
        "note": ("Regenerate when a configuration changes. Notebook 8 reads these and "
                 "does not run the searches itself: grand takes about 25 minutes, "
                 "tambo about 1 (its targets are 2-5 km away against grand's 10-40, "
                 "so the profile walks are far shorter)."),
    }
    path = os.path.join(paths["store"], "manifest.json")
    with open(path, "w") as f:
        json.dump(manifest, f, indent=4)
    print(f"\nmanifest: {os.path.relpath(path, REPO)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--region", default="arequipa", choices=sorted(REGIONS),
                        help="which full DEM to search (default: arequipa)")
    parser.add_argument("--only", choices=sorted(EXPERIMENTS),
                        help="run just one of the searches")
    parser.add_argument("--out", default=os.path.join(REPO, "output"),
                        help="where the full outputs go (default: output/)")
    parser.add_argument("--max-memory-gb", type=float, default=None,
                        help="Address-space ceiling in GiB, passed to the pipeline. The "
                             "default (None) uses 80%% of what the system reports "
                             "available, which is not enough on a machine whose desktop "
                             "already holds half of RAM: this DEM needs about 6 GiB at "
                             "the scoring stage. 0 disables the cap, which risks the OOM "
                             "killer rather than a clean MemoryError -- prefer a number.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what the real run would cost, then stop without "
                             "starting a search, writing a file or touching the store. "
                             "Prints the DEM it would search, the pre-flight memory "
                             "estimate against what is free, which searches would run "
                             "(honouring --only), the wall time to expect per search, "
                             "and where the artefacts would land. No memory cap is "
                             "applied, since nothing is allocated.")
    args = parser.parse_args()

    paths = region_paths(args.region)
    dem = paths["dem"]
    if not os.path.exists(dem):
        raise SystemExit(
            f"the full {args.region} DEM is not here: {dem}\n"
            f"fetch it with: oroscope-fetch-dem --region {args.region} "
            f"--open_topography_api_key YOUR_KEY")
    missing = [c for c in paths["configs"].values() if not os.path.exists(c)]
    if missing:
        raise SystemExit("missing configuration(s): "
                         + ", ".join(os.path.relpath(m, REPO) for m in missing))

    labels = [args.only] if args.only else list(EXPERIMENTS)

    report = ss.preflight_memory(dem, downsample_factor=4, candidate_stride=5,
                                 max_memory_gb=0, quiet=args.dry_run)
    if args.dry_run:
        print(f"region:    {args.region}")
        print(f"DEM:       {os.path.relpath(dem, REPO)}")
        print(f"estimate:  {report['estimate_gb']:.2f} GiB at downsample_factor 4")
        print(f"available: {report['available_gb']:.1f} GiB"
              if report["available_gb"] else "available: unknown")
        print(f"would run: {', '.join(labels)}, then combine")
        print(f"expected:  {paths['expect']}")
        print(f"store:     {os.path.relpath(paths['store'], REPO)}")
        return

    roads, places = default_map_context(dem)
    for what, path in (("roads", roads), ("places", places)):
        print(f"map {what}: {os.path.relpath(path, REPO) if path else 'none found'}")

    kept = []
    for label in labels:
        _, out_dir = run_one(paths, label, args.out, args.max_memory_gb,
                             roads=roads, places=places)
        kept += store(paths, label, out_dir)

    if not args.only:
        kept += combine(paths, args.out, roads=roads, places=places) or []

    write_manifest(paths, args.region, kept)
    print(f"stored {len(kept)} files in {os.path.relpath(paths['store'], REPO)}")
    print(f"The {args.region} notebook reads these. Re-execute it to refresh its "
          f"outputs.")


if __name__ == "__main__":
    main()
