#!/usr/bin/env python
"""
Runs the full Arequipa DEM locally and stores the small artefacts for the notebook.

Notebook 8 reads results rather than producing them, for one reason: GRAND over this
DEM takes about 25 minutes and CI executes every notebook on every push. A tutorial
that costs half an hour of compute per commit is not a tutorial, it is a bill. So the
expensive part runs here, on a machine that already has the DEM, and what it leaves
behind is a few hundred kilobytes of JSON that the notebook opens instantly.

(TAMBO is the cheap one at about a minute: its targets are 2-5 km away against GRAND's
10-40, so the profile walks are a fraction as long. The cost is GRAND's alone.)

**Run this when the configuration changes, not otherwise.** The stored results record
which commit and which parameters produced them, so a stale store is detectable rather
than merely suspected.

    python tools/run_arequipa_full.py --dry-run       # what it would cost, then stop
    python tools/run_arequipa_full.py                 # all three, for real
    python tools/run_arequipa_full.py --only grand    # one of them

``--dry-run`` starts nothing, writes nothing and touches no store. It reports the DEM
it would search, the pre-flight memory estimate against what the system reports free,
which searches would run, the wall time to expect for each, and where the artefacts
would land -- the five things worth knowing before committing an hour of a machine.
No memory cap is applied in a dry run, because nothing is allocated.

Three searches, all over the same DEM at the same ``downsample_factor`` so that their
masks are pixel-aligned and can be overlaid:

1. **GRAND alone** -- ``config/grand_arequipa_full.json``
2. **TAMBO alone** -- ``config/tambo_arequipa_full.json``
3. **The combination** -- ``combine_experiments`` over the two, giving joint, union and
   co-location.

What is stored, per run, in ``results/arequipa_full/``: the results JSON, the
provenance record and the explanation. Not the rasters -- a GeoTIFF of a 129 Mpx mask
is far too large for a repository, and the notebook does not need it. The full outputs
stay in ``output/``, which is gitignored.
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

DEM = os.path.join(REPO, "input", "dem", "arequipa_SRTMGL1.tif")
STORE = os.path.join(REPO, "results", "arequipa_full")

RUNS = {
    "grand": os.path.join(REPO, "config", "grand_arequipa_full.json"),
    "tambo": os.path.join(REPO, "config", "tambo_arequipa_full.json"),
}


def run_one(label, config_path, out_root, max_memory_gb=None):
    """
    Runs one configuration over the full DEM and returns its results dictionary.

    The configuration-to-pipeline translation used to be written out here as well, a
    third copy alongside ``main()`` and the sensitivity child. It is
    ``ss.run_from_config`` now, so a parameter added to the pipeline reaches this
    runner without anyone having to remember that this file exists.

    ``dem_path`` is overridden because the configs spell it relative to ``src/``.
    """
    overrides = {"dem_path": DEM}
    if max_memory_gb is not None:
        overrides["max_memory_gb"] = max_memory_gb

    out_dir = os.path.join(out_root, f"arequipa_full_{label}")
    print(f"\n=== {label.upper()} ===")
    print(f"config: {os.path.relpath(config_path, REPO)}")
    print(f"output: {os.path.relpath(out_dir, REPO)}")

    started = time.time()
    results = ss.run_from_config(config_path, run_output_dir=out_dir, **overrides)
    print(f"\n{label}: finished in {(time.time() - started) / 60:.1f} minutes")
    return results, out_dir


def store(label, out_dir):
    """Copies the small, readable artefacts into the committed results store."""
    os.makedirs(STORE, exist_ok=True)
    kept = []

    found = ss.find_results_json(out_dir)
    if found:
        shutil.copy(found, os.path.join(STORE, f"{label}_results.json"))
        kept.append(f"{label}_results.json")
    for name, target in (("provenance.json", f"{label}_provenance.json"),
                         ("explanation.txt", f"{label}_explanation.txt")):
        path = os.path.join(out_dir, name)
        if os.path.exists(path):
            shutil.copy(path, os.path.join(STORE, target))
            kept.append(target)
    return kept


def combine(out_root):
    """Overlays the two runs and stores the combined report."""
    dirs = [os.path.join(out_root, f"arequipa_full_{label}") for label in RUNS]
    missing = [d for d in dirs if not os.path.exists(d)]
    if missing:
        print(f"\ncannot combine: {missing} not present. Run both searches first.")
        return None

    out_dir = os.path.join(out_root, "arequipa_full_combined")
    print("\n=== COMBINED ===")
    argv = sys.argv
    sys.argv = ["combine_experiments.py", *dirs,
                "--labels", "GRAND", "TAMBO", "--out", out_dir]
    try:
        ce.main()
    finally:
        sys.argv = argv

    os.makedirs(STORE, exist_ok=True)
    report = os.path.join(out_dir, "combined_report.json")
    if os.path.exists(report):
        shutil.copy(report, os.path.join(STORE, "combined_report.json"))
        return ["combined_report.json"]
    return []


def write_manifest(kept):
    """
    Records what the store holds and what produced it.

    Without this a stale store is indistinguishable from a current one, and the whole
    point of storing rather than recomputing is that nobody looks again.
    """
    manifest = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "generated_by": "tools/run_arequipa_full.py",
        "dem": os.path.relpath(DEM, REPO),
        "configs": {k: os.path.relpath(v, REPO) for k, v in RUNS.items()},
        "files": sorted(kept),
        "note": ("Regenerate when a configuration changes. Notebook 8 reads these and "
                 "does not run the searches itself: grand takes about 25 minutes, "
                 "tambo about 1 (its targets are 2-5 km away against grand's 10-40, "
                 "so the profile walks are far shorter)."),
    }
    path = os.path.join(STORE, "manifest.json")
    with open(path, "w") as f:
        json.dump(manifest, f, indent=4)
    print(f"\nmanifest: {os.path.relpath(path, REPO)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--only", choices=sorted(RUNS),
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

    if not os.path.exists(DEM):
        raise SystemExit(
            f"the full Arequipa DEM is not here: {DEM}\n"
            f"fetch it with: oroscope-fetch-dem --open_topography_api_key YOUR_KEY")

    labels = [args.only] if args.only else list(RUNS)

    report = ss.preflight_memory(DEM, downsample_factor=4, candidate_stride=5,
                                 max_memory_gb=0, quiet=args.dry_run)
    if args.dry_run:
        print(f"DEM:       {os.path.relpath(DEM, REPO)}")
        print(f"estimate:  {report['estimate_gb']:.2f} GiB at downsample_factor 4")
        print(f"available: {report['available_gb']:.1f} GiB"
              if report["available_gb"] else "available: unknown")
        print(f"would run: {', '.join(labels)}, then combine")
        print("expected:  ~25 min for grand, ~1 min for tambo")
        print(f"store:     {os.path.relpath(STORE, REPO)}")
        return

    kept = []
    for label in labels:
        _, out_dir = run_one(label, RUNS[label], args.out, args.max_memory_gb)
        kept += store(label, out_dir)

    if not args.only:
        kept += combine(args.out) or []

    write_manifest(kept)
    print(f"stored {len(kept)} files in {os.path.relpath(STORE, REPO)}")
    print("Notebook 8 reads these. Re-execute it to refresh its outputs.")


if __name__ == "__main__":
    main()
