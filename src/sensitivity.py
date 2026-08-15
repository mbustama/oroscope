#!/usr/bin/env python
"""
One-at-a-time sensitivity sweep over a search's parameters.

A single site-search result is only as firm as the assumptions behind it, and several
of the criteria are choices rather than measurements: what fraction of shower maximum
still counts as a usable shower, how steep the far wall must be, which energy stands in
for a spectrum in the tau-decay term, where the score cut sits. This runs the pipeline
repeatedly, varying one parameter at a time about a baseline, and tabulates how much
each one moves the answer.

One-at-a-time rather than a full grid on purpose: the question here is "which
assumption is this result most sensitive to", which OAT answers directly and cheaply.
It does not capture interactions between parameters, and a result that hinges on an
interaction will not show up -- for that, sweep the pair explicitly.

    python sensitivity.py ../config/tambo_colca_config.json \\
        --sweep decay_energy_pev 3 10 55 100 1000 \\
        --sweep min_score 0.0 0.2 0.35 0.5
"""

import argparse
import contextlib
import glob
import io
import json
import os
import shutil
import sys
import tempfile
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("MPLBACKEND", "Agg")

import site_searcher as ss                       # noqa: E402


def _quiet():
    """Silences the pipeline's console output, including tqdm's stderr bars."""
    return contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO())


def run_once(config, out_dir, verbose=False):
    """
    Runs the pipeline once for a fully-resolved configuration.

    Parameters
    ----------
    config : dict
        Parameters as they appear in a config file. Keys beginning with ``_`` are
        treated as comments and dropped, matching the config files in ``config/``.
    out_dir : str
        Directory to write this run's outputs into.
    verbose : bool
        Leave the pipeline's own output on stdout rather than swallowing it.

    Returns
    -------
    dict
        The parsed results JSON, or None when the run produced nothing.
    """
    params = {k: v for k, v in config.items() if not k.startswith("_")}
    for drop in ("dem_path", "origin_lat", "origin_lon", "print_info",
                 "output_directory_base_with_given_json", "output_image_format"):
        params.pop(drop, None)
    params["generate_kml"] = False
    if params.get("score_weights") is not None:
        params["score_weights"] = ss.parse_score_weights(params["score_weights"])
    for tup in ("depth_band_gcm2", "grammage_band_gcm2", "distance_band_m"):
        if params.get(tup) is not None:
            params[tup] = tuple(params[tup])

    os.makedirs(out_dir, exist_ok=True)
    ctx = _quiet() if not verbose else (contextlib.nullcontext(), contextlib.nullcontext())
    with ctx[0], ctx[1]:
        ss.find_grand_regions_interactive(
            dem_path=config["dem_path"],
            origin_lat=config["origin_lat"], origin_lon=config["origin_lon"],
            run_output_dir=out_dir, **params)

    matches = glob.glob(os.path.join(out_dir, "grand_search_results_*.json"))
    if not matches:
        return None
    with open(matches[0]) as f:
        return json.load(f)


def summarise(results):
    """Reduces a results JSON to the few numbers a sensitivity table compares."""
    if results is None:
        return dict(sites=0, capacity=0, area_km2=0.0, accepted=0, kept=0, acceptance=0.0)
    r = results["results"]
    funnel = results.get("funnel", {})
    kept = next((v for k, v in funnel.items() if k.startswith("kept by stride")), 0)
    accepted = funnel.get("directions accepted", 0)
    return dict(
        sites=r["total_sites"],
        capacity=r["total_capacity"],
        area_km2=round(sum(s["area_km2"] for s in r["sites"]), 1),
        accepted=accepted,
        kept=kept,
        acceptance=(accepted / kept) if kept else 0.0,
    )


def _coerce(text):
    """Config values arrive as strings on the command line; keep JSON's types."""
    lowered = text.strip().lower()
    if lowered in ("none", "null"):
        return None
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        value = float(text)
    except ValueError:
        return text
    return int(value) if value.is_integer() and "." not in text and "e" not in lowered else value


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="baseline configuration to vary about")
    ap.add_argument("--sweep", nargs="+", action="append", metavar=("PARAM", "VALUE"),
                    required=True,
                    help="a parameter and the values to try, repeatable")
    ap.add_argument("--out", default=None, help="where to write the report (default: alongside output/)")
    ap.add_argument("--keep_runs", action="store_true",
                    help="keep each run's output directory instead of discarding it")
    args = ap.parse_args()

    with open(args.config) as f:
        base = json.load(f)

    tmp_root = tempfile.mkdtemp(prefix="sensitivity_")
    report = {"config": os.path.abspath(args.config), "baseline": {}, "sweeps": {}}

    try:
        t0 = time.perf_counter()
        print(f"baseline: {os.path.basename(args.config)}", flush=True)
        baseline = summarise(run_once(base, os.path.join(tmp_root, "baseline")))
        report["baseline"] = baseline
        print(f"   {baseline['sites']} sites, {baseline['capacity']:,} DUs, "
              f"{baseline['area_km2']:,} km², acceptance {baseline['acceptance']*100:.1f}%"
              f"   [{time.perf_counter()-t0:.0f}s]\n", flush=True)

        for spec in args.sweep:
            param, raw_values = spec[0], spec[1:]
            if param not in base:
                print(f"   WARNING: '{param}' is not in the baseline config; adding it")
            values = [_coerce(v) for v in raw_values]
            rows = []
            print(f"{param}  (baseline {base.get(param)!r})", flush=True)
            print(f"   {'value':>12} {'sites':>7} {'capacity':>10} {'area km²':>11}"
                  f" {'accept':>8} {'vs base':>9}")
            for value in values:
                cfg = dict(base)
                cfg[param] = value
                out = os.path.join(tmp_root, f"{param}_{value}")
                got = summarise(run_once(cfg, out))
                ratio = (got["capacity"] / baseline["capacity"]) if baseline["capacity"] else 0.0
                rows.append(dict(value=value, **got, capacity_ratio=ratio))
                print(f"   {str(value):>12} {got['sites']:>7} {got['capacity']:>10,}"
                      f" {got['area_km2']:>11,} {got['acceptance']*100:>7.1f}%"
                      f" {ratio:>8.2f}x", flush=True)
                if not args.keep_runs:
                    shutil.rmtree(out, ignore_errors=True)
            report["sweeps"][param] = rows
            print(flush=True)

        out_path = args.out or os.path.join(REPO_ROOT, "output", "sensitivity.json")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(report, f, indent=2)
        print(f"written to {out_path}")
    finally:
        if not args.keep_runs:
            shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
