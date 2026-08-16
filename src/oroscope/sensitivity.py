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

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

# Three levels up now: src/oroscope/sensitivity.py -> src/oroscope -> src -> repo
# Three levels up now: src/oroscope/sensitivity.py -> src/oroscope -> src -> repo
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oroscope import site_searcher as ss         # noqa: E402

__all__ = ["run_once", "summarise", "main"]


# The child process: reads a JSON of parameters, runs one search, exits. Kept as a
# string rather than a separate file so the tool stays self-contained.
_CHILD = r"""
import json, os, sys
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, sys.argv[3])
from oroscope import site_searcher as ss

with open(sys.argv[1]) as f:
    payload = json.load(f)

# run_from_config does the configuration-to-pipeline translation, the same one main()
# uses. This child used to splat the payload straight in, which meant rfi_zones reached
# the pipeline as the preset *name*: iterated character by character, every character
# failing the 'circle' test, and the point silently searching with no exclusion zones.
ss.run_from_config(payload, run_output_dir=sys.argv[2])
"""


def run_once(config, out_dir, verbose=False, max_memory_gb=None, timeout=3600):
    """
    Runs the pipeline once, in a subprocess.

    A subprocess rather than a function call, for two reasons that a sweep makes
    unavoidable. Memory is reclaimed completely between points: running ten searches in
    one process took 6.9 GB and was killed by the kernel, because matplotlib retains
    every figure it is not explicitly asked to close and the leak compounded. And one
    point that fails -- an impossible parameter, or a genuine out-of-memory -- reports
    a failed row rather than ending the sweep.

    The cost is a few seconds of Numba compilation per point, which is the right trade
    for a sweep that otherwise cannot finish.

    Parameters
    ----------
    config : dict
        Parameters as they appear in a config file. Keys beginning with ``_`` are
        treated as comments and dropped.
    out_dir : str
        Directory for this run's outputs.
    verbose : bool, optional
        Let the child's output through rather than capturing it.
    max_memory_gb : float, optional
        Address-space ceiling for the child.
    timeout : float, optional
        Seconds before the child is killed and the point reported as failed.

    Returns
    -------
    dict or None
        The parsed results JSON, or ``None`` when the run produced nothing. Read back
        from disk rather than returned directly: the pipeline hands its results to its
        own caller, and that caller is in another process here -- which is the whole
        point of running each point in one.
    """
    # Only the sweep's own overrides here; the translation itself belongs to
    # run_from_config, which the child calls.
    params = {k: v for k, v in config.items() if not k.startswith("_")}
    params.pop("output_image_format", None)
    params["generate_kml"] = False
    # A sweep point is a row in a table, not something anyone reads: skip the
    # per-run summary rather than writing eighty of them nobody opens.
    params["explain"] = False
    if max_memory_gb:
        # The pipeline applies the cap itself now, so this is an ordinary parameter
        # rather than something the child had to remember to do first.
        params["max_memory_gb"] = max_memory_gb

    os.makedirs(out_dir, exist_ok=True)
    payload = os.path.join(out_dir, "_params.json")
    with open(payload, "w") as f:
        json.dump(params, f)

    src_dir = os.path.join(REPO_ROOT, "src")
    proc = subprocess.run(
        [sys.executable, "-c", _CHILD, payload, out_dir, src_dir],
        cwd=src_dir, timeout=timeout,
        capture_output=not verbose, text=True, check=False)
    if proc.returncode != 0:
        detail = (proc.stderr or "").strip().splitlines()[-1:] if not verbose else []
        print(f"      run failed (exit {proc.returncode})"
              + (f": {detail[0][:120]}" if detail else ""))
        return None

    found = ss.find_results_json(out_dir)
    if not found:
        return None
    with open(found) as f:
        return json.load(f)


def summarise(results: dict | None) -> dict:
    """
    Reduces a results JSON to the few numbers a sensitivity table compares.

    Parameters
    ----------
    results : dict or None
        A parsed results JSON, or ``None`` when the run produced nothing.

    Returns
    -------
    dict
        ``sites``, ``capacity``, ``area_km2``, ``accepted``, ``kept`` and
        ``acceptance``. All zero for a run that produced nothing, so a sweep row still
        appears rather than vanishing.
    """
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
    ap.add_argument("--max_memory_gb", type=float, default=None,
                    help="address-space ceiling for each child, in GiB "
                         "(default: 70%% of what the system reports available)")
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
        cap = args.max_memory_gb
        if cap is None:
            have = ss.available_memory_gb()
            cap = 0.7 * have if have else None
        if cap:
            print(f"   each run capped at {cap:.1f} GiB of address space\n")
        baseline = summarise(run_once(base, os.path.join(tmp_root, "baseline"),
                                      max_memory_gb=cap))
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
                got = summarise(run_once(cfg, out, max_memory_gb=cap))
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
