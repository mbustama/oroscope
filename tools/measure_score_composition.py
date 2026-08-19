#!/usr/bin/env python
"""
Regenerates the measured constants behind ``figures.score_composition_measured``.

The figure is executed by the documentation build, so it cannot call the pipeline:
a builder that needs a DEM fails on every machine that has never run this project.
The constants are therefore embedded in ``figures.py``, and nothing in that file can
notice on its own when scoring has changed underneath them. **Run this and paste the
output over the ``_MEASURED_*`` block whenever the score components, their
composition, or the Colca configuration change.**

The run is not instrumented and ``src/`` is not modified: ``run_arrival_scan`` already
returns the per-candidate arrays for per-site aggregation, and this wraps it on the way
past to keep them.

Colca is the cheap region to do this on --- a 6M-pixel crop whose scan stages take a
few seconds, against the twenty-five minutes a department costs.

Usage
-----
    python tools/measure_score_composition.py [--keep-npz]
"""
import argparse
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))

import oroscope.site_searcher as ss  # noqa: E402

# Least restrictive first is imposed later by sorting on the median, so that the
# collapse is attributable to a named component rather than to how many there are.
COMPONENTS = ["score_depth", "score_distance", "score_shower",
              "score_decay", "score_footprint", "score_solid_angle"]
BINS = 60
CUT = 0.35


def capture(out_dir):
    """Runs the Colca TAMBO search and returns its per-candidate arrays."""
    captured = []
    real = ss.run_arrival_scan

    def capturing(*args, **kwargs):
        n_hits, obs = real(*args, **kwargs)
        captured.append({k: np.asarray(v).copy() for k, v in obs.items()
                         if isinstance(v, np.ndarray) and v.ndim == 1})
        return n_hits, obs

    ss.run_arrival_scan = capturing
    try:
        ss.run_from_config(os.path.join(REPO, "config", "tambo_colca_config.json"),
                           run_output_dir=out_dir,
                           dem_path=os.path.join(REPO, "input", "dem", "colca.tif"))
    finally:
        ss.run_arrival_scan = real

    if not captured:
        raise SystemExit("no scan calls captured; has run_arrival_scan been renamed?")
    keys = sorted(set().union(*[set(c) for c in captured]))
    return {k: np.concatenate([c[k] for c in captured if k in c]) for k in keys}


def summarise(arrays):
    """Bins the running products over the candidates the score cut actually acts on."""
    # `accepted` in the observables is the *post*-cut flag, so it would give a censored
    # sample in which everything is above the cut by construction. The population the
    # cut acts on is the geometrically accepted one, which is where the score is
    # non-zero -- 360,939 at Colca, matching the funnel's `directions accepted`.
    geo = arrays["score"] > 0
    values = {c: arrays[c][geo] for c in COMPONENTS}
    order = sorted(COMPONENTS, key=lambda c: -np.median(values[c]))

    running = np.ones(int(geo.sum()))
    counts, above, comp_med, run_med = [], [], [], []
    for name in order:
        comp_med.append(round(float(np.median(values[name])), 4))
        running = running * values[name]
        counts.append(np.histogram(running, bins=BINS, range=(0, 1))[0].tolist())
        above.append(round(100.0 * float((running >= CUT).mean()), 2))
        run_med.append(round(float(np.median(running)), 4))

    product = np.ones(int(geo.sum()))
    for name in COMPONENTS:
        product = product * values[name]
    drift = float(np.abs(product - arrays["score"][geo]).max())
    return order, counts, above, comp_med, run_med, int(geo.sum()), drift


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--keep-npz", action="store_true",
                    help="also write the raw per-candidate arrays beside the run")
    args = ap.parse_args()

    out_dir = os.path.join(REPO, "output", "colca_realscores")
    arrays = capture(out_dir)
    order, counts, above, comp_med, run_med, n, drift = summarise(arrays)

    if drift > 1e-12:
        raise SystemExit(f"components no longer reproduce 'score' (max drift {drift:.2e}). "
                         "COMPONENTS is out of date -- fix it before pasting anything.")
    if args.keep_npz:
        np.savez_compressed(os.path.join(out_dir, "components.npz"), **arrays)

    names = [c.replace("score_", "") for c in order]
    print(f"# components reproduce the composed score to {drift:.1e}")
    print(f"_MEASURED_N = {n}")
    print(f"_MEASURED_NAMES = {names!r}")
    print(f"_MEASURED_ABOVE_CUT = {above!r}")
    print(f"_MEASURED_COMPONENT_MEDIAN = {comp_med!r}")
    print(f"_MEASURED_RUNNING_MEDIAN = {run_med!r}")
    print("_MEASURED_COUNTS = (")
    for name, row in zip(names, counts):
        print(f"    # {name}")
        print("    (" + ", ".join(str(v) for v in row[:20]) + ",")
        print("     " + ", ".join(str(v) for v in row[20:40]) + ",")
        print("     " + ", ".join(str(v) for v in row[40:]) + "),")
    print(")")


if __name__ == "__main__":
    main()
