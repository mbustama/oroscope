#!/usr/bin/env python
"""
Benchmark harness: per-stage wall time and peak memory on fixed inputs.

Establishes the baseline that phase 3's optimisation work is measured against. Runs
are cold — the cached .npy is removed first — so numbers are comparable between
invocations rather than depending on what a previous run left behind.

    python bench/benchmark.py              # run and compare against the baseline
    python bench/benchmark.py --update     # rewrite the baseline
    python bench/benchmark.py --quick      # skip the large cases
"""

import argparse
import json
import os
import resource
import shutil
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The searcher is not yet an installed package (phase 4); reuse the test scaffolding
# so the benchmark drives the pipeline exactly the way the tests do.
sys.path.insert(0, os.path.join(REPO_ROOT, "tests"))

import _support                                    # noqa: E402
from _support import run_pipeline                  # noqa: E402
import synthetic                                   # noqa: E402

BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "baseline.json")
ORIGIN_LAT, ORIGIN_LON = -15.6, -72.3

# Stage timings above this fraction of a slowdown are called out as regressions
REGRESSION_FACTOR = 1.30


def warm_up_jit(tmp, n_jobs):
    """
    Compiles the numba kernels before timing starts, in the processes that will use them.

    The ray-caster runs inside joblib workers, so compiling it in the parent does not
    help: the first timed case would still absorb ~0.9 s of worker-side compilation and
    report it as ray-tracing. Driving one tiny job through the real parallel path
    compiles it where it is needed. Note that a real invocation always pays this cost
    once — joblib's pool only persists within a single process.
    """
    import numpy as np
    from _support import ss
    elevation = np.zeros((64, 64), dtype=np.float32)
    buf = np.lib.format.open_memmap(os.path.join(tmp, "warm.npy"), mode="w+",
                                    shape=(64, 64), dtype=bool)
    with _support.quiet():
        ss.run_ray_tracing_parallel(np.array([[32.0, 32.0, 90.0]]), elevation,
                                    30.7, 29.7, 64, 64, 200.0, 1.0, 2.0, n_jobs, buf)
    del buf
    ss.count_grid_capacity(np.ones((8, 8), dtype=bool), 2, 2, 1)
    ss.apply_poly_mask_numba(np.zeros(1), np.zeros(1),
                             np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]),
                             np.ones(1, dtype=bool))


def check_machine_is_quiet():
    """
    Warns when other work is competing for the CPU.

    Timings are only comparable against a baseline measured under similar conditions;
    an unrelated job at full tilt can double a single-threaded stage and read as a
    regression in code that did not change.
    """
    try:
        load1 = os.getloadavg()[0]
    except (OSError, AttributeError):               # pragma: no cover
        return
    if load1 > 1.0:
        print(f"   WARNING: 1-minute load average is {load1:.2f}. Other work is competing "
              f"for the CPU, so these timings are not comparable with a baseline taken "
              f"on an idle machine.\n")


def peak_rss_mb():
    """High-water memory for this process and its children, in MiB."""
    me = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    kids = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return max(me, kids) / 1024.0        # ru_maxrss is KiB on Linux


def make_synthetic(tmp, n):
    cell_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
    z = synthetic.ridge_and_slope(n, cell_x)
    path = os.path.join(tmp, f"ridge_{n}.tif")
    return synthetic.write_geotiff(path, z, ORIGIN_LAT, ORIGIN_LON), ORIGIN_LAT, ORIGIN_LON


def make_real_crop(tmp, n, r0=2500, c0=4000):
    import tifffile as tiff
    with tiff.TiffFile(_support.REAL_DEM) as tf:
        page = tf.pages[0]
        scale = page.tags["ModelPixelScaleTag"].value
        tie = page.tags["ModelTiepointTag"].value
        data = page.asarray()
    crop = data[r0:r0 + n, c0:c0 + n]
    lat = tie[4] - r0 * scale[1]
    lon = tie[3] + c0 * scale[0]
    path = os.path.join(tmp, f"crop_{n}.tif")
    return synthetic.write_geotiff(path, crop, lat, lon, cell_size_deg=scale[1]), lat, lon


CASES = [
    ("synthetic_900",   lambda t: make_synthetic(t, 900),   dict(tile_size=256), False),
    ("synthetic_1800",  lambda t: make_synthetic(t, 1800),  dict(tile_size=512), True),
    ("arequipa_900",    lambda t: make_real_crop(t, 900),   dict(tile_size=256), False),
    ("arequipa_2500",   lambda t: make_real_crop(t, 2500),  dict(tile_size=1024), True),
]


def run_case(name, builder, params, tmp):
    dem, lat, lon = builder(tmp)
    # Cold run: drop the cached .npy so the load stage is measured, not skipped
    cache = dem.replace(".tif", ".npy")
    if os.path.exists(cache):
        os.remove(cache)

    out_dir = os.path.join(tmp, f"out_{name}")
    before = peak_rss_mb()
    results = run_pipeline(dem, out_dir, lat, lon, **params)
    return {
        "timings_sec": {k: round(v, 4) for k, v in results["timings_sec"].items()},
        "funnel": results["funnel"],
        "total_sites": results["results"]["total_sites"],
        "total_capacity": results["results"]["total_capacity"],
        "peak_rss_mb": round(max(peak_rss_mb(), before), 1),
    }


def report(current, baseline):
    """Prints a stage-by-stage table, flagging slowdowns against the baseline."""
    regressed = []
    for case, data in current.items():
        print(f"\n{case}   ({data['total_sites']} sites, {data['total_capacity']} DUs,"
              f" peak RSS {data['peak_rss_mb']} MiB)")
        print(f"   {'stage':<22} | {'seconds':>9} | {'baseline':>9} | {'change':>9}")
        print("   " + "-" * 58)
        base = (baseline or {}).get(case, {}).get("timings_sec", {})
        for stage, secs in data["timings_sec"].items():
            was = base.get(stage)
            if was:
                ratio = secs / was
                change = f"{(ratio - 1) * 100:+8.1f}%"
                if ratio > REGRESSION_FACTOR and secs > 0.05:
                    regressed.append(f"{case}/{stage}: {was:.2f}s -> {secs:.2f}s")
            else:
                change = "        -"
            print(f"   {stage:<22} | {secs:>9.3f} | {(f'{was:.3f}' if was else '-'):>9} | {change:>9}")
    return regressed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true", help="rewrite baseline.json")
    ap.add_argument("--quick", action="store_true", help="skip the large cases")
    args = ap.parse_args()

    baseline = None
    if os.path.exists(BASELINE) and not args.update:
        with open(BASELINE) as f:
            baseline = json.load(f).get("cases")

    check_machine_is_quiet()

    tmp = tempfile.mkdtemp(prefix="sitesearch_bench_")
    current = {}
    try:
        print("warming up JIT ...", flush=True)
        warm_up_jit(tmp, n_jobs=2)
        for name, builder, params, is_large in CASES:
            if args.quick and is_large:
                continue
            if name.startswith("arequipa") and not _support.have_real_dem():
                print(f"skipping {name}: real DEM not present under input/dem/")
                continue
            print(f"running {name} ...", flush=True)
            current[name] = run_case(name, builder, params, tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    regressed = report(current, baseline)

    if args.update:
        import platform
        with open(BASELINE, "w") as f:
            json.dump({
                "note": "Cold-cache runs. Regenerate with: python bench/benchmark.py --update",
                "host": {"platform": platform.platform(), "python": platform.python_version(),
                         "cpu_count": os.cpu_count()},
                "cases": current,
            }, f, indent=2, sort_keys=True)
        print(f"\nbaseline written to {BASELINE}")
        return 0

    if regressed:
        print("\nSLOWER THAN BASELINE:")
        for line in regressed:
            print(f"   {line}")
        return 1
    if baseline:
        print("\nno stage regressed beyond "
              f"{int((REGRESSION_FACTOR - 1) * 100)}% of baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
