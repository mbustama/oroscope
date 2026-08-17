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

# Cap threads so benchmarking does not saturate a shared workstation
BENCH_CORES = 8

# Stage timings above this fraction of a slowdown are called out as regressions
REGRESSION_FACTOR = 1.30


def warm_up_jit(tmp):
    """
    Compiles the numba kernels before timing starts.

    Without this the first timed case absorbs the compilation and reports it as scan
    time, which is how an earlier measurement came to look like a large memory-locality
    effect. Everything now runs in-process under numba's own threads, so compiling here
    is enough; a real invocation still pays the cost once per run.
    """
    import numpy as np
    from _support import ss
    elevation = np.zeros((64, 64), dtype=np.float32)
    from oroscope import arrival_scan
    grid = ss.resolve_grid_geometry("nonexistent.tif", -15.6, cell_size_deg=1 / 3600)
    arrival_scan.scan(np.array([[32.0, 32.0, 90.0]]), elevation, grid,
                      n_azimuths=1, half_width_deg=0.0, max_range_m=500.0,
                      min_dist_km=0.0, max_dist_km=1.0)
    ss.count_grid_capacity(np.ones((8, 8), dtype=bool), 30.0, 30.0, 60.0, 1)
    ss.apply_poly_mask_numba(np.zeros(1), np.zeros(1),
                             np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]),
                             np.ones(1, dtype=bool))


def check_machine_is_quiet():
    """
    Warns when other work is competing for the CPU, and returns the load it saw.

    Timings are only comparable against a baseline measured under similar conditions;
    an unrelated job at full tilt can double a single-threaded stage and read as a
    regression in code that did not change.

    The value is returned so that ``--update`` can record the load *before* the run.
    Sampling it afterwards instead recorded the benchmark's own load — always around 3
    on this box — which made the field useless for the one thing it is there for.
    """
    try:
        load1 = os.getloadavg()[0]
    except (OSError, AttributeError):               # pragma: no cover
        return None
    if load1 > 1.0:
        print(f"   WARNING: 1-minute load average is {load1:.2f}. Other work is competing "
              f"for the CPU, so these timings are not comparable with a baseline taken "
              f"on an idle machine.\n")
    return load1


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
    results = run_pipeline(dem, out_dir, lat, lon, num_cores=BENCH_CORES, **params)
    return {
        "timings_sec": {k: round(v, 4) for k, v in results["timings_sec"].items()},
        "funnel": results["funnel"],
        "total_sites": results["results"]["total_sites"],
        "total_capacity": results["results"]["total_capacity"],
        "peak_rss_mb": round(max(peak_rss_mb(), before), 1),
    }


def best_of(name, builder, params, tmp, repeat):
    """
    Runs a case ``repeat`` times and keeps the *fastest* time for each stage.

    The minimum, not the mean, and for a measured reason. Two consecutive passes over
    identical code on this host disagreed by **48%** on ``arequipa_900/ray_tracing``
    (1.23 s then 1.82 s) and 38% on ``synthetic_1800/ray_tracing`` -- both far outside
    the 30% regression gate, which makes a single-pass baseline worse than useless: it
    bakes one sample of the noise in as if it were the cost.

    The noise is one-sided. Nothing makes a stage run faster than its true cost, while
    a great many things make it run slower -- a scheduler moving the thread from a
    4.6 GHz P-core to a 3.4 GHz E-core, another process waking, a cache evicted. So the
    minimum over several passes is the closest estimate of the real cost available, and
    it converges as ``repeat`` rises rather than wandering as the mean does.

    Note which cases actually needed it: ``arequipa_2500``, at 17 seconds, was stable to
    6.5% across the same two passes. It is the short stages that are unmeasurable
    singly, because scheduler placement is a fixed cost amortised over a longer run.

    Returns the fastest run's structure, with per-stage minima substituted and the
    observed spread recorded so a later reader can see what the machine could resolve.
    """
    runs = []
    for i in range(repeat):
        if repeat > 1:
            print(f"   pass {i + 1}/{repeat} ...", flush=True)
        runs.append(run_case(name, builder, params, tmp))

    best = dict(runs[0])
    stages = runs[0]["timings_sec"]
    best["timings_sec"] = {
        stage: round(min(r["timings_sec"][stage] for r in runs), 4)
        for stage in stages
    }
    if repeat > 1:
        best["spread_pct"] = {}
        for stage in stages:
            times = [r["timings_sec"][stage] for r in runs]
            lo, hi = min(times), max(times)
            best["spread_pct"][stage] = round((hi / lo - 1) * 100, 1) if lo > 0 else None
        best["repeat"] = repeat
    best["peak_rss_mb"] = round(max(r["peak_rss_mb"] for r in runs), 1)
    return best


# A stage may only be gated when the machine can resolve a change smaller than the
# gate. Half the gate is the working rule: if repeated passes over identical code
# already spread by that much, a "regression" at the gate is as likely to be the
# scheduler as the code.
RESOLVABLE_SPREAD = (REGRESSION_FACTOR - 1) * 100 / 2.0


def report(current, baseline):
    """
    Prints a stage-by-stage table, flagging slowdowns the machine can actually resolve.

    A stage whose baseline records a ``spread_pct`` at or above half the regression
    gate is reported but never failed on. On this host that silences
    ``synthetic_900/ray_tracing``, which spread **149.6%** over five identical passes --
    gating on it would fail builds at random while telling nobody anything.
    """
    regressed = []
    unresolvable = []
    for case, data in current.items():
        print(f"\n{case}   ({data['total_sites']} sites, {data['total_capacity']} DUs,"
              f" peak RSS {data['peak_rss_mb']} MiB)")
        print(f"   {'stage':<22} | {'seconds':>9} | {'baseline':>9} | {'change':>9}")
        print("   " + "-" * 58)
        base_case = (baseline or {}).get(case, {})
        base = base_case.get("timings_sec", {})
        spreads = base_case.get("spread_pct", {})
        for stage, secs in data["timings_sec"].items():
            was = base.get(stage)
            noisy = (spreads.get(stage) or 0) >= RESOLVABLE_SPREAD
            if was:
                ratio = secs / was
                change = f"{(ratio - 1) * 100:+8.1f}%"
                if noisy:
                    change += " ~"
                if ratio > REGRESSION_FACTOR and secs > 0.05:
                    line = f"{case}/{stage}: {was:.2f}s -> {secs:.2f}s"
                    if noisy:
                        unresolvable.append(
                            f"{line}  (baseline spread {spreads[stage]:.0f}%, "
                            f"not gated)")
                    else:
                        regressed.append(line)
            else:
                change = "        -"
            print(f"   {stage:<22} | {secs:>9.3f} | {(f'{was:.3f}' if was else '-'):>9} | {change:>9}")

    if unresolvable:
        print("\nSLOWER, BUT BELOW THIS MACHINE'S RESOLUTION (~ marks these stages):")
        for line in unresolvable:
            print(f"   {line}")

    # The baseline has always stored the funnel, the site count and the capacity, and
    # has never compared them -- so the *results* could drift while the timing table
    # went on looking healthy. They did: the fan-tiling change moved arequipa_2500 from
    # 6301 detector positions to 6294, and nothing here said so. A timing baseline whose
    # answers no longer match is measuring two different programs against each other,
    # which is worse than having no baseline, because it looks like one.
    #
    # Reported, never gated. This is a benchmark, and a changed answer is usually an
    # intended change of code rather than a fault -- it just has to be *seen*, and the
    # baseline refreshed deliberately with --update on a quiet machine.
    drifted = []
    for case, data in current.items():
        base_case = (baseline or {}).get(case)
        if not base_case:
            continue
        for key in ("total_sites", "total_capacity"):
            was, now = base_case.get(key), data.get(key)
            if was is not None and now is not None and was != now:
                drifted.append(f"{case}/{key}: {was:,} -> {now:,}")
        was_f, now_f = base_case.get("funnel") or {}, data.get("funnel") or {}
        for stage in sorted(set(was_f) | set(now_f)):
            if was_f.get(stage) != now_f.get(stage):
                drifted.append(f"{case}/funnel[{stage}]: "
                               f"{was_f.get(stage, '-')} -> {now_f.get(stage, '-')}")
    if drifted:
        print("\nRESULTS DIFFER FROM THE BASELINE (not a timing matter):")
        for line in drifted:
            print(f"   {line}")
        print("   The baseline was measured against different behaviour. Refresh it with"
              "\n   `--update --repeat 5` on an idle machine once the change is intended.")
    return regressed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true", help="rewrite baseline.json")
    ap.add_argument("--quick", action="store_true", help="skip the large cases")
    ap.add_argument("--repeat", type=int, default=1, metavar="N",
                    help="Run each case N times and keep the fastest time per stage. "
                         "The minimum is the honest estimate: timing noise is one-sided, "
                         "since nothing runs faster than its true cost. Two identical "
                         "passes on this host disagreed by 48%% on a short stage, so "
                         "--repeat 5 is the minimum worth using with --update.")
    args = ap.parse_args()

    if args.update and args.repeat < 3:
        print(f"   WARNING: --update with --repeat {args.repeat}. Identical code has "
              f"varied by 48% between consecutive passes on this host, which is beyond "
              f"the {int((REGRESSION_FACTOR - 1) * 100)}% regression gate -- a baseline "
              f"from one pass records the noise as if it were the cost. Use --repeat 5.\n")

    baseline = None
    if os.path.exists(BASELINE) and not args.update:
        with open(BASELINE) as f:
            baseline = json.load(f).get("cases")

    load_at_start = check_machine_is_quiet()

    tmp = tempfile.mkdtemp(prefix="oroscope_bench_")
    current = {}
    try:
        print("warming up JIT ...", flush=True)
        warm_up_jit(tmp)
        for name, builder, params, is_large in CASES:
            if args.quick and is_large:
                continue
            if name.startswith("arequipa") and not _support.have_real_dem():
                print(f"skipping {name}: real DEM not present under input/dem/")
                continue
            print(f"running {name} ...", flush=True)
            current[name] = best_of(name, builder, params, tmp, args.repeat)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    regressed = report(current, baseline)

    if args.update:
        import platform
        with open(BASELINE, "w") as f:
            json.dump({
                "note": ("Cold-cache runs, per-stage MINIMUM over host.repeat passes. "
                         "Regenerate with: python bench/benchmark.py --update --repeat 5. "
                         "The minimum, not the mean: timing noise here is one-sided, so "
                         "the fastest pass is the closest estimate of the real cost. "
                         "cases.*.spread_pct records what the machine could actually "
                         "resolve for each stage -- a stage whose spread approaches the "
                         "30% regression gate cannot be gated on this host, however many "
                         "passes are taken. Check host.load_average_1min_at_start too. "
                         "See docs/ROADMAP.md 6.12 and 6.37."),
                "host": {"platform": platform.platform(), "python": platform.python_version(),
                         "cpu_count": os.cpu_count(), "bench_cores": BENCH_CORES,
                         "repeat": args.repeat,
                         "load_average_1min_at_start": (round(load_at_start, 2)
                                                        if load_at_start is not None else None)},
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
