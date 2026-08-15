"""
Shared test scaffolding.

The searcher is still a standalone script rather than an installed package
(packaging is phase 4), so tests put ``src/`` on the path themselves. Matplotlib is
forced to a headless backend before the import, since importing the searcher pulls in
pyplot.
"""

import contextlib
import io
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(REPO_ROOT, "src")
GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

import site_searcher as ss  # noqa: E402

# Real DEMs live under input/, which is gitignored. Tests that need them skip when absent.
REAL_DEM = os.path.join(REPO_ROOT, "input", "dem", "arequipa_SRTMGL1.tif")


def have_real_dem():
    return os.path.exists(REAL_DEM)


def updating_golden():
    """True when golden files should be rewritten rather than compared."""
    return os.environ.get("UPDATE_GOLDEN", "").lower() in ("1", "true", "yes")


@contextlib.contextmanager
def quiet():
    """Silences the pipeline's console output, including tqdm's stderr bars."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        yield out


def run_pipeline(dem_path, out_dir, origin_lat, origin_lon, **overrides):
    """
    Runs the full search and returns the parsed results JSON.

    Exercising the real entry point (rather than the internals) means the golden
    files also cover the output-writing stages.
    """
    params = dict(
        target_antennas=100,
        min_width_km=1.0,
        antenna_spacing_km=1.0,
        min_dist_km=3.0,
        max_dist_km=20.0,
        grid_type="hex",
        search_mode="distributed",
        min_sub_array_size=20,
        downsample_factor=2,
        tile_size=256,
        candidate_stride=5,
        num_cores=2,
        generate_kml=False,
        rfi_zones=None,
        region_name="test",
    )
    params.update(overrides)

    os.makedirs(out_dir, exist_ok=True)
    with quiet():
        ss.find_grand_regions_interactive(
            dem_path=dem_path, origin_lat=origin_lat, origin_lon=origin_lon,
            run_output_dir=out_dir, **params
        )

    import json
    found = ss.find_results_json(out_dir)
    if not found:
        raise AssertionError(f"pipeline produced no results JSON in {out_dir}")
    with open(found) as f:
        return json.load(f)


def summarize(results):
    """Reduces a results JSON to the stable quantities a golden file should pin."""
    return {
        "total_sites": results["results"]["total_sites"],
        "total_capacity": results["results"]["total_capacity"],
        "sites": [
            {
                "area_km2": s["area_km2"],
                "capacity_exact": s["capacity_exact"],
                "facing_direction": s["facing_direction"],
            }
            for s in results["results"]["sites"]
        ],
        "funnel": results["funnel"],
        "regions": results["regions"],
    }
