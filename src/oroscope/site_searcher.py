"""
The search itself: from a digital elevation model to sites with detector capacity.

Six stages, each streamed so that a DEM larger than memory is not a special case.
Terrain is **screened** by slope, aspect, altitude and exclusion zones; the survivors
are **scanned** over a fan of arrival directions (:mod:`arrival_scan`) and **scored**
against per-experiment criteria (:mod:`scoring`); the accepted mask is **cleaned**
morphologically, **labelled** into sites and packed with a detector lattice; and the
result is **written** as GeoTIFF, world file, KML, PNG and JSON, with a selection
funnel, a provenance record and a plain-language summary (:mod:`explain`).

Three things are worth knowing before reading further.

**GRAND and TAMBO are configurations, not code paths.** Adding an experiment means
writing a JSON file. Nothing that shapes a result is hard-coded; every criterion is a
parameter of :func:`find_grand_regions_interactive`, and the command line, the
configuration file and the library all reach the same function.

**The funnel is the diagnostic.** Every filter records how many pixels survived it, so
a search that returns little or nothing names the constraint responsible rather than
leaving it to be guessed. It is printed, stored in the results JSON, and read back by
:func:`explain.binding_constraint`.

**Geometry comes from the file.** Pixel size and the north-west corner are read from
the DEM's own GeoTIFF tags, because an origin typed by hand does not fail when it is
wrong -- it silently georeferences every output to the wrong ground. A supplied origin
that disagrees with the file is reported rather than honoured in silence.

The pipeline returns its results dictionary, so a caller does not have to find and
re-read the file it was just handed the path to.
"""

from __future__ import annotations

import argparse
import inspect
import sys
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon as MplPolygon
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
from tqdm import tqdm
import multiprocessing
from scipy.ndimage import (binary_dilation, binary_erosion, label,
                           sum as ndi_sum, mean as ndi_mean, find_objects,
                           uniform_filter)
import os
import shutil
import math
import time
import json
import xml.etree.ElementTree as ET
from collections import namedtuple
from datetime import datetime
import re

from oroscope import arrival_scan
from oroscope import aperture as aperture_mod
from oroscope import explain as explain_mod
from oroscope import physics
from oroscope import scoring

# The public surface. Without this, autodoc documents every symbol the module imports
# -- Ellipse, FuncFormatter, tqdm, namedtuple -- as though they were ours.
__all__ = [
    "find_grand_regions_interactive", "main",
    "resolve_grid_geometry", "read_dem_geometry", "read_dem_origin",
    "resolve_origin", "build_elevation_cache",
    "load_dem_and_init_buffers",
    "terrain_gradients", "terrain_derivatives", "slope_band_gradient_sq",
    "slope_baseline_pixels", "get_candidates_chunked",
    "run_arrival_scan", "summarize_observables_by_site",
    "clean_shape_artifacts", "apply_morphology_pingpong",
    "separable_closing", "separable_opening",
    "analyze_sites_and_capacity", "count_grid_capacity",
    "create_world_file", "generate_kml_file", "generate_visualizations_and_outputs",
    "collect_provenance", "validate_parameters", "parse_score_weights",
    "explicitly_passed", "is_point_in_poly", "apply_poly_mask_numba",
    "Funnel", "MapGrid", "RESULTS_PREFIX", "LEGACY_RESULTS_PREFIX",
    "find_results_json",
    # Everything the command line can do, the library can do too: configuration
    # files, the memory pre-flight, and the run summary.
    "default_config", "generate_config", "load_config", "CONFIG_PRESETS",
    "estimate_peak_memory_gb", "apply_memory_cap", "available_memory_gb",
    "preflight_memory", "emit_explanation", "resolve_config_paths",
    "stride_gap_m", "closing_element_m", "warn_stride_outruns_closing", "add_scale_bar",
    "altitude_limits", "add_north_arrow", "SEA_LEVEL_M", "WATER_COLOUR", "NODATA_COLOUR",
]

# Try to import psutil for RAM stats
try:
    import psutil
except ImportError:
    psutil = None

# ==========================================
#          UI THEME & FORMATTING
# ==========================================
# Enables ANSI escape sequences in Windows 10+ terminals
if sys.platform == 'win32':
    os.system('') 

def supports_color():
    """Checks if the terminal supports ANSI colors."""
    supported_platform = sys.platform != 'win32' or 'ANSICON' in os.environ or 'WT_SESSION' in os.environ or os.environ.get('TERM') == 'xterm-256color'
    is_a_tty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
    return supported_platform and is_a_tty

def supports_emoji():
    """Checks if the terminal supports UTF-8 for emojis."""
    if sys.stdout.encoding:
        return sys.stdout.encoding.lower() == 'utf-8'
    return False

USE_COLOR = supports_color()
USE_EMOJI = supports_emoji()

class C:
    HEADER = '\033[96m' if USE_COLOR else ''  # Cyan
    OK = '\033[92m' if USE_COLOR else ''      # Green
    WARN = '\033[93m' if USE_COLOR else ''    # Yellow
    FAIL = '\033[91m' if USE_COLOR else ''    # Red
    BOLD = '\033[1m' if USE_COLOR else ''     # Bold
    MAGENTA = '\033[95m' if USE_COLOR else '' # Magenta
    RESET = '\033[0m' if USE_COLOR else ''    # Reset

class Icon:
    MAP = '🗺️  ' if USE_EMOJI else '[*] '
    GEAR = '⚙️  ' if USE_EMOJI else '[~] '
    BROOM = '🧹 ' if USE_EMOJI else '[C] '
    DISK = '💾 ' if USE_EMOJI else '[S] '
    WARN = '⚠️  ' if USE_EMOJI else '[!] '
    CHECK = '✅ ' if USE_EMOJI else '[✓] '
    INFO = 'ℹ️  ' if USE_EMOJI else '[i] '
    CROSS = '❌ ' if USE_EMOJI else '[x] '

# ==========================================
#               CONFIGURATION
# ==========================================
# Pre-defined Radio Frequency Interference (RFI) exclusion zones
AREQUIPA_RFI_ZONES = [
    ('circle', -16.409, -71.537, 25.0, "Arequipa"),
    ('circle', -16.264, -71.956, 10.0, "Majes"),
    ('circle', -16.533, -71.658, 15.0, "Cerro Verde"),
    ('circle', -16.480, -71.930, 8.0, "La Joya"),
    ('circle', -17.015, -72.015, 10.0, "Mollendo"),
]
LIMA_RFI_ZONES = [
    # 1. Metropolitan Area (Consolidated, Coastal)
    ('circle', -12.080, -77.010, 40.0, "Greater Lima Urban Area"),
    # 2. Inland Valleys (Rímac/Sierra Hubs)
    ('circle', -11.950, -76.680, 8.0, "Chosica"),
    ('circle', -11.850, -76.360, 6.0, "Matucana"),
    ('circle', -11.750, -76.220, 5.0, "San Mateo"),
    ('circle', -11.470, -76.630, 7.0, "Canta"),
    # 3. Provincial Coastal Hubs (North/South)
    ('circle', -11.100, -77.600, 12.0, "Huacho"),
    ('circle', -13.060, -76.380, 10.0, "Cañete"),
    ('circle', -11.080, -77.560, 6.0, "Huaura / Carquin"),
    ('circle', -10.750, -77.750, 7.0, "Barranca"),
    ('circle', -12.480, -76.650, 8.0, "Chilca"),
    ('circle', -11.480, -77.200, 8.0, "Huaral"),
    ('circle', -12.670, -76.620, 5.0, "Mala"),
    ('circle', -13.000, -76.350, 4.0, "Asia"),
    ('circle', -11.560, -77.270, 6.0, "Chancay"),
]

ORIGIN_LAT_AREQUIPA = -14.555380967667489
ORIGIN_LON_AREQUIPA = -73.58612537384033

ORIGIN_LAT_LIMA = -10.228479499469358
ORIGIN_LON_LIMA = -78.07665824890137

# Geodetic constants used to translate angular map units into physical distances.
# A degree of latitude is very nearly constant; 110.6 km is a good mid-latitude value.
KM_PER_DEG_LAT = 110.6
# A degree of longitude spans this at the equator, shrinking as cos(latitude).
KM_PER_DEG_LON_EQUATOR = 111.32
# Fallback pixel size (1 arc-second, i.e. SRTMGL1 / AW3D30) used only when the DEM
# carries no georeferencing tags and the user supplied no explicit value.
DEFAULT_CELL_SIZE_DEG = 1.0 / 3600.0

# ==========================================
#           NUMBA & PHYSICS KERNELS
# ==========================================
try:
    import numba
    from numba import jit, prange
    HAS_NUMBA = True
except ImportError:
    numba = None
    HAS_NUMBA = False
    def jit(*args, **kwargs):
        def decorator(func): return func
        return decorator

@jit(nopython=True, fastmath=True)
def is_point_in_poly(x, y, poly_verts):
    """
    Determines if a given 2D point lies inside a polygon using the Ray-Casting algorithm.
    Optimized for Numba execution.
    
    Parameters
    ----------
    x, y : float
        Coordinates of the test point.
    poly_verts : ndarray
        ``(M, 2)`` array of polygon vertices as ``(x, y)``.

    Returns
    -------
    bool
        ``True`` if the point is inside the polygon.
    """
    n = len(poly_verts)
    inside = False
    p1x, p1y = poly_verts[0]
    
    for i in range(n + 1):
        p2x, p2y = poly_verts[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
        
    return inside

@jit(nopython=True, parallel=True)
def apply_poly_mask_numba(valid_rows, valid_cols, poly_verts, mask_out):
    """
    Parallelized application of the polygon ray-casting check across an array of coordinates.
    Used for excluding regions defined by arbitrary polygonal RFI zones.
    
    Parameters
    ----------
    valid_rows, valid_cols : ndarray
        Row and column coordinates of the points to check.
    poly_verts : ndarray
        ``(M, 2)`` array of polygon vertices.
    mask_out : ndarray
        Boolean array modified in place; entries inside the polygon are cleared. Only
        ever clears bits, so several polygons can be applied in sequence.
    """
    n = len(valid_rows)
    for i in prange(n):
        if is_point_in_poly(valid_cols[i], valid_rows[i], poly_verts):
            mask_out[i] = False

# sin(60 deg): the row pitch of a triangular lattice, as a fraction of the spacing
_SIN_60 = 0.8660254037844386


@jit(nopython=True, fastmath=True)
def count_grid_capacity(mask_chunk, cell_size_y, cell_size_x, spacing_m, grid_type_code):
    """
    Counts the detectors that fit on the validated terrain at a given ground spacing.

    Detector positions are laid out in **metres on the ground** and only then looked up
    in the pixel grid. The earlier version did the reverse: it converted the spacing to
    an integer number of pixels and stepped the array by that stride, which truncated
    three separate times --- ``int()`` on the row stride, on the column stride, and on
    the hexagonal row pitch ``int(spacing_r * sin60)``. Every truncation shortens the
    spacing, so detectors ended up closer together than asked for and the count came out
    high: +7.4% at GRAND's 1 km, and **+58% at TAMBO's 100 m**, where only about three
    pixels span one separation on a 30 m DEM and the hex pitch collapsed from 2.6 to 2
    pixels. Placing points in continuous coordinates has no stride to truncate, so the
    count follows the requested geometry at any spacing.

    The layout is anchored at the chunk's own corner rather than fitted to it, so this
    is a capacity estimate for an arbitrarily-placed array, not the best packing
    achievable by sliding the grid around. That is the same convention as before.

    A spacing finer than the DEM's own pixels is permitted and yields several detectors
    per pixel. That is the honest continuum limit --- capacity is usable area divided by
    area per detector --- but note the terrain mask cannot resolve whether those
    sub-pixel positions really are usable.

    Parameters
    ----------
    mask_chunk : ndarray
        2D boolean array; ``True`` marks valid terrain.
    cell_size_y, cell_size_x : float
        Ground size of one pixel, in metres. They differ on a geographic grid, which is
        why an equal ground spacing is a different number of pixels on each axis.
    spacing_m : float
        Distance between neighbouring detectors, in metres. Zero or less returns 0.
    grid_type_code : int
        0 for a square grid, 1 for a hexagonal (triangular) one.

    Returns
    -------
    int
        Detectors fitting inside the valid terrain.

    Examples
    --------
    >>> import numpy as np
    >>> from oroscope import site_searcher as ss
    >>> mask = np.ones((100, 100), dtype=bool)          # 3 km square of 30 m pixels
    >>> ss.count_grid_capacity(mask, 30.0, 30.0, 1000.0, 1)
    12
    >>> ss.count_grid_capacity(mask, 30.0, 30.0, 0.0, 1)   # degenerate spacing
    0
    """
    h, w = mask_chunk.shape
    if spacing_m <= 0.0 or cell_size_y <= 0.0 or cell_size_x <= 0.0:
        return 0

    height_m = h * cell_size_y
    width_m = w * cell_size_x
    # A triangular lattice puts its rows sin(60) apart and offsets alternate rows by
    # half a spacing; that is what makes every neighbour distance equal to spacing_m.
    row_pitch = spacing_m * _SIN_60 if grid_type_code == 1 else spacing_m

    count = 0
    k = 0
    while True:
        # Indexed rather than accumulated, so the position of the thousandth row does
        # not depend on rounding in the nine hundred and ninety-nine before it
        y = k * row_pitch
        if y >= height_m:
            break
        r = int(y / cell_size_y)
        if r >= h:
            break
        x0 = 0.5 * spacing_m if (grid_type_code == 1 and k % 2 == 1) else 0.0
        j = 0
        while True:
            x = x0 + j * spacing_m
            if x >= width_m:
                break
            c = int(x / cell_size_x)
            if c < w and mask_chunk[r, c]:
                count += 1
            j += 1
        k += 1
    return count

# ==========================================
#            RUN ACCOUNTING
# ==========================================

class Funnel:
    """
    Records how many pixels survive each successive stage of the search.

    A search that returns nothing gives the user no clue which constraint was
    responsible. Every stage reports the count of pixels that passed it *and all
    preceding stages*, so the table reads as a funnel from the raw DEM down to the
    selected sites. Counts accumulate across tiles.
    """

    def __init__(self):
        self.stages = []          # ordered [name, count] pairs
        self._index = {}

    def add(self, name, count):
        """
        Adds to a stage's running total, creating the stage on first use.

        Parameters
        ----------
        name : str
            Stage label, as it will appear in the funnel table.
        count : int
            Survivors to add. Stages accumulate across tiles, so this is called once
            per tile per filter.
        """
        count = int(count)
        if name in self._index:
            self.stages[self._index[name]][1] += count
        else:
            self._index[name] = len(self.stages)
            self.stages.append([name, count])

    def get(self, name, default=0):
        return self.stages[self._index[name]][1] if name in self._index else default

    def as_dict(self):
        return {name: count for name, count in self.stages}

    def render(self):
        """Formats the funnel as a table: count, share of the DEM, share of the previous stage."""
        if not self.stages:
            return ""
        total = max(self.stages[0][1], 1)
        width = max(len(name) for name, _ in self.stages)
        lines = [f"   {'stage'.ljust(width)} | {'pixels':>15} | {'of DEM':>8} | {'of prev':>8}",
                 "   " + "-" * (width + 40)]
        prev = None
        for name, count in self.stages:
            of_dem = f"{100.0*count/total:7.3f}%"
            of_prev = "        -" if prev is None else f"{100.0*count/prev:7.3f}%" if prev else "      n/a"
            lines.append(f"   {name.ljust(width)} | {count:>15,} | {of_dem:>8} | {of_prev:>8}")
            prev = count
        return "\n".join(lines)


def _package_versions():
    """Best-effort version lookup for the third-party stack, for provenance."""
    import importlib
    versions = {}
    for name in ("numpy", "scipy", "numba", "tifffile", "matplotlib", "tqdm", "psutil", "imagecodecs"):
        try:
            mod = importlib.import_module(name)
            versions[name] = getattr(mod, "__version__", "unknown")
        except Exception:
            versions[name] = None
    return versions


def _git_state():
    """Returns the repository commit/branch and whether the tree was dirty."""
    import subprocess
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    def run(*args):
        return subprocess.run(args, cwd=repo_dir, capture_output=True, text=True,
                              timeout=10, check=True).stdout.strip()
    try:
        return {
            "commit": run("git", "rev-parse", "HEAD"),
            "branch": run("git", "rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(run("git", "status", "--porcelain")),
        }
    except Exception:
        return {"commit": None, "branch": None, "dirty": None}


def _file_digest(path, block=1 << 20):
    """SHA-256 of a file, streamed so multi-gigabyte DEMs do not land in RAM."""
    import hashlib
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(block), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def collect_provenance(dem_path, map_grid):
    """
    Captures everything needed to reproduce a run: code version, input identity,
    environment and invocation. Written alongside the scientific outputs.

    Parameters
    ----------
    dem_path : str
        Path to the DEM, whose sha256 is recorded.
    map_grid : MapGrid
        Resolved grid geometry, recorded along with where its resolution came from.

    Returns
    -------
    dict
        Git commit and dirty flag, DEM path, size and checksum, resolved grid
        geometry, and third-party package versions -- enough to say what produced a
        result months later.
    """
    dem_abs = os.path.abspath(dem_path)
    return {
        "timestamp": datetime.now().isoformat(),
        "git": _git_state(),
        "dem": {
            "path": dem_abs,
            "sha256": _file_digest(dem_abs),
            "bytes": os.path.getsize(dem_abs) if os.path.exists(dem_abs) else None,
            "cell_size_deg": map_grid.cell_size_deg,
            "cell_size_y_m": map_grid.cell_size_y,
            "cell_size_x_m": map_grid.cell_size_x,
            "cell_size_source": map_grid.source,
        },
        "platform": {
            "python": sys.version.split()[0],
            "system": f"{os.uname().sysname} {os.uname().release}" if hasattr(os, "uname") else sys.platform,
            "cpu_count": multiprocessing.cpu_count(),
        },
        "packages": _package_versions(),
        "command": " ".join(sys.argv),
    }


# ==========================================
#           CORE PIPELINE HELPERS
# ==========================================

def slope_baseline_pixels(map_grid, slope_baseline_m: float | None) -> tuple[int, int]:
    """
    Converts a slope measurement baseline in metres to a per-axis window in pixels.

    Slope is scale-dependent: on real Andean terrain the median slope falls from
    ~17.8 deg measured over the DEM's native ~61 m to ~10.8 deg over 1 km, and the
    fraction passing a 3-25 deg band rises from 60% to 78%. Which of those is
    "the" slope depends on the footprint being deployed, so the baseline is an
    explicit parameter rather than an accident of the DEM's resolution.

    Parameters
    ----------
    map_grid : MapGrid
        Angular and metric pixel sizes of the DEM.
    slope_baseline_m : float or None
        Ground distance over which slope is measured, in metres. ``None`` or 0 uses the
        DEM's native resolution, which on 30 m data is dominated by DEM noise.

    Returns
    -------
    tuple of int
        Smoothing window as ``(rows, columns)`` in pixels. ``(0, 0)`` when no baseline
        is requested, meaning the native gradient.
    """
    if not slope_baseline_m:
        return 0, 0
    ny = max(1, int(round(slope_baseline_m / map_grid.cell_size_y)))
    nx = max(1, int(round(slope_baseline_m / map_grid.cell_size_x)))
    return ny, nx


def terrain_gradients(elevation_block: np.ndarray, cell_size_y: float,
                      cell_size_x: float, smooth_y: int = 0,
                      smooth_x: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """
    Smoothed partial derivatives of the surface, the raw material for slope and aspect.

    Smoothing before differentiating gives the average gradient over the window,
    which is what "slope at 1 km scale" means physically. Callers must supply a
    block with a halo of at least max(smooth)//2 + 1 and crop the result, otherwise
    the window reaches past the block edge.

    Kept separate from :func:`terrain_derivatives` because the screening stage wants
    the gradients themselves: a slope *band* can be tested without ever forming the
    angle (see :func:`slope_band_gradient_sq`), and aspect is needed only at the few
    pixels that survive.

    Parameters
    ----------
    elevation_block : ndarray
        Elevation tile, including a halo of at least ``max(smooth)//2 + 1``.
    cell_size_y, cell_size_x : float
        Ground size of one pixel on each axis, in metres. They differ on a geographic
        grid, which is why they are separate.
    smooth_y, smooth_x : int, optional
        Smoothing window in pixels, from :func:`slope_baseline_pixels`.

    Returns
    -------
    tuple of ndarray
        ``(d/dy, d/dx)``, in metres per metre.
    """
    block = elevation_block
    if smooth_y > 1 or smooth_x > 1:
        block = uniform_filter(block, size=(max(1, smooth_y), max(1, smooth_x)), mode="nearest")
    return np.gradient(block, cell_size_y, cell_size_x)


def slope_band_gradient_sq(min_slope_deg: float | None,
                           max_slope_deg: float | None
                           ) -> tuple[float | None, float | None]:
    """
    The slope band restated as bounds on the squared gradient magnitude.

    ``slope = atan(|grad|)`` rises monotonically with the gradient magnitude, so

        min <= atan(sqrt(g)) <= max   <=>   tan(min)^2 <= g <= tan(max)^2

    which tests the same pixels without a sqrt or an arctan. Bounds at or beyond the
    vertical, and non-positive lower bounds, are returned as None meaning "unbounded":
    tan is singular at 90 degrees and every real gradient satisfies them anyway.

    Parameters
    ----------
    min_slope_deg, max_slope_deg : float or None
        Edges of the accepted slope band, in degrees.

    Returns
    -------
    tuple
        Lower and upper bounds on ``dx^2 + dy^2``. Either may be ``None``, meaning
        unbounded on that side.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> lo, hi = ss.slope_band_gradient_sq(3.0, 25.0)
    >>> f"{lo:.4f} {hi:.4f}"
    '0.0027 0.2174'
    >>> ss.slope_band_gradient_sq(0.0, 90.0)      # both edges degenerate
    (None, None)
    """
    lo = None
    if min_slope_deg is not None and min_slope_deg > 0.0:
        if min_slope_deg >= 90.0:
            lo = np.inf                      # nothing on a real surface qualifies
        else:
            lo = math.tan(math.radians(min_slope_deg)) ** 2
    hi = None
    if max_slope_deg is not None and max_slope_deg < 90.0:
        hi = math.tan(math.radians(max_slope_deg)) ** 2
    return lo, hi


def terrain_derivatives(elevation_block: np.ndarray, cell_size_y: float,
                        cell_size_x: float, smooth_y: int = 0,
                        smooth_x: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """
    Slope and aspect over a stated measurement baseline.

    Parameters
    ----------
    elevation_block : ndarray
        Elevation tile, including a halo. See :func:`terrain_gradients`.
    cell_size_y, cell_size_x : float
        Ground size of one pixel on each axis, in metres.
    smooth_y, smooth_x : int, optional
        Smoothing window in pixels.

    Returns
    -------
    tuple of ndarray
        Slope in degrees, and aspect in degrees clockwise from north.

    See Also
    --------
    slope_band_gradient_sq : tests a slope band without forming the angle at all,
        which is what the screening stage uses.
    """
    dy, dx = terrain_gradients(elevation_block, cell_size_y, cell_size_x, smooth_y, smooth_x)
    slope = np.degrees(np.arctan(np.sqrt(dx ** 2 + dy ** 2)))
    aspect = np.degrees(np.arctan2(-dx, dy)) % 360
    return slope, aspect


def find_results_json(run_dir):
    """
    Locates a run's results JSON, under either the current or the legacy prefix.

    Outputs used to be named ``grand_search_results_*`` whatever the experiment. The
    prefix is now ``oroscope_results_``, and both are accepted so that runs made before
    the rename still load -- a reader that could not open last week's output would make
    the rename cost more than it saves.

    Parameters
    ----------
    run_dir : str
        A run's output directory.

    Returns
    -------
    str or None
        Path to the results JSON, or ``None`` if the directory holds none.
    """
    import glob as _glob
    for prefix in (RESULTS_PREFIX, LEGACY_RESULTS_PREFIX):
        found = sorted(_glob.glob(os.path.join(run_dir, prefix + "*.json")))
        if found:
            return found[0]
    return None


def read_dem_origin(dem_path):
    """
    North-west corner of a GeoTIFF, from its ``ModelTiepointTag``.

    Standard geographic DEMs carry their own corner, so asking a user to type it is
    asking for a mistake that nothing catches: an origin that disagrees with the file
    does not fail, it silently georeferences every output to the wrong ground. Reading
    it removes the most error-prone input the tool has.

    Parameters
    ----------
    dem_path : str
        Path to the input elevation GeoTIFF.

    Returns
    -------
    tuple
        ``(latitude, longitude)`` of the north-west corner in degrees, or
        ``(None, None)`` when the file or the tag cannot be read -- which is not an
        error, since a caller may supply the origin explicitly.
    """
    try:
        with tiff.TiffFile(dem_path) as tf:
            tie = tf.pages[0].tags["ModelTiepointTag"].value
        return float(tie[4]), float(tie[3])
    except Exception:
        return None, None


def resolve_origin(dem_path, origin_lat=None, origin_lon=None, tolerance_deg=1e-3):
    """
    Settles the DEM's origin, preferring the file and checking anything supplied.

    Two failure modes, and the second is the dangerous one. An origin nobody supplied
    used to be a fatal error even though the file knows it. And an origin supplied
    *wrongly* was accepted in silence, mis-georeferencing every output -- the GeoTIFF,
    the world file, the KML and every coordinate in the results -- while the search
    itself ran perfectly and looked right.

    So the tag wins when nothing is given, and disagreement past ``tolerance_deg`` is
    reported loudly rather than resolved quietly. 1e-3 degrees is about 100 m, which is
    a few pixels: closer than that is rounding in a config file, further is a mistake.

    Parameters
    ----------
    dem_path : str
        Path to the input elevation GeoTIFF.
    origin_lat, origin_lon : float, optional
        Origin as supplied by the user, if any.
    tolerance_deg : float, optional
        Disagreement beyond which the supplied value is called out.

    Returns
    -------
    tuple
        ``(latitude, longitude, source)``, where source describes where the value came
        from and is recorded in the run's provenance.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> lat, lon, source = ss.resolve_origin("nonexistent.tif", -15.3, -72.4)
    >>> (lat, lon, source)
    (-15.3, -72.4, 'supplied (DEM carries no tiepoint)')
    """
    file_lat, file_lon = read_dem_origin(dem_path)

    if origin_lat is None or origin_lon is None:
        if file_lat is None:
            return None, None, "missing"
        return file_lat, file_lon, "auto-detected from the GeoTIFF tiepoint"

    if file_lat is None:
        return float(origin_lat), float(origin_lon), "supplied (DEM carries no tiepoint)"

    off = max(abs(file_lat - origin_lat), abs(file_lon - origin_lon))
    if off > tolerance_deg:
        print(f"{C.FAIL}{Icon.WARN}The supplied origin disagrees with the DEM's own "
              f"tiepoint by {off:.4f} deg ({off * 111:.0f} km).{C.RESET}")
        print(f"{C.WARN}   supplied: {origin_lat:.6f}, {origin_lon:.6f}{C.RESET}")
        print(f"{C.WARN}   GeoTIFF:  {file_lat:.6f}, {file_lon:.6f}{C.RESET}")
        print(f"{C.WARN}   Using the supplied value, but every output -- the GeoTIFF, "
              f"the world file, the KML and every coordinate in the results -- will be "
              f"georeferenced to that. Omit origin_lat/origin_lon to use the file's "
              f"own.{C.RESET}")
        return float(origin_lat), float(origin_lon), "supplied (DISAGREES with the tiepoint)"

    return float(origin_lat), float(origin_lon), "supplied (agrees with the tiepoint)"


def read_dem_geometry(dem_path):
    """
    Reads the angular pixel size and row count of a GeoTIFF DEM from its header.

    Standard geographic (EPSG:4326) DEMs such as SRTMGL1 or AW3D30 store the pixel
    size in degrees, which is what the georeferenced outputs (.tfw, .kml) require.
    Only the header is touched, so this stays cheap on multi-gigabyte files.

    Parameters
    ----------
    dem_path : str
        Path to the input elevation GeoTIFF.

    Returns
    -------
    tuple
        Pixel size in degrees and the number of rows. Either is ``None`` when the file
        or the tag cannot be read, which is not an error: the caller falls back to an
        explicit value or to 1 arc-second.
    """
    try:
        with tiff.TiffFile(dem_path) as tf:
            page = tf.pages[0]
            n_rows = int(page.shape[0])
            tag = page.tags.get('ModelPixelScaleTag')
            if tag is None:
                return None, n_rows
            scale_x, scale_y = float(tag.value[0]), float(tag.value[1])
            if scale_y <= 0:
                return None, n_rows
            # A geographic grid is square in degrees; warn if the DEM says otherwise.
            if scale_x > 0 and abs(scale_x - scale_y) / scale_y > 1e-6:
                print(f"      {C.WARN}{Icon.WARN}WARNING: Non-square DEM pixels "
                      f"({scale_x:.8f} x {scale_y:.8f} deg). Using the latitude scale.{C.RESET}")
            return scale_y, n_rows
    except Exception:
        return None, None

# Describes the sampling grid of the DEM. Angular pixel size is identical on both
# axes (that is what "geographic" means), but the two metric sizes are not.
# Stem of every output file. Was "grand_search_results_", which a TAMBO run also wrote
# and which was plainly wrong once one engine served more than one experiment. Readers
# accept the old prefix too, so runs made before the rename still load.
RESULTS_PREFIX = "oroscope_results_"
LEGACY_RESULTS_PREFIX = "grand_search_results_"

MapGrid = namedtuple("MapGrid", "cell_size_deg cell_size_y cell_size_x center_lat source")
MapGrid.__doc__ = """
Resolved pixel geometry of a DEM.

Parameters
----------
cell_size_deg : float
    Angular pixel size, in degrees. A geographic DEM steps by the same angle on both
    axes, which is why this is a single number while the metric sizes are two.
cell_size_y : float
    North-south ground size of one pixel, in metres.
cell_size_x : float
    East-west ground size of one pixel, in metres. Smaller than ``cell_size_y`` away
    from the equator, by the cosine of the latitude.
center_lat : float
    Latitude at which ``cell_size_x`` was evaluated, in degrees.
source : str
    Where the resolution came from -- detected from the GeoTIFF, supplied explicitly,
    or defaulted. Recorded so a run's provenance says which.
"""

def resolve_grid_geometry(dem_path, origin_lat, cell_size_deg=None):
    """
    Determines the sampling geometry used by the whole pipeline.

    Resolution priority: explicit user value > GeoTIFF ModelPixelScaleTag > 1 arc-second.

    A geographic raster has pixels that are square in degrees but *not* in metres: a
    degree of longitude shrinks with the cosine of the latitude, so a 1 arc-second
    pixel spans roughly 30.7 m north-south but only ~29.5 m east-west at 17 degrees
    south. The longitude scale is evaluated at the DEM's centre latitude so the
    residual error from ignoring its north-south variation is spread evenly over the
    map rather than accumulating towards one edge.

    Parameters
    ----------
    dem_path : str
        Path to the input elevation GeoTIFF.
    origin_lat : float
        Latitude of the DEM's northern edge, in degrees, for the latitude-dependent
        east-west scaling.
    cell_size_deg : float, optional
        Explicit pixel size in degrees, overriding whatever the file says.

    Returns
    -------
    MapGrid
        Angular pixel size, both metric pixel sizes, the centre latitude used for the
        longitude scaling, and where the resolution value came from -- recorded so a
        run's provenance says whether the resolution was detected or asserted.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> grid = ss.resolve_grid_geometry("nonexistent.tif", -15.6, cell_size_deg=1/3600)
    >>> f"{grid.cell_size_y:.1f} m x {grid.cell_size_x:.1f} m"
    '30.7 m x 29.8 m'
    """
    detected_deg, n_rows = read_dem_geometry(dem_path)

    if cell_size_deg is not None:
        source = "user-specified"
    elif detected_deg is not None:
        cell_size_deg = detected_deg
        source = "auto-detected from GeoTIFF"
    else:
        cell_size_deg = DEFAULT_CELL_SIZE_DEG
        source = "assumed 1 arc-second (DEM carries no georeferencing tags)"

    center_lat = origin_lat
    if n_rows:
        center_lat = origin_lat - (n_rows / 2.0) * cell_size_deg

    cell_size_y = cell_size_deg * KM_PER_DEG_LAT * 1000.0
    cell_size_x = cell_size_deg * KM_PER_DEG_LON_EQUATOR * math.cos(math.radians(center_lat)) * 1000.0

    return MapGrid(float(cell_size_deg), float(cell_size_y), float(cell_size_x),
                   float(center_lat), source)

# Values below this are ocean or void in the DEMs this tool reads
NODATA_BELOW_M = -100.0


def build_elevation_cache(dem_path, npy_path, block_rows=2048):
    """
    Converts a DEM to the memory-mapped float32 cache, without ever holding it in RAM.

    The obvious ``tiff.imread(path).astype(np.float32)`` materialises the whole DEM and
    then a second full copy of it — which defeats the point of the out-of-core design
    the rest of the pipeline is built around, and fails outright on the multi-gigabyte
    DEMs this tool is meant to handle. Instead the page is decoded straight into a
    native-dtype file, then converted a block of rows at a time. Peak memory is one
    block, whatever the size of the DEM.

    float32 with NaN is kept rather than the DEM's own integer dtype: NaN propagates
    through the gradient and comparison chain in the screening stage, so nodata is
    excluded without a sentinel test in every kernel. That costs twice the disk of an
    int16 cache and buys correctness that would otherwise have to be re-established in
    half a dozen places.

    Parameters
    ----------
    dem_path : str
        Path to the input GeoTIFF.
    npy_path : str
        Path to write the float32 memory-mapped cache to.
    block_rows : int, optional
        Rows converted at a time. Peak memory is one block, whatever the DEM's size.
    """
    raw_path = npy_path + ".raw"
    try:
        with tiff.TiffFile(dem_path) as tf:
            page = tf.pages[0]
            rows, cols = int(page.shape[0]), int(page.shape[1])
            raw = np.lib.format.open_memmap(raw_path, mode="w+", shape=(rows, cols),
                                            dtype=page.dtype)
            page.asarray(out=raw)
            raw.flush()
            del raw

        raw = np.load(raw_path, mmap_mode="r")
        out = np.lib.format.open_memmap(npy_path, mode="w+", shape=(rows, cols),
                                        dtype=np.float32)
        for r in range(0, rows, block_rows):
            block = raw[r:r + block_rows].astype(np.float32)
            block[block < NODATA_BELOW_M] = np.nan      # ocean and void
            out[r:r + block_rows] = block
        out.flush()
        del out, raw
    finally:
        if os.path.exists(raw_path):
            os.remove(raw_path)


def load_dem_and_init_buffers(dem_path, temp_dir, resume=False, resume_dir=None):
    """
    Step 1 Pipeline: Converts TIF to memory-mapped NPY for rapid random access
    and initializes the ping-pong buffers for later morphology steps.
    If resume is True and resume_dir is provided, it attempts to load an existing ray-tracing buffer.
    
    Parameters
    ----------
    dem_path : str
        Path to the input elevation GeoTIFF.
    temp_dir : str
        Directory for the working buffers.
    resume : bool, optional
        Reuse a previous run's scan buffer instead of recomputing it.
    resume_dir : str, optional
        Directory holding that buffer.

    Returns
    -------
    elevation : ndarray
        Memory-mapped DEM.
    rows, cols : int
        Array dimensions.
    path_A, path_B : str
        Paths to the boolean ping-pong buffers.
    buf_a : ndarray
        Open memory map of buffer A.
    is_resuming : bool
        ``True`` if a previous scan buffer was loaded successfully.
    """
    npy_path = dem_path.replace(".tif", ".npy")
    if not os.path.exists(npy_path):
        build_elevation_cache(dem_path, npy_path)
    elevation = np.load(npy_path, mmap_mode='r')
    rows, cols = elevation.shape
    
    path_A = os.path.join(temp_dir, "buffer_A.npy")
    path_B = os.path.join(temp_dir, "buffer_B.npy")
    
    is_resuming = False
    if resume and resume_dir and os.path.exists(os.path.join(resume_dir, "buffer_A.npy")):
        src_buffer = os.path.join(resume_dir, "buffer_A.npy")
        if os.path.abspath(src_buffer) != os.path.abspath(path_A):
            print(f"      {Icon.INFO}{C.WARN}Resuming: Copying existing physics buffer from {src_buffer}...{C.RESET}")
            shutil.copy(src_buffer, path_A)
        else:
            print(f"      {Icon.INFO}{C.WARN}Resuming: Using existing physics buffer at {src_buffer}...{C.RESET}")
        buf_a = np.lib.format.open_memmap(path_A, mode='r+', shape=(rows, cols), dtype=bool)
        is_resuming = True
    else:
        buf_a = np.lib.format.open_memmap(path_A, mode='w+', shape=(rows, cols), dtype=bool)
        
    buf_b = np.lib.format.open_memmap(path_B, mode='w+', shape=(rows, cols), dtype=bool)
    del buf_b # Release lock immediately; only keep file on disk
    
    return elevation, rows, cols, path_A, path_B, buf_a, is_resuming

def get_candidates_chunked(elevation, map_grid, rfi_zones, origin_lat, origin_lon,
                           min_alt=None, max_alt=None, min_aspect_deg=None, max_aspect_deg=None,
                           road_map_path=None, max_road_dist_km=None,
                           min_slope_deg=3.0, max_slope_deg=25.0,
                           tile_size=2048, candidate_stride=5, slope_baseline_m=None, funnel=None):
    """
    Step 2 Pipeline: Memory-efficient topographic screening. Iterates over the large DEM in chunks (tiles) 
    to find pixels that meet the primary geometrical criteria (slope, aspect, altitude) 
    and logistics constraints (RFI distance, road distance) prior to running ray-tracing.
    
    Parameters:

    - elevation (ndarray): Full DEM array (usually memory-mapped).
    - map_grid (MapGrid): Angular and metric pixel sizes of the DEM.
    - rfi_zones (list): List of configured exclusion zones (circles/polygons).
    - origin_lat, origin_lon (float): Reference coordinates for converting km to pixels.
    - min_alt, max_alt (float): Elevation restrictions.
    - min_aspect_deg, max_aspect_deg (float): Required facing directions for slopes.
    - road_map_path (str): Path to an aligned TIFF containing distance-to-road values.
    - max_road_dist_km (float): Maximum allowed distance from a road.
    - min_slope_deg, max_slope_deg (float): Required steepness limits for detector slopes.
    - tile_size (int): Size of the square chunk to process in RAM at one time.
    - candidate_stride (int): Keeps every Nth surviving pixel before ray-tracing. Higher
      values trade spatial sampling density for speed; 1 keeps every candidate.
    - slope_baseline_m (float): Ground distance over which slope is measured. None uses
      the DEM's native resolution, which on 30 m data is dominated by DEM noise.
    - funnel (Funnel): Optional accounting object recording per-filter survivor counts.

    Parameters
    ----------
    elevation : ndarray
        Full DEM, usually memory-mapped.
    map_grid : MapGrid
        Angular and metric pixel sizes.
    rfi_zones : list or None
        Exclusion zones as ``('circle', lat, lon, radius_km, name)`` or
        ``('poly', [(lat, lon), ...], name)``.
    origin_lat, origin_lon : float
        North-west corner of the DEM, in degrees, for converting zones to pixels.
    min_alt, max_alt : float, optional
        Altitude bounds, in metres.
    min_aspect_deg, max_aspect_deg : float, optional
        Required facing directions, in degrees clockwise from north. Wraps through 360
        when the lower bound exceeds the upper.
    road_map_path : str, optional
        Aligned GeoTIFF of distance-to-road values.
    max_road_dist_km : float, optional
        Maximum allowed distance from a road, in km.
    min_slope_deg, max_slope_deg : float, optional
        Slope band, in degrees. Tested on the squared gradient, so neither a square
        root nor an arctangent is formed over the tile.
    tile_size : int, optional
        Side of the square chunk processed in RAM at once.
    candidate_stride : int, optional
        Keeps every Nth surviving pixel. Measured to be unbiased: acceptance is
        identical at strides 1 and 5, and the stride-corrected area matches the
        stride-1 truth to 0.05%.
    slope_baseline_m : float, optional
        Ground distance over which slope is measured, in metres. ``None`` uses the
        DEM's native resolution, which on 30 m data is dominated by DEM noise.
    funnel : Funnel, optional
        Accounting object recording per-filter survivor counts.

    Returns
    -------
    ndarray
        ``(N, 3)`` array of surviving pixels as ``[row, col, aspect_deg]``, ready for
        :func:`run_arrival_scan`.
    """
    rows, cols = elevation.shape
    candidates_list = []

    # A geographic DEM steps by the same angle on both axes, so one pixel is one
    # cell_size_deg of latitude *and* of longitude. The two metric sizes differ, and
    # are applied wherever a real ground distance (rather than an angle) is meant.
    cell_size_deg = map_grid.cell_size_deg
    cell_size_y = map_grid.cell_size_y
    cell_size_x = map_grid.cell_size_x

    # Slope measurement baseline, and the halo the derivative window needs
    smooth_y, smooth_x = slope_baseline_pixels(map_grid, slope_baseline_m)
    halo_y = max(1, smooth_y // 2 + 1)
    halo_x = max(1, smooth_x // 2 + 1)

    # The slope band, once, as bounds on the squared gradient
    grad_lo, grad_hi = slope_band_gradient_sq(min_slope_deg, max_slope_deg)

    # Load Logistics Road map if provided
    road_dist_map = None
    if road_map_path and max_road_dist_km:
        if os.path.exists(road_map_path):
            try:
                road_dist_map = tiff.imread(road_map_path, out='memmap')
                print(f"      {Icon.INFO}Logistics: Loaded Road Distance Map ({road_map_path})")
            except Exception:
                print(f"      {C.WARN}{Icon.WARN}WARNING: Could not load road map.{C.RESET}")
        else:
            print(f"      {C.WARN}{Icon.WARN}WARNING: Road map file not found.{C.RESET}")

    # Convert geographic RFI definitions into pixel coordinates for local checking
    rfi_circles = [] 
    rfi_polys = []   
    if rfi_zones:
        for item in rfi_zones:
            type_tag = item[0]
            if type_tag == 'circle':
                _, zlat, zlon, zrad_km, _ = item
                z_r = (origin_lat - zlat) / cell_size_deg
                z_c = (zlon - origin_lon) / cell_size_deg
                # Keep the radius in metres: the distance test below is done on the
                # ground, so the zone stays a true circle instead of a pixel-space one
                rfi_circles.append((z_r, z_c, (zrad_km * 1000.0)**2))
            elif type_tag == 'poly':
                _, coords, _ = item
                pixel_verts = []
                for (plat, plon) in coords:
                    pr = (origin_lat - plat) / cell_size_deg
                    pc = (plon - origin_lon) / cell_size_deg
                    pixel_verts.append((pc, pr))
                rfi_polys.append(np.array(pixel_verts, dtype=np.float64))

    r_steps = range(0, rows, tile_size)
    c_steps = range(0, cols, tile_size)
    
    # Process the map in chunks to avoid blowing out system RAM
    with tqdm(total=len(r_steps)*len(c_steps), desc="   Scanning Topography", unit="tile", colour='magenta' if USE_COLOR else None) as pbar:
        for r in r_steps:
            for c in c_steps:
                r_end = min(r + tile_size, rows)
                c_end = min(c + tile_size, cols)
                chunk = elevation[r:r_end, c:c_end]
                if funnel is not None:
                    funnel.add("DEM pixels", chunk.size)
                    funnel.add("finite elevation", np.count_nonzero(~np.isnan(chunk)))

                # Read a haloed block so the derivative window never reaches past the
                # tile edge; without it, slope at tile boundaries is a tiling artefact
                r_lo, r_hi = max(0, r - halo_y), min(rows, r_end + halo_y)
                c_lo, c_hi = max(0, c - halo_x), min(cols, c_end + halo_x)
                block = elevation[r_lo:r_hi, c_lo:c_hi]

                dy_block, dx_block = terrain_gradients(
                    block, cell_size_y, cell_size_x, smooth_y, smooth_x)

                core = (slice(r - r_lo, r - r_lo + (r_end - r)),
                        slice(c - c_lo, c - c_lo + (c_end - c)))
                dy, dx = dy_block[core], dx_block[core]

                # Filter 1: Fundamental detector slope requirement, tested on the
                # squared gradient so the whole tile needs neither sqrt nor arctan.
                # NaN elevations propagate into the gradient and compare False, which
                # is how they were excluded before.
                grad_sq = dx * dx + dy * dy
                if grad_lo is None and grad_hi is None:
                    mask = np.ones(grad_sq.shape, dtype=bool)
                elif grad_lo is None:
                    mask = grad_sq <= grad_hi
                elif grad_hi is None:
                    mask = grad_sq >= grad_lo
                else:
                    mask = (grad_sq >= grad_lo) & (grad_sq <= grad_hi)
                if funnel is not None:
                    funnel.add(f"slope {min_slope_deg}-{max_slope_deg} deg", np.count_nonzero(mask))

                # Filter 2: Altitude bounds
                if min_alt is not None: mask &= (chunk >= min_alt)
                if max_alt is not None: mask &= (chunk <= max_alt)
                if funnel is not None and (min_alt is not None or max_alt is not None):
                    funnel.add("altitude bounds", np.count_nonzero(mask))

                # From here on the tile is represented by its survivors alone. The
                # remaining filters each touched the whole tile before; on real terrain
                # the subset is a fraction of it, and aspect -- the one transcendental
                # left -- is now evaluated only where it is actually read.
                #
                # An empty tile is deliberately not short-circuited: every filter still
                # reports its (zero) count, so a search that rejects everything still
                # produces a complete funnel, which is the case the funnel is for.
                cr, cc = np.where(mask)
                aspect_vals = np.degrees(np.arctan2(-dx[cr, cc], dy[cr, cc])) % 360
                keep = np.ones(len(cr), dtype=bool)

                # Filter 3: Aspect bounds (handle wrapping around 360 degrees)
                if min_aspect_deg is not None and max_aspect_deg is not None:
                    min_a, max_a = min_aspect_deg, max_aspect_deg
                    if min_a > max_a:
                        keep &= (aspect_vals >= min_a) | (aspect_vals <= max_a)
                    else:
                        keep &= (aspect_vals >= min_a) & (aspect_vals <= max_a)
                    if funnel is not None:
                        funnel.add("aspect bounds", int(np.count_nonzero(keep)))

                # Filter 4: Road Logistics
                if road_dist_map is not None:
                    road_chunk = road_dist_map[r:r_end, c:c_end]
                    keep &= (road_chunk[cr, cc] <= (max_road_dist_km * 1000))
                    if funnel is not None:
                        funnel.add("road distance", int(np.count_nonzero(keep)))

                # Filter 5: Dynamic Exclusion Zones (RFI)
                if rfi_zones:
                    abs_r = (r + cr).astype(np.float64)
                    abs_c = (c + cc).astype(np.float64)

                    # Circular exclusion zones, measured in metres on the ground
                    for (zr, zc, zrad_m_sq) in rfi_circles:
                        dist_sq = ((abs_r - zr) * cell_size_y)**2 + ((abs_c - zc) * cell_size_x)**2
                        keep &= dist_sq >= zrad_m_sq

                    # Polygonal exclusion zones. apply_poly_mask_numba only ever clears
                    # bits, so it can be handed `keep` directly and accumulate in place.
                    for poly in rfi_polys:
                        apply_poly_mask_numba(abs_r, abs_c, poly, keep)
                    if funnel is not None:
                        funnel.add("outside RFI zones", int(np.count_nonzero(keep)))

                # Extract surviving pixels for the physics simulation step
                cr, cc, aspect_vals = cr[keep], cc[keep], aspect_vals[keep]
                if funnel is not None:
                    funnel.add(f"kept by stride {candidate_stride}", len(cr[::candidate_stride]))
                if len(cr) > 0:
                    chunk_cands = np.column_stack((cr + r, cc + c, aspect_vals))
                    # Thin the candidates to speed up ray tracing; assumption is terrain is continuous
                    kept = chunk_cands[::candidate_stride]
                    candidates_list.append(kept)
                    pbar.set_postfix(candidates=f"{len(kept):,}")
                pbar.update(1)

    if not candidates_list: return np.zeros((0, 3))
    return np.vstack(candidates_list)

def run_arrival_scan(candidates_arr, elevation, map_grid, buf_a, scan_params,
                     score_config=None, min_score=0.0, rfi_zones_px=None,
                     score_percentile=None):
    """
    Step 3 alternative: scan arrival directions instead of casting one ray per pixel.

    Marks a candidate as valid when at least one accepted (azimuth, elevation)
    direction strikes rock within the decay-baseline window with enough column depth.
    See arrival_scan.py for the geometry.

    Parameters
    ----------
    candidates_arr : ndarray
        ``(N, 3)`` array of ``[row, col, aspect_deg]``.
    elevation : ndarray
        The DEM.
    map_grid : MapGrid
        Angular and metric pixel sizes.
    buf_a : ndarray
        Open memory map to mark accepted pixels in.
    scan_params : dict
        Keyword arguments for :func:`arrival_scan.scan`.
    score_config : dict, optional
        Overrides for :data:`scoring.DEFAULT_SCORE_CONFIG`.
    min_score : float, optional
        Absolute score a candidate must reach. Used only when ``score_percentile`` is
        not given. The default composition is a product, whose distribution piles up
        near zero, so any threshold in the middle sits on a cliff.
    rfi_zones_px : sequence, optional
        Radio-noise sources in pixel coordinates, enabling the exposure observable.
    score_percentile : float, optional
        Keep this percentage of viable candidates, by score. Rank-based and so
        scale-free: preferred over ``min_score`` for exactly the reason above.

    Returns
    -------
    n_hits : int
        Number of accepted candidates.
    observables : dict
        Per-candidate arrays, including the scores and their named components, kept
        for per-site aggregation.
    """
    observables = arrival_scan.scan(candidates_arr, elevation, map_grid, **scan_params)

    # Site altitude enters the footprint term: a higher site has a narrower Cherenkov
    # cone, so the same array spacing samples the footprint less well
    observables["altitude_m"] = elevation[candidates_arr[:, 0].astype(np.int64),
                                          candidates_arr[:, 1].astype(np.int64)]
    if rfi_zones_px:
        observables["rfi_exposure"] = arrival_scan.rfi_exposure(
            candidates_arr, elevation, map_grid, rfi_zones_px)

    # Score every candidate, then keep those clearing the floor. Scores travel with the
    # observables so per-site records can report their distribution.
    window = (scan_params.get("min_dist_km", 0.0) * 1000.0,
              scan_params.get("max_dist_km", 0.0) * 1000.0)
    total, components = scoring.score_candidates(observables, score_config, window)
    observables["score"] = total
    for name, values in components.items():
        observables[f"score_{name}"] = values

    # A rank-based cut where one is asked for, and an absolute one otherwise.
    #
    # The default score is a *product* of several components each in [0, 1], so its
    # distribution piles up near zero and an absolute threshold sits on a cliff:
    # measured on one search, min_score 0.0, 0.35 and 0.5 gave 45928, 2056 and zero
    # detector positions. A percentile asks the question that was meant all along --
    # keep the best fraction of what this terrain offers -- and is scale-free, so it
    # does not move when the composition or the number of components changes.
    viable = observables["cells"] > 0
    if score_percentile is not None and np.any(viable):
        floor = float(np.percentile(total[viable], 100.0 - float(score_percentile)))
        accepted = viable & (total >= floor)
    else:
        accepted = viable & (total >= min_score)
    n_hits = int(np.count_nonzero(accepted))
    if n_hits:
        buf_a[candidates_arr[accepted, 0].astype(np.int64),
              candidates_arr[accepted, 1].astype(np.int64)] = True
    buf_a.flush()
    return n_hits, observables


def summarize_observables_by_site(labeled, downsample_factor, candidates_arr, observables,
                                  site_ids):
    """
    Aggregates per-candidate scan observables over each labelled site.

    Storing the distributions rather than a single score is deliberate: absolute
    apertures can then be obtained later by folding these against an acceptance table,
    without re-running the terrain analysis (roadmap 4.10).

    Parameters
    ----------
    labeled : ndarray
        Downsampled labelled site map.
    downsample_factor : int
        Factor relating full-resolution candidate coordinates to ``labeled``.
    candidates_arr : ndarray
        ``(N, 3)`` array of ``[row, col, aspect_deg]``.
    observables : dict
        Per-candidate arrays from :func:`run_arrival_scan`.
    site_ids : sequence of int
        Sites to summarise.

    Returns
    -------
    dict
        Site id to summary statistics -- mean, median and 90th percentile of each
        observable -- over that site's accepted candidates. Empty when there are no
        accepted candidates at all.
    """
    if observables is None or candidates_arr is None or len(site_ids) == 0:
        return {}

    accepted = observables["cells"] > 0
    if not np.any(accepted):
        return {}

    rows = (candidates_arr[accepted, 0].astype(np.int64)) // downsample_factor
    cols = (candidates_arr[accepted, 1].astype(np.int64)) // downsample_factor
    inside = (rows >= 0) & (rows < labeled.shape[0]) & (cols >= 0) & (cols < labeled.shape[1])
    rows, cols = rows[inside], cols[inside]
    site_of = labeled[rows, cols]

    fields = ["solid_angle_sr", "mean_distance_m", "max_depth_gcm2", "horizon_deg"]
    for extra in ("score", "best_clearance_ratio", "geomag_solid_angle_sr",
                  "path_grammage_gcm2", "earth_chord_gcm2", "altitude_m",
                  "target_slope_deg", "rfi_exposure"):
        if extra in observables:
            fields.append(extra)
    # The named score components, each in [0, 1]. Without them the record carries the
    # total and nothing else, so a weak site cannot be attributed to the criterion that
    # weakened it -- which is the whole reason the components are named (scoring.py).
    fields.extend(sorted(k for k in observables
                         if k.startswith("score_") and k not in fields))
    fields = tuple(fields)
    values = {f: observables[f][accepted][inside] for f in fields}

    # Group the candidates by site with one sort rather than one full-array comparison
    # per site. The sort is stable, so each site's values stay in the order a boolean
    # mask would have produced them and the statistics are unchanged to the last bit.
    order = np.argsort(site_of, kind="stable")
    sorted_ids = site_of[order]
    starts = np.searchsorted(sorted_ids, site_ids, side="left")
    stops = np.searchsorted(sorted_ids, site_ids, side="right")

    summary = {}
    for site_id, start, stop in zip(site_ids, starts, stops):
        if stop == start:
            continue
        idx = order[start:stop]
        entry = {"scanned_pixels": int(stop - start)}
        for f in fields:
            v = values[f][idx]
            entry[f"{f}_mean"] = float(np.mean(v))
            entry[f"{f}_p50"] = float(np.median(v))
            entry[f"{f}_p90"] = float(np.percentile(v, 90))
        summary[int(site_id)] = entry
    return summary


def apply_morphology_pingpong(source_path, dest_path, shape, dtype, operation_func, structure, desc="Processing", tile_size=2048):
    """
    Applies image morphology operations (closing/opening) on a massive memory-mapped array
    without loading the whole array into RAM. It reads from one file and writes to another ("ping-pong").

    Parameters
    ----------
    source_path, dest_path : str
        Paths to the two ``.npy`` buffers, read and written respectively.
    shape : tuple of int
        Shape of the arrays.
    dtype : dtype
        Element type of the destination.
    operation_func : callable
        Morphological operation, applied tile by tile.
    structure : ndarray
        Structuring element.
    desc : str, optional
        Label for the progress bar.
    tile_size : int, optional
        Side of the square tile held in RAM at once.

    Returns
    -------
    int
        Set pixels in the result, counted while writing so the funnel accounting costs
        nothing extra.
    """
    surviving = 0
    source = np.lib.format.open_memmap(source_path, mode='r')
    dest = np.lib.format.open_memmap(dest_path, mode='r+', shape=shape, dtype=dtype)
    rows, cols = shape
    
    with tqdm(total=(rows//tile_size + 1)*(cols//tile_size + 1), desc=f"   {desc}", unit="tile", colour='magenta' if USE_COLOR else None) as pbar:
        # Pad the chunk by half the structure size to prevent edge artifacts between chunks
        pad = max(structure.shape) // 2
        for r in range(0, rows, tile_size):
            for c in range(0, cols, tile_size):
                r_end = min(r + tile_size, rows)
                c_end = min(c + tile_size, cols)
                r_start = max(0, r - pad)
                c_start = max(0, c - pad)
                
                chunk = source[r_start:min(rows, r_end+pad), c_start:min(cols, c_end+pad)]
                processed = operation_func(chunk, structure)
                
                loc_r_start = r - r_start
                loc_c_start = c - c_start
                # Write back only the non-padded, processed core of the chunk
                core = processed[loc_r_start:loc_r_start + (r_end - r), loc_c_start:loc_c_start + (c_end - c)]
                dest[r:r_end, c:c_end] = core
                surviving += int(np.count_nonzero(core))
                pbar.update(1)
    dest.flush()
    return surviving

def separable_closing(chunk, structure):
    """
    Binary closing with a rectangular structuring element, done separably.

    A rectangle of ones factorises into a column and a row, so dilation or erosion by
    (h, w) is dilation by (h, 1) followed by (1, w). That turns an O(N h w) operation
    into O(N (h + w)) -- about 10x for the 33x33 element a 1 km antenna spacing implies
    -- and the result is bit-identical, not an approximation.

    Parameters
    ----------
    chunk : ndarray
        Boolean tile to operate on.
    structure : ndarray
        Rectangle of ones. Its two side lengths are what the operation factorises into.

    Returns
    -------
    ndarray
        The closed tile.
    """
    h, w = structure.shape
    col = np.ones((h, 1), dtype=bool)
    row = np.ones((1, w), dtype=bool)
    grown = binary_dilation(binary_dilation(chunk, col), row)
    return binary_erosion(binary_erosion(grown, col), row)


def separable_opening(chunk, structure):
    """
    Binary opening with a rectangular element, separably.

    Prunes features narrower than the element. See :func:`separable_closing` for why
    the factorisation is exact rather than an approximation.

    Parameters
    ----------
    chunk : ndarray
        Boolean tile to operate on.
    structure : ndarray
        Rectangle of ones.

    Returns
    -------
    ndarray
        The opened tile.
    """
    h, w = structure.shape
    col = np.ones((h, 1), dtype=bool)
    row = np.ones((1, w), dtype=bool)
    shrunk = binary_erosion(binary_erosion(chunk, col), row)
    return binary_dilation(binary_dilation(shrunk, col), row)


def clean_shape_artifacts(path_A, path_B, rows, cols, cell_size_y, cell_size_x, antenna_spacing_km, min_width_km, tile_size, gap_close_km=None):
    """
    Step 4 Pipeline: Prunes spatial artifacts to ensure solid, block-like arrays.
    Applies closing to fill gaps and opening to prune unusable tendrils.

    The structuring elements are sized per axis so that they cover the requested
    ground distance in both directions rather than only north-south.

    ``min_width_km = 0`` degenerates the opening to a 1x1 element, i.e. an identity
    that only carries the closed map back into ``path_A``. That is deliberate: a
    "block-like array" is a GRAND assumption, and an experiment deployed along a
    canyon wall is a strip a few hundred metres wide and tens of kilometres long,
    which the opening would delete outright.

    Parameters
    ----------
    path_A, path_B : str
        The two ping-pong buffers. The result is left in ``path_A``.
    rows, cols : int
        Array dimensions.
    cell_size_y, cell_size_x : float
        Ground size of one pixel on each axis, in metres.
    antenna_spacing_km : float
        Detector spacing, used as the default closing scale.
    min_width_km : float
        Narrowest feature to keep. 0 disables pruning, which is what a strip-shaped
        array needs.
    tile_size : int
        Side of the square tile held in RAM at once.
    gap_close_km : float, optional
        Size of the closing element, in km. Defaults to ``antenna_spacing_km``.

    Returns
    -------
    tuple of int
        Set-pixel counts after closing and after pruning.
    """
    # Gap closing is its own criterion, not a consequence of detector spacing. It used
    # to be tied to antenna_spacing_km, which coupled two unrelated things and hid how
    # much of the reported area it creates: measured at Colca, closing with a 1 km
    # element more than doubles the accepted area (2.29x, §6.17). 0 disables it.
    close_km = antenna_spacing_km if gap_close_km is None else gap_close_km
    close_r = max(1, int(close_km * 1000 / cell_size_y))
    close_c = max(1, int(close_km * 1000 / cell_size_x))
    tendril_r = max(1, int((min_width_km * 0.5 * 1000) / cell_size_y))
    tendril_c = max(1, int((min_width_km * 0.5 * 1000) / cell_size_x))
    n_closed = apply_morphology_pingpong(path_A, path_B, (rows, cols), bool, separable_closing, np.ones((close_r, close_c)), desc="Closing", tile_size=tile_size)
    n_pruned = apply_morphology_pingpong(path_B, path_A, (rows, cols), bool, separable_opening, np.ones((tendril_r, tendril_c)), desc="Pruning", tile_size=tile_size)
    return n_closed, n_pruned

def analyze_sites_and_capacity(path_A, elevation, rows, cols, cell_size_y, cell_size_x, downsample_factor, search_mode,
                               target_antennas, min_sub_array_size, antenna_spacing_km, grid_type, funnel=None,
                               origin_lat=None, origin_lon=None, cell_size_deg=None,
                               candidates_arr=None, observables=None, stop_at_target=False):
    """
    Step 5 Pipeline: Isolates unique sites and measures their capacity mathematically.
    Uses SciPy labeling to find continuous regions and simulates physical grid placement.
    
    Returns:
    - small_final (ndarray): Downsampled binary mask of the validated sites.
    - labeled_viz (ndarray): Multi-integer labeled array for color coding visualizations.
    - site_details (list): Dictionaries containing metadata about each valid site found.
    - cumulative_capacity (int): Sum of all antennas fitting in valid sites.
    - count (int): Total number of independent valid sites found.
    - region_stats (dict): Region-level accounting for the funnel report.

    Parameters
    ----------
    path_A : str
        Buffer holding the cleaned, full-resolution site mask.
    elevation : ndarray
        The DEM, for the per-site mean aspect.
    rows, cols : int
        Full-resolution dimensions.
    cell_size_y, cell_size_x : float
        Ground size of one pixel on each axis, in metres.
    downsample_factor : int
        Factor at which labelling and area are computed. Note area is measured on the
        downsampled map while capacity is measured at full resolution, so a feature
        only a few pixels wide loses area it keeps detectors on.
    search_mode : str
        ``single`` or ``distributed``, deciding which capacity threshold applies.
    target_antennas : int
        Capacity wanted from a single site.
    min_sub_array_size : int
        Capacity a sub-array must reach in distributed mode.
    antenna_spacing_km : float
        Detector spacing, in km.
    grid_type : str
        ``square`` or ``hex``.
    funnel : Funnel, optional
        Accounting object recording survivor counts.
    origin_lat, origin_lon : float, optional
        The DEM's north-west corner. With ``cell_size_deg``, each site record gains its
        centre coordinates and bounding box, so a reader can find the ground without
        opening the raster.
    cell_size_deg : float, optional
        Pixel size in degrees, at full resolution.
    candidates_arr : ndarray, optional
        Candidates, for folding scan observables into each site's record.
    observables : dict, optional
        Their per-candidate observables.
    stop_at_target : bool, optional
        In distributed mode, stop selecting sites once ``target_antennas`` is reached.
        Sites are sorted by capacity, so this takes the best ones and reports the array
        actually wanted rather than every patch of qualifying ground.

    Returns
    -------
    small_final : ndarray
        Downsampled binary mask of the validated sites.
    labeled_viz : ndarray
        Site labels for colour coding, sized from the label count so that selecting
        more than 255 sites does not overflow.
    site_details : list of dict
        Per-site metadata for every site that cleared the thresholds, sorted by
        capacity, each carrying a ``selected`` flag. With ``stop_at_target`` the list
        is longer than the selection: ``cumulative_capacity``, ``count`` and the
        exported mask cover the selected sites only, so anything totalling this list
        must filter on ``selected`` or it will over-report both area and site count.
    cumulative_capacity : int
        Total capacity across the selected sites.
    count : int
        Number of sites selected.
    region_stats : dict
        Region-level accounting for the funnel report.
    """
    final_map_disk = np.lib.format.open_memmap(path_A, mode='r')
    small_map = final_map_disk[::downsample_factor, ::downsample_factor]
    labeled, num = label(small_map) # Give unique integer IDs to disconnected array zones
    
    eff_cell_y = cell_size_y * downsample_factor
    eff_cell_x = cell_size_x * downsample_factor
    px_area_km2 = (eff_cell_y / 1000.0) * (eff_cell_x / 1000.0)
    
    if search_mode == 'single':
        threshold_antennas = target_antennas
    else:
        threshold_antennas = min_sub_array_size
        
    req_pixels = int((threshold_antennas * antenna_spacing_km**2) / px_area_km2)
    small_final = np.zeros_like(labeled, dtype=np.uint8)
    # Wide enough for the labels actually present. This was uint8, which raised
    # OverflowError as soon as a search selected a 256th site -- reachable in
    # distributed mode with a small min_sub_array_size, and the normal case for a
    # layout made of many small sub-arrays rather than one blob.
    viz_dtype = np.min_scalar_type(max(1, num))
    labeled_viz = np.zeros_like(labeled, dtype=viz_dtype)
    cumulative_capacity = 0
    site_details = []
    count = 0
    
    if num > 0:
        sizes = ndi_sum(small_map, labeled, index=np.arange(1, num+1))
        potential_ids = np.where(sizes >= req_pixels)[0] + 1
        valid_ids_final = []
        
        if len(potential_ids) > 0:
            dy_ds, dx_ds = np.gradient(elevation[::downsample_factor, ::downsample_factor], eff_cell_y, eff_cell_x)
            aspect_ds = np.degrees(np.arctan2(-dx_ds, dy_ds)) % 360

            # Passed in metres: converting to a pixel stride here was one of the three
            # truncations that inflated capacity (see count_grid_capacity)
            spacing_m = antenna_spacing_km * 1000.0
            grid_code = 1 if grid_type == 'hex' else 0
            all_slices = find_objects(labeled)

            # Mean aspect for every candidate region in one labelled pass. Doing it as
            # `aspect_ds[labeled == site_id]` inside the loop re-scanned the whole
            # downsampled map once per site, which is O(sites x pixels).
            mean_aspects = np.atleast_1d(ndi_mean(aspect_ds, labeled, index=potential_ids))
            aspect_by_id = {int(sid): float(m)
                            for sid, m in zip(potential_ids, mean_aspects)}

            # Iterate through found blobs to calculate physical internal placement of DUs
            for site_id in potential_ids:
                loc = all_slices[site_id - 1]
                r_start = loc[0].start * downsample_factor
                r_stop = loc[0].stop * downsample_factor
                c_start = loc[1].start * downsample_factor
                c_stop = loc[1].stop * downsample_factor
                
                r_stop = min(r_stop, rows)
                c_stop = min(c_stop, cols)
                
                mask_chunk = final_map_disk[r_start:r_stop, c_start:c_stop]

                # Restrict to this region's own pixels. A bounding box is not the
                # region: it also contains whatever else happens to fall inside it,
                # including other sites and regions that failed the area threshold.
                # Counting those attributes their detectors to this site as well, and
                # since the totals are summed over sites the same ground is then sold
                # twice. One compact site barely notices; a canyon network of thirty
                # interleaved strips inflated its total by about 38%.
                own = labeled[loc] == site_id
                if downsample_factor > 1:
                    own = np.repeat(np.repeat(own, downsample_factor, axis=0),
                                    downsample_factor, axis=1)
                h = min(own.shape[0], mask_chunk.shape[0])
                w = min(own.shape[1], mask_chunk.shape[1])
                site_chunk = np.logical_and(mask_chunk[:h, :w], own[:h, :w])
                antennas_fit = count_grid_capacity(site_chunk, cell_size_y, cell_size_x,
                                                   spacing_m, grid_code)
                
                if antennas_fit >= threshold_antennas:
                    valid_ids_final.append(site_id)
                    mean_aspect = aspect_by_id[int(site_id)]
                    dirs = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
                    aspect_str = dirs[round(mean_aspect / 45) % 8]
                    area_km2 = sizes[site_id-1] * px_area_km2

                    # Where the site actually is. The records carried area, capacity
                    # and facing but no position, so answering "which ground is this?"
                    # meant opening the raster in a GIS. The bounding box is already
                    # in hand, and the centroid costs one pass over a region that has
                    # just been scanned anyway.
                    where = {}
                    if origin_lat is not None and cell_size_deg:
                        step = cell_size_deg * downsample_factor
                        own_ds = labeled[loc] == site_id
                        rr, cc = np.nonzero(own_ds)
                        centre_row = loc[0].start + (rr.mean() if rr.size else 0.0)
                        centre_col = loc[1].start + (cc.mean() if cc.size else 0.0)
                        where = {
                            "center_lat": float(f"{origin_lat - centre_row * step:.6f}"),
                            "center_lon": float(f"{origin_lon + centre_col * step:.6f}"),
                            "bounds": {
                                "north": float(f"{origin_lat - loc[0].start * step:.6f}"),
                                "south": float(f"{origin_lat - loc[0].stop * step:.6f}"),
                                "west": float(f"{origin_lon + loc[1].start * step:.6f}"),
                                "east": float(f"{origin_lon + loc[1].stop * step:.6f}"),
                            },
                        }
                    
                    site_details.append({
                        "site_id": int(site_id),
                        "area_km2": float(f"{area_km2:.2f}"),
                        "capacity_exact": int(antennas_fit),
                        "grid_type": grid_type,
                        "mean_aspect_deg": float(f"{mean_aspect:.1f}"),
                        "facing_direction": aspect_str,
                        **where,
                        # Set below, once selection has run. Every site that clears the
                        # thresholds is listed, because the ones just below the cut are
                        # worth seeing -- but only the selected ones are in the exported
                        # mask, in total_sites and in total_capacity, so the record has
                        # to say which it is rather than leaving a reader to infer it.
                        "selected": False,
                    })

        site_details.sort(key=lambda x: x['capacity_exact'], reverse=True)
        final_selection_ids = []
        
        if search_mode == 'distributed':
            # Sites are already sorted by capacity, so this takes the best ones. With
            # stop_at_target the run reports the array actually wanted rather than
            # every patch of qualifying ground -- "the best sites for 5000 detectors"
            # is a different and usually more useful question than "all terrain that
            # passes", and it does not depend on where a score threshold was put.
            for site in site_details:
                if stop_at_target and cumulative_capacity >= target_antennas:
                    break
                final_selection_ids.append(site['site_id'])
                cumulative_capacity += site['capacity_exact']
        else:
            final_selection_ids = [s['site_id'] for s in site_details]

        chosen = set(final_selection_ids)
        for site in site_details:
            site['selected'] = site['site_id'] in chosen

        if len(final_selection_ids) > 0:
            # Recolour by table lookup: one pass over the labelled map instead of one
            # pass per selected site (plus a second for np.isin). The table is indexed
            # by original label, so `labeled` is read exactly once.
            lut = np.zeros(num + 1, dtype=viz_dtype)
            for current_viz_id, original_id in enumerate(final_selection_ids, start=1):
                lut[original_id] = current_viz_id
            labeled_viz = lut[labeled]
            small_final = (labeled_viz > 0).astype(np.uint8)
            count = len(final_selection_ids)

            # Fold the scan observables of each site's candidates into its record
            per_site = summarize_observables_by_site(
                labeled, downsample_factor, candidates_arr, observables, final_selection_ids)
            for site in site_details:
                if site['site_id'] in per_site:
                    site['arrival_scan'] = per_site[site['site_id']]

    region_stats = {
        "labelled_regions": int(num),
        "passed_area_threshold": int(len(potential_ids)) if num > 0 else 0,
        "passed_capacity_threshold": int(len(site_details)),
        "selected": int(count),
        "required_pixels_per_region": int(req_pixels),
        "capacity_threshold_antennas": int(threshold_antennas),
    }
    if funnel is not None:
        # small_final is downsampled; scale back to full-resolution pixels
        funnel.add("pixels in selected sites (est.)",
                   int(small_final.sum()) * downsample_factor * downsample_factor)

    return small_final, labeled_viz, site_details, cumulative_capacity, count, region_stats

def create_world_file(tif_filename, top_left_lat, top_left_lon, cell_size_deg):
    """
    Creates an ESRI World File (.tfw) which accompanies a standard TIFF image, 
    allowing GIS software (like QGIS or ArcGIS) to project it correctly on a map.

    Parameters
    ----------
    tif_filename : str
        Path to the raster the world file accompanies. The ``.tfw`` is written beside
        it with a matching stem.
    top_left_lat, top_left_lon : float
        Coordinates of the raster's north-west corner, in degrees.
    cell_size_deg : float
        Pixel size in degrees, after any downsampling.
    """
    tfw_name = os.path.splitext(tif_filename)[0] + ".tfw"
    try:
        with open(tfw_name, "w") as f:
            # Format: Pixel X size, Rotation, Rotation, Negative Pixel Y size, Top-Left X, Top-Left Y
            f.write(f"{cell_size_deg:.10f}\n0.0\n0.0\n-{cell_size_deg:.10f}\n{top_left_lon:.10f}\n{top_left_lat:.10f}\n") 
    except Exception: pass

def generate_kml_file(mask, elevation, filename, origin_lat, origin_lon, cell_size_deg, downsample=1):
    """
    Generates a Google Earth compatible KML file representing the valid site polygons.
    It extracts polygon contours from the binary mask using Matplotlib's contour tool.
    
    Parameters:
    - mask (ndarray): Binary mask indicating valid deployment sites.
    - filename (str): Output path for the KML file.
    - origin_lat, origin_lon, cell_size_deg: Used to convert array pixel indices to GPS coordinates.

    Parameters
    ----------
    mask : ndarray
        Boolean site mask.
    elevation : ndarray
        The DEM, used to place the contours in height.
    filename : str
        Path to write the ``.kml`` to.
    origin_lat, origin_lon : float
        North-west corner of the mask, in degrees.
    cell_size_deg : float
        Pixel size in degrees.
    downsample : int, optional
        Factor by which ``mask`` is already downsampled relative to the DEM.
    """
    print(f"      {Icon.INFO}Generating KML: {os.path.basename(filename)} ...")
    
    try:
        root = ET.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
        doc = ET.SubElement(root, "Document")
        
        # Define visual styles for the KML polygons (Yellow, semi-transparent)
        style = ET.SubElement(doc, "Style", id="grand_site")
        lstyle = ET.SubElement(style, "LineStyle")
        ET.SubElement(lstyle, "color").text = "ffffff00" 
        ET.SubElement(lstyle, "width").text = "3"
        pstyle = ET.SubElement(style, "PolyStyle")
        ET.SubElement(pstyle, "color").text = "40ffff00" 
        
        # Use Matplotlib to trace the boundaries of the mask areas
        fig = plt.figure()
        ax = fig.add_subplot(111)
        contours = ax.contour(mask, levels=[0.5])
        plt.close(fig)
        
        # Iterate over traced contours to build KML geometry blocks
        site_idx = 1
        for path in contours.get_paths():
            placemark = ET.SubElement(doc, "Placemark")
            ET.SubElement(placemark, "name").text = f"Site {site_idx}"
            ET.SubElement(placemark, "styleUrl").text = "#grand_site"
            
            poly = ET.SubElement(placemark, "Polygon")
            outer = ET.SubElement(poly, "outerBoundaryIs")
            ring = ET.SubElement(outer, "LinearRing")
            coords_str = ""
            
            # Map Matplotlib vertices back to real-world Long/Lat coordinates
            for (c, r) in path.vertices:
                r_full = r * downsample
                c_full = c * downsample
                lat = origin_lat - (r_full * cell_size_deg)
                lon = origin_lon + (c_full * cell_size_deg)
                coords_str += f"{lon},{lat},0 "
            
            ET.SubElement(ring, "coordinates").text = coords_str
            site_idx += 1
                
        ET.ElementTree(root).write(filename, encoding='UTF-8', xml_declaration=True)
    except Exception as e:
        print(f"      {C.WARN}{Icon.WARN}WARNING: KML generation failed (Skipping KML). Error: {e}{C.RESET}")

def generate_visualizations_and_outputs(dem_path, elevation, small_final, labeled_viz, site_details, count, cumulative_capacity,
                                        origin_lat, origin_lon, map_grid, downsample_factor, generate_kml, run_output_dir,
                                        output_image_format, rfi_zones, search_mode, grid_type, antenna_spacing_km, 
                                        min_altitude, max_altitude, region_name, final_params, run_info=None):
    """
    Step 6 Pipeline: Formats and exports all scientific products including geo-registered TIFs, KML models, 
    an annotated map graphic, and a serialized JSON summary of the run parameters and results 
    to the designated unified output directory.

    Parameters
    ----------
    dem_path : str
        Path to the DEM, used for the output stem.
    elevation : ndarray
        The DEM, as the map background.
    small_final : ndarray
        Downsampled binary mask of the selected sites.
    labeled_viz : ndarray
        Site labels for colour coding.
    site_details : list of dict
        Per-site records.
    count : int
        Number of selected sites.
    cumulative_capacity : int
        Total detector capacity across them.
    origin_lat, origin_lon : float
        North-west corner of the DEM, in degrees.
    map_grid : MapGrid
        Resolved grid geometry.
    downsample_factor : int
        Factor relating ``small_final`` to the DEM.
    generate_kml : bool
        Also write a Google Earth ``.kml``.
    run_output_dir : str
        Directory to write into.
    output_image_format : str
        Extension for the overview map, such as ``png`` or ``pdf``.
    rfi_zones : sequence
        Exclusion zones, drawn on the map.
    search_mode : str
        ``single`` or ``distributed``.
    grid_type : str
        ``square`` or ``hex``.
    antenna_spacing_km : float
        Detector spacing, in km.
    min_altitude, max_altitude : float or None
        Altitude bounds applied, for the annotation.
    region_name : str
        Human-readable region label.
    final_params : dict
        Resolved parameters, serialised into the results JSON.
    run_info : dict, optional
        Timings, funnel and provenance to record alongside the results.

    Returns
    -------
    generated_files : list of str
        Absolute paths of everything written.
    out_data : dict
        The results as serialised into the JSON. Returned as well as written so the
        caller does not have to find and re-read the file it was just handed the path
        to, which is what every caller was doing.
    """
    generated_files = []
    cell_size_deg = map_grid.cell_size_deg
    base_filename = RESULTS_PREFIX + os.path.splitext(os.path.basename(dem_path))[0]
    
    # Save TIF
    out_tif = os.path.join(run_output_dir, base_filename + ".tif")
    tiff.imwrite(out_tif, small_final)
    generated_files.append(os.path.abspath(out_tif))
    
    # Save TFW (the exported raster is downsampled, so its pixels are that much larger)
    new_res_deg = cell_size_deg * downsample_factor
    create_world_file(out_tif, origin_lat, origin_lon, new_res_deg)
    generated_files.append(os.path.abspath(os.path.splitext(out_tif)[0] + ".tfw"))
    
    # Save KML
    if generate_kml:
        kml_name = os.path.join(run_output_dir, base_filename + ".kml")
        generate_kml_file(small_final, elevation, kml_name, origin_lat, origin_lon, new_res_deg)
        generated_files.append(os.path.abspath(kml_name))

    # Save Custom Visualization
    try:
        fig, ax = plt.subplots(figsize=(14, 12))
        viz_ds = downsample_factor * 2 
        elev_viz = elevation[::viz_ds, ::viz_ds]
        mask_viz = small_final[::2, ::2] 
        
        mask_viz_labeled = labeled_viz[::2, ::2]
        
        mr = min(elev_viz.shape[0], mask_viz.shape[0])
        mc = min(elev_viz.shape[1], mask_viz.shape[1])
        elev_viz = elev_viz[:mr, :mc]
        mask_viz_labeled = mask_viz_labeled[:mr, :mc]
        
        # Scaled to this DEM rather than a fixed 0-6000 m, so the whole colour range
        # describes ground that is actually in the picture. Water is drawn as water
        # instead of as the bottom of the terrain ramp, where it reads as low land.
        vmin, vmax = altitude_limits(elev_viz)
        terrain = plt.get_cmap('terrain').with_extremes(under=WATER_COLOUR,
                                                       bad=NODATA_COLOUR)
        # NaN survives this comparison as NaN -- `NaN <= x` is False -- so nodata
        # reaches the colormap as "bad" and water as "under", and the two stay
        # distinguishable. Collapsing both to water claimed sea where the DEM is
        # simply silent.
        land_only = np.where(elev_viz <= SEA_LEVEL_M, -np.inf, elev_viz)
        im = ax.imshow(land_only, cmap=terrain, vmin=vmin, vmax=vmax)

        legend_handles = []
        legend_labels = []

        if count > 0:
            cmap = plt.get_cmap('tab10')
            for i in range(1, count + 1):
                color = cmap((i - 1) % 10)
                ax.contour((mask_viz_labeled == i), levels=[0.5], colors=[color], linewidths=2.5)

                site_data = site_details[i - 1]
                label_str = f"Site {site_data['site_id']}: {site_data['capacity_exact']} DUs ({site_data['area_km2']} km²)"
                legend_handles.append(Line2D([0], [0], color=color, lw=2.5))
                legend_labels.append(label_str)

                # The number on the map, so a site in the legend can be found in the
                # picture. At its own centroid, which for a ring-shaped region may not
                # be inside it -- acceptable for a label, and better than a corner.
                where = np.argwhere(mask_viz_labeled == i)
                if where.size:
                    cy, cx = where.mean(axis=0)
                    tag = ax.text(cx, cy, str(site_data['site_id']), color=color,
                                  fontsize=11, fontweight='bold', ha='center',
                                  va='center', zorder=12, clip_on=True)
                    tag.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),
                                          path_effects.Normal()])
        
        if rfi_zones:
            deg_viz = cell_size_deg * viz_ds
            viz_rows, viz_cols = mask_viz_labeled.shape[:2]

            def on_map(x0, x1, y0, y1):
                """Whether a zone's bounding box touches the raster at all."""
                return not (x1 < 0 or x0 > viz_cols or y1 < 0 or y0 > viz_rows)

            drawn, off_map = 0, []
            for item in rfi_zones:
                type_tag = item[0]
                if type_tag == 'circle':
                    _, lat, lon, radius_km, name = item
                    px_x = (lon - origin_lon) / deg_viz
                    px_y = (origin_lat - lat) / deg_viz
                    # A circle on the ground is an ellipse on an angular pixel grid, since
                    # a pixel covers less ground east-west than north-south
                    w_px = 2.0 * (radius_km * 1000.0 / map_grid.cell_size_x) / viz_ds
                    h_px = 2.0 * (radius_km * 1000.0 / map_grid.cell_size_y) / viz_ds
                    if not on_map(px_x - w_px / 2, px_x + w_px / 2,
                                  px_y - h_px / 2, px_y + h_px / 2):
                        off_map.append(name)
                        continue
                    ax.add_patch(Ellipse((px_x, px_y), w_px, h_px, edgecolor='red',
                                         facecolor='none', ls='--', lw=2, clip_on=True))
                    text = ax.text(px_x, px_y-h_px/4, name, color='red', fontsize=12, ha='center', clip_on=True)
                    text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'), path_effects.Normal()])
                    drawn += 1
                elif type_tag == 'poly':
                    _, coords, name = item
                    verts = []
                    for (plat, plon) in coords:
                        px = (plon - origin_lon) / deg_viz
                        py = (origin_lat - plat) / deg_viz
                        verts.append((px, py))
                    xs = [p[0] for p in verts]
                    ys = [p[1] for p in verts]
                    if not on_map(min(xs), max(xs), min(ys), max(ys)):
                        off_map.append(name)
                        continue
                    ax.add_patch(MplPolygon(verts, closed=True, edgecolor='red',
                                            facecolor='none', ls='--', lw=2, clip_on=True))
                    cx = sum(xs)/len(xs)
                    cy = sum(ys)/len(ys)
                    text = ax.text(cx, cy, name, color='red', fontsize=8, ha='center', clip_on=True)
                    text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'), path_effects.Normal()])
                    drawn += 1

            # Zones outside the DEM are dropped rather than drawn. They still exclude
            # nothing inside it, and drawing them wrecked the figure: an artist beyond
            # the image expands the axes, and `bbox_inches='tight'` then grew the saved
            # PNG to reach a label 150 km off the south edge -- leaving the map itself
            # in the top fifth of a mostly empty page.
            if drawn:
                legend_handles.append(Line2D([0], [0], color='red', linestyle='--', lw=2))
                legend_labels.append("RFI exclusion zone")
            if off_map:
                print(f"      {Icon.INFO}{len(off_map)} RFI zone(s) lie outside the DEM "
                      f"and are not drawn: {', '.join(off_map)}")

        # Hold the view to the raster. Every RFI ellipse is an artist, and an artist
        # outside the image expands the axes to contain it -- so a zone like Mollendo,
        # 150 km off the south edge of the Colca crop, pushed the map into the top
        # fifth of the frame and filled the rest with white. The zones are still drawn;
        # the ones off the map are simply off the map, which is what they are.
        viz_rows, viz_cols = mask_viz_labeled.shape[:2]
        ax.set_xlim(-0.5, viz_cols - 0.5)
        ax.set_ylim(viz_rows - 0.5, -0.5)

        deg_viz = cell_size_deg * viz_ds
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p: f"{origin_lon + x*deg_viz:.2f}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y,p: f"{origin_lat - y*deg_viz:.2f}"))
        ax.set_xlabel("Longitude (°)")
        ax.set_ylabel("Latitude (°)")
        # The axes are pixels, so the bar needs the metric pixel size, not a latitude.
        add_scale_bar(ax, map_grid.cell_size_x * viz_ds / 1000.0)
        cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02, extend='min')
        cbar.set_label('Altitude (m)', rotation=270, labelpad=15)
        add_north_arrow(ax)

        # No title. Everything it carried is either on the figure already (the sites,
        # in the legend), in the run's own summary, or in the caption of whatever this
        # ends up in -- and a two-line title over a map is a caption in the wrong place.
        #
        # The legend goes outside the axes for the same reason: it was a filled box
        # sitting on top of the data it describes, and on a dense map it hid whichever
        # sites happened to be in the top right.
        fs = 'small' if len(legend_labels) > 8 else 'medium'
        if legend_labels:
            ax.legend(legend_handles, legend_labels, fontsize=fs, framealpha=0.9,
                      loc='upper left', bbox_to_anchor=(0.0, -0.06),
                      ncol=max(1, min(3, len(legend_labels))), borderaxespad=0.0)
        
        img_name = os.path.join(run_output_dir, base_filename + "." + output_image_format.strip('.'))
        
        plt.savefig(img_name, format=output_image_format.strip('.'), dpi=150, bbox_inches='tight')
        generated_files.append(os.path.abspath(img_name))
        print(f"      {Icon.CHECK}Map saved.")

    except Exception as e:
        print(f"      {C.FAIL}{Icon.CROSS}Viz Error: {e}{C.RESET}")
    finally:
        # pyplot holds a global reference to every figure it creates, so one that is
        # never closed can never be collected. A single search does not notice; a
        # process that runs several -- a parameter sweep, a notebook, a service --
        # accumulates the whole figure each time, artists and canvas included. That is
        # what took a 10-point sensitivity sweep to 6.9 GB and into the OOM killer.
        #
        # In the `finally` because the failure path leaks just as readily as the happy
        # one, and an exception here is caught and reported rather than fatal.
        plt.close('all')

    # Save JSON output log
    run_info = run_info or {}
    out_data = {
        "timestamp": datetime.now().isoformat(),
        "mode": search_mode,
        "parameters": final_params,
        "results": {
            "total_sites": count,
            "total_capacity": cumulative_capacity if search_mode=='distributed' else 'N/A',
            "sites": site_details
        },
        "funnel": run_info.get("funnel", {}),
        "regions": run_info.get("regions", {}),
        "timings_sec": run_info.get("timings_sec", {}),
        "aperture": run_info.get("aperture", {}),
    }
    json_name = os.path.join(run_output_dir, base_filename + ".json")
    with open(json_name, "w") as f:
        json.dump(out_data, f, indent=4)
    generated_files.append(os.path.abspath(json_name))
    print(f"      {Icon.CHECK}JSON Data Summary saved.")

    # Provenance goes in its own file so it stays readable next to the science outputs
    if "provenance" in run_info:
        prov_name = os.path.join(run_output_dir, "provenance.json")
        with open(prov_name, "w") as f:
            json.dump(run_info["provenance"], f, indent=4)
        generated_files.append(os.path.abspath(prov_name))
        print(f"      {Icon.CHECK}Provenance saved.")

    return generated_files, out_data

def print_tool_explanation():
    """
    Prints what the tool is about to do, before it does it.

    Suppressed with ``--no_print_info``. Kept current deliberately: this described a
    single ray cast to a target mountain, and a clearance buffer over intervening
    terrain, for some time after both had been replaced by the arrival-direction scan.
    """
    print(f"""
{C.HEADER}================================================================================
{C.BOLD}OROSCOPE - TERRAIN SITE SEARCH FOR PARTICLE-ASTROPHYSICS OBSERVATORIES{C.RESET}{C.HEADER}
================================================================================{C.RESET}
Searches a digital elevation model for ground that can host an observatory, against
one experiment's criteria. GRAND and TAMBO are configurations of the same engine, not
separate code paths: adding an experiment means writing a JSON file.

{C.BOLD}The question it answers, for every patch of ground:{C.RESET}
  Is there a target surface at the right range, in the right direction, at the right
  relative orientation, with the right matter behind it?

{C.BOLD}Core workflow:{C.RESET}
1. Topographic screen: slope (default 3-25 degrees), altitude, facing direction,
   distance to roads, and radio-quiet exclusion zones. Survivors are thinned by
   --candidate_stride, which was measured unbiased against a stride-1 control.
2. Arrival scan: from each survivor, walks terrain profiles along a fan of azimuths
   and reports what each arrival direction meets -- accepted solid angle, distance to
   the exit point, column depth of rock, horizon, atmospheric depth, Earth chord and
   the slope of the terrain struck. One walk serves every elevation bin, so the
   azimuth count is what sets the cost. Earth curvature throughout; the radio path
   alone uses the 4/3 refraction convention.
3. Scoring: each observable becomes a named component in [0, 1] -- depth, distance,
   solid angle, shower, decay, geomagnetic, clearance -- combined by --score_composition.
   Naming them is what lets a weak site be attributed rather than merely reported.
4. Spatial pruning: morphological closing fills gaps between accepted pixels, opening
   removes tendrils narrower than --min_width_km. Closing inflates reported area, so
   --gap_close_km is worth setting deliberately.
5. Sites and capacity: connected regions are labelled and a detector lattice is packed
   into each, in 'hex' or 'square'. Regions below the area or capacity threshold are
   dropped; --stop_at_target reports the best sites for the array actually wanted.
6. Outputs: georeferenced GeoTIFF and world file, KML, annotated map, a results JSON,
   a selection funnel, a provenance record (git commit, DEM checksum, versions), the
   full log, and a plain-language summary of what was found and why.

{C.BOLD}Two things worth knowing before reading the numbers:{C.RESET}
- The funnel is the diagnostic. When a search returns little or nothing, the stage
  where the survivor count collapses IS the constraint responsible.
- Reported area is not physics-accepted area. Morphological closing inflated it 2.29x
  at Colca, measured against a stride-1 control. The run summary reports the factor
  for the run in front of you.

{C.BOLD}Memory:{C.RESET} the DEM is memory-mapped and processed in tiles (--tile_size), the run
estimates its own peak against available RAM and caps its address space
(--max_memory_gb), and --resume skips the scan on a failed run.
{C.HEADER}================================================================================{C.RESET}
    """)

def explicitly_passed(parser, argv=None):
    """
    The set of options the user actually typed, as opposed to argparse's defaults.

    argparse gives no way to distinguish ``--candidate_stride 5`` from the default of
    5, which is why the configuration merge used to prefer a config file over the
    command line: with no way to tell a typed flag from an untyped one, honouring the
    command line would have let every default silently overwrite the config.

    Re-parsing with every default suppressed answers the question directly --- with
    ``SUPPRESS``, argparse only sets an attribute for an option that actually appeared.
    The defaults are restored afterwards, so the original ``args`` is untouched.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser to interrogate. Left exactly as it was found.
    argv : list of str, optional
        Arguments to parse. Defaults to ``sys.argv``.

    Returns
    -------
    set of str
        Destinations of the options that actually appeared on the command line.
    """
    saved = [(action, action.default) for action in parser._actions]
    try:
        for action, _ in saved:
            action.default = argparse.SUPPRESS
        return set(vars(parser.parse_args(argv)))
    finally:
        for action, default in saved:
            action.default = default


def _one_or_pair(value):
    """
    Normalises a parameter that may be a single number or a (low, high) range.

    The command line gives ``nargs="+"``, so one value arrives as a one-element list; a
    config file may give a bare number or a two-element list. All three mean the same
    things, and the physics accepts either form, so this only has to unwrap the
    one-element case.

    Parameters
    ----------
    value : float, sequence of float, or None
        As supplied by the configuration or the command line.

    Returns
    -------
    float, tuple of float, or None
        A scalar for one value, a ``(low, high)`` tuple for two, ``None`` for nothing.

    Raises
    ------
    SystemExit
        If more than two values are given, which no parameter here means.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> ss._one_or_pair([2.0])
    2.0
    >>> ss._one_or_pair([1.5, 2.7])
    (1.5, 2.7)
    >>> ss._one_or_pair(None) is None
    True
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    values = [float(v) for v in value]
    if len(values) == 1:
        return values[0]
    if len(values) == 2:
        return (values[0], values[1])
    raise SystemExit(f"expected one value or a (low, high) pair, got {len(values)}")


def parse_score_weights(value):
    """
    Normalises per-component score weights from either input form.

    A config file is JSON, so it can carry a mapping directly. The command line
    cannot, so it takes ``shower=2,solid_angle=1`` instead. Both end up as a dict, and
    anything unnamed keeps weight 1.

    Parameters
    ----------
    value : str, dict or None
        Either ``name=value`` pairs separated by commas, or a mapping, or ``None``.

    Returns
    -------
    dict or None
        Component name to weight, or ``None`` when nothing was supplied, which leaves
        the composition unweighted.

    Raises
    ------
    SystemExit
        If a pair lacks ``=`` or its value is not a number.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> ss.parse_score_weights("shower=2,depth=0.5") == {"shower": 2.0, "depth": 0.5}
    True
    >>> ss.parse_score_weights(None) is None
    True
    """
    if value is None or value == "":
        return None
    if isinstance(value, dict):
        return {str(k): float(v) for k, v in value.items()}
    weights = {}
    for pair in str(value).split(","):
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            raise SystemExit(f"--score_weights expects name=value pairs, got {pair!r}")
        name, _, raw = pair.partition("=")
        try:
            weights[name.strip()] = float(raw)
        except ValueError:
            raise SystemExit(f"--score_weights value for {name.strip()!r} is not a number: {raw!r}")
    return weights or None


def estimate_peak_memory_gb(rows, cols, downsample_factor=1, candidate_stride=5,
                            survival_fraction=0.6, n_observables=12,
                            n_scoring_arrays=24):
    """
    Rough estimate of the anonymous memory one search will need, in GiB.

    Only the allocations that can exhaust RAM are counted. The DEM itself is
    memory-mapped and file-backed, so the kernel can evict it under pressure and it is
    excluded deliberately -- counting it would make every large search look impossible
    when the streaming design exists precisely so that it is not.

    This is an estimate and says so. ``survival_fraction`` in particular is the fraction
    of pixels passing the topographic screen, which is terrain-dependent and not known
    until the screen has run; 0.6 is typical of Andean terrain at a 3-25 degree band.
    It is meant to catch the order-of-magnitude mistake -- a full DEM at
    ``downsample_factor: 1`` -- rather than to predict a number.

    Notes
    -----
    **The peak is in the scoring, not in the scan.** This counted only the arrays
    ``arrival_scan.scan`` returns, which is not where the high-water mark is: by the
    time :func:`scoring.compose` runs, the scan's arrays are still live, a score
    component has been built alongside them for each criterion, ``compose`` has clipped
    a float64 *copy* of every component, and the composition and the scoring
    intermediates need several more. About three times the scan's own count is live at
    once, all of it ``n_cand`` long.

    Under-counting that term is not academic: it advertised 2.32 GiB for the full
    Arequipa DEM, which then peaked at **5.68 GiB measured RSS** and died against its
    own cap 23 minutes in. ``n_scoring_arrays`` is calibrated on that run -- 15.1M
    candidates, 7 components -- where the anonymous share of the peak implies ~36
    live per-candidate arrays against the 12 this modelled.

    Note also which knob moves it. ``downsample_factor`` scales only the labelling and
    gradient terms, because candidates are taken on the *native* grid; at full-DEM
    scale the per-candidate terms dominate, so going from 1 to 4 cuts the estimate by
    about 1.4x rather than the 16x the inverse-square scaling suggests. To move the
    dominant term, raise ``candidate_stride`` or crop the DEM.

    Parameters
    ----------
    rows, cols : int
        DEM dimensions in pixels.
    downsample_factor : int, optional
        Factor at which sites are labelled and areas measured. Scales the labelling
        arrays as its inverse square, and nothing else -- see the note above.
    candidate_stride : int, optional
        Keeps every Nth screened pixel. Scales the dominant term directly.
    survival_fraction : float, optional
        Fraction of pixels expected to pass the topographic screen.
    n_observables : int, optional
        Per-candidate arrays the scan returns.
    n_scoring_arrays : int, optional
        Further per-candidate arrays live at the peak, inside ``compose``: the score
        components, the float64 copy made of each, and the temporaries.

    Returns
    -------
    float
        Estimated peak anonymous memory, in GiB.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> round(ss.estimate_peak_memory_gb(1981, 3061, downsample_factor=1), 2)
    0.77
    >>> round(ss.estimate_peak_memory_gb(10204, 12603, downsample_factor=4), 2)
    5.08

    Downsampling is the weaker of the two levers at this scale, and striding the
    stronger, because the candidates are taken on the native grid either way:

    >>> round(ss.estimate_peak_memory_gb(10204, 12603, downsample_factor=1), 2)
    7.21
    >>> round(ss.estimate_peak_memory_gb(10204, 12603, downsample_factor=4,
    ...                                  candidate_stride=10), 2)
    2.83
    """
    n_pixels = float(rows) * float(cols)
    n_small = n_pixels / float(max(1, downsample_factor) ** 2)

    # Labelling and per-site geometry, all on the downsampled map
    labelling = n_small * (1 + 4 + 2)            # mask, int32 labels, viz
    gradients = n_small * 4 * 3                  # d/dy, d/dx, aspect, float32

    # Candidates and their observables, at full resolution. Candidates are taken on the
    # native grid -- the stride subsamples the surviving-pixel list, not the map -- so
    # downsample_factor does not touch these two terms at all.
    n_cand = n_pixels * survival_fraction / float(max(1, candidate_stride))
    candidates = n_cand * 3 * 8
    observables = n_cand * n_observables * 8

    # The high-water mark, inside compose(): components, their clipped copies, and the
    # temporaries, all still holding the scan's arrays above. See the notes.
    scoring = n_cand * n_scoring_arrays * 8

    # Interpreter, numba, matplotlib and the tiled screening buffers
    baseline = 0.45 * 1024 ** 3

    total = labelling + gradients + candidates + observables + scoring + baseline
    return total / 1024 ** 3


SEA_LEVEL_M = 0.5          # at or below this is water, not ground at zero altitude
WATER_COLOUR = "#5A7FA6"
# Nodata is neither water nor low ground, and painting it as either is a claim the DEM
# does not support. A warm neutral, so it cannot be mistaken for the grey ramp or the sea.
NODATA_COLOUR = "#E0DACE"


def altitude_limits(elevation, low_percentile=0.5, high_percentile=99.8):
    """
    Altitude range for a colour scale, from the DEM rather than from a constant.

    A fixed 0-6000 m scale spends most of its range on altitudes a given DEM does not
    contain: the Colca crop runs 1500-6300 m, so half the colour bar described ground
    that is not in the picture and the relief that is there got half the contrast it
    could have had.

    Percentiles rather than the extremes, because one spurious pixel -- a nodata
    sentinel, a spike -- otherwise sets the whole scale. Water is excluded from the
    upper end and pinned to the bottom, so a coastal DEM does not spend a third of its
    range on ocean.

    Parameters
    ----------
    elevation : ndarray
        Elevation in metres. NaN is ignored.
    low_percentile, high_percentile : float, optional
        Percentiles of the land pixels to clip to.

    Returns
    -------
    tuple of float
        ``(vmin, vmax)`` in metres, rounded outward to a round number.

    Examples
    --------
    >>> import numpy as np
    >>> from oroscope import site_searcher as ss
    >>> z = np.linspace(1500.0, 6300.0, 1000)
    >>> ss.altitude_limits(z)
    (1500.0, 6300.0)

    Ocean does not drag the floor down with it:

    >>> z = np.concatenate([np.zeros(500), np.linspace(2000.0, 5000.0, 500)])
    >>> ss.altitude_limits(z)
    (0.0, 5000.0)
    """
    values = np.asarray(elevation, dtype=np.float64).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:                             # pragma: no cover - empty DEM
        return (0.0, 6000.0)
    land = values[values > SEA_LEVEL_M]
    if land.size == 0:                               # pragma: no cover - all water
        return (0.0, 1.0)
    lo = float(np.percentile(land, low_percentile))
    hi = float(np.percentile(land, high_percentile))
    if (values <= SEA_LEVEL_M).any():
        lo = 0.0
    step = 100.0
    return (math.floor(lo / step) * step, math.ceil(hi / step) * step)


def add_north_arrow(ax, x=0.965, y=0.955, size=0.055):
    """
    Draws a north arrow in axes coordinates.

    Both maps this project writes are north-up, so the arrow is a convention rather
    than information -- but a map without one asks the reader to assume, and a map in a
    talk gets read by people who did not make it.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    x, y : float, optional
        Position of the arrow's tip, in axes coordinates.
    size : float, optional
        Length of the arrow, as a fraction of the axes height.
    """
    ax.annotate("N", xy=(x, y), xytext=(x, y - size), xycoords="axes fraction",
                textcoords="axes fraction", ha="center", va="top",
                fontsize=11, fontweight="bold", color="black", zorder=25,
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6,
                                shrinkA=0, shrinkB=0),
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.7))


def add_scale_bar(ax, km_per_x_unit, fraction=0.22, colour="black"):
    """
    Draws a kilometre scale bar on a map, and returns the length it chose.

    A map axis labelled in degrees or in pixels does not tell a reader how far anything
    is, and neither unit converts to distance without knowing where on the Earth it
    sits: a degree of longitude at Arequipa is 4% shorter than a degree of latitude,
    and a pixel is whatever the DEM says it is.

    Taking ``km_per_x_unit`` rather than a latitude keeps one function usable by both
    maps this project writes -- the search map, whose axes are pixels, and the
    combination overlay, whose axes are degrees.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on. Its limits must already be final; the bar is placed relative
        to them.
    km_per_x_unit : float
        Kilometres per unit of the x axis. For degrees of longitude that is
        ``111.32 * cos(latitude)``; for pixels it is the metric pixel size over 1000.
    fraction : float, optional
        Roughly what fraction of the map width the bar should span, before rounding to
        a human number.
    colour : str, optional
        Bar colour. The default reads on both terrain and shaded relief.

    Returns
    -------
    float
        Length of the bar drawn, in km. Always 1, 2 or 5 times a power of ten.

    Examples
    --------
    A two-degree map at 16 degrees south is about 214 km wide, so it gets a 50 km bar:

    >>> import matplotlib; matplotlib.use("Agg")
    >>> import matplotlib.pyplot as plt, numpy as np
    >>> from oroscope import site_searcher as ss
    >>> fig, ax = plt.subplots()
    >>> _ = ax.set_xlim(-73.0, -71.0); _ = ax.set_ylim(-17.0, -15.0)
    >>> ss.add_scale_bar(ax, 111.32 * np.cos(np.radians(-16.0)))
    50.0

    The same function on a pixel axis, 3061 pixels of 30 m:

    >>> _ = ax.set_xlim(0, 3061); _ = ax.set_ylim(1981, 0)
    >>> ss.add_scale_bar(ax, 0.030)
    20.0
    >>> plt.close(fig)
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    span_km = abs(x1 - x0) * km_per_x_unit
    if not np.isfinite(span_km) or span_km <= 0:     # pragma: no cover - defensive
        return 0.0

    target = span_km * fraction
    power = 10.0 ** math.floor(math.log10(max(target, 1e-9)))
    length_km = min((1.0, 2.0, 5.0, 10.0), key=lambda m: abs(m * power - target)) * power
    length_units = length_km / km_per_x_unit

    pad_x, pad_y = 0.04 * (x1 - x0), 0.05 * (y1 - y0)
    bx, by = x0 + pad_x, y0 + pad_y
    ax.plot([bx, bx + length_units], [by, by], color=colour, lw=3.4,
            solid_capstyle="butt", zorder=20)
    ax.plot([bx, bx + length_units], [by, by], color="white", lw=1.4,
            solid_capstyle="butt", zorder=21)
    ax.text(bx + 0.5 * length_units, by + 0.30 * pad_y, f"{length_km:g} km",
            ha="center", va="bottom", fontsize=9, color=colour, zorder=22,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.75))
    return float(length_km)


def stride_gap_m(candidate_stride, cell_size_y_m):
    """
    Distance between kept candidates, in metres, after striding.

    ``candidate_stride`` subsamples the *list* of surviving pixels rather than the map,
    so the gap it leaves is a stride's worth of pixels along a scanline.

    Parameters
    ----------
    candidate_stride : int
        Keeps every Nth surviving pixel.
    cell_size_y_m : float
        Metric pixel size, N-S.

    Returns
    -------
    float
        Gap between kept candidates, in metres.
    """
    return float(max(1, int(candidate_stride))) * float(cell_size_y_m)


def closing_element_m(gap_close_km, antenna_spacing_km):
    """
    Size of the morphological closing element in metres, defaulting as the pipeline does.

    Parameters
    ----------
    gap_close_km : float or None
        Closing element in km. ``None`` defaults to ``antenna_spacing_km``.
    antenna_spacing_km : float
        Detector spacing, which the closing element defaults to.

    Returns
    -------
    float
        Closing element size, in metres.
    """
    km = antenna_spacing_km if gap_close_km is None else gap_close_km
    return float(km) * 1000.0


def warn_stride_outruns_closing(candidate_stride, cell_size_y_m,
                                gap_close_km, antenna_spacing_km, quiet=False):
    """
    Warns when the closing element is too small to bridge the gaps striding leaves.

    Striding is unbiased in *acceptance* -- measured at both scales, 60.1% against
    60.1% for GRAND and 17.491% against 17.494% for TAMBO -- so it is tempting to treat
    it as free. It is not. Accepted pixels are marked one in ``candidate_stride``, and
    the mask is then closed morphologically before areas are measured. If the closing
    element is smaller than the gap the stride leaves, the mask never reconnects: it
    stays a scatter of isolated pixels, small regions fall below the size and capacity
    thresholds, and the reported area collapses.

    Measured at Colca with a 100 m element against a 154 m stride-5 gap: **83.6 km²
    reported against 396.9 km² at stride 1, a 4.75x under-report**, with acceptance
    identical to three decimal places. The same run at GRAND's 1 km element -- 32 px,
    against the same 154 m gap -- is unaffected, which is why this went unnoticed.

    The rule is simply that the element must outrun the gap. Raise ``gap_close_km``,
    lower ``candidate_stride``, or accept the area as a lower bound and say so.

    Parameters
    ----------
    candidate_stride : int
        Keeps every Nth surviving pixel.
    cell_size_y_m : float
        Metric pixel size, N-S.
    gap_close_km : float or None
        Closing element in km. ``None`` defaults to ``antenna_spacing_km``.
    antenna_spacing_km : float
        Detector spacing, which the closing element defaults to.
    quiet : bool, optional
        Suppress the printed warning, keeping the returned verdict.

    Returns
    -------
    dict or None
        ``{"gap_m", "element_m", "ratio"}`` when the element cannot bridge the gap,
        and ``None`` when it can.

    Examples
    --------
    GRAND's 1 km element easily bridges a stride-5 gap at 30 m pixels:

    >>> from oroscope import site_searcher as ss
    >>> ss.warn_stride_outruns_closing(5, 30.72, None, 1.0, quiet=True) is None
    True

    TAMBO's 100 m element does not, and that is the 4.75x:

    >>> r = ss.warn_stride_outruns_closing(5, 30.72, None, 0.1, quiet=True)
    >>> round(r["gap_m"]), round(r["element_m"]), round(r["ratio"], 2)
    (154, 100, 1.54)

    Closing disabled entirely is not this failure, so it does not warn:

    >>> ss.warn_stride_outruns_closing(5, 30.72, 0.0, 0.1, quiet=True) is None
    True
    """
    element = closing_element_m(gap_close_km, antenna_spacing_km)
    if element <= 0.0:                    # closing switched off deliberately
        return None
    gap = stride_gap_m(candidate_stride, cell_size_y_m)
    if element >= gap:
        return None

    if not quiet:
        print(f"{C.WARN}{Icon.WARN}The closing element ({element:.0f} m) is smaller "
              f"than the gap candidate_stride {int(candidate_stride)} leaves "
              f"({gap:.0f} m). Accepted pixels will not reconnect, so the reported "
              f"AREA will be an under-report while acceptance stays unbiased -- "
              f"measured 4.75x at Colca on TAMBO's settings. Raise gap_close_km, "
              f"lower candidate_stride, or read the area as a lower bound.{C.RESET}")
    return {"gap_m": gap, "element_m": element, "ratio": gap / element}


def apply_memory_cap(max_memory_gb):
    """
    Caps this process's address space, so a runaway fails instead of taking the machine.

    Without a cap, an over-large search does not fail: it grows until the kernel's OOM
    killer chooses a victim, which may well be something else the user cares about. A
    ``MemoryError`` inside this process is a far better outcome than a dead editor, and
    it names the run that caused it.

    Parameters
    ----------
    max_memory_gb : float or None
        Ceiling in GiB. ``None`` or non-positive leaves the limit alone.

    Returns
    -------
    bool
        Whether a cap was applied. False on platforms without ``RLIMIT_AS``.

    Notes
    -----
    ``RLIMIT_AS`` limits *virtual* address space, which is larger than resident memory:
    numba and BLAS reserve address space they never touch. Set it generously -- a
    little above physical RAM is usually right -- or it will refuse runs that would
    have fit.
    """
    if not max_memory_gb or max_memory_gb <= 0:
        return False
    try:
        import resource
    except ImportError:                                  # pragma: no cover - Windows
        return False
    limit = int(max_memory_gb * 1024 ** 3)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    if hard != resource.RLIM_INFINITY:
        limit = min(limit, hard)
    resource.setrlimit(resource.RLIMIT_AS, (limit, hard))
    return True


def available_memory_gb():
    """
    Memory the system can give us right now, in GiB, or None if it cannot be told.

    Reads ``MemAvailable`` from ``/proc/meminfo``, which accounts for reclaimable page
    cache and so is the figure that matters; ``free`` alone understates it badly on a
    machine that has been running a while.

    Returns
    -------
    float or None
        Available memory in GiB, or ``None`` on a platform without ``/proc/meminfo``.
    """
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return float(line.split()[1]) / 1024 ** 2
    except OSError:                                      # pragma: no cover - non-Linux
        return None
    return None                                          # pragma: no cover


def preflight_memory(dem_path, downsample_factor=1, candidate_stride=5,
                     max_memory_gb=None, quiet=False):
    """
    Estimates what a search will need, says so, and caps the address space.

    This ran only inside ``main()``, so a library user -- a sweep, a notebook, a
    service -- got neither the warning nor the cap unless they knew to ask for both.
    That is exactly the caller most likely to need them: a ten-point sweep once reached
    6.9 GB and was killed by the OOM killer, taking other work with it.

    Parameters
    ----------
    dem_path : str
        DEM whose dimensions set the estimate. An unreadable file skips the estimate
        but not the cap.
    downsample_factor : int, optional
        As the pipeline's. Dominates the estimate: the labelling arrays scale as its
        inverse square.
    candidate_stride : int, optional
        As the pipeline's.
    max_memory_gb : float, optional
        Ceiling in GiB. ``None`` uses 80% of what the system reports available, so the
        cap bites before the kernel does; 0 disables capping entirely.
    quiet : bool, optional
        Suppress the printed report, keeping the cap and the returned numbers.

    Returns
    -------
    dict
        ``{"estimate_gb", "available_gb", "cap_gb", "capped"}``. ``estimate_gb`` is
        ``None`` when the DEM could not be measured.

    Notes
    -----
    The estimate is rough by construction -- it assumes a survival fraction the
    topographic screen has not yet measured -- so an over-large search is warned about
    rather than refused.
    """
    rows = cols = None
    try:
        with tiff.TiffFile(dem_path) as handle:
            rows, cols = handle.pages[0].shape[:2]
    except Exception:
        pass

    have = available_memory_gb()
    need = None
    if rows and cols:
        need = estimate_peak_memory_gb(rows, cols,
                                       downsample_factor=int(downsample_factor or 1),
                                       candidate_stride=int(candidate_stride or 5))
        if not quiet:
            print(f"   {Icon.GEAR}Estimated peak memory: {C.MAGENTA}{need:.1f} GiB{C.RESET}"
                  + (f", available {have:.1f} GiB" if have else ""))
        if have and need > 0.8 * have and not quiet:
            print(f"{C.WARN}{Icon.WARN}This search is estimated to need {need:.1f} GiB "
                  f"against {have:.1f} GiB available. Raise downsample_factor "
                  f"(memory scales as its inverse square) or crop the DEM. The estimate "
                  f"is rough, so this is a warning rather than a refusal.{C.RESET}")

    cap = max_memory_gb
    if cap is None:
        # 80% of what is available, so the cap bites before the kernel does while
        # leaving room for the address space numba and BLAS reserve without touching.
        cap = 0.8 * have if have else None
    capped = bool(cap and apply_memory_cap(cap))
    if capped and not quiet:
        print(f"   {Icon.GEAR}Address space capped at {C.MAGENTA}{cap:.1f} GiB{C.RESET}"
              f" (max_memory_gb=0 disables)")
    return {"estimate_gb": need, "available_gb": have,
            "cap_gb": cap if capped else None, "capped": capped}


def validate_parameters(params):
    """
    Pre-flight validation checks to enforce 'Fail Fast' mechanisms. 
    Verifies the existence of critical files and the physical logic of search bounds 
    before engaging the memory-heavy processing loops.

    Parameters
    ----------
    params : dict
        Fully resolved parameters, after the config, fallback and command line have
        been reconciled.

    Raises
    ------
    SystemExit
        If any check fails. Every problem is collected and reported at once rather
        than one per run, since the expensive stages come afterwards.
    """
    errors = []
    
    # 1. Check DEM path existence
    if not os.path.exists(params['dem_path']):
        errors.append(f"DEM file not found: {params['dem_path']}")
    
    # 2. Check physical layout impossibilities
    # 0 is legitimate and means "do not prune": an experiment whose array is a long
    # thin strip rather than a compact blob -- TAMBO along a canyon wall -- is exactly
    # what the opening would erase, so it has to be possible to turn off.
    if params['min_width_km'] < 0:
        errors.append("min_width_km cannot be negative (0 disables tendril pruning).")
        
    if params['target_antennas'] <= 0:
        errors.append("target_antennas must be strictly positive (> 0).")
        
    if params['antenna_spacing_km'] <= 0:
        errors.append("antenna_spacing_km must be strictly positive (> 0).")
        
    # 3. Verify Road Map existence if specified
    if params['road_map_path'] is not None and not os.path.exists(params['road_map_path']):
        errors.append(f"Road map file not found: {params['road_map_path']}")
        
    # 4. Verify Ray-Tracing bounds logic
    if params['min_dist_km'] >= params['max_dist_km']:
        errors.append("min_dist_km must be strictly less than max_dist_km.")
        
    # 5. Verify Altitude bounds logic
    if params['min_altitude'] is not None and params['max_altitude'] is not None:
        if params['min_altitude'] >= params['max_altitude']:
            errors.append("min_altitude must be strictly less than max_altitude.")
            
    # 6. Verify Slope bounds logic
    if params['min_slope_deg'] < 0:
        errors.append("min_slope_deg cannot be negative.")
    if params['min_slope_deg'] >= params['max_slope_deg']:
        errors.append("min_slope_deg must be strictly less than max_slope_deg.")

    # 7. Verify memory constraints
    if params['tile_size'] <= 0:
        errors.append("tile_size must be strictly positive (> 0).")

    if params.get('candidate_stride', 5) <= 0:
        errors.append("candidate_stride must be strictly positive (> 0).")

    if params.get('cell_size_deg') is not None and params['cell_size_deg'] <= 0:
        errors.append("cell_size_deg must be strictly positive (> 0) when specified.")

    # 8. Verify Resume Directory
    if params.get('resume'):
        res_dir = params.get('resume_dir')
        if not res_dir or not os.path.exists(res_dir):
            errors.append(f"Resume directory not found: {res_dir}")
        elif not os.path.exists(os.path.join(res_dir, 'buffer_A.npy')):
            errors.append(f"Cannot resume: 'buffer_A.npy' not found inside {res_dir}")

    # 9. Verify CPU Cores
    sys_cores = multiprocessing.cpu_count()
    if params.get('num_cores', -1) != -1:
        if params['num_cores'] <= 0:
            errors.append(f"num_cores must be a positive integer or -1. Received: {params['num_cores']}")
        elif params['num_cores'] > sys_cores:
            errors.append(f"num_cores requested ({params['num_cores']}) exceeds available system cores ({sys_cores}).")

    # Execute Fail-Fast
    if errors:
        print(f"\n{C.FAIL}================================================================================{C.RESET}")
        print(f"{C.FAIL}{C.BOLD}PRE-FLIGHT VALIDATION FAILED:{C.RESET}")
        for error in errors:
            print(f"  {C.FAIL}{Icon.CROSS}{error}{C.RESET}")
        print(f"{C.WARN}Please correct the parameters in your config file or CLI and try again.{C.RESET}")
        print(f"{C.FAIL}================================================================================{C.RESET}\n")
        sys.exit(1)


# ==========================================
#             CONFIGURATION FILES
# ==========================================
# Config handling used to live inside main(), which meant a library user could load
# neither a config file nor a template without reimplementing both. They are ordinary
# functions on the module now, and main() calls them like anyone else.

CONFIG_PRESETS = ("default", "lima", "arequipa")


def default_config(preset="default"):
    """
    The full set of knobs with their default values, as a plain dictionary.

    Every key the pipeline understands appears here, which is the point: a template
    with holes in it silently falls back for whatever it omits, and the fallbacks are
    the least visible input the tool has.

    Parameters
    ----------
    preset : str, optional
        ``"default"``, or ``"lima"``/``"arequipa"`` to fill in that region's origin,
        RFI zones, name and DEM filename.

    Returns
    -------
    dict
        Configuration, ready to serialise or to pass to the pipeline.

    Raises
    ------
    ValueError
        If the preset is not one of :data:`CONFIG_PRESETS`.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> cfg = ss.default_config("arequipa")
    >>> cfg["rfi_zones"], cfg["min_slope_deg"], cfg["explain"]
    ('arequipa', 3.0, True)
    """
    if preset not in CONFIG_PRESETS:
        raise ValueError(f"unknown preset {preset!r}; expected one of {CONFIG_PRESETS}")

    config = {
        "dem_path": "path_to_your_dem.tif",
        # Deliberately null rather than 0.0. Zero is a *valid* coordinate -- it is in
        # the Gulf of Guinea -- so a placeholder someone forgets to edit produces a
        # run georeferenced to the wrong continent rather than an error. Null means
        # "read the DEM's own ModelTiepointTag", which is the recommended use.
        "origin_lat": None,
        "origin_lon": None,
        "target_antennas": 10000,
        "min_width_km": 2.0,
        "min_altitude": None,
        "max_altitude": None,
        "min_slope_deg": 3.0,
        "max_slope_deg": 25.0,
        "antenna_spacing_km": 1.0,
        "min_dist_km": 10.0,
        "max_dist_km": 80.0,
        "grid_type": "hex",
        "downsample_factor": 4,
        "cell_size_deg": None,
        "candidate_stride": 5,
        "slope_baseline_m": None,
        "energy_min_pev": None,
        "energy_max_pev": None,
        "n_azimuths": 9,
        "azimuth_half_width_deg": 60.0,
        "elev_min_deg": -3.0,
        "elev_max_deg": 3.0,
        "n_elev_bins": 12,
        "min_column_depth_gcm2": 0.0,
        # The CLI's --require_sky, in its config-file spelling. The pipeline takes the
        # positive form (require_terrain); main() inverts it.
        "require_sky": False,
        "decay_energy_pev": None,
        "max_range_km": None,
        # What weights the spectrum-folded decay probability: "flux" (the default, and
        # what every published number used), "acceptance", or "flux_times_acceptance".
        # The latter two need decay_response_csv, a two-column A(E) table.
        "decay_weight_by": "flux",
        "decay_response_csv": None,
        "score_percentile": None,
        "stop_at_target": False,
        "max_memory_gb": None,
        "decay_energy_min_pev": None,
        "decay_energy_max_pev": None,
        "decay_spectral_index": None,
        "shower_development_m": 3000.0,
        "gap_close_km": None,
        "min_target_slope_deg": None,
        "max_target_slope_deg": None,
        "fresnel_frequency_mhz": None,
        "antenna_height_m": 2.0,
        "exclude_near_field": True,
        "fresnel_near_field_m": 500.0,
        "refraction_k": None,
        "depth_band_gcm2": None,
        "geomag_declination_deg": None,
        "geomag_inclination_deg": None,
        "muon_shielding_km": None,
        "bilinear_sampling": True,
        "use_geomagnetic": True,
        "grammage_mode": "radio",
        "grammage_band_gcm2": None,
        "grammage_maturity_gcm2": None,
        "grammage_band_fraction": None,
        "shower_elongation_rate_gcm2": None,
        "shower_lambda_gcm2": None,
        "solid_angle_half_sr": None,
        "distance_band_m": None,
        "clearance_full_at": None,
        "score_weights": None,
        "nu_interaction_length_gcm2": None,
        "score_composition": "product",
        "min_score": 0.0,
        "tile_size": 2048,
        "num_cores": -1,
        "rfi_zones": "none",
        "road_map_path": None,
        "max_road_dist_km": 20.0,
        "search_mode": "distributed",
        "min_sub_array_size": 500,
        "min_aspect_deg": None,
        "max_aspect_deg": None,
        "region_name": "Custom Region",
        "generate_kml": False,
        "print_info": True,
        "explain": True,
        "output_directory_base_with_given_json": "../output/",
        "output_image_format": "png",
        "resume": False,
        "resume_dir": None
    }

    if preset == 'lima':
        config['origin_lat'] = ORIGIN_LAT_LIMA
        config['origin_lon'] = ORIGIN_LON_LIMA
        config['rfi_zones'] = 'lima'
        config['region_name'] = 'Lima, Peru'
        config['dem_path'] = 'lima_AW3D30.tif'
    elif preset == 'arequipa':
        config['origin_lat'] = ORIGIN_LAT_AREQUIPA
        config['origin_lon'] = ORIGIN_LON_AREQUIPA
        config['rfi_zones'] = 'arequipa'
        config['region_name'] = 'Arequipa, Peru'
        config['dem_path'] = 'arequipa_SRTMGL1.tif'
    return config


def generate_config(path, preset="default"):
    """
    Writes a configuration template to ``path``, creating its directory if needed.

    What ``--generate_config`` does, available to anyone driving the pipeline in a
    loop or generating a family of runs.

    Parameters
    ----------
    path : str
        Destination JSON file.
    preset : str, optional
        As :func:`default_config`.

    Returns
    -------
    dict
        The configuration that was written.

    Notes
    -----
    A generated config names every key, so it also overrides every fallback. That is
    intended -- but note the command line still wins over both, which it did not
    always do.
    """
    config = default_config(preset)
    dir_name = os.path.dirname(os.path.abspath(path))
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)
    with open(path, 'w') as f:
        json.dump(config, f, indent=4)
    return config


# Configuration keys naming a file or directory the run needs to find. Resolved
# relative to the configuration file, which is the only fixed point a config can
# reason about; the working directory is not.
_PATH_KEYS = ("dem_path", "road_map_path", "resume_dir")


def resolve_config_paths(config, config_dir, quiet=False):
    """
    Makes a configuration's relative paths absolute, against the config's own directory.

    A configuration that says ``"dem_path": "../input/dem/colca.tif"`` is describing
    where the DEM sits relative to *itself*, which is the only thing it can know. The
    pipeline resolved it against the working directory instead, so the bundled configs
    ran only from ``src/`` and produced a `FileNotFoundError` anywhere else -- the
    long-standing "you must ``cd src`` first" wart.

    A path that does not resolve against the configuration's directory but does resolve
    against the working directory is left alone, with a warning: that is the old
    behaviour, and silently breaking someone's working setup to fix a wart is a poor
    trade. Absolute paths are untouched.

    Parameters
    ----------
    config : dict
        A configuration mapping. Not modified; a copy is returned.
    config_dir : str
        Directory holding the configuration file.
    quiet : bool, optional
        Suppress the warning about a working-directory-relative fallback.

    Returns
    -------
    dict
        A copy with the path keys resolved.

    Examples
    --------
    The repository's own layout is why this is safe to change: ``config/`` and ``src/``
    are both one level below the root, so ``../input/dem/colca.tif`` names the same file
    read either way, and no shipped configuration has to change.

    >>> import os
    >>> from oroscope import site_searcher as ss
    >>> a = os.path.normpath(os.path.join("/repo/config", "../input/dem/colca.tif"))
    >>> b = os.path.normpath(os.path.join("/repo/src", "../input/dem/colca.tif"))
    >>> a == b == "/repo/input/dem/colca.tif"
    True

    An absolute path is left as it is:

    >>> cfg = ss.resolve_config_paths({"dem_path": "/data/x.tif"}, "/repo/config")
    >>> cfg["dem_path"]
    '/data/x.tif'
    """
    resolved = dict(config)
    for key in _PATH_KEYS:
        value = resolved.get(key)
        if not value or not isinstance(value, str) or os.path.isabs(value):
            continue
        candidate = os.path.normpath(os.path.join(config_dir, value))
        if os.path.exists(candidate):
            resolved[key] = candidate
            continue
        if os.path.exists(value):
            if not quiet:
                print(f"{C.WARN}{Icon.WARN}{key}={value!r} was found relative to the "
                      f"working directory but not to the configuration file. That is "
                      f"the old behaviour and still works; making it relative to the "
                      f"config (or absolute) will keep working from anywhere.{C.RESET}")
            continue
        # Neither exists. Resolve against the config anyway, so the error names the
        # path the configuration actually asked for.
        resolved[key] = candidate
    return resolved


def load_config(path):
    """
    Reads a configuration JSON, resolving its relative paths against its own directory.

    Parameters
    ----------
    path : str
        Path to the file.

    Returns
    -------
    dict
        Its contents, or an empty dictionary when the file does not exist -- which is
        how the command line has always treated a missing ``--config_path``, and
        matching it here keeps one behaviour rather than two.

        Path-valued keys (``dem_path``, ``road_map_path``, ``resume_dir``) come back
        absolute, resolved against the directory holding the configuration rather than
        the working directory. See :func:`resolve_config_paths`.

    Raises
    ------
    json.JSONDecodeError
        If the file exists but is not valid JSON. Unlike a missing file, that is a
        mistake worth failing on.
    """
    if not path or not os.path.exists(path):
        return {}
    with open(path, 'r') as f:
        config = json.load(f)
    return resolve_config_paths(config, os.path.dirname(os.path.abspath(path)))


def emit_explanation(results, run_output_dir=None, print_it=True):
    """
    Composes the run's plain-language summary, prints it, and saves it beside the run.

    A thin wrapper over :func:`explain.explain_results`: the words themselves are that
    function's business, so they can be regenerated from an old results file without
    this one. What is added here is the placement -- last on the console, so it is what
    a reader is left with, and in ``explanation.txt`` so the run can be handed on
    without the terminal it was run in.

    Failures are reported and swallowed. A summary that cannot be written is not a
    reason to lose a search that already succeeded.

    Parameters
    ----------
    results : dict
        The run's results. Gains an ``"explanation"`` key.
    run_output_dir : str, optional
        Directory to write ``explanation.txt`` into. Omitted writes no file.
    print_it : bool, optional
        Whether to print. False still composes and stores the text.

    Returns
    -------
    str or None
        The summary, or ``None`` if it could not be composed.
    """
    try:
        text = explain_mod.explain_results(results, results.get("provenance"))
    except Exception as e:                               # pragma: no cover - defensive
        print(f"   {C.WARN}{Icon.WARN}Could not compose the run summary: {e}{C.RESET}")
        return None

    results["explanation"] = text
    if print_it:
        print("\n" + text)
    if run_output_dir:
        try:
            path = os.path.join(run_output_dir, "explanation.txt")
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
            if print_it:
                print(f"   {Icon.INFO}{os.path.abspath(path)}")
        except OSError as e:                             # pragma: no cover - defensive
            print(f"   {C.WARN}{Icon.WARN}Could not write explanation.txt: {e}{C.RESET}")
    return text


# ==========================================
#        CONFIGURATION -> PIPELINE
# ==========================================
# Keys a configuration file carries that are not pipeline parameters: they steer the
# command line, not the search.
_NOT_PIPELINE_KEYS = frozenset({
    "print_info", "output_directory_base_with_given_json", "require_sky",
})


def resolve_rfi_zones(value, quiet=False):
    """
    Turns the many spellings of ``rfi_zones`` into the list the pipeline wants.

    A configuration may carry a preset name, a JSON string from the command line, an
    explicit list of zones, or nothing. The pipeline accepts only the list, and it
    iterates whatever it is given -- so a preset name reaching it unresolved is iterated
    *character by character*, each character failing the ``item[0] == 'circle'`` test.
    That is silent: no exception, no warning, and a search that believes it is excluding
    five radio-noise zones runs with none.

    Parameters
    ----------
    value : str, list, or None
        ``'lima'`` or ``'arequipa'`` for a preset, ``'none'`` or ``None`` for no zones,
        a JSON string, or an explicit list of zone tuples.
    quiet : bool, optional
        Suppress the warning printed when a string cannot be parsed.

    Returns
    -------
    list or None
        Zones for the pipeline, or ``None``.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> len(ss.resolve_rfi_zones("arequipa"))
    5
    >>> ss.resolve_rfi_zones("none") is None
    True
    >>> ss.resolve_rfi_zones(None) is None
    True

    An explicit list is passed through untouched:

    >>> ss.resolve_rfi_zones([("circle", -16.4, -71.5, 25.0, "Arequipa")])
    [('circle', -16.4, -71.5, 25.0, 'Arequipa')]
    """
    if value is None:
        return None
    if isinstance(value, str):
        key = value.lower()
        if key == "lima":
            return LIMA_RFI_ZONES
        if key == "arequipa":
            return AREQUIPA_RFI_ZONES
        if key == "none":
            return None
        try:
            return json.loads(value)
        except Exception as e:                           # pragma: no cover - defensive
            if not quiet:
                print(f"{C.WARN}{Icon.WARN}Could not parse rfi_zones {value!r}; "
                      f"proceeding with none. {e}{C.RESET}")
            return None
    return list(value) if value else None


def config_to_pipeline_kwargs(config, quiet=False, **overrides):
    """
    Translates a configuration mapping into :func:`find_grand_regions_interactive` kwargs.

    This translation was written out three times -- in ``main()``, in the child process
    ``sensitivity`` spawns, and in ``tools/run_arequipa_full.py`` -- and the copies had
    already drifted. The sweep child passed ``rfi_zones`` through as the raw preset
    name, which the pipeline then iterated character by character and silently resolved
    to no zones at all; it never inverted ``require_sky`` either. Having one function
    do it means a new parameter is added once, and means the drift cannot recur.

    What the translation actually consists of:

    - comment keys (``_``-prefixed) and command-line-only keys are dropped;
    - ``require_sky`` is inverted into the pipeline's positive ``require_terrain``;
    - ``rfi_zones`` presets are resolved to lists (:func:`resolve_rfi_zones`);
    - the bands are made tuples, since JSON gives lists;
    - ``score_weights`` is parsed from its ``"name=value"`` spelling;
    - ``decay_spectral_index`` accepts a scalar or a ``(low, high)`` pair;
    - a negative ``azimuth_half_width_deg`` means "unbounded", i.e. ``None``;
    - anything the pipeline does not accept is dropped, with a warning naming it, so a
      misspelled key is reported rather than silently ignored.

    Parameters
    ----------
    config : dict
        A configuration mapping, as :func:`load_config` returns.
    quiet : bool, optional
        Suppress the warnings about unknown keys and unparseable zones.
    **overrides
        Applied after the translation, so a caller can set ``run_output_dir`` or force
        ``generate_kml=False`` without having to edit the mapping.

    Returns
    -------
    dict
        Keyword arguments for :func:`find_grand_regions_interactive`.

    Examples
    --------
    >>> from oroscope import site_searcher as ss
    >>> kw = ss.config_to_pipeline_kwargs(
    ...     {"_comment": "ignored", "dem_path": "d.tif", "require_sky": True,
    ...      "rfi_zones": "arequipa", "grammage_band_gcm2": [236.0, 1287.0],
    ...      "print_info": True})
    >>> kw["require_terrain"], kw["grammage_band_gcm2"], len(kw["rfi_zones"])
    (False, (236.0, 1287.0), 5)

    The command-line-only keys do not reach the pipeline:

    >>> "print_info" in kw or "require_sky" in kw
    False
    """
    params = {k: v for k, v in config.items()
              if not k.startswith("_") and k not in _NOT_PIPELINE_KEYS}

    # The pipeline asks the positive question; the config and the CLI ask the negative.
    params["require_terrain"] = not config.get("require_sky", False)

    params["rfi_zones"] = resolve_rfi_zones(config.get("rfi_zones"), quiet=quiet)

    # JSON has no tuples, and these are compared and unpacked as pairs downstream.
    for key in ("depth_band_gcm2", "grammage_band_gcm2", "distance_band_m"):
        if params.get(key) is not None:
            params[key] = tuple(params[key])

    if params.get("score_weights") is not None:
        params["score_weights"] = parse_score_weights(params["score_weights"])

    if params.get("decay_spectral_index") is not None:
        params["decay_spectral_index"] = _one_or_pair(params["decay_spectral_index"])

    # Negative means unbounded: scan every azimuth rather than a wedge around aspect.
    half_width = params.get("azimuth_half_width_deg")
    if half_width is not None and half_width < 0:
        params["azimuth_half_width_deg"] = None

    params.update(overrides)

    unknown = sorted(k for k in params if k not in _PIPELINE_PARAMS)
    for key in unknown:
        del params[key]
    if unknown and not quiet:
        print(f"{C.WARN}{Icon.WARN}Ignoring {len(unknown)} key(s) the pipeline does not "
              f"take: {', '.join(unknown)}. Check the spelling against "
              f"default_config().{C.RESET}")
    return params


def run_from_config(config, run_output_dir=".", quiet=False, **overrides):
    """
    Runs one search from a configuration mapping or a path to one.

    The library counterpart of ``oroscope --config_path ...``: one call, from a file a
    run can be reproduced from, returning the results dictionary.

    Parameters
    ----------
    config : dict or str
        A configuration mapping, or a path to a JSON configuration file.
    run_output_dir : str, optional
        Where the run writes its outputs.
    quiet : bool, optional
        Passed to :func:`config_to_pipeline_kwargs`.
    **overrides
        Applied after the translation.

    Returns
    -------
    dict
        The results dictionary, as :func:`find_grand_regions_interactive` returns.

    Examples
    --------
    Translation and execution are separable, which is what makes the mapping testable
    without running a search:

    >>> from oroscope import site_searcher as ss
    >>> kw = ss.config_to_pipeline_kwargs(ss.default_config(), quiet=True)
    >>> kw["require_terrain"], kw["rfi_zones"]
    (True, None)
    """
    if isinstance(config, str):
        config = load_config(config)
    kwargs = config_to_pipeline_kwargs(config, quiet=quiet, **overrides)
    return find_grand_regions_interactive(run_output_dir=run_output_dir, **kwargs)


# ==========================================
#             MAIN EXECUTION ORCHESTRATOR
# ==========================================
def find_grand_regions_interactive(dem_path, cell_size_deg=None, target_antennas=10000,
                            rfi_zones=None, origin_lat=None, origin_lon=None,
                            min_width_km=2.0, min_altitude=None, max_altitude=None,
                            antenna_spacing_km=1.0, min_dist_km=10.0, max_dist_km=80.0,
                            road_map_path=None, max_road_dist_km=20.0,
                            grid_type='hex', generate_kml=False,
                            search_mode='distributed', min_sub_array_size=500,
                            min_aspect_deg=None, max_aspect_deg=None,
                            min_slope_deg=3.0, max_slope_deg=25.0,
                            region_name=None,
                            downsample_factor=4, run_output_dir=".", 
                            output_image_format='png', tile_size=2048,
                            resume=False, resume_dir=None, num_cores=-1,
                            candidate_stride=5, slope_baseline_m=None,
                            energy_min_pev=None, energy_max_pev=None,
                            n_azimuths=9, azimuth_half_width_deg=60.0,
                            elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
                            min_column_depth_gcm2=0.0, require_terrain=True,
                            min_target_slope_deg=None, max_target_slope_deg=None,
                            max_range_km=None, score_percentile=None,
                            decay_weight_by='flux', decay_response_csv=None,
                            stop_at_target=False,
                            decay_energy_pev=None,
                            decay_energy_min_pev=None, decay_energy_max_pev=None,
                            decay_spectral_index=None,
                            shower_development_m=3000.0,
                            gap_close_km=None,
                            fresnel_frequency_mhz=None, refraction_k=None,
                            antenna_height_m=2.0, fresnel_near_field_m=500.0,
                            exclude_near_field=True,
                            depth_band_gcm2=None, score_composition='product',
                            score_weights=None, distance_band_m=None,
                            solid_angle_half_sr=None, clearance_full_at=None,
                            min_score=0.0,
                            geomag_declination_deg=None, geomag_inclination_deg=None,
                            use_geomagnetic=True, grammage_mode='radio',
                            grammage_band_gcm2=None, grammage_maturity_gcm2=None,
                            grammage_band_fraction=None,
                            shower_elongation_rate_gcm2=None, shower_lambda_gcm2=None,
                            muon_shielding_km=None, bilinear_sampling=True,
                            nu_interaction_length_gcm2=None,
                            max_memory_gb=None, explain=True):
    """
    The main orchestrator. Now decoupled from logic, it sets up the environment,
    calls the pipeline helpers in sequence, and manages memory cleanup and checkpointing.

    The map resolution (cell_size_deg) is read from the DEM's own georeferencing tags
    unless the caller overrides it; every metric conversion downstream derives from it.

    Parameters
    ----------
    dem_path : str
        Path to the input elevation GeoTIFF.

    cell_size_deg : float, optional
        Pixel size in degrees, overriding the GeoTIFF's own tag.

    target_antennas : int, optional
        Capacity wanted from a single site.

    rfi_zones : list or str, optional
        Exclusion zones, or the name of a bundled set.

    origin_lat, origin_lon : float, optional
        Coordinates of the DEM's north-west corner, in degrees. Read from the file's
        own ``ModelTiepointTag`` when omitted, which is the recommended use; a supplied
        value that disagrees with the tag by more than ~100 m is reported rather than
        silently honoured.

    min_width_km : float, optional
        Narrowest feature to keep, in km. 0 disables pruning, which is what a
        strip-shaped array along a canyon wall needs.

    min_altitude, max_altitude : float, optional
        Altitude bounds, in metres.

    antenna_spacing_km : float, optional
        Detector spacing, in km.

    min_dist_km, max_dist_km : float, optional
        Accepted range to the first terrain intersection, in km -- the decay-baseline
        window. Also sets how far the profile is walked.

    road_map_path : str, optional
        Aligned GeoTIFF of distance-to-road values.

    max_road_dist_km : float, optional
        Maximum allowed distance from a road, in km.

    grid_type : {'square', 'hex'}, optional
        Lattice the detectors are placed on.

    generate_kml : bool, optional
        Also write a Google Earth ``.kml``.

    search_mode : {'single', 'distributed'}, optional
        Whether one site must hold the whole array.

    min_sub_array_size : int, optional
        Capacity a sub-array must reach in distributed mode.

    min_aspect_deg, max_aspect_deg : float, optional
        Required facing directions, in degrees clockwise from north.

    min_slope_deg, max_slope_deg : float, optional
        Slope band the detector site must fall in, in degrees. This is the *near* wall,
        the ground the array stands on.

    region_name : str, optional
        Human-readable label for the outputs.

    downsample_factor : int, optional
        Factor at which sites are labelled and areas measured. Above 1, a feature only
        a few pixels wide loses area it keeps detectors on.

    run_output_dir : str, optional
        Directory to write into.

    output_image_format : str, optional
        Extension for the overview map, such as ``png`` or ``pdf``.

    tile_size : int, optional
        Side of the square tile held in RAM at once.

    resume : bool, optional
        Reuse a previous run's scan buffer instead of recomputing it.

    resume_dir : str, optional
        Directory holding that buffer.

    num_cores : int, optional
        Threads for the scan. ``-1`` uses all of them.

    candidate_stride : int, optional
        Keeps every Nth screened pixel. Measured to be unbiased; see
        :func:`get_candidates_chunked`.

    slope_baseline_m : float, optional
        Ground distance over which slope is measured, in metres. Slope is
        scale-dependent, so this is an explicit choice rather than an accident of the
        DEM's resolution.

    energy_min_pev, energy_max_pev : float, optional
        Tau energy range. When given it *overrides* the distance window with one
        derived from the decay length, and in particle mode also sets the shower band.

    n_azimuths : int, optional
        Azimuths scanned per candidate. This is what sets the cost of a run.

    azimuth_half_width_deg : float, optional
        Half-width of the fan about each candidate's aspect, in degrees.

    elev_min_deg, elev_max_deg : float, optional
        Edges of the accepted arrival window, in degrees.

    n_elev_bins : int, optional
        Bins across that window. Nearly free: one walk serves them all.

    min_column_depth_gcm2 : float, optional
        Column depth a direction must have to count, in g/cm^2.

    require_terrain : bool, optional
        ``True`` selects directions striking rock; ``False`` selects directions
        escaping to the sky, which is the cosmic-ray channel.

    min_target_slope_deg, max_target_slope_deg : float, optional
        Bounds on the struck terrain's slope along the arrival azimuth. This is the
        *far* wall. Unset by default, which asks only that rock is present -- true
        almost everywhere in mountainous terrain.

    max_range_km : float, optional
        How far each profile is walked, in km. Defaults to ``max_dist_km``. Worth
        setting larger for a short-range search: column depth accumulates over the
        whole walk, so tying the two makes the reported depth a property of where the
        walk stopped rather than of the target's thickness.
    score_percentile : float, optional
        Keep this percentage of viable candidates, ranked by score, instead of cutting
        at an absolute ``min_score``. Preferred, and for the same reason: a rank is
        scale-free, so it does not move when the composition or the number of
        components changes.

    decay_weight_by : {'flux', 'acceptance', 'flux_times_acceptance'}, optional
        What weights the spectrum-folded decay probability. ``'flux'`` asks what
        fraction of arriving neutrinos decays usefully and is the default;
        ``'acceptance'`` asks the same over the energies the detector responds to, with
        no assumed spectrum; ``'flux_times_acceptance'`` is the event-rate integrand.
        The latter two require ``decay_response_csv``.
    decay_response_csv : str or callable, optional
        Detection response ``A(E)`` for the acceptance weightings: a path to a
        two-column CSV of energy in PeV against relative response, or a callable.
        :func:`aperture.infer_response` recovers one from a published integral curve.

    stop_at_target : bool, optional
        In distributed mode, stop selecting sites once ``target_antennas`` is reached.
        Sites are ranked by capacity, so this reports the best sites for the array
        actually wanted rather than every patch of qualifying ground.

    decay_energy_pev : float, optional
        Single energy at which to score the probability the tau decays in the gap.
        Superseded by the range below and kept for asking what one energy would have
        said: measured, the answer ran from 10878 detector positions at 3 PeV to zero
        at 100, so one number does not stand in for a spectrum.
    decay_energy_min_pev, decay_energy_max_pev : float, optional
        Tau energy range over which to fold the decay probability against the flux.
        The defensible form, and what makes the result a property of the terrain rather
        than of the energy someone picked.
    decay_spectral_index : float, optional
        Gamma in dN/dE ~ E^-gamma for that folding (default 2.0). A softer spectrum
        weights low energies, where the tau decays readily, so it drives the term
        toward 1 -- an assumption deserving the same scrutiny as any other.

    shower_development_m : float, optional
        Path the shower needs after the tau decays, in metres.

    gap_close_km : float, optional
        Size of the morphological closing element, in km. Defaults to
        ``antenna_spacing_km``. Closing more than doubles reported area on real
        terrain, so it is worth setting deliberately; 0 disables it.

    fresnel_frequency_mhz : float, optional
        Radio band for the Fresnel clearance measurement. ``None`` skips that pass, as
        a particle experiment wants.

    refraction_k : float, optional
        Refraction factor for the radio path only. The particle geometry always uses
        the true Earth radius.

    antenna_height_m : float, optional
        Receiver height above ground, in metres.

    fresnel_near_field_m : float, optional
        Stretch of path the clearance measurement skips, in metres.

    exclude_near_field : bool, optional
        Apply that near-field cut-off.

    depth_band_gcm2 : tuple of float, optional
        Column-depth band scoring 1, in g/cm^2.

    score_composition : {'product', 'mean', 'min'}, optional
        How components combine.

    score_weights : dict, optional
        Per-component weights.

    distance_band_m : tuple of float, optional
        Exit-distance band scoring 1, in metres. Defaults to the decay window.

    solid_angle_half_sr : float, optional
        Accepted solid angle scoring 0.5, in steradians. The 0.05 default is
        GRAND-scale and saturates against a canyon's much larger acceptance.

    clearance_full_at : float, optional
        Fresnel clearance ratio scoring 1.

    min_score : float, optional
        Score a candidate must reach. Note a product composition concentrates near
        zero, so any threshold in the middle sits on a cliff.

    geomag_declination_deg, geomag_inclination_deg : float, optional
        Field direction, in degrees. Supply the IGRF values for the site.

    use_geomagnetic : bool, optional
        Weight directions by ``|v x B|``. Radio only; particles do not care.

    grammage_mode : {'radio', 'particle'}, optional
        Whether atmospheric depth is scored as a maturity threshold or as a band.

    grammage_band_gcm2 : tuple of float, optional
        Explicit shower band in particle mode. Setting it disables
        ``grammage_band_fraction``.

    grammage_maturity_gcm2 : float, optional
        Depth at which the radio maturity ramp reaches 1, in g/cm^2.

    grammage_band_fraction : float, optional
        Fraction of peak particle content still counted as a usable shower, when the
        band is derived from an energy range.

    shower_elongation_rate_gcm2, shower_lambda_gcm2 : float, optional
        Shower-profile parameters: how much deeper maximum sits per decade of energy,
        and the Gaisser-Hillas interaction length.

    muon_shielding_km : float, optional
        Rock overburden required for muon rejection, in km.

    bilinear_sampling : bool, optional
        Interpolate the terrain profile between pixel centres. Costs about 1.44x and
        removes an asymmetric half-pixel bias.

    nu_interaction_length_gcm2 : float, optional
        Neutrino interaction length, enabling the Earth-chord attenuation term.

    max_memory_gb : float, optional
        Ceiling on this process's address space, in GiB. ``None`` uses 80% of what the
        system reports available, so a search that outgrows the machine fails with
        ``MemoryError`` instead of inviting the OOM killer to choose a victim; 0
        disables the cap. See :func:`preflight_memory`.

    explain : bool, optional
        Print the plain-language account of the run -- what was found, which
        constraint set the size of the answer, and which numbers are assumptions --
        and save it as ``explanation.txt`` beside the results. On by default. The
        text is also in the returned dictionary under ``"explanation"``, and can be
        regenerated from any results file with :func:`explain.explain_results`.

    Returns
    -------
    dict
        The run's results: ``parameters``, ``results`` (sites, capacity), ``funnel``,
        ``regions``, ``timings_sec``, ``aperture``, ``provenance``, ``explanation``
        and the paths of the files written. The same content as the results JSON, so
        a caller no longer has to find and re-read the file this just wrote. A run
        that finds no candidate at all still returns its funnel, which is the case
        where the funnel matters most.

    Notes
    -----
    Writes GeoTIFF, world file, PNG, optional KML, a results JSON and a provenance
    record into ``run_output_dir``, and prints a selection funnel. When a search
    returns nothing, the funnel is the first place to look: the constraint responsible
    is the line where the survivor count collapses.
    """
    # The DEM knows its own corner. Prefer it, and say so loudly when a supplied one
    # disagrees -- a wrong origin does not fail, it mis-georeferences everything.
    origin_lat, origin_lon, origin_source = resolve_origin(dem_path, origin_lat, origin_lon)
    if origin_lat is None:
        raise ValueError(
            "origin_lat/origin_lon were not given and the DEM carries no "
            "ModelTiepointTag to read them from")
    print(f"      Origin: {C.MAGENTA}{origin_lat:.6f}, {origin_lon:.6f}{C.RESET}"
          f"  ({origin_source})")

    # The run's own directory, which main() created and the library did not -- so a
    # caller who passed a path that did not exist got a FileNotFoundError from inside
    # numpy's open_memmap, naming a scratch buffer rather than the directory.
    os.makedirs(run_output_dir, exist_ok=True)

    # Estimate, warn and cap before anything expensive is allocated. In the pipeline
    # rather than in main() so that every caller is protected, not only the CLI one.
    preflight_memory(dem_path, downsample_factor=downsample_factor,
                     candidate_stride=candidate_stride, max_memory_gb=max_memory_gb)

    # Cast safety to ensure slice logic doesn't fail if passed as float via JSON
    downsample_factor = int(downsample_factor)
    tile_size = int(tile_size)
    candidate_stride = int(candidate_stride)

    # Establish the sampling geometry before anything else depends on it
    map_grid = resolve_grid_geometry(dem_path, origin_lat, cell_size_deg)
    cell_size_deg = map_grid.cell_size_deg
    cell_size_y, cell_size_x = map_grid.cell_size_y, map_grid.cell_size_x

    # An energy range, when given, sets the decay-baseline window (roadmap 4.7)
    energy_note = ""
    if energy_min_pev and energy_max_pev:
        lo_m, hi_m = arrival_scan.distance_window_from_energy(energy_min_pev, energy_max_pev)
        min_dist_km, max_dist_km = lo_m / 1000.0, hi_m / 1000.0
        energy_note = f" (from {energy_min_pev:g}-{energy_max_pev:g} PeV)"

        # The same energies also fix a particle array's shower band, but through the
        # shower profile rather than the tau decay length. Only derived when no band
        # was given explicitly, and only in particle mode -- for radio the criterion
        # is a maturity threshold, not a band.
        if grammage_mode == 'particle' and grammage_band_gcm2 is None:
            band_kw = {}
            if grammage_band_fraction is not None:
                band_kw["fraction"] = grammage_band_fraction
            if shower_lambda_gcm2 is not None:
                band_kw["lambda_gcm2"] = shower_lambda_gcm2
            if shower_elongation_rate_gcm2 is not None:
                band_kw["elongation_rate"] = shower_elongation_rate_gcm2
            grammage_band_gcm2 = physics.grammage_band_from_energy(
                energy_min_pev, energy_max_pev, **band_kw)
            frac = (grammage_band_fraction
                    if grammage_band_fraction is not None else 0.1)
            print(f"      Shower band: {C.MAGENTA}{grammage_band_gcm2[0]:.0f}"
                  f"-{grammage_band_gcm2[1]:.0f} g/cm²{C.RESET}"
                  f" (from {energy_min_pev:g}-{energy_max_pev:g} PeV"
                  f" at {frac:g} of shower maximum)")

    observables = None

    # Field for this site: inclination follows from its own coordinates, so moving
    # the search elsewhere gets that right automatically (roadmap 4.12b)
    # A(E) for the acceptance weightings. A config carries a path rather than a
    # callable, since JSON cannot hold a function; a library caller may pass either.
    _decay_response = decay_response_csv
    if isinstance(_decay_response, str):
        _decay_response = aperture_mod.TabulatedResponse.from_csv(_decay_response)
    if decay_weight_by != "flux" and _decay_response is None:
        raise ValueError(
            f"decay_weight_by={decay_weight_by!r} needs decay_response_csv: a two-column "
            f"CSV of energy in PeV against relative response. See data/ for the "
            f"published curves and aperture.infer_response() for recovering A(E) from "
            f"one of them.")

    geomag_declination_deg, geomag_inclination_deg = physics.default_field_for_site(
        map_grid.center_lat, origin_lon, geomag_declination_deg, geomag_inclination_deg)

    # RFI sources as pixel coordinates, weighted by radius as a crude strength proxy.
    # The scan can then drop the ones terrain hides, which a circular exclusion cannot.
    rfi_zones_px = []
    if rfi_zones:
        for item in rfi_zones:
            if item[0] == 'circle':
                _, zlat, zlon, zrad_km, _ = item
                rfi_zones_px.append((
                    (origin_lat - zlat) / map_grid.cell_size_deg,
                    (zlon - origin_lon) / map_grid.cell_size_deg,
                    float(zrad_km)))

    scan_params = dict(
        elev_min_deg=elev_min_deg, elev_max_deg=elev_max_deg, n_elev_bins=int(n_elev_bins),
        n_azimuths=int(n_azimuths), half_width_deg=azimuth_half_width_deg,
        # How far the profile is walked, which is *not* the same question as which
        # intersections are accepted. Column depth accumulates over the whole walk, so
        # tying the two meant a short-range search reported a depth set by where the
        # walk stopped rather than by the target's thickness -- at TAMBO's 5 km the
        # depth term scored ~1 for everything. Defaults to max_dist_km, preserving the
        # old behaviour where it is not set.
        max_range_m=(max_range_km if max_range_km else max_dist_km) * 1000.0,
        min_dist_km=min_dist_km, max_dist_km=max_dist_km,
        min_depth_gcm2=min_column_depth_gcm2, require_terrain=require_terrain,
        min_target_slope_deg=min_target_slope_deg,
        max_target_slope_deg=max_target_slope_deg,
        geomag_declination_deg=(geomag_declination_deg if use_geomagnetic else None),
        geomag_inclination_deg=(geomag_inclination_deg if use_geomagnetic else None),
        frequency_mhz=fresnel_frequency_mhz, bilinear=bilinear_sampling,
        shower_offset_m=shower_development_m,
        antenna_height_m=antenna_height_m,
        near_field_m=(fresnel_near_field_m if exclude_near_field else 0.0),
        # Particle geometry is not refracted; only the radio path is
        earth_radius_m=arrival_scan.TRUE_EARTH_RADIUS_M,
        radio_earth_radius_m=(arrival_scan.earth_radius_for_k(refraction_k) if refraction_k
                              else arrival_scan.RADIO_EARTH_RADIUS_M),
    )


    # Store explicit params for final JSON export
    export_params = {
        "dem": dem_path, "origin": [origin_lat, origin_lon],
        "origin_source": origin_source,
        "cell_size_deg": cell_size_deg, "cell_size_y_m": cell_size_y,
        "cell_size_x_m": cell_size_x, "cell_size_center_lat": map_grid.center_lat,
        "cell_size_source": map_grid.source, "candidate_stride": candidate_stride,
        "slope_baseline_m": slope_baseline_m,
        "energy_min_pev": energy_min_pev, "energy_max_pev": energy_max_pev,
        "scan": scan_params,
        "refraction_k": refraction_k, "fresnel_frequency_mhz": fresnel_frequency_mhz,
        "antenna_height_m": antenna_height_m, "fresnel_near_field_m": fresnel_near_field_m,
        "exclude_near_field": exclude_near_field,
        "depth_band_gcm2": depth_band_gcm2, "score_composition": score_composition,
        "score_weights": score_weights, "distance_band_m": distance_band_m,
        "solid_angle_half_sr": solid_angle_half_sr, "clearance_full_at": clearance_full_at,
        "grammage_band_fraction": grammage_band_fraction,
        "min_score": min_score,
        "geomag_declination_deg": geomag_declination_deg,
        "geomag_inclination_deg": geomag_inclination_deg,
        "nu_interaction_length_gcm2": nu_interaction_length_gcm2,
        "decay_energy_pev": decay_energy_pev,
        "decay_energy_min_pev": decay_energy_min_pev,
        "decay_energy_max_pev": decay_energy_max_pev,
        "decay_spectral_index": decay_spectral_index,
        "shower_development_m": shower_development_m, "gap_close_km": gap_close_km,
        "max_range_km": max_range_km, "score_percentile": score_percentile,
        "decay_weight_by": decay_weight_by, "decay_response_csv": decay_response_csv,
        "stop_at_target": stop_at_target,
        "min_target_slope_deg": min_target_slope_deg,
        "max_target_slope_deg": max_target_slope_deg,
        "use_geomagnetic": use_geomagnetic, "grammage_mode": grammage_mode,
        "grammage_band_gcm2": grammage_band_gcm2,
        "grammage_maturity_gcm2": grammage_maturity_gcm2,
        "muon_shielding_km": muon_shielding_km,
        "bilinear_sampling": bilinear_sampling,
        "target": target_antennas, "spacing_km": antenna_spacing_km,
        "min_dist_km": min_dist_km, "max_dist_km": max_dist_km,
        "min_sub_array": min_sub_array_size,
        "grid_type": grid_type, "road_map": road_map_path,
        "downsample_factor": downsample_factor,
        "min_altitude": min_altitude, "max_altitude": max_altitude,
        "min_slope_deg": min_slope_deg, "max_slope_deg": max_slope_deg,
        # Screening and shaping knobs that were resolved but never recorded, so a
        # summary written from the results file could not name the filter that did
        # the work. They cost nothing to carry and the funnel is unreadable without
        # them: "after pruning (< 2.0 km wide)" means little if min_width_km is absent.
        "min_width_km": min_width_km,
        "min_aspect_deg": min_aspect_deg, "max_aspect_deg": max_aspect_deg,
        "max_road_dist_km": max_road_dist_km,
        "search_mode": search_mode, "region_name": region_name,
        "rfi_zone_count": len(rfi_zones) if rfi_zones else 0,
        "tile_size": tile_size, "resume": resume, "resume_dir": resume_dir,
        "num_cores": num_cores
    }
    
    print(f"\n{C.HEADER}============================================={C.RESET}")
    print(f"   {C.BOLD}OROSCOPE SITE SEARCH: RUN PARAMETERS{C.RESET}")
    print(f"{C.HEADER}============================================={C.RESET}")
    print(f"   -> DEM File: {C.MAGENTA}{dem_path}{C.RESET}")
    print(f"   -> Origin: {C.MAGENTA}{origin_lat}, {origin_lon}{C.RESET}")
    print(f"   -> Target: {C.MAGENTA}{target_antennas} antennas{C.RESET}")
    print(f"   -> Spacing: {C.MAGENTA}{antenna_spacing_km} km{C.RESET} ({grid_type} grid)")
    print(f"   -> Min Width: {C.MAGENTA}{min_width_km} km{C.RESET}")
    print(f"   -> Slope Range: {C.MAGENTA}{min_slope_deg}° to {max_slope_deg}°{C.RESET}")
    print(f"   -> Target Dist: {C.MAGENTA}{min_dist_km:g} - {max_dist_km:g} km{C.RESET}{energy_note}")
    print(f"      Arrival window: {C.MAGENTA}{elev_min_deg:g}° to {elev_max_deg:g}°{C.RESET}"
          f" in {n_elev_bins} bins, {n_azimuths} azimuths"
          f"{f' within ±{azimuth_half_width_deg:g}° of aspect' if azimuth_half_width_deg is not None else ' (full sweep)'}")
    print(f"      Requires: {C.MAGENTA}{'rock' if require_terrain else 'clear sky'}{C.RESET}"
          f", min column depth {min_column_depth_gcm2:,.0f} g/cm²")
    _lo = arrival_scan.energy_pev_for_decay_length(min_dist_km * 1000.0)
    _hi = arrival_scan.energy_pev_for_decay_length(max_dist_km * 1000.0)
    print(f"      Baseline implies tau energies {C.MAGENTA}{_lo:.3g} - {_hi:.3g} PeV{C.RESET}")
    print(f"   -> Downsample Factor: {C.MAGENTA}{downsample_factor}{C.RESET}")
    print(f"   -> Resolution: {C.MAGENTA}{cell_size_deg:.8f} deg/px{C.RESET} [{map_grid.source}]")
    print(f"   -> Pixel Size: {C.MAGENTA}{cell_size_y:.2f} m N-S x {cell_size_x:.2f} m E-W{C.RESET} (at lat {map_grid.center_lat:.3f})")
    print(f"   -> Candidate Stride: {C.MAGENTA}every {candidate_stride} px{C.RESET}")
    warn_stride_outruns_closing(candidate_stride, cell_size_y,
                                gap_close_km, antenna_spacing_km)
    _sb = f"{slope_baseline_m:.0f} m" if slope_baseline_m else "native DEM resolution"
    print(f"   -> Slope Baseline: {C.MAGENTA}{_sb}{C.RESET}")
    print(f"   -> Memory: Tile Size {C.MAGENTA}{tile_size}x{tile_size} px{C.RESET}")
    if road_map_path:
        print(f"   -> Logistics: Require road within {C.MAGENTA}{max_road_dist_km} km{C.RESET}")
    if min_altitude or max_altitude:
        min_s = f"{min_altitude}m" if min_altitude else "0m"
        max_s = f"{max_altitude}m" if max_altitude else "Inf"
        print(f"   -> Altitude: {C.MAGENTA}{min_s} < h < {max_s}{C.RESET}")
    if min_aspect_deg is not None and max_aspect_deg is not None:
        print(f"   -> Aspect Range: {C.MAGENTA}{min_aspect_deg}° to {max_aspect_deg}°{C.RESET}") 
    print(f"   -> RFI Zones: {C.MAGENTA}{len(rfi_zones) if rfi_zones else 0} active{C.RESET} (Numba Optimized)")

    print(f"\n{C.HEADER}============================================={C.RESET}")
    print(f"   {C.BOLD}SYSTEM & RESOURCE REPORT{C.RESET}")
    print(f"{C.HEADER}============================================={C.RESET}")
    sys_cores = multiprocessing.cpu_count()
    active_cores = sys_cores if num_cores == -1 else int(num_cores)
    # Cap Numba's own thread pool too, otherwise the parallel kernels ignore the request
    if HAS_NUMBA:
        try:
            numba.set_num_threads(active_cores)
        except Exception as e:
            print(f"   {C.WARN}{Icon.WARN}Could not limit Numba threads: {e}{C.RESET}")
    print(f"   -> CPU Cores: {C.MAGENTA}{active_cores} allocated{C.RESET} (System Max: {sys_cores})")
    print(f"   -> Numba JIT: {C.OK}ENABLED{C.RESET}" if HAS_NUMBA else f"   -> Numba JIT: {C.FAIL}DISABLED{C.RESET}")
    if psutil:
        mem = psutil.virtual_memory()
        print(f"   -> System RAM: {C.MAGENTA}{mem.total/1024**3:.1f} GB{C.RESET} (Free: {mem.available/1024**3:.1f} GB)")
    
    print(f"   -> Working Dir: {C.MAGENTA}{run_output_dir}{C.RESET}")
    if resume:
        print(f"   -> Resuming From: {C.WARN}{resume_dir}{C.RESET}")
    print(f"{C.HEADER}============================================={C.RESET}\n")

    t_start_total = time.time()

    # Tracking variables for safe cleanup
    path_A = None
    path_B = None
    success_flag = False
    results = None
    generated_files = []

    # Run accounting: per-stage wall time, per-filter survivor counts, and provenance
    timings = {}
    funnel = Funnel()
    region_stats = {}
    provenance = collect_provenance(dem_path, map_grid)

    try:
        # Step 1: Disk Setup
        print(f"{C.BOLD}[1/6]{C.RESET} {Icon.MAP}Loading Map Data...")
        t0 = time.time()
        elevation, rows, cols, path_A, path_B, buf_a, is_resuming = load_dem_and_init_buffers(dem_path, run_output_dir, resume, resume_dir)
        est_disk_gb = (rows * cols * 2) / (1024**3)
        print(f"      Map: {C.MAGENTA}{rows} x {cols}{C.RESET} pixels")
        print(f"      Estimated Temp Disk Usage: {C.MAGENTA}~{est_disk_gb:.2f} GB{C.RESET}")
        timings["load_dem"] = time.time() - t0
        print(f"      Time: {timings['load_dem']:.2f}s")

        if not is_resuming:
            # Step 2: Topographic Screen
            print(f"\n{C.BOLD}[2/6]{C.RESET} {Icon.GEAR}Identifying Candidates...")
            t0 = time.time()
            candidates_arr = get_candidates_chunked(
                elevation, map_grid, rfi_zones, origin_lat, origin_lon,
                min_alt=min_altitude, max_alt=max_altitude,
                road_map_path=road_map_path, max_road_dist_km=max_road_dist_km,
                min_aspect_deg=min_aspect_deg, max_aspect_deg=max_aspect_deg,
                min_slope_deg=min_slope_deg, max_slope_deg=max_slope_deg,
                tile_size=tile_size, candidate_stride=candidate_stride,
                slope_baseline_m=slope_baseline_m, funnel=funnel
            )
            total = candidates_arr.shape[0]
            timings["topographic_screen"] = time.time() - t0
            print(f"      Time: {timings['topographic_screen']:.2f}s")

            if total == 0:
                print(f"      {Icon.CROSS}{C.WARN}No viable candidates found in topographic pass.{C.RESET}")
                success_flag = True
                # Still a result, and the one where the funnel earns its keep: the
                # stage where the count collapsed is the whole answer. Returning it
                # (rather than None) means the explanation below can say which
                # screening filter emptied the map.
                results = {
                    "timestamp": datetime.now().isoformat(),
                    "mode": search_mode,
                    "parameters": export_params,
                    "results": {"total_sites": 0, "total_capacity": 0, "sites": []},
                    "funnel": funnel.as_dict(),
                    "regions": {},
                    "timings_sec": timings,
                    "aperture": {},
                    "provenance": provenance,
                }
                return results
    
            # Step 3: Physics Simulation
            print(f"\n{C.BOLD}[3/6]{C.RESET} {Icon.GEAR}Ray Tracing ({C.MAGENTA}{total}{C.RESET} candidates)...")
            t0 = time.time()
            n_hits, observables = run_arrival_scan(
                    candidates_arr, elevation, map_grid, buf_a, scan_params,
                    score_config={"depth_band_gcm2": depth_band_gcm2,
                                  "composition": score_composition,
                                  "weights": score_weights,
                                  "nu_interaction_length_gcm2": nu_interaction_length_gcm2,
                                  "spacing_m": antenna_spacing_km * 1000.0,
                                  "grammage_mode": grammage_mode,
                                  "grammage_band_gcm2": grammage_band_gcm2,
                                  "grammage_maturity_gcm2": grammage_maturity_gcm2,
                                  "distance_band_m": distance_band_m,
                                  "decay_energy_pev": decay_energy_pev,
                                  "decay_energy_min_pev": decay_energy_min_pev,
                                  "decay_energy_max_pev": decay_energy_max_pev,
                                  "decay_spectral_index": decay_spectral_index,
                                  "decay_weight_by": decay_weight_by,
                                  "decay_response": _decay_response,
                                  "shower_development_m": shower_development_m,
                                  "solid_angle_half_sr": solid_angle_half_sr,
                                  "clearance_full_at": clearance_full_at,
                                  "muon_shielding_km": muon_shielding_km},
                min_score=min_score, rfi_zones_px=rfi_zones_px,
                score_percentile=score_percentile)
            funnel.add("directions accepted", n_hits)
            if min_score > 0:
                funnel.add(f"score >= {min_score:g}", n_hits)
            timings["ray_tracing"] = time.time() - t0
            print(f"      Time: {timings['ray_tracing']:.2f}s")
        else:
            print(f"\n{C.BOLD}[2/6 & 3/6]{C.RESET} {Icon.GEAR}Resuming: Skipping candidate search and ray tracing...")

        del buf_a  # Release memory map write lock to flush to disk securely

        # Step 4: Spatial Pruning
        print(f"\n{C.BOLD}[4/6]{C.RESET} {Icon.BROOM}Cleaning Shapes...")
        t0 = time.time()
        n_closed, n_pruned = clean_shape_artifacts(path_A, path_B, rows, cols, cell_size_y, cell_size_x, antenna_spacing_km, min_width_km, tile_size, gap_close_km)
        funnel.add("after gap closing", n_closed)
        funnel.add(f"after pruning (< {min_width_km} km wide)", n_pruned)
        timings["morphology"] = time.time() - t0
        print(f"      Time: {timings['morphology']:.2f}s")

        # Step 5: Capacity Analysis
        print(f"\n{C.BOLD}[5/6]{C.RESET} {Icon.INFO}Final Analysis...")
        t0 = time.time()
        small_final, labeled_viz, site_details, cumulative_capacity, count, region_stats = analyze_sites_and_capacity(
            path_A, elevation, rows, cols, cell_size_y, cell_size_x, downsample_factor, search_mode,
            target_antennas, min_sub_array_size, antenna_spacing_km, grid_type, funnel=funnel,
            origin_lat=origin_lat, origin_lon=origin_lon, cell_size_deg=cell_size_deg,
            candidates_arr=candidates_arr, observables=observables,
            stop_at_target=stop_at_target
        )
        if search_mode == 'distributed':
            print(f"      Distributed: {C.MAGENTA}{count}{C.RESET} sites found.")
            print(f"      Total Cap: {C.MAGENTA}{cumulative_capacity}{C.RESET} (Target: {target_antennas})")
        else:
            print(f"      Single: {C.MAGENTA}{count}{C.RESET} valid sites found.")
        timings["capacity_analysis"] = time.time() - t0
        print(f"      Time: {timings['capacity_analysis']:.2f}s")

        # Step 6: Create Outputs
        print(f"\n{C.BOLD}[6/6]{C.RESET} {Icon.DISK}Saving & Visualization...")
        t0 = time.time()
        generated_files, results = generate_visualizations_and_outputs(
            dem_path, elevation, small_final, labeled_viz, site_details, count, cumulative_capacity,
            origin_lat, origin_lon, map_grid, downsample_factor, generate_kml, run_output_dir,
            output_image_format, rfi_zones, search_mode, grid_type, antenna_spacing_km,
            min_altitude, max_altitude, region_name, export_params,
            run_info={"funnel": funnel.as_dict(), "regions": region_stats,
                      "timings_sec": timings, "provenance": provenance,
                      "aperture": aperture_mod.summarize_sites(
                          site_details, min_dist_km * 1000.0, max_dist_km * 1000.0,
                          np.logspace(0, 5, 26))}
        )
        timings["outputs"] = time.time() - t0
        print(f"      Time Elapsed: {timings['outputs']:.2f}s")

        # Funnel Report: where the candidate pixels went
        print(f"\n{C.HEADER}============================================={C.RESET}")
        print(f"   {C.BOLD}SELECTION FUNNEL{C.RESET}")
        print(f"{C.HEADER}============================================={C.RESET}")
        print(funnel.render())
        if region_stats:
            print(f"\n   Regions: {C.MAGENTA}{region_stats['labelled_regions']}{C.RESET} labelled"
                  f" -> {C.MAGENTA}{region_stats['passed_area_threshold']}{C.RESET} above area threshold"
                  f" ({region_stats['required_pixels_per_region']:,} px)"
                  f" -> {C.MAGENTA}{region_stats['passed_capacity_threshold']}{C.RESET} above capacity threshold"
                  f" ({region_stats['capacity_threshold_antennas']:,} antennas)")

        # Results Table Printout
        print(f"\n{C.HEADER}============================================={C.RESET}")
        print(f"   {Icon.CHECK}{C.BOLD}RESULTS SUMMARY: {count} Sites Found{C.RESET}")
        print(f"{C.HEADER}============================================={C.RESET}")
        
        if count > 0:
            print(f"   {C.BOLD}{'ID':>4} | {'Area (km²)':>12} | {'Capacity':>10} | {'Grid':>6} | {'Facing'}{C.RESET}")
            print("   " + "-" * 50)
            for site in site_details:
                print(f"   {site['site_id']:>4} | {site['area_km2']:>12.2f} | {site['capacity_exact']:>10} | {site['grid_type']:>6} | {site['facing_direction']}")
        else:
            print(f"   {C.WARN}{Icon.CROSS}No valid sites met all constraints.{C.RESET}")

        # Print outputs generated block for log
        print(f"\n{C.HEADER}============================================={C.RESET}")
        print(f"   {C.BOLD}OUTPUTS GENERATED{C.RESET}")
        print(f"{C.HEADER}============================================={C.RESET}")
        for fpath in generated_files:
            print(f"   {Icon.INFO}{fpath}")

        # Mark as cleanly finished so finally block knows to purge temporary arrays
        success_flag = True
        return results

    finally:
        # Smart Cleanup
        if success_flag:
            try:
                if path_A and os.path.exists(path_A): os.remove(path_A)
                if path_B and os.path.exists(path_B): os.remove(path_B)
            except Exception: 
                pass
        else:
            print(f"\n   {C.FAIL}{Icon.CROSS}[!] Run did not complete successfully.{C.RESET}")
            print(f"   {C.WARN}Buffer files have been retained in the workspace.{C.RESET}")
            print(f"   {C.WARN}Resume this exact run later using:{C.RESET} --resume --resume_dir {os.path.abspath(run_output_dir)}")
            
        timings["total"] = time.time() - t_start_total
        print(f"\n{C.OK}Total Execution Time: {timings['total']:.2f} seconds{C.RESET}")
        print(f"{C.BOLD}Done.{C.RESET}")

        # Last, so it is what a reader is left with. In the `finally` so that the
        # run that found nothing is explained too -- it is the one most in need of it.
        if results is not None:
            results.setdefault("provenance", provenance)
            results["output_files"] = list(generated_files or [])
            emit_explanation(results, run_output_dir if explain else None,
                             print_it=explain)

# Custom Logger Interceptor
class TeeLogger:
    """Duplicates stream writes to both the original terminal and an attached log file.
    Strips ANSI color codes from the log file output to maintain readability."""
    def __init__(self, terminal, log_file):
        self.terminal = terminal
        self.log_file = log_file
        # Regex to match standard ANSI escape sequences (colors, bold, etc.)
        self.ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
        
    def write(self, message):
        self.terminal.write(message)
        # Strip colors before writing to the text log
        clean_message = self.ansi_escape.sub('', message)
        self.log_file.write(clean_message)
        
    def flush(self):
        self.terminal.flush()
        self.log_file.flush()


# Bound once, to the real function, rather than read inside
# config_to_pipeline_kwargs on every call. Reading it at call time meant the filter
# followed whatever `find_grand_regions_interactive` happened to be at that moment --
# so a test double, or any decorator, presented a bare (*args, **kwargs) signature and
# the translation quietly dropped every parameter it was supposed to pass.
_PIPELINE_PARAMS = frozenset(
    inspect.signature(find_grand_regions_interactive).parameters)


def main():
    """
    Command-line entry point: parses arguments, reconciles them against the config
    file and the fallbacks, validates, and runs one search.

    Kept as a function rather than a bare ``__main__`` block so the console script
    declared in pyproject.toml has something to call, and so the argument handling
    can be exercised from a test without spawning a subprocess.
    """
    parser = argparse.ArgumentParser(
        description="Oroscope - terrain site search for particle-astrophysics "
                    "observatories. GRAND and TAMBO are configurations of the same "
                    "engine, not separate code paths.")
    
    # Made DEM and Origin optional here so they can be exclusively provided via config or fallbacks
    parser.add_argument("--dem_path", type=str, help="Path to the Digital Elevation Model (.tif) file.")
    parser.add_argument("--origin_lat", type=float, help="Reference origin latitude (e.g., -10.228).")
    parser.add_argument("--origin_lon", type=float, help="Reference origin longitude (e.g., -78.076).")
    
    # Configuration and Layout Arguments
    parser.add_argument("--target_antennas", type=int, default=10000, help="Total target capacity for the array (default: 10000).")
    parser.add_argument("--min_width_km", type=float, default=2.0, help="Minimum acceptable width of the array site in km (default: 2.0).")
    parser.add_argument("--min_altitude", type=float, default=None, help="Minimum allowable altitude in meters (optional).")
    parser.add_argument("--max_altitude", type=float, default=None, help="Maximum allowable altitude in meters (optional).")
    parser.add_argument("--antenna_spacing_km", type=float, default=1.0, help="Distance between antennas in km (default: 1.0).")
    parser.add_argument("--min_dist_km", type=float, default=10.0, help="Minimum required distance to target mountain in km (default: 10.0).")
    parser.add_argument("--max_dist_km", type=float, default=80.0, help="Maximum required distance to target mountain in km (default: 80.0).")
    parser.add_argument("--grid_type", type=str, choices=['square', 'hex'], default='hex', help="Antenna layout grid type (default: 'hex').")
    
    # Internal Math & Physics Parameters
    parser.add_argument("--min_slope_deg", type=float, default=3.0, help="Minimum terrain steepness in degrees (default: 3.0).")
    parser.add_argument("--max_slope_deg", type=float, default=25.0, help="Maximum terrain steepness in degrees (default: 25.0).")
    parser.add_argument("--downsample_factor", type=int, default=4, help="Internal capacity mask downsampling factor for processing speed (default: 4).")
    parser.add_argument("--cell_size_deg", type=float, default=None, help="Map resolution in degrees per pixel. Defaults to reading the DEM's GeoTIFF tags.")
    parser.add_argument("--slope_baseline_m", type=float, default=None, help="Ground distance in metres over which slope is measured. Default: the DEM's native resolution.")
    parser.add_argument("--candidate_stride", type=int, default=5, help="Keep every Nth candidate pixel before ray tracing (default: 5). Use 1 for no thinning.")
    parser.add_argument("--tile_size", type=int, default=2048, help="Size of the square memory chunk for RAM management (default: 2048).")
    parser.add_argument("--num_cores", type=int, default=-1, help="Number of CPU cores to use. Set to -1 to use all available cores (default: -1).")
    
    # Arrival-direction scan (roadmap phase 1)
    parser.add_argument("--energy_min_pev", type=float, default=None, help="Lower tau energy in PeV. With --energy_max_pev, derives the decay-baseline distance window.")
    parser.add_argument("--energy_max_pev", type=float, default=None, help="Upper tau energy in PeV.")
    parser.add_argument("--n_azimuths", type=int, default=9, help="Azimuths scanned per candidate in scan mode (default: 9).")
    parser.add_argument("--azimuth_half_width_deg", type=float, default=60.0, help="Half-width of the azimuth fan about the aspect. Use -1 for a full 360 sweep (default: 60).")
    parser.add_argument("--elev_min_deg", type=float, default=-3.0, help="Lower edge of the accepted arrival elevation window (default: -3).")
    parser.add_argument("--elev_max_deg", type=float, default=3.0, help="Upper edge of the accepted arrival elevation window (default: +3).")
    parser.add_argument("--n_elev_bins", type=int, default=12, help="Elevation bins across the window (default: 12). Nearly free: cost scales with azimuths.")
    parser.add_argument("--min_column_depth_gcm2", type=float, default=0.0, help="Column depth a direction must have to count, in g/cm2 (default: 0).")
    parser.add_argument("--require_sky", action="store_true", dest="require_sky", help="Invert the test: accept directions that reach clear sky, for cosmic-ray style channels.")

    parser.add_argument("--fresnel_frequency_mhz", type=float, default=None, help="Radio band for the Fresnel clearance measurement, e.g. 50. Omitted skips the second pass.")
    parser.add_argument("--antenna_height_m", type=float, default=2.0, help="Antenna height above ground, for the Fresnel measurement (default: 2).")
    parser.add_argument("--include_near_field", action="store_false", dest="exclude_near_field", help="Measure Fresnel clearance from the antenna outward instead of skipping the near field. Included for study: the result is then dominated by ground beside the antenna rather than by intervening terrain.")
    parser.add_argument("--fresnel_near_field_m", type=float, default=500.0, help="Skip this much of the path when measuring Fresnel clearance (default: 500). Below ~500 m the measure is dominated by ground beside the antenna rather than by intervening terrain.")
    parser.add_argument("--nearest_sampling", action="store_false", dest="bilinear_sampling", help="Sample terrain profiles at pixel centres instead of interpolating. Faster, but treats terrain as blocky, which over-estimates how much it blocks a ray.")
    parser.add_argument("--muon_shielding_km", type=float, default=None, help="Rock overburden required along the arrival direction to reject atmospheric muons, in km (TAMBO quotes >4). A floor on column depth, not a band.")
    parser.add_argument("--geomag_declination_deg", type=float, default=None, help="Geomagnetic declination, degrees east of north. Defaults to the Arequipa IGRF 2026 value (-6.9); supply the IGRF value for other regions.")
    parser.add_argument("--geomag_inclination_deg", type=float, default=None, help="Geomagnetic inclination, degrees, positive downward. Defaults to a centered-dipole estimate at the DEM's own centre, so it follows the site automatically.")
    parser.add_argument("--no_geomagnetic", action="store_false", dest="use_geomagnetic", help="Ignore the geomagnetic angle and weight all directions equally.")
    parser.add_argument("--grammage_mode", type=str, choices=['radio', 'particle'], default='radio', help="How atmospheric depth is scored. 'radio' is a maturity threshold, since emission comes from shower maximum and then propagates through transparent air. 'particle' is a band, since particle content dies after maximum (default: radio).")
    parser.add_argument("--grammage_band_gcm2", type=float, nargs=2, default=None, metavar=("LO", "HI"), help="Atmospheric depth band scoring 1 in 'particle' mode, in g/cm2. Defaults to (X_max, 4*X_max) = (700, 2800), which suits a long path to a distant target. A short crossing gives far less: Colca supplies about 170 g/cm2, so a detector there sees a shower that is still developing and this band must be lowered or nothing scores.")
    parser.add_argument("--grammage_maturity_gcm2", type=float, default=None, help="Atmospheric depth at which the 'radio' maturity ramp reaches 1, in g/cm2 (default: X_max = 700).")
    parser.add_argument("--decay_energy_pev", type=float, default=None, help="Tau energy, in PeV, at which to score the probability that it decays in the gap with room left for a shower. Left out by default because the probability is strongly energy-dependent and one number cannot stand in for a spectrum. Matters most across a canyon: at 1 EeV the decay length is ~49 km against a ~3 km crossing.")
    parser.add_argument("--decay_weight_by", type=str, default="flux",
                        choices=list(physics.DECAY_WEIGHTINGS),
                        help="What weights the spectrum-folded decay probability. "
                             "'flux' (default) asks what fraction of arriving neutrinos "
                             "decay usefully, and is what every published number here "
                             "used. 'acceptance' asks the same over the energies the "
                             "detector responds to, with no assumed spectrum -- useful "
                             "precisely because the spectral index is an assumption. "
                             "'flux_times_acceptance' is the event-rate integrand. The "
                             "latter two need --decay_response_csv.")
    parser.add_argument("--decay_response_csv", type=str, default=None,
                        help="Two-column CSV of energy in PeV against relative detector "
                             "response A(E), for the acceptance weightings. data/ holds "
                             "the published integral curves; aperture.infer_response() "
                             "recovers A(E) from one by dividing out the geometric "
                             "model.")
    parser.add_argument("--max_range_km", type=float, default=None, help="How far to walk each profile, in km. Defaults to max_dist_km. Worth setting larger for a short-range search: column depth accumulates over the whole walk, so tying the two makes the reported depth a property of where the walk stopped rather than of the target's thickness.")
    parser.add_argument("--score_percentile", type=float, default=None, help="Keep this percentage of viable candidates, ranked by score, instead of cutting at an absolute --min_score. Preferred: the default score is a product whose distribution piles up near zero, so an absolute threshold sits on a cliff, while a percentile is scale-free.")
    parser.add_argument("--stop_at_target", action="store_true", help="In distributed mode, stop selecting sites once target_antennas is reached. Sites are ranked by capacity, so this reports the best sites for the array actually wanted rather than every patch of qualifying ground.")
    parser.add_argument("--max_memory_gb", type=float, default=None, help="Ceiling on this process's address space, in GiB. Defaults to 80%% of what the system reports available, so a search that outgrows the machine fails with MemoryError instead of inviting the OOM killer to choose a victim. 0 disables the cap.")
    parser.add_argument("--decay_energy_min_pev", type=float, default=None, help="Lower end of the tau energy range for the decay term. With --decay_energy_max_pev this folds the decay probability over a power-law spectrum, which is the defensible form: the probability runs over three decades across one experiment's reach, so a single energy chooses the answer rather than approximating it.")
    parser.add_argument("--decay_energy_max_pev", type=float, default=None, help="Upper end of that range, in PeV.")
    parser.add_argument("--decay_spectral_index", type=float, nargs="+", default=None, metavar="GAMMA", help="Spectral index gamma in dN/dE ~ E^-gamma for the folded decay term (default: 2.0). Give one value to pin the spectrum, or two to marginalise uniformly over that range when the index is not known -- a flat prior says so rather than pretending to a value. A softer spectrum weights low energies, where the tau decays readily, so it drives the term toward 1.")
    parser.add_argument("--shower_development_m", type=float, default=3000.0, help="Path the shower needs after the tau decays, in metres (default: 3000). Used both by the decay term and as the far endpoint of the Fresnel clearance measurement.")
    parser.add_argument("--gap_close_km", type=float, default=None, help="Size of the morphological closing element that fills gaps between accepted pixels, in km. Defaults to antenna_spacing_km, which couples two unrelated things. Closing more than doubles the reported area on real terrain (measured 2.29x at Colca), so this is worth setting deliberately; 0 disables it.")
    parser.add_argument("--min_target_slope_deg", type=float, default=None, help="Require the terrain a ray strikes to be at least this steep, measured along the arrival azimuth. Unset by default, which asks only that rock is present -- true almost everywhere in the Andes. TAMBO's tau exits a canyon *wall*, so this is what separates a canyon from a hillside.")
    parser.add_argument("--max_target_slope_deg", type=float, default=None, help="Upper bound on the struck terrain's slope along the arrival azimuth. Unset by default. Note a ceiling does not empty the result: a flat valley floor passes any ceiling, so this removes walls rather than everything.")
    parser.add_argument("--grammage_band_fraction", type=float, default=None, help="When the shower band is derived from an energy range, the fraction of peak particle content that still counts as a usable shower (default: 0.1). Lower admits younger and older showers, so it widens the band and accepts narrower canyons.")
    parser.add_argument("--shower_elongation_rate_gcm2", type=float, default=None, help="How much deeper shower maximum sits per decade of primary energy, in g/cm2 (default: 55, the usual hadronic value; a purely electromagnetic cascade is nearer 85).")
    parser.add_argument("--shower_lambda_gcm2", type=float, default=None, help="Gaisser-Hillas interaction length setting how fast the shower profile rises and falls, in g/cm2 (default: 70).")
    parser.add_argument("--solid_angle_half_sr", type=float, default=None, help="Accepted solid angle scoring 0.5, in steradians (default: 0.05). This is a GRAND-scale value: an experiment looking across a canyon sees far more sky, and leaving it at 0.05 saturates the term so it stops discriminating.")
    parser.add_argument("--distance_band_m", type=float, nargs=2, default=None, metavar=("LO", "HI"), help="Exit-point distance band scoring 1, in metres. Defaults to the configured decay-baseline window.")
    parser.add_argument("--clearance_full_at", type=float, default=None, help="Fresnel clearance ratio, in first-Fresnel radii, that scores 1 (default: 1.0).")
    parser.add_argument("--score_weights", type=str, default=None, help="Per-component weights for --score_composition weighted, as name=value pairs, e.g. 'shower=2,solid_angle=1,depth=0.5'. Components not named default to weight 1.")
    parser.add_argument("--nu_interaction_length_gcm2", type=float, default=None, help="Neutrino interaction length for the Earth-chord attenuation term, g/cm2 (order 1e8 near an EeV). Omitted reports the chord without weighting by it.")
    parser.add_argument("--refraction_k", type=float, default=None, help="Refraction k-factor for the RADIO path only (default: 4/3). Particle trajectories always use the true Earth radius, since neutrinos and taus are not refracted.")
    parser.add_argument("--depth_band_gcm2", type=float, nargs=2, default=None, metavar=("LO", "HI"), help="Column depth band scoring 1, in g/cm2. The tau must be produced and must escape, so this is a band, not a floor.")
    parser.add_argument("--score_composition", type=str, choices=['product', 'mean', 'min'], default='product', help="How component scores combine (default: product).")
    parser.add_argument("--min_score", type=float, default=0.0, help="Discard candidates scoring below this (default: 0, keep all).")

    # Logistics and Geography Arguments
    parser.add_argument("--rfi_zones", type=str, default='none', help="Can be preset ('lima', 'arequipa') or a valid JSON string outlining custom exclusion zones.")
    parser.add_argument("--road_map_path", type=str, default=None, help="Path to a raster mapping distance-to-roads (optional).")
    parser.add_argument("--max_road_dist_km", type=float, default=20.0, help="Maximum distance allowed from a road in km (default: 20.0).")
    
    # Execution modes and constraints
    parser.add_argument("--search_mode", type=str, choices=['single', 'distributed'], default='distributed', help="'single' finds one monolithic site, 'distributed' allows sub-arrays.")
    parser.add_argument("--min_sub_array_size", type=int, default=500, help="Minimum required capacity for a sub-array to be considered valid (default: 500).")
    parser.add_argument("--min_aspect_deg", type=float, default=None, help="Minimum bound for site facing direction in degrees (0-360).")
    parser.add_argument("--max_aspect_deg", type=float, default=None, help="Maximum bound for site facing direction in degrees (0-360).")
    
    # Metadata and output flags
    parser.add_argument("--region_name", type=str, default=None, help="Cosmetic region name to print on the final visualization chart.")
    parser.add_argument("--generate_kml", action="store_true", help="Include this flag to generate a Google Earth KML file of the findings.")
    parser.add_argument("--no_print_info", action="store_false", dest="print_info", help="Include this flag to skip printing the detailed explanatory text.")
    parser.add_argument("--no_explain", action="store_false", dest="explain", help="Skip the plain-language summary of the run. It is printed by default, and saved as explanation.txt beside the results: what was found, which constraint set the size of the answer, what held the surviving sites back, and which of the numbers are assumptions rather than measurements. A results file can be re-explained at any time with explain.explain_results().")
    
    # IO / Configs mapping & Tools
    parser.add_argument("--config_path", type=str, default=None, help="Path to external JSON configuration file.")
    parser.add_argument("--output_directory_base_with_given_json", type=str, default="../output/", help="Base directory for outputs when a JSON config is supplied (default: ../output/).")
    parser.add_argument("--output_image_format", type=str, default="png", help="Format of the saved map visual, e.g., png, pdf, svg (default: png).")
    parser.add_argument("--resume", action="store_true", help="Include this flag to resume a previous run from the ray-tracing checkpoint.")
    parser.add_argument("--resume_dir", type=str, default=None, help="Path to an output folder from a previously failed run to resume from the ray-tracing checkpoint.")
    
    # Tool Generation Arguments
    parser.add_argument("--generate_config", type=str, default=None, help="Supply a filepath to generate a default JSON config template and exit.")
    parser.add_argument("--config_preset", type=str, choices=list(CONFIG_PRESETS), default='default', help="Optional presets to inject when using --generate_config.")

    args = parser.parse_args()
    explicit_cli = explicitly_passed(parser)

    # --- Tool Execution: Generate Config Template ---
    if args.generate_config:
        generate_config(args.generate_config, args.config_preset)
        print(f"Configuration file generated successfully at: {args.generate_config}"
              f" (Preset: {args.config_preset})")
        sys.exit(0)


    # 1. Initialize Configuration Maps
    config_params = load_config(args.config_path)

    # 2. Retrieve Fallbacks
    fallback_path = os.path.join("..", "config", "fallbacks.json")
    fallback_params = load_config(fallback_path)

    # 3. Determine Unified Logging and Output Directory Hierarchically
    # Same precedence as every other parameter: an explicitly typed option wins, then
    # the config file, then the fallbacks. This one resolved before the merge loop and
    # so kept the old rule -- the config beat the command line -- which meant
    # --output_directory_base_with_given_json was silently ignored whenever a config
    # set it, and every other flag on that command line was honoured.
    key = "output_directory_base_with_given_json"
    if key in explicit_cli:
        base_dir = getattr(args, key)
    elif key in config_params:
        base_dir = config_params[key]
    elif key in fallback_params:
        base_dir = fallback_params[key]
    else:
        base_dir = getattr(args, key)

    if args.config_path and os.path.exists(args.config_path):
        config_basename = os.path.splitext(os.path.basename(args.config_path))[0]
        # Relative to the configuration, not the working directory -- the same rule
        # load_config applies to dem_path, and for the same reason. Without it, fixing
        # the inputs to run from anywhere would leave the *outputs* still landing
        # wherever the caller happened to stand: from the repository root the default
        # "../output/" writes a sibling of the repository. A base typed on the command
        # line is the caller's own instruction and is left relative to where they are.
        if key not in explicit_cli and not os.path.isabs(base_dir):
            config_dir = os.path.dirname(os.path.abspath(args.config_path))
            base_dir = os.path.normpath(os.path.join(config_dir, base_dir))
        run_output_dir = os.path.join(base_dir, config_basename)
    else:
        timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_output_dir = os.path.join("..", "output", timestamp_str)

    os.makedirs(run_output_dir, exist_ok=True)
    
    # 4. Apply Custom Standard-Out / Standard-Error interceptors for the log file
    log_path = os.path.join(run_output_dir, "log.txt")
    log_file = open(log_path, "a", encoding="utf-8")

    saved_stdout, saved_stderr = sys.stdout, sys.stderr
    sys.stdout = TeeLogger(sys.stdout, log_file)
    sys.stderr = TeeLogger(sys.stderr, log_file)
    try:
        return _run_from_arguments(parser, args, explicit_cli, config_params,
                                   fallback_params, fallback_path, run_output_dir)
    finally:
        # Both the stream swap and the open log used to outlive the call. A process
        # that ran main() twice -- a test, a sweep, anything driving the CLI in a loop
        # -- stacked a TeeLogger on the previous one and leaked a file handle each
        # time, and the streams never came back.
        sys.stdout, sys.stderr = saved_stdout, saved_stderr
        log_file.close()


def _run_from_arguments(parser, args, explicit_cli, config_params, fallback_params,
                        fallback_path, run_output_dir):
    """
    The body of :func:`main`, once its logging is installed.

    Split out only so that the log file and the redirected streams can be restored in
    a ``finally`` without indenting two hundred lines. Not part of the public surface:
    call :func:`main`, or the pipeline directly.
    """

    # Ensure log captures initiation context
    print(f"\n{C.HEADER}================================================================================{C.RESET}")
    print(f"Execution started at: {datetime.now().isoformat()}")
    if args.config_path:
        print(f"Using config file: {os.path.abspath(args.config_path)}")
    else:
        print("No config file provided. Relying on CLI arguments and fallbacks.")
    print(f"Using fallbacks file: {os.path.abspath(fallback_path)}")
    print(f"Unified output directory initialized at: {os.path.abspath(run_output_dir)}")
    print(f"{C.HEADER}================================================================================{C.RESET}\n")

    # 5. Reconcile Configuration Strategy (Config > Fallback > CLI / Standard defaults)
    final_params = {}
    
    # Collect all available arguments parsed from CLI framework
    param_names = [action.dest for action in parser._actions if action.dest not in ('help', 'config_path', 'generate_config', 'config_preset')]
    
    for param in param_names:
        # An explicitly typed command-line option wins over everything. It used to lose
        # to the config file, silently: since --generate_config writes all 67 keys, a
        # generated config made every flag on the command line a no-op with no warning.
        if param in explicit_cli:
            final_params[param] = getattr(args, param)
            if param in config_params and config_params[param] != final_params[param]:
                print(f"{C.WARN}{Icon.WARN}Command line overrides config for '{param}': "
                      f"{config_params[param]!r} -> {final_params[param]!r}{C.RESET}")
        elif param in config_params:
            final_params[param] = config_params[param]
        elif param in fallback_params:
            final_params[param] = fallback_params[param]
            print(f"{C.WARN}{Icon.WARN}WARNING: Parameter '{param}' not explicitly specified in config. Using fallback value: {fallback_params[param]}{C.RESET}")
        else:
            final_params[param] = getattr(args, param)

    # Post-validation of absolutely required parameters to prevent early crashes during Fail-Fast
    if final_params.get('dem_path') is None:
        print(f"{C.FAIL}{Icon.CROSS}ERROR: Critical parameter 'dem_path' must be provided via config file, fallback, or CLI.{C.RESET}")
        sys.exit(1)
    if final_params.get('origin_lat') is None or final_params.get('origin_lon') is None:
        # Not an error any more: a standard geographic DEM carries its own corner, and
        # reading it removes the most error-prone input the tool had.
        probe_lat, probe_lon = read_dem_origin(final_params['dem_path'])
        if probe_lat is None:
            print(f"{C.FAIL}{Icon.CROSS}ERROR: 'origin_lat'/'origin_lon' were not given "
                  f"and the DEM carries no tiepoint to read them from.{C.RESET}")
            sys.exit(1)

    # Default resume_dir to run_output_dir if resume is True and resume_dir not provided
    if final_params.get('resume') and not final_params.get('resume_dir'):
        final_params['resume_dir'] = run_output_dir

    # 6. Run Pre-Flight Validation (Fail-Fast Mechanism)
    validate_parameters(final_params)

    # The memory estimate, the warning and the address-space cap now live in the
    # pipeline itself (preflight_memory), so a library caller gets all three without
    # having to know they exist. main() only passes the ceiling through.

    # Output explanations if requested (either via flag or configured true)
    if final_params.get('print_info', True):
        print_tool_explanation()

    # Execute the search. The configuration-to-pipeline translation lives in
    # config_to_pipeline_kwargs, so main(), the sensitivity child and the full-DEM
    # runner all use the same one and a new parameter is added in a single place.
    # run_output_dir is an override because main() derives it from the config's name.
    return run_from_config(final_params, run_output_dir=run_output_dir)


if __name__ == "__main__":
    main()
