import argparse
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

import arrival_scan
import aperture as aperture_mod
import physics
import scoring

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
    
    Parameters:
    - x, y (float): Coordinates of the test point.
    - poly_verts (ndarray): A list of (x, y) coordinates defining the polygon vertices.
    
    Returns:
    - bool: True if the point is inside the polygon, False otherwise.
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
    
    Parameters:
    - valid_rows, valid_cols (ndarray): Arrays containing the row and col coordinates to check.
    - poly_verts (ndarray): Polygon vertices.
    - mask_out (ndarray): A boolean array modified in-place. Sets to False if point is inside polygon.
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

    Parameters:
    - mask_chunk (ndarray): A 2D boolean array where True indicates valid terrain.
    - cell_size_y, cell_size_x (float): Ground size of one pixel, in metres. They differ
      on a geographic grid, which is why an equal ground spacing is a different number
      of pixels on each axis.
    - spacing_m (float): Distance between neighbouring detectors, in metres.
    - grid_type_code (int): 0 for 'square' grid, 1 for 'hexagonal' (triangular) grid.

    Returns:
    - int: The total number of detectors that fit inside the valid terrain.
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
        """Adds to a stage's running total, creating the stage on first use."""
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

def slope_baseline_pixels(map_grid, slope_baseline_m):
    """
    Converts a slope measurement baseline in metres to a per-axis window in pixels.

    Slope is scale-dependent: on real Andean terrain the median slope falls from
    ~17.8 deg measured over the DEM's native ~61 m to ~10.8 deg over 1 km, and the
    fraction passing a 3-25 deg band rises from 60% to 78%. Which of those is
    "the" slope depends on the footprint being deployed, so the baseline is an
    explicit parameter rather than an accident of the DEM's resolution.

    Returns (0, 0) when no baseline is requested, meaning the native gradient.
    """
    if not slope_baseline_m:
        return 0, 0
    ny = max(1, int(round(slope_baseline_m / map_grid.cell_size_y)))
    nx = max(1, int(round(slope_baseline_m / map_grid.cell_size_x)))
    return ny, nx


def terrain_gradients(elevation_block, cell_size_y, cell_size_x, smooth_y=0, smooth_x=0):
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

    Returns:
    - tuple(ndarray, ndarray): d/dy and d/dx, in metres per metre.
    """
    block = elevation_block
    if smooth_y > 1 or smooth_x > 1:
        block = uniform_filter(block, size=(max(1, smooth_y), max(1, smooth_x)), mode="nearest")
    return np.gradient(block, cell_size_y, cell_size_x)


def slope_band_gradient_sq(min_slope_deg, max_slope_deg):
    """
    The slope band restated as bounds on the squared gradient magnitude.

    slope = atan(|grad|) rises monotonically with |grad|, so

        min <= atan(sqrt(g)) <= max   <=>   tan(min)^2 <= g <= tan(max)^2

    which tests the same pixels without a sqrt or an arctan. Bounds at or beyond the
    vertical, and non-positive lower bounds, are returned as None meaning "unbounded":
    tan is singular at 90 degrees and every real gradient satisfies them anyway.

    Returns:
    - tuple(float or None, float or None): lower and upper bounds on dx^2 + dy^2.
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


def terrain_derivatives(elevation_block, cell_size_y, cell_size_x, smooth_y=0, smooth_x=0):
    """
    Slope and aspect over a stated measurement baseline.

    Returns:
    - tuple(ndarray, ndarray): slope in degrees, aspect in degrees clockwise from north.
    """
    dy, dx = terrain_gradients(elevation_block, cell_size_y, cell_size_x, smooth_y, smooth_x)
    slope = np.degrees(np.arctan(np.sqrt(dx ** 2 + dy ** 2)))
    aspect = np.degrees(np.arctan2(-dx, dy)) % 360
    return slope, aspect


def read_dem_geometry(dem_path):
    """
    Reads the angular pixel size and row count of a GeoTIFF DEM from its header.

    Standard geographic (EPSG:4326) DEMs such as SRTMGL1 or AW3D30 store the pixel
    size in degrees, which is what the georeferenced outputs (.tfw, .kml) require.
    Only the header is touched, so this stays cheap on multi-gigabyte files.

    Parameters:
    - dem_path (str): Path to the input elevation .tif file.

    Returns:
    - tuple(float or None, int or None): Pixel size in degrees and the number of
      rows, either of which is None when the file or the tag cannot be read.
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
MapGrid = namedtuple("MapGrid", "cell_size_deg cell_size_y cell_size_x center_lat source")

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

    Parameters:
    - dem_path (str): Path to the input elevation .tif file.
    - origin_lat (float): Latitude of the DEM's top-left corner.
    - cell_size_deg (float or None): Explicit override in degrees per pixel.

    Returns:
    - MapGrid: Angular pixel size, both metric pixel sizes, the centre latitude used
      for the longitude scaling, and where the resolution value came from.
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
    
    Returns:
    - elevation (ndarray): Memory mapped DEM.
    - rows, cols (int): Array dimensions.
    - path_A, path_B (str): Paths to the initialized boolean buffer arrays.
    - buf_a (ndarray): Open memory map of Buffer A.
    - is_resuming (bool): True if successfully loaded a previous physics buffer.
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

    Returns:
    - ndarray: Nx3 array of valid candidate pixels formatted as [row, col, aspect_degrees].
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
            except:
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
                     score_config=None, min_score=0.0, rfi_zones_px=None):
    """
    Step 3 alternative: scan arrival directions instead of casting one ray per pixel.

    Marks a candidate as valid when at least one accepted (azimuth, elevation)
    direction strikes rock within the decay-baseline window with enough column depth.
    See arrival_scan.py for the geometry.

    Returns:
    - tuple(int, dict): number of accepted candidates, and the per-candidate
      observables kept for per-site aggregation.
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

    accepted = (observables["cells"] > 0) & (total >= min_score)
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

    Returns:
    - dict mapping site id to summary statistics of that site's accepted candidates.
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

    Returns:
    - int: Number of set pixels in the result, counted while writing so the funnel
      accounting costs nothing extra.
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
    """
    h, w = structure.shape
    col = np.ones((h, 1), dtype=bool)
    row = np.ones((1, w), dtype=bool)
    grown = binary_dilation(binary_dilation(chunk, col), row)
    return binary_erosion(binary_erosion(grown, col), row)


def separable_opening(chunk, structure):
    """Binary opening with a rectangular element, separably. See :func:`separable_closing`."""
    h, w = structure.shape
    col = np.ones((h, 1), dtype=bool)
    row = np.ones((1, w), dtype=bool)
    shrunk = binary_erosion(binary_erosion(chunk, col), row)
    return binary_dilation(binary_dilation(shrunk, col), row)


def clean_shape_artifacts(path_A, path_B, rows, cols, cell_size_y, cell_size_x, antenna_spacing_km, min_width_km, tile_size):
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

    Returns:
    - tuple(int, int): Set-pixel counts after closing and after pruning.
    """
    close_r = max(1, int(antenna_spacing_km * 1000 / cell_size_y))
    close_c = max(1, int(antenna_spacing_km * 1000 / cell_size_x))
    tendril_r = max(1, int((min_width_km * 0.5 * 1000) / cell_size_y))
    tendril_c = max(1, int((min_width_km * 0.5 * 1000) / cell_size_x))
    n_closed = apply_morphology_pingpong(path_A, path_B, (rows, cols), bool, separable_closing, np.ones((close_r, close_c)), desc="Closing", tile_size=tile_size)
    n_pruned = apply_morphology_pingpong(path_B, path_A, (rows, cols), bool, separable_opening, np.ones((tendril_r, tendril_c)), desc="Pruning", tile_size=tile_size)
    return n_closed, n_pruned

def analyze_sites_and_capacity(path_A, elevation, rows, cols, cell_size_y, cell_size_x, downsample_factor, search_mode,
                               target_antennas, min_sub_array_size, antenna_spacing_km, grid_type, funnel=None,
                               candidates_arr=None, observables=None):
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
                    
                    site_details.append({
                        "site_id": int(site_id),
                        "area_km2": float(f"{area_km2:.2f}"),
                        "capacity_exact": int(antennas_fit),
                        "grid_type": grid_type,
                        "mean_aspect_deg": float(f"{mean_aspect:.1f}"),
                        "facing_direction": aspect_str
                    })

        site_details.sort(key=lambda x: x['capacity_exact'], reverse=True)
        final_selection_ids = []
        
        if search_mode == 'distributed':
            for site in site_details:
                final_selection_ids.append(site['site_id'])
                cumulative_capacity += site['capacity_exact']
        else:
            final_selection_ids = [s['site_id'] for s in site_details]

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
    """
    tfw_name = os.path.splitext(tif_filename)[0] + ".tfw"
    try:
        with open(tfw_name, "w") as f:
            # Format: Pixel X size, Rotation, Rotation, Negative Pixel Y size, Top-Left X, Top-Left Y
            f.write(f"{cell_size_deg:.10f}\n0.0\n0.0\n-{cell_size_deg:.10f}\n{top_left_lon:.10f}\n{top_left_lat:.10f}\n") 
    except: pass

def generate_kml_file(mask, elevation, filename, origin_lat, origin_lon, cell_size_deg, downsample=1):
    """
    Generates a Google Earth compatible KML file representing the valid site polygons.
    It extracts polygon contours from the binary mask using Matplotlib's contour tool.
    
    Parameters:
    - mask (ndarray): Binary mask indicating valid deployment sites.
    - filename (str): Output path for the KML file.
    - origin_lat, origin_lon, cell_size_deg: Used to convert array pixel indices to GPS coordinates.
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
            ET.SubElement(placemark, "name").text = f"GRAND Site {site_idx}"
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
    """
    generated_files = []
    cell_size_deg = map_grid.cell_size_deg
    base_filename = "grand_search_results_" + os.path.splitext(os.path.basename(dem_path))[0]
    
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
        
        im = ax.imshow(elev_viz, cmap='terrain', vmin=0, vmax=6000)
        
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
        
        if rfi_zones:
            deg_viz = cell_size_deg * viz_ds
            legend_handles.append(Line2D([0], [0], color='red', linestyle='--', lw=2))
            legend_labels.append("RFI exclusion zone")
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
                    ax.add_patch(Ellipse((px_x, px_y), w_px, h_px, edgecolor='red', facecolor='none', ls='--', lw=2))
                    text = ax.text(px_x, px_y-h_px/4, name, color='red', fontsize=12, ha='center')
                    text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'), path_effects.Normal()])
                elif type_tag == 'poly':
                    _, coords, name = item
                    verts = []
                    for (plat, plon) in coords:
                        px = (plon - origin_lon) / deg_viz
                        py = (origin_lat - plat) / deg_viz
                        verts.append((px, py))
                    ax.add_patch(MplPolygon(verts, closed=True, edgecolor='red', facecolor='none', ls='--', lw=2))
                    cx = sum(p[0] for p in verts)/len(verts)
                    cy = sum(p[1] for p in verts)/len(verts)
                    text = ax.text(cx, cy, name, color='red', fontsize=8, ha='center')
                    text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'), path_effects.Normal()])

        deg_viz = cell_size_deg * viz_ds
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p: f"{origin_lon + x*deg_viz:.2f}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y,p: f"{origin_lat - y*deg_viz:.2f}"))
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
        cbar = plt.colorbar(im, fraction=0.035, pad=0.04)
        cbar.set_label('Altitude (m)', rotation=270, labelpad=15)
        ax.set_title(f"GRAND site search | {region_name if region_name is not None else ''} {'|' if region_name is not None else ''} {search_mode.title()} mode\nFound {count} sites | Total capacity: {cumulative_capacity if search_mode=='distributed' else 'N/A'} DUs | Grid: {grid_type} | Spacing: {antenna_spacing_km} km | Altitude restriction: {min_altitude}-{max_altitude} m")
        
        fs = 'small' if len(legend_labels) > 8 else 'medium'
        ax.legend(legend_handles, legend_labels, loc='upper right', fontsize=fs, framealpha=0.8)
        
        img_name = os.path.join(run_output_dir, base_filename + "." + output_image_format.strip('.'))
        
        plt.savefig(img_name, format=output_image_format.strip('.'), dpi=150, bbox_inches='tight')
        generated_files.append(os.path.abspath(img_name))
        print(f"      {Icon.CHECK}Map saved.")
        
    except Exception as e:
        print(f"      {C.FAIL}{Icon.CROSS}Viz Error: {e}{C.RESET}")

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

    return generated_files

def print_tool_explanation():
    """Outputs a formatted explanation of the tool's capabilities and logic to the console."""
    print(f"""
{C.HEADER}================================================================================
{C.BOLD}GRAND NEUTRINO OBSERVATORY - AUTOMATED SITE SEARCH TOOL{C.RESET}{C.HEADER}
================================================================================{C.RESET}
This tool performs a high-performance topographic and physics simulation to 
identify suitable deployment sites for the GRAND array.

{C.BOLD}Core Workflow:{C.RESET}
1. Topographic Filtering: Scans the DEM for terrain with suitable slopes (configurable, default 3-25 degrees),
   enforcing altitude limits and specified facing directions (Aspect).
2. Logistics & RFI: Masks out areas overlapping populated centers and, optionally,
   areas situated too far from road infrastructure.
3. Ray-Tracing (Physics): Simulates line-of-sight from candidates to target 
   mountain ranges. It actively accounts for Earth's curvature and maintains a 
   clearance buffer over intermediate terrain.
4. Spatial Pruning: Implements morphological math (Closing/Opening) to remove 
   isolated "tendril" ridges that are unsuitable for wide array deployments.
5. Grid Packing: Simulates placing antennas in 'hex' or 'square' grids to calculate 
   the true physical capacity of the resulting sites.
   
{C.BOLD}Customizable Constraints & Processing Parameters:{C.RESET}
- Slope Bounds: Customizable minimum and maximum terrain steepness in degrees.
- RFI Zones: Accept pre-defined sets ('lima', 'arequipa') or custom geometry lists via JSON config.
- Map Resolution: Read from the DEM's own GeoTIFF tags, or forced via `--cell_size_deg`. Pixels
  are square in degrees but not in metres, so each axis carries its own ground scale.
- Candidate Stride: Thins the candidate set before ray-tracing via `--candidate_stride`.
- Downsample Factor: Modifies the internal resolution of the capacity masking, speeding up processing.
- Tile Size: Configures the size of memory-mapped square chunks for RAM management.
- Core Scaling: Control CPU thread allocation via `--num_cores` (defaults to all available cores).
- Checkpointing: The tool automatically saves progress. Use `--resume` to bypass ray-tracing on a failed run.
- Unified Output Generation: Exports georeferenced TIFFs, a KML file, an annotated graphical map, 
  a JSON run-summary, and the execution log into a single dynamically named output directory.
{C.HEADER}================================================================================{C.RESET}
    """)

def parse_score_weights(value):
    """
    Normalises per-component score weights from either input form.

    A config file is JSON, so it can carry a mapping directly. The command line
    cannot, so it takes ``shower=2,solid_angle=1`` instead. Both end up as a dict, and
    anything unnamed keeps weight 1.

    Returns None when nothing was supplied, which leaves the composition unweighted.
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


def validate_parameters(params):
    """
    Pre-flight validation checks to enforce 'Fail Fast' mechanisms. 
    Verifies the existence of critical files and the physical logic of search bounds 
    before engaging the memory-heavy processing loops.
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
#             MAIN EXECUTION ORCHESTRATOR
# ==========================================
def find_grand_regions_interactive(dem_path, cell_size_deg=None, target_antennas=1000,
                            rfi_zones=None, origin_lat=-15.0, origin_lon=-73.0,
                            min_width_km=2.0, min_altitude=None, max_altitude=None,
                            antenna_spacing_km=1.0, min_dist_km=30.0, max_dist_km=80.0,
                            road_map_path=None, max_road_dist_km=None,
                            grid_type='square', generate_kml=False,
                            search_mode='single', min_sub_array_size=100,
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
                            nu_interaction_length_gcm2=None):
    """
    The main orchestrator. Now decoupled from logic, it sets up the environment,
    calls the pipeline helpers in sequence, and manages memory cleanup and checkpointing.

    The map resolution (cell_size_deg) is read from the DEM's own georeferencing tags
    unless the caller overrides it; every metric conversion downstream derives from it.
    """

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
        max_range_m=max_dist_km * 1000.0,
        min_dist_km=min_dist_km, max_dist_km=max_dist_km,
        min_depth_gcm2=min_column_depth_gcm2, require_terrain=require_terrain,
        min_target_slope_deg=min_target_slope_deg,
        max_target_slope_deg=max_target_slope_deg,
        geomag_declination_deg=(geomag_declination_deg if use_geomagnetic else None),
        geomag_inclination_deg=(geomag_inclination_deg if use_geomagnetic else None),
        frequency_mhz=fresnel_frequency_mhz, bilinear=bilinear_sampling,
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
        "tile_size": tile_size, "resume": resume, "resume_dir": resume_dir,
        "num_cores": num_cores
    }
    
    print(f"\n{C.HEADER}============================================={C.RESET}")
    print(f"   {C.BOLD}GRAND SITE SEARCH: RUN PARAMETERS{C.RESET}")
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
                return
    
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
                                  "solid_angle_half_sr": solid_angle_half_sr,
                                  "clearance_full_at": clearance_full_at,
                                  "muon_shielding_km": muon_shielding_km},
                min_score=min_score, rfi_zones_px=rfi_zones_px)
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
        n_closed, n_pruned = clean_shape_artifacts(path_A, path_B, rows, cols, cell_size_y, cell_size_x, antenna_spacing_km, min_width_km, tile_size)
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
            candidates_arr=candidates_arr, observables=observables
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
        generated_files = generate_visualizations_and_outputs(
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
            print(f"   " + "-" * 50)
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

    finally:
        # Smart Cleanup
        if success_flag:
            try:
                if path_A and os.path.exists(path_A): os.remove(path_A)
                if path_B and os.path.exists(path_B): os.remove(path_B)
            except: 
                pass
        else:
            print(f"\n   {C.FAIL}{Icon.CROSS}[!] Run did not complete successfully.{C.RESET}")
            print(f"   {C.WARN}Buffer files have been retained in the workspace.{C.RESET}")
            print(f"   {C.WARN}Resume this exact run later using:{C.RESET} --resume --resume_dir {os.path.abspath(run_output_dir)}")
            
        timings["total"] = time.time() - t_start_total
        print(f"\n{C.OK}Total Execution Time: {timings['total']:.2f} seconds{C.RESET}")
        print(f"{C.BOLD}Done.{C.RESET}")

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GRAND Neutrino Array - Automated Site Search Tool")
    
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
    
    # IO / Configs mapping & Tools
    parser.add_argument("--config_path", type=str, default=None, help="Path to external JSON configuration file.")
    parser.add_argument("--output_directory_base_with_given_json", type=str, default="../output/", help="Base directory for outputs when a JSON config is supplied (default: ../output/).")
    parser.add_argument("--output_image_format", type=str, default="png", help="Format of the saved map visual, e.g., png, pdf, svg (default: png).")
    parser.add_argument("--resume", action="store_true", help="Include this flag to resume a previous run from the ray-tracing checkpoint.")
    parser.add_argument("--resume_dir", type=str, default=None, help="Path to an output folder from a previously failed run to resume from the ray-tracing checkpoint.")
    
    # Tool Generation Arguments
    parser.add_argument("--generate_config", type=str, default=None, help="Supply a filepath to generate a default JSON config template and exit.")
    parser.add_argument("--config_preset", type=str, choices=['default', 'lima', 'arequipa'], default='default', help="Optional presets to inject when using --generate_config.")

    args = parser.parse_args()

    # --- Tool Execution: Generate Config Template ---
    if args.generate_config:
        preset = args.config_preset
        default_config = {
            "dem_path": "path_to_your_dem.tif",
            "origin_lat": 0.0,
            "origin_lon": 0.0,
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
            "generate_kml": True,
            "print_info": True,
            "output_directory_base_with_given_json": "../output/",
            "output_image_format": "png",
            "resume": False,
            "resume_dir": None
        }
        
        # Inject presets if specifically requested
        if preset == 'lima':
            default_config['origin_lat'] = ORIGIN_LAT_LIMA
            default_config['origin_lon'] = ORIGIN_LON_LIMA
            default_config['rfi_zones'] = 'lima'
            default_config['region_name'] = 'Lima, Peru'
            default_config['dem_path'] = 'lima_AW3D30.tif'
        elif preset == 'arequipa':
            default_config['origin_lat'] = ORIGIN_LAT_AREQUIPA
            default_config['origin_lon'] = ORIGIN_LON_AREQUIPA
            default_config['rfi_zones'] = 'arequipa'
            default_config['region_name'] = 'Arequipa, Peru'
            default_config['dem_path'] = 'arequipa_SRTMGL1.tif'

        # Safely create directory structure if necessary and save
        dir_name = os.path.dirname(os.path.abspath(args.generate_config))
        if dir_name:
            os.makedirs(dir_name, exist_ok=True)
            
        with open(args.generate_config, 'w') as f:
            json.dump(default_config, f, indent=4)
        print(f"Configuration file generated successfully at: {args.generate_config} (Preset: {preset})")
        sys.exit(0)


    # 1. Initialize Configuration Maps
    config_params = {}
    if args.config_path and os.path.exists(args.config_path):
        with open(args.config_path, 'r') as f:
            config_params = json.load(f)

    # 2. Retrieve Fallbacks
    fallback_path = os.path.join("..", "config", "fallbacks.json")
    fallback_params = {}
    if os.path.exists(fallback_path):
        with open(fallback_path, 'r') as f:
            fallback_params = json.load(f)

    # 3. Determine Unified Logging and Output Directory Hierarchically
    base_dir = args.output_directory_base_with_given_json
    if "output_directory_base_with_given_json" in config_params:
        base_dir = config_params["output_directory_base_with_given_json"]
    elif "output_directory_base_with_given_json" in fallback_params:
        base_dir = fallback_params["output_directory_base_with_given_json"]

    if args.config_path and os.path.exists(args.config_path):
        config_basename = os.path.splitext(os.path.basename(args.config_path))[0]
        run_output_dir = os.path.join(base_dir, config_basename)
    else:
        timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_output_dir = os.path.join("..", "output", timestamp_str)

    os.makedirs(run_output_dir, exist_ok=True)
    
    # 4. Apply Custom Standard-Out / Standard-Error interceptors for the log file
    log_path = os.path.join(run_output_dir, "log.txt")
    log_file = open(log_path, "a", encoding="utf-8")
    
    sys.stdout = TeeLogger(sys.stdout, log_file)
    sys.stderr = TeeLogger(sys.stderr, log_file)

    # Ensure log captures initiation context
    print(f"\n{C.HEADER}================================================================================{C.RESET}")
    print(f"Execution started at: {datetime.now().isoformat()}")
    if args.config_path:
        print(f"Using config file: {os.path.abspath(args.config_path)}")
    else:
        print(f"No config file provided. Relying on CLI arguments and fallbacks.")
    print(f"Using fallbacks file: {os.path.abspath(fallback_path)}")
    print(f"Unified output directory initialized at: {os.path.abspath(run_output_dir)}")
    print(f"{C.HEADER}================================================================================{C.RESET}\n")

    # 5. Reconcile Configuration Strategy (Config > Fallback > CLI / Standard defaults)
    final_params = {}
    
    # Collect all available arguments parsed from CLI framework
    param_names = [action.dest for action in parser._actions if action.dest not in ('help', 'config_path', 'generate_config', 'config_preset')]
    
    for param in param_names:
        if param in config_params:
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
        print(f"{C.FAIL}{Icon.CROSS}ERROR: Critical parameters 'origin_lat' and 'origin_lon' must be provided via config file, fallback, or CLI.{C.RESET}")
        sys.exit(1)

    # Default resume_dir to run_output_dir if resume is True and resume_dir not provided
    if final_params.get('resume') and not final_params.get('resume_dir'):
        final_params['resume_dir'] = run_output_dir

    # 6. Run Pre-Flight Validation (Fail-Fast Mechanism)
    validate_parameters(final_params)

    # Handle RFI Zone selection mapping (Checks Config-passed custom lists, or matches string presets)
    rfi_input = final_params.get('rfi_zones', 'none')
    selected_rfi = None
    
    if isinstance(rfi_input, str):
        if rfi_input.lower() == 'lima':
            selected_rfi = LIMA_RFI_ZONES
        elif rfi_input.lower() == 'arequipa':
            selected_rfi = AREQUIPA_RFI_ZONES
        elif rfi_input.lower() != 'none':
            # Attempt to parse as JSON if a raw string array was passed via CLI
            try:
                selected_rfi = json.loads(rfi_input)
            except Exception as e:
                print(f"{C.WARN}{Icon.WARN}WARNING: Could not parse custom rfi_zones string. Proceeding with 'none'. Error: {e}{C.RESET}")
    elif isinstance(rfi_input, list):
        # Naturally supports custom RFI arrays loaded cleanly from the JSON config file
        selected_rfi = rfi_input

    # Output explanations if requested (either via flag or configured true)
    if final_params.get('print_info', True):
        print_tool_explanation()

    # Execute main search pipeline with our integrated parameters
    find_grand_regions_interactive(
        dem_path=final_params['dem_path'],
        target_antennas=final_params['target_antennas'], 
        rfi_zones=selected_rfi,
        min_width_km=final_params['min_width_km'],
        origin_lat=final_params['origin_lat'],
        origin_lon=final_params['origin_lon'],
        min_altitude=final_params['min_altitude'], 
        max_altitude=final_params['max_altitude'],
        antenna_spacing_km=final_params['antenna_spacing_km'],
        min_dist_km=final_params['min_dist_km'],
        max_dist_km=final_params['max_dist_km'],
        grid_type=final_params['grid_type'],       
        generate_kml=final_params['generate_kml'],     
        road_map_path=final_params['road_map_path'],    
        max_road_dist_km=final_params['max_road_dist_km'],
        search_mode=final_params['search_mode'],
        min_sub_array_size=final_params['min_sub_array_size'],
        min_aspect_deg=final_params['min_aspect_deg'], 
        max_aspect_deg=final_params['max_aspect_deg'],
        min_slope_deg=final_params['min_slope_deg'],
        max_slope_deg=final_params['max_slope_deg'],
        region_name=final_params['region_name'],
        downsample_factor=final_params['downsample_factor'],
        run_output_dir=run_output_dir,
        output_image_format=final_params['output_image_format'],
        tile_size=final_params['tile_size'],
        resume=final_params.get('resume', False),
        resume_dir=final_params.get('resume_dir'),
        num_cores=final_params.get('num_cores', -1),
        cell_size_deg=final_params.get('cell_size_deg'),
        candidate_stride=final_params.get('candidate_stride', 5),
        slope_baseline_m=final_params.get('slope_baseline_m'),
        energy_min_pev=final_params.get('energy_min_pev'),
        energy_max_pev=final_params.get('energy_max_pev'),
        n_azimuths=final_params.get('n_azimuths', 9),
        azimuth_half_width_deg=(None if (final_params.get('azimuth_half_width_deg', 60.0) or 0) < 0
                                else final_params.get('azimuth_half_width_deg', 60.0)),
        elev_min_deg=final_params.get('elev_min_deg', -3.0),
        elev_max_deg=final_params.get('elev_max_deg', 3.0),
        n_elev_bins=final_params.get('n_elev_bins', 12),
        min_column_depth_gcm2=final_params.get('min_column_depth_gcm2', 0.0),
        require_terrain=not final_params.get('require_sky', False),
        fresnel_frequency_mhz=final_params.get('fresnel_frequency_mhz'),
        refraction_k=final_params.get('refraction_k'),
        antenna_height_m=final_params.get('antenna_height_m', 2.0),
        exclude_near_field=final_params.get('exclude_near_field', True),
        fresnel_near_field_m=final_params.get('fresnel_near_field_m', 500.0),
        depth_band_gcm2=(tuple(final_params['depth_band_gcm2'])
                         if final_params.get('depth_band_gcm2') else None),
        score_composition=final_params.get('score_composition', 'product'),
        min_score=final_params.get('min_score', 0.0),
        geomag_declination_deg=final_params.get('geomag_declination_deg'),
        geomag_inclination_deg=final_params.get('geomag_inclination_deg'),
        muon_shielding_km=final_params.get('muon_shielding_km'),
        bilinear_sampling=final_params.get('bilinear_sampling', True),
        use_geomagnetic=final_params.get('use_geomagnetic', True),
        grammage_mode=final_params.get('grammage_mode', 'radio'),
        grammage_band_gcm2=(tuple(final_params['grammage_band_gcm2'])
                            if final_params.get('grammage_band_gcm2') else None),
        grammage_maturity_gcm2=final_params.get('grammage_maturity_gcm2'),
        min_target_slope_deg=final_params.get('min_target_slope_deg'),
        max_target_slope_deg=final_params.get('max_target_slope_deg'),
        grammage_band_fraction=final_params.get('grammage_band_fraction'),
        shower_elongation_rate_gcm2=final_params.get('shower_elongation_rate_gcm2'),
        shower_lambda_gcm2=final_params.get('shower_lambda_gcm2'),
        solid_angle_half_sr=final_params.get('solid_angle_half_sr'),
        distance_band_m=(tuple(final_params['distance_band_m'])
                         if final_params.get('distance_band_m') else None),
        clearance_full_at=final_params.get('clearance_full_at'),
        score_weights=parse_score_weights(final_params.get('score_weights')),
        nu_interaction_length_gcm2=final_params.get('nu_interaction_length_gcm2')
    )