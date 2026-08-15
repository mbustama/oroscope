"""
Synthetic terrain with analytically known properties.

Every generator here produces terrain whose slope, aspect, target distance or canyon
geometry is known in closed form, so tests can assert against arithmetic rather than
against a previous run. The fixtures themselves are verified in
``test_fixtures.py`` — a fixture that is wrong verifies nothing.

Shared by the test suite and the benchmark harness.
"""

import numpy as np
import tifffile as tiff

# 1 arc-second, matching SRTMGL1 / AW3D30
CELL_DEG = 1.0 / 3600.0


def cell_sizes(latitude_deg):
    """Metric pixel sizes at a latitude, matching the searcher's own convention."""
    cell_y = CELL_DEG * 110.6 * 1000.0
    cell_x = CELL_DEG * 111.32 * np.cos(np.radians(latitude_deg)) * 1000.0
    return cell_y, cell_x


def planar(n, slope_deg, aspect_deg, cell_y, cell_x, base=2000.0):
    """
    A perfect inclined plane with exactly the requested slope and aspect.

    Aspect is the downhill direction, measured clockwise from north, matching the
    searcher's ``degrees(arctan2(-dx, dy)) % 360``. Deriving elevation as

        z(r, c) = tan(S) * (r*cell_y*cos(A) - c*cell_x*sin(A))

    makes dz/d(east) = -tan(S)sin(A) and dz/d(south) = tan(S)cos(A), which recovers
    slope S and aspect A exactly.
    """
    rr, cc = np.mgrid[0:n, 0:n].astype(np.float64)
    t = np.tan(np.radians(slope_deg))
    a = np.radians(aspect_deg)
    z = t * (rr * cell_y * np.cos(a) - cc * cell_x * np.sin(a))
    return (z - z.min() + base).astype(np.float32)


def flat_with_peak(n, peak_r, peak_c, height, half_width=3, base=0.0):
    """A flat plain carrying one rectangular block, for ray-targeting tests."""
    z = np.full((n, n), base, dtype=np.float32)
    z[peak_r - half_width:peak_r + half_width + 1,
      peak_c - half_width:peak_c + half_width + 1] = base + height
    return z


def ridge_and_slope(n, cell_x, ridge_col=120, ridge_height=3500.0, ridge_sigma_px=60.0,
                    valley_col=220, rise_per_px=2.0, base=500.0):
    """
    A tall ridge in the west with terrain rising eastward away from it.

    Slopes east of the valley face west (downhill toward the ridge), so candidates
    there look across the valley at the ridge — the GRAND viewing geometry. Fully
    deterministic: no noise, so results are reproducible bit-for-bit.
    """
    cols = np.arange(n, dtype=np.float64)[None, :].repeat(n, 0)
    z = base + np.clip(cols - valley_col, 0, None) * rise_per_px
    z = z + ridge_height * np.exp(-((cols - ridge_col) ** 2) / (2 * ridge_sigma_px ** 2))
    return z.astype(np.float32)


def canyon(n, cell_x, floor_width_m=1000.0, depth_m=1500.0, wall_slope_deg=35.0,
           rim_elevation=3500.0):
    """
    A north-south canyon with two opposing walls of known geometry.

    Rim-to-rim separation is ``floor_width_m + 2*depth_m/tan(wall_slope)``, which is
    what a TAMBO-style search must recover. Colca's published figures — ~1.5 km deep,
    ~4.5 km median between valley sides — correspond to a ~35 degree wall.
    """
    cols = np.arange(n, dtype=np.float64)[None, :].repeat(n, 0)
    x = cols * cell_x
    center = (n * cell_x) / 2.0
    from_edge = np.abs(x - center) - floor_width_m / 2.0
    rise = np.clip(from_edge, 0.0, None) * np.tan(np.radians(wall_slope_deg))
    z = rim_elevation - depth_m + np.clip(rise, 0.0, depth_m)
    return z.astype(np.float32)


def canyon_rim_separation_m(floor_width_m, depth_m, wall_slope_deg):
    """Closed-form rim-to-rim separation of :func:`canyon`."""
    return floor_width_m + 2.0 * depth_m / np.tan(np.radians(wall_slope_deg))


# Colca Canyon as published in ref. [2]: ~1.5 km deep, ~4.5 km median between valley
# sides. The wall slope follows from those two, given a ~1 km floor:
#     4500 = 1000 + 2*1500/tan(s)  ->  s = 40.6 deg
# Note the walls are far steeper than GRAND's 3-25 deg deployable band, which is one
# reason the slope criterion has to be per-experiment.
COLCA = dict(floor_width_m=1000.0, depth_m=1500.0, wall_slope_deg=40.6, rim_elevation=3500.0)


def colca_like(n, cell_x):
    """A canyon with Colca's published depth and rim-to-rim separation."""
    return canyon(n, cell_x, **COLCA)


def write_geotiff(path, array, origin_lat, origin_lon, cell_size_deg=CELL_DEG):
    """Writes a geographic GeoTIFF carrying the tags the searcher reads."""
    tiff.imwrite(
        path, array,
        extratags=[
            (33550, "d", 3, (cell_size_deg, cell_size_deg, 0.0)),           # ModelPixelScale
            (33922, "d", 6, (0.0, 0.0, 0.0, origin_lon, origin_lat, 0.0)),  # ModelTiepoint
        ],
    )
    return path
