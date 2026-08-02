"""
Arrival-direction scanning: what a candidate site can actually see, and what lies behind it.

The question a site search must answer is not "is there a tall mountain out there" but
"from which arrival directions does a backward ray from this pixel enter rock, and how
much rock does it cross". This module answers that directly.

Tracing backward from a candidate along an arrival direction (azimuth phi, elevation
angle theta measured from horizontal):

  * rays above the local horizon escape to the sky and see no matter;
  * rays below it strike terrain. That first intersection is the tau exit point, the
    distance to it is the decay baseline, and the path length beyond it that runs
    under the surface is the column depth.

Different experiments accept different arrival directions, so the elevation window is
a parameter rather than a constant. GRAND accepts neutrinos within roughly +/-3 deg of
the horizon and cosmic rays from above it; TAMBO looks across a canyon. The cosmic-ray
case inverts the test -- terrain in the accepted directions is an obstruction rather
than a target -- which ``require_terrain=False`` expresses without a second code path.

Algorithm
---------
For a fixed azimuth, one walk outward yields every elevation bin at once. Writing the
elevation angle of the terrain at ground distance d as

    theta_terrain(d) = atan( (z(d) - d^2/2R - z0) / d )

then a ray at angle theta first meets terrain at the smallest d where
theta_terrain(d) >= theta. Since the running maximum of theta_terrain only increases,
each new maximum claims a contiguous band of elevation bins, so first-intersection
distances for all bins are filled in a single pass.

Column depth follows from the same samples: the ray at angle theta is underground
wherever theta_terrain(d) > theta, so binning theta_terrain and taking an inclusive
suffix sum gives the underground path length for every bin at once. Rays crossing
several ridges accumulate all of the rock they traverse, not just the first chord.

This costs one profile walk per (candidate, azimuth) regardless of how finely the
elevation window is sampled.
"""

import numpy as np

try:
    from numba import jit, prange
    HAS_NUMBA = True
except ImportError:                                     # pragma: no cover
    HAS_NUMBA = False
    prange = range

    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


# Standard rock, the usual convention for tau propagation through the crust
STANDARD_ROCK_DENSITY = 2650.0          # kg/m^3
# Effective radius for radio line-of-sight (the 4/3-Earth refraction convention)
DEFAULT_EARTH_RADIUS_M = 8.5e6
# Column depth is reported in g/cm^2: 1 kg/m^2 = 0.1 g/cm^2
KGM2_TO_GCM2 = 0.1


@jit(nopython=True, nogil=True, fastmath=True)
def _scan_one_direction(elevation, r0, c0, z0, azimuth_deg,
                        cell_size_y, cell_size_x, rows, cols,
                        elev_min_deg, elev_bin_deg, n_bins,
                        step_m, max_range_m, inv_2R,
                        first_dist, angle_hist):
    """
    Walks one azimuth and fills, for every elevation bin:

      * ``first_dist`` -- ground distance to the first terrain intersection, or -1
        where the ray escapes to the sky;
      * ``angle_hist`` -- path length binned by the terrain's elevation angle, which
        an inclusive suffix sum turns into underground path length per bin.

    Both arrays are caller-supplied and overwritten. Returns the horizon angle.
    """
    for b in range(n_bins):
        first_dist[b] = -1.0
        angle_hist[b] = 0.0
    angle_hist[n_bins] = 0.0            # overflow bin: terrain above the window

    look = np.radians(azimuth_deg)
    sin_look = np.sin(look)
    cos_look = np.cos(look)

    running_max = -1.0e30
    fill_bin = 0                        # lowest bin not yet assigned a distance
    horizon = -1.0e30

    d = step_m
    while d <= max_range_m:
        # Ground displacement to pixel offsets, one scale per axis
        tc = c0 + int(d * sin_look / cell_size_x)
        tr = r0 - int(d * cos_look / cell_size_y)
        if tr < 0 or tr >= rows or tc < 0 or tc >= cols:
            break

        z = elevation[tr, tc]
        if not np.isnan(z):
            # Earth curvature lowers distant terrain: apparent height drops as d^2/2R
            apparent = z - (d * d) * inv_2R - z0
            theta = np.degrees(np.arctan(apparent / d))

            if theta > horizon:
                horizon = theta

            # Path length attributed to this terrain angle, for the suffix sum
            k = int((theta - elev_min_deg) / elev_bin_deg)
            if k >= n_bins:
                k = n_bins              # above the window: counts for every bin
            if k >= 0:
                angle_hist[k] += step_m

            # Each new running maximum claims the elevation bins it has just risen past
            if theta > running_max:
                running_max = theta
                while fill_bin < n_bins and (elev_min_deg + fill_bin * elev_bin_deg) <= theta:
                    first_dist[fill_bin] = d
                    fill_bin += 1
        d += step_m

    return horizon


@jit(nopython=True, nogil=True, parallel=True)
def scan_candidates(candidates, elevation, cell_size_y, cell_size_x, rows, cols,
                    azimuth_offsets_deg, use_aspect,
                    elev_min_deg, elev_max_deg, n_bins,
                    step_m, max_range_m,
                    min_dist_m, max_dist_m,
                    min_depth_gcm2, require_terrain,
                    rock_density, earth_radius_m,
                    out_cells, out_solid_angle, out_mean_dist,
                    out_max_depth, out_mean_depth, out_horizon):
    """
    Scans arrival directions for every candidate and reports what each one sees.

    Parameters
    ----------
    candidates : (N, 3) float array of [row, col, aspect_deg].
    azimuth_offsets_deg : (A,) azimuths. Offsets from each candidate's aspect when
        ``use_aspect`` is true, otherwise absolute compass bearings.
    elev_min_deg, elev_max_deg, n_bins : the accepted elevation window and its sampling.
    step_m, max_range_m : profile sampling step and how far to look.
    min_dist_m, max_dist_m : accepted range to the first intersection, i.e. the decay
        baseline window.
    min_depth_gcm2 : column depth a direction must have to count.
    require_terrain : true selects directions that strike rock (neutrino channels);
        false selects directions that escape to the sky (cosmic-ray channels), where
        terrain is an obstruction and the depth and distance criteria do not apply.
    rock_density : kg/m^3 for converting path length to column depth.

    Results are written into the ``out_*`` arrays, one entry per candidate:
    accepted-direction count, accepted solid angle (sr), mean distance to the exit
    point (m), maximum and mean column depth (g/cm^2), and the horizon angle (deg).
    """
    n = candidates.shape[0]
    n_az = azimuth_offsets_deg.shape[0]
    elev_bin_deg = (elev_max_deg - elev_min_deg) / n_bins
    inv_2R = 1.0 / (2.0 * earth_radius_m)

    # Solid angle of one (azimuth, elevation) cell: dOmega = cos(theta) dtheta dphi
    d_phi = 2.0 * np.pi / n_az if n_az > 0 else 0.0
    d_theta = np.radians(elev_bin_deg)
    depth_scale = rock_density * KGM2_TO_GCM2

    for i in prange(n):
        r0 = int(candidates[i, 0])
        c0 = int(candidates[i, 1])
        aspect = candidates[i, 2]
        z0 = elevation[r0, c0]

        first_dist = np.empty(n_bins, dtype=np.float64)
        angle_hist = np.empty(n_bins + 1, dtype=np.float64)

        cells = 0
        solid_angle = 0.0
        dist_sum = 0.0
        depth_sum = 0.0
        depth_max = 0.0
        horizon_max = -1.0e30

        if np.isnan(z0):
            out_cells[i] = 0
            out_solid_angle[i] = 0.0
            out_mean_dist[i] = 0.0
            out_max_depth[i] = 0.0
            out_mean_depth[i] = 0.0
            out_horizon[i] = 0.0
            continue

        for a in range(n_az):
            azimuth = azimuth_offsets_deg[a]
            if use_aspect:
                azimuth = aspect + azimuth

            horizon = _scan_one_direction(
                elevation, r0, c0, z0, azimuth,
                cell_size_y, cell_size_x, rows, cols,
                elev_min_deg, elev_bin_deg, n_bins,
                step_m, max_range_m, inv_2R,
                first_dist, angle_hist)
            if horizon > horizon_max:
                horizon_max = horizon

            # Inclusive suffix sum: underground path length for each elevation bin
            running = angle_hist[n_bins]
            for b in range(n_bins - 1, -1, -1):
                running += angle_hist[b]
                theta = elev_min_deg + (b + 0.5) * elev_bin_deg
                cos_theta = np.cos(np.radians(theta))

                if require_terrain:
                    d_hit = first_dist[b]
                    if d_hit < 0.0 or d_hit < min_dist_m or d_hit > max_dist_m:
                        continue
                    # Slant path through rock, beyond the exit point
                    depth = (running / cos_theta) * depth_scale
                    if depth < min_depth_gcm2:
                        continue
                    cells += 1
                    solid_angle += cos_theta * d_theta * d_phi
                    dist_sum += d_hit
                    depth_sum += depth
                    if depth > depth_max:
                        depth_max = depth
                else:
                    # Cosmic-ray style: the direction counts only if nothing blocks it
                    if first_dist[b] >= 0.0:
                        continue
                    cells += 1
                    solid_angle += cos_theta * d_theta * d_phi

        out_cells[i] = cells
        out_solid_angle[i] = solid_angle
        out_mean_dist[i] = dist_sum / cells if cells > 0 else 0.0
        out_max_depth[i] = depth_max
        out_mean_depth[i] = depth_sum / cells if cells > 0 else 0.0
        out_horizon[i] = horizon_max if horizon_max > -1.0e29 else 0.0


def azimuth_fan(n_azimuths, half_width_deg=None):
    """
    Azimuths to scan, as offsets from each candidate's aspect.

    ``half_width_deg`` restricts the fan to a forward arc, which is what a slope-mounted
    array sees; leave it None for a full 360 degree sweep, appropriate when the array
    orientation does not constrain the acceptance.
    """
    if half_width_deg is None:
        return np.linspace(0.0, 360.0, n_azimuths, endpoint=False)
    if n_azimuths == 1:
        return np.zeros(1)
    return np.linspace(-half_width_deg, half_width_deg, n_azimuths)


def scan(candidates, elevation, map_grid, *,
         elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
         n_azimuths=9, half_width_deg=60.0, use_aspect=True,
         step_m=None, max_range_m=80000.0,
         min_dist_km=0.0, max_dist_km=80.0,
         min_depth_gcm2=0.0, require_terrain=True,
         rock_density=STANDARD_ROCK_DENSITY,
         earth_radius_m=DEFAULT_EARTH_RADIUS_M):
    """
    Convenience wrapper over :func:`scan_candidates` with defaults for GRAND neutrinos.

    Sampling defaults to one DEM pixel, since a coarser step can miss a ridge entirely.

    Returns:
    - dict of per-candidate arrays: cells, solid_angle_sr, mean_distance_m,
      max_depth_gcm2, mean_depth_gcm2, horizon_deg.
    """
    candidates = np.ascontiguousarray(candidates, dtype=np.float64)
    n = candidates.shape[0]
    rows, cols = elevation.shape
    if step_m is None:
        step_m = min(map_grid.cell_size_y, map_grid.cell_size_x)

    offsets = azimuth_fan(n_azimuths, half_width_deg)

    out = {
        "cells": np.zeros(n, dtype=np.int64),
        "solid_angle_sr": np.zeros(n, dtype=np.float64),
        "mean_distance_m": np.zeros(n, dtype=np.float64),
        "max_depth_gcm2": np.zeros(n, dtype=np.float64),
        "mean_depth_gcm2": np.zeros(n, dtype=np.float64),
        "horizon_deg": np.zeros(n, dtype=np.float64),
    }
    if n == 0:
        return out

    scan_candidates(
        candidates, elevation, map_grid.cell_size_y, map_grid.cell_size_x, rows, cols,
        offsets, use_aspect,
        elev_min_deg, elev_max_deg, n_elev_bins,
        float(step_m), float(max_range_m),
        min_dist_km * 1000.0, max_dist_km * 1000.0,
        float(min_depth_gcm2), bool(require_terrain),
        float(rock_density), float(earth_radius_m),
        out["cells"], out["solid_angle_sr"], out["mean_distance_m"],
        out["max_depth_gcm2"], out["mean_depth_gcm2"], out["horizon_deg"],
    )
    return out
