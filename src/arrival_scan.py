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

import math

import numpy as np

import physics

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
# Tau lepton, for the boosted decay length
TAU_MASS_GEV = 1.77686
TAU_CTAU_M = 87.03e-6
# Two radii, because two different things propagate.
#
# Particles travel in straight lines: the neutrino and the tau are not refracted, so
# the geometry deciding where the tau exits uses the true Earth radius. That is not a
# modelling choice, it is what the trajectory is.
TRUE_EARTH_RADIUS_M = 6.371e6
# The radio signal *is* refracted by the tropospheric density gradient, and the 4/3
# convention makes a refracted ray straight again. It applies to the Fresnel clearance
# of the signal path, and to nothing else.
RADIO_EARTH_RADIUS_M = 8.5e6
DEFAULT_EARTH_RADIUS_M = TRUE_EARTH_RADIUS_M      # retained for older callers
# Column depth is reported in g/cm^2: 1 kg/m^2 = 0.1 g/cm^2
KGM2_TO_GCM2 = 0.1
SPEED_OF_LIGHT = 2.99792458e8            # m/s, for the Fresnel wavelength


def earth_radius_for_k(k_factor):
    """
    Effective Earth radius for a refraction k-factor.

    k = 1 is true geometry; k = 4/3 is the standard radio convention and gives the
    8500 km the searcher has always used. The choice is not negligible: over an 80 km
    path the apparent drop is 376 m at k = 4/3 against 502 m at k = 1, a difference
    comparable to the Fresnel clearance itself.
    """
    return 6371000.0 * float(k_factor)


@jit(nopython=True, nogil=True, fastmath=True)
def _scan_one_direction(elevation, r0, c0, z0, azimuth_deg,
                        cell_size_y, cell_size_x, rows, cols,
                        tan_edges, n_bins,
                        step_m, max_range_m, inv_2R, bilinear,
                        first_dist, angle_hist):
    """
    Walks one azimuth and fills, for every elevation bin:

      * ``first_dist`` -- ground distance to the first terrain intersection, or -1
        where the ray escapes to the sky;
      * ``angle_hist`` -- path length binned by the terrain's elevation angle, which
        an inclusive suffix sum turns into underground path length per bin.

    Works in **slope** rather than angle. Every comparison the walk makes is monotonic
    in the elevation angle, so comparing ``apparent/d`` against pre-computed tangents of
    the bin edges gives identical results without an arctangent per sample -- which at
    roughly 15 ns per sample was about half the cost of the whole scan. The per-axis
    pixel steps are hoisted out of the loop for the same reason.

    ``tan_edges`` holds tan of the n_bins+1 bin edges. Both output arrays are
    caller-supplied and overwritten. Returns the horizon angle in degrees.
    """
    for b in range(n_bins):
        first_dist[b] = -1.0
        angle_hist[b] = 0.0
    angle_hist[n_bins] = 0.0            # overflow bin: terrain above the window

    look = np.radians(azimuth_deg)
    # Pixels per metre of ground distance, one scale per axis, hoisted out of the loop
    dc_per_m = np.sin(look) / cell_size_x
    dr_per_m = np.cos(look) / cell_size_y

    running_max = -1.0e30               # in slope, not angle
    fill_bin = 0                        # lowest bin not yet assigned a distance
    horizon_slope = -1.0e30

    d = step_m
    while d <= max_range_m:
        tc = c0 + int(d * dc_per_m)
        tr = r0 - int(d * dr_per_m)
        if tr < 0 or tr >= rows or tc < 0 or tc >= cols:
            break
        # Kept as a per-sample test on purpose. Solving for the exit distance up front
        # removes four comparisons but has to reproduce int()'s truncation exactly at
        # every edge; getting it wrong silently truncates rays, and an attempt at it
        # altered 405 of 40,000 candidates for about 5%, which the branch predictor
        # handles nearly as cheaply anyway.

        if bilinear:
            # Sub-pixel sampling. Nearest-neighbour quantises the profile to pixel
            # centres, and int() truncates toward zero, which biases the sampled point
            # back toward the candidate by up to half a pixel -- asymmetrically, since
            # the sign of the offset depends on the azimuth.
            fc = c0 + d * dc_per_m
            fr = r0 - d * dr_per_m
            ci = int(fc) if fc >= 0.0 else int(fc) - 1
            ri = int(fr) if fr >= 0.0 else int(fr) - 1
            if ri < 0 or ri + 1 >= rows or ci < 0 or ci + 1 >= cols:
                z = elevation[tr, tc]
            else:
                u = fc - ci
                v = fr - ri
                z00 = elevation[ri, ci]
                z01 = elevation[ri, ci + 1]
                z10 = elevation[ri + 1, ci]
                z11 = elevation[ri + 1, ci + 1]
                if np.isnan(z00) or np.isnan(z01) or np.isnan(z10) or np.isnan(z11):
                    z = elevation[tr, tc]       # a nodata neighbour: fall back
                else:
                    z = ((1.0 - v) * ((1.0 - u) * z00 + u * z01)
                         + v * ((1.0 - u) * z10 + u * z11))
        else:
            z = elevation[tr, tc]

        if not np.isnan(z):
            # Earth curvature lowers distant terrain: apparent height drops as d^2/2R
            slope = (z - (d * d) * inv_2R - z0) / d

            if slope > horizon_slope:
                horizon_slope = slope

            # Path length attributed to this terrain angle, for the suffix sum. The
            # sample counts for every bin whose lower edge lies below it.
            if slope >= tan_edges[0]:
                k = n_bins                      # above the window: counts for every bin
                for b in range(1, n_bins + 1):
                    if tan_edges[b] > slope:
                        k = b - 1
                        break
                angle_hist[k] += step_m

            # Each new running maximum claims the elevation bins it has just risen past
            if slope > running_max:
                running_max = slope
                while fill_bin < n_bins and tan_edges[fill_bin] <= slope:
                    first_dist[fill_bin] = d
                    fill_bin += 1
        d += step_m

    if horizon_slope < -1.0e29:
        return -1.0e30
    return np.degrees(np.arctan(horizon_slope))


@jit(nopython=True, nogil=True, fastmath=True)
def _min_clearance_ratio(elevation, r0, c0, z0, azimuth_deg,
                         cell_size_y, cell_size_x, rows, cols,
                         theta_deg, d_hit, step_m, inv_2R, wavelength_m, shower_offset_m,
                         antenna_height_m, near_field_m, radio_inv_2R):
    """
    Worst Fresnel clearance along the path to an intersection, in units of r1.

    The scan already guarantees nothing blocks the line of sight -- the intersection is
    by construction the first terrain met -- so this is a refinement, not a gate: a ray
    that merely grazes a ridge on the way suffers diffraction loss even though the
    geometric path is clear. The first Fresnel radius at distance d along a path of
    length D is

        r1 = sqrt(lambda * d * (D - d) / D)

    which is why this needs the hit distance and so runs as a second pass, over
    accepted directions only.

    The far endpoint is the *shower*, not the exit point. The radio source is the air
    shower developing over some kilometres after the tau decays, and taking the exit
    point instead makes the measure degenerate: approaching the target both the
    clearance and r1 go to zero, so their ratio collapses for every path regardless of
    whether anything actually obstructs it.

    The antenna sits ``antenna_height_m`` above the ground. Without that the measure is
    meaningless: a receiver at ground level always has terrain inside the first Fresnel
    zone immediately beside it, so every path scores near zero whatever the terrain
    beyond does.

    ``near_field_m`` skips the first stretch of the path. The criterion is meant to
    catch an intervening ridge, but within a few hundred metres the first Fresnel zone
    is narrow and the ground the antenna stands on fills it, so including that stretch
    measures roughness at a scale a 30 m DEM cannot resolve.

    Returns the minimum of clearance/r1, or a large value when there is no path to
    measure.
    """
    d_end = d_hit - shower_offset_m
    if d_end <= 2.0 * step_m or wavelength_m <= 0.0:
        return 1.0e30

    look = np.radians(azimuth_deg)
    sin_look = np.sin(look)
    cos_look = np.cos(look)
    tan_theta = np.tan(np.radians(theta_deg))

    worst = 1.0e30
    d = max(step_m, near_field_m)
    while d <= d_end:
        tc = c0 + int(d * sin_look / cell_size_x)
        tr = r0 - int(d * cos_look / cell_size_y)
        if tr < 0 or tr >= rows or tc < 0 or tc >= cols:
            break
        # Kept as a per-sample test on purpose. Solving for the exit distance up front
        # removes four comparisons but has to reproduce int()'s truncation exactly at
        # every edge; getting it wrong silently truncates rays, and an attempt at it
        # altered 405 of 40,000 candidates for about 5%, which the branch predictor
        # handles nearly as cheaply anyway.
        z = elevation[tr, tc]
        if not np.isnan(z):
            ray_z = z0 + antenna_height_m + d * tan_theta
            # The signal path is refracted; it uses the radio radius, not the true
            # one that the particle geometry uses
            terrain_z = z - (d * d) * radio_inv_2R
            clearance = ray_z - terrain_z
            r1 = np.sqrt(wavelength_m * d * (d_end - d) / d_end)
            if r1 > 0.0:
                ratio = clearance / r1
                if ratio < worst:
                    worst = ratio
        d += step_m
    return worst


@jit(nopython=True, nogil=True, parallel=True)
def scan_candidates(candidates, elevation, cell_size_y, cell_size_x, rows, cols,
                    azimuth_offsets_deg, use_aspect,
                    elev_min_deg, elev_max_deg, n_bins,
                    step_m, max_range_m,
                    min_dist_m, max_dist_m,
                    min_depth_gcm2, require_terrain,
                    rock_density, earth_radius_m, wavelength_m, shower_offset_m,
                    antenna_height_m, near_field_m, radio_earth_radius_m,
                    bx, by, bz, use_geomag,
                    sea_level_density, scale_height_m, crust_density, bilinear,
                    out_cells, out_solid_angle, out_mean_dist,
                    out_max_depth, out_mean_depth, out_horizon, out_clearance,
                    out_geomag_omega, out_grammage, out_earth_chord):
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
    # Tangents of the bin edges, so the walk needs no arctangent. There are n_bins+1
    # of them: the lower edges plus the top of the last bin, without which samples in
    # the highest bin fall through into the overflow bin.
    tan_edges = np.empty(n_bins + 1, dtype=np.float64)
    for b in range(n_bins + 1):
        tan_edges[b] = np.tan(np.radians(elev_min_deg + b * elev_bin_deg))
    radio_inv_2R = 1.0 / (2.0 * radio_earth_radius_m)

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
        clearance_best = -1.0e30
        geomag_omega = 0.0
        grammage_sum = 0.0
        chord_sum = 0.0

        if np.isnan(z0):
            out_cells[i] = 0
            out_solid_angle[i] = 0.0
            out_mean_dist[i] = 0.0
            out_max_depth[i] = 0.0
            out_mean_depth[i] = 0.0
            out_horizon[i] = 0.0
            out_clearance[i] = 0.0
            out_geomag_omega[i] = 0.0
            out_grammage[i] = 0.0
            out_earth_chord[i] = 0.0
            continue

        for a in range(n_az):
            azimuth = azimuth_offsets_deg[a]
            if use_aspect:
                azimuth = aspect + azimuth

            horizon = _scan_one_direction(
                elevation, r0, c0, z0, azimuth,
                cell_size_y, cell_size_x, rows, cols,
                tan_edges, n_bins,
                step_m, max_range_m, inv_2R, bilinear,
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
                    if wavelength_m > 0.0:
                        # Worst clearance *along* this path, but the best *across*
                        # directions: a site can use its best-cleared direction, so
                        # one grazing direction must not condemn the whole pixel.
                        ratio = _min_clearance_ratio(
                            elevation, r0, c0, z0, azimuth,
                            cell_size_y, cell_size_x, rows, cols,
                            theta, d_hit, step_m, inv_2R, wavelength_m, shower_offset_m,
                            antenna_height_m, near_field_m, radio_inv_2R)
                        if ratio > clearance_best:
                            clearance_best = ratio
                    cell_omega = cos_theta * d_theta * d_phi

                    # Geomagnetic emission goes as |v x B|: a shower travelling along
                    # the field radiates almost none of it, so the azimuth of a target
                    # matters and not merely its existence.
                    if use_geomag:
                        sin_t = np.sin(np.radians(theta))
                        vx = np.sin(np.radians(azimuth)) * cos_theta
                        vy = np.cos(np.radians(azimuth)) * cos_theta
                        dot = vx * bx + vy * by + sin_t * bz
                        sin_alpha = np.sqrt(max(0.0, 1.0 - dot * dot))
                        geomag_omega += cell_omega * sin_alpha
                    else:
                        geomag_omega += cell_omega

                    # Atmospheric depth over the path: a shower develops through
                    # g/cm^2, and air at 4000 m is a third thinner than at sea level
                    slant = d_hit / cos_theta
                    sin_t2 = np.sin(np.radians(theta))
                    base = sea_level_density * np.exp(-z0 / scale_height_m)
                    if abs(sin_t2) < 1.0e-9:
                        gram = base * slant
                    else:
                        rise = slant * sin_t2
                        gram = base * (scale_height_m / sin_t2) * (
                            1.0 - np.exp(-rise / scale_height_m))
                    grammage_sum += gram * 0.1          # kg/m^2 -> g/cm^2

                    # Earth chord for downgoing directions, 2R sin(theta): it dwarfs
                    # local topography and governs neutrino attenuation
                    if theta < 0.0:
                        chord_sum += (2.0 * earth_radius_m
                                      * np.sin(np.radians(-theta))) * 100.0 * crust_density

                    cells += 1
                    solid_angle += cell_omega
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
                    geomag_omega += cos_theta * d_theta * d_phi

        out_cells[i] = cells
        out_solid_angle[i] = solid_angle
        out_mean_dist[i] = dist_sum / cells if cells > 0 else 0.0
        out_max_depth[i] = depth_max
        out_mean_depth[i] = depth_sum / cells if cells > 0 else 0.0
        out_horizon[i] = horizon_max if horizon_max > -1.0e29 else 0.0
        out_clearance[i] = clearance_best if clearance_best > -1.0e29 else 0.0
        out_geomag_omega[i] = geomag_omega
        out_grammage[i] = grammage_sum / cells if cells > 0 else 0.0
        out_earth_chord[i] = chord_sum / cells if cells > 0 else 0.0


def tau_decay_length_m(energy_pev):
    """
    Lorentz-boosted tau decay length, ``(E/m) * c*tau``.

    This is what sets the scale of the useful detector-to-exit-point distance, and it
    is exactly analytic — no simulation input required. Worth noting that it reproduces
    the published numbers on both sides: 1-100 PeV gives 49 m to 4.9 km, matching
    TAMBO's quoted 50 m - 5 km range, while the searcher's inherited 10-80 km GRAND
    window corresponds to 0.2-1.6 EeV.
    """
    return (energy_pev * 1.0e6 / TAU_MASS_GEV) * TAU_CTAU_M


def energy_pev_for_decay_length(distance_m):
    """Inverse of :func:`tau_decay_length_m`, for reporting what a distance implies."""
    return distance_m / TAU_CTAU_M * TAU_MASS_GEV / 1.0e6


def decay_probability(min_dist_m, max_dist_m, energy_pev):
    """
    Probability the tau decays inside the accepted baseline window.

    ``exp(-d_min/L) - exp(-d_max/L)``. One of the few factors that needs no acceptance
    table, and the one that couples most strongly to site geometry.
    """
    length = tau_decay_length_m(energy_pev)
    if length <= 0:
        return 0.0
    return math.exp(-min_dist_m / length) - math.exp(-max_dist_m / length)


def distance_window_from_energy(energy_min_pev, energy_max_pev, shower_development_m=3000.0):
    """
    A decay-baseline window implied by an energy range.

    The tau must decay before reaching the detector and the shower then needs room to
    develop, so the window runs from about one decay length at the low end to one at
    the high end plus the shower length. Ref. [2] quotes 3-10 km of shower development.

    This is a stated convention rather than a derivation: it fixes the *scale*
    correctly, but the useful window also depends on acceptance details this tool does
    not model. Callers can always set the distances directly.

    Returns (min_dist_m, max_dist_m).
    """
    return (tau_decay_length_m(energy_min_pev),
            tau_decay_length_m(energy_max_pev) + shower_development_m)


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


def balanced_order(n_candidates, n_threads, block=256):
    """
    Candidate ordering that balances thread load without destroying locality.

    Numba's ``prange`` schedules statically, giving each thread one contiguous slice of
    the index range. Candidates leave the topographic screen in spatial order, so that
    slice is a contiguous patch of map -- and walk cost varies enormously across the
    map, since rays near an edge terminate early while interior ones run the full
    range. The result is that some threads finish long before others: measured scaling
    was 2.4x on 12 cores against 4-5x for randomly scattered candidates.

    Shuffling fixes the balance but destroys cache locality, and measured barely better
    overall. Dealing *blocks* of neighbouring candidates round-robin keeps locality
    inside a block while spreading each thread's slice across the whole map.

    Returns an index array, or None when reordering cannot help.
    """
    if n_threads <= 1 or n_candidates < block * n_threads * 2:
        return None
    n_blocks = (n_candidates + block - 1) // block
    block_ids = np.arange(n_blocks)
    dealt = np.concatenate([block_ids[r::n_threads] for r in range(n_threads)])
    order = np.concatenate([np.arange(b * block, min((b + 1) * block, n_candidates))
                            for b in dealt])
    return order


def scan(candidates, elevation, map_grid, *,
         elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
         n_azimuths=9, half_width_deg=60.0, use_aspect=True,
         step_m=None, max_range_m=80000.0,
         min_dist_km=0.0, max_dist_km=80.0,
         min_depth_gcm2=0.0, require_terrain=True,
         rock_density=STANDARD_ROCK_DENSITY,
         earth_radius_m=TRUE_EARTH_RADIUS_M,
         radio_earth_radius_m=RADIO_EARTH_RADIUS_M,
         frequency_mhz=None, shower_offset_m=3000.0, antenna_height_m=2.0,
         near_field_m=500.0, bilinear=True,
         geomag_declination_deg=None, geomag_inclination_deg=None):
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
    # Fresnel clearance is measured only when a band is given; 0 disables the second pass
    wavelength_m = 0.0 if not frequency_mhz else SPEED_OF_LIGHT / (frequency_mhz * 1.0e6)

    # Geomagnetic weighting only when a field is supplied: guessing a field vector
    # would be worse than declining to weight at all
    use_geomag = (geomag_declination_deg is not None and geomag_inclination_deg is not None)
    bx, by, bz = (physics.geomagnetic_unit_vector(geomag_declination_deg,
                                                  geomag_inclination_deg)
                  if use_geomag else (0.0, 0.0, 0.0))

    out = {
        "cells": np.zeros(n, dtype=np.int64),
        "solid_angle_sr": np.zeros(n, dtype=np.float64),
        "mean_distance_m": np.zeros(n, dtype=np.float64),
        "max_depth_gcm2": np.zeros(n, dtype=np.float64),
        "mean_depth_gcm2": np.zeros(n, dtype=np.float64),
        "horizon_deg": np.zeros(n, dtype=np.float64),
        "best_clearance_ratio": np.zeros(n, dtype=np.float64),
        "geomag_solid_angle_sr": np.zeros(n, dtype=np.float64),
        "path_grammage_gcm2": np.zeros(n, dtype=np.float64),
        "earth_chord_gcm2": np.zeros(n, dtype=np.float64),
    }
    if n == 0:
        return out

    if elevation.dtype != np.float32 and elevation.dtype != np.float64:
        elevation = elevation.astype(np.float32)

    # Reorder for thread balance, then undo it so callers see their own ordering
    try:
        import numba as _numba
        order = balanced_order(n, _numba.get_num_threads())
    except Exception:
        order = None
    if order is not None:
        candidates = np.ascontiguousarray(candidates[order])

    scan_candidates(
        candidates, elevation, map_grid.cell_size_y, map_grid.cell_size_x, rows, cols,
        offsets, use_aspect,
        elev_min_deg, elev_max_deg, n_elev_bins,
        float(step_m), float(max_range_m),
        min_dist_km * 1000.0, max_dist_km * 1000.0,
        float(min_depth_gcm2), bool(require_terrain),
        float(rock_density), float(earth_radius_m), wavelength_m, float(shower_offset_m),
        float(antenna_height_m), float(near_field_m), float(radio_earth_radius_m),
        float(bx), float(by), float(bz), bool(use_geomag),
        physics.SEA_LEVEL_DENSITY_KGM3, physics.DENSITY_SCALE_HEIGHT_M,
        physics.CRUST_DENSITY_GCM3, bool(bilinear),
        out["cells"], out["solid_angle_sr"], out["mean_distance_m"],
        out["max_depth_gcm2"], out["mean_depth_gcm2"], out["horizon_deg"],
        out["best_clearance_ratio"],
        out["geomag_solid_angle_sr"], out["path_grammage_gcm2"], out["earth_chord_gcm2"],
    )

    if order is not None:
        inverse = np.empty_like(order)
        inverse[order] = np.arange(len(order))
        for key in out:
            out[key] = out[key][inverse]
    return out


@jit(nopython=True, nogil=True, parallel=True)
def _rfi_exposure(candidates, elevation, zone_rows, zone_cols, zone_weight,
                  cell_size_y, cell_size_x, rows, cols, step_m, inv_2R, out_exposure):
    """
    Radio-noise exposure, counting only sources a candidate can actually see.

    A circular exclusion zone treats a town behind a ridge exactly like one in plain
    view, which the terrain says is wrong. This walks the straight line to each source
    and drops the ones the terrain occludes; survivors contribute as 1/d^2.

    The horizon machinery this reuses is already computed for the arrival scan, so the
    only new cost is one short walk per candidate per source.
    """
    n = candidates.shape[0]
    n_zones = zone_rows.shape[0]

    for i in prange(n):
        r0 = int(candidates[i, 0])
        c0 = int(candidates[i, 1])
        z0 = elevation[r0, c0]
        exposure = 0.0
        if np.isnan(z0):
            out_exposure[i] = 0.0
            continue

        for k in range(n_zones):
            dr = (zone_rows[k] - r0) * cell_size_y
            dc = (zone_cols[k] - c0) * cell_size_x
            dist = np.sqrt(dr * dr + dc * dc)
            if dist < 1.0:
                out_exposure[i] += zone_weight[k] * 1.0e6
                continue

            zr = int(zone_rows[k])
            zc = int(zone_cols[k])
            if zr < 0 or zr >= rows or zc < 0 or zc >= cols:
                continue
            z_zone = elevation[zr, zc]
            if np.isnan(z_zone):
                continue

            # Walk the straight line and stop at the first terrain above it
            blocked = False
            steps = int(dist / step_m)
            for sidx in range(1, steps):
                frac = sidx / steps
                d = dist * frac
                tr = int(r0 + (zone_rows[k] - r0) * frac)
                tc = int(c0 + (zone_cols[k] - c0) * frac)
                if tr < 0 or tr >= rows or tc < 0 or tc >= cols:
                    continue
                z = elevation[tr, tc]
                if np.isnan(z):
                    continue
                line_z = z0 + (z_zone - z0) * frac
                if z - (d * d) * inv_2R > line_z:
                    blocked = True
                    break

            if not blocked:
                exposure += zone_weight[k] / (dist * dist)

        out_exposure[i] = exposure


def rfi_exposure(candidates, elevation, map_grid, zones_rowcol_weight, step_m=None,
                 earth_radius_m=RADIO_EARTH_RADIUS_M):
    """
    Line-of-sight-weighted radio noise exposure for each candidate.

    ``zones_rowcol_weight`` is an iterable of (row, col, weight); weight is normally the
    zone's radius or population proxy. Sources hidden behind terrain contribute nothing.

    Returns an array of exposure in weight per metre squared; smaller is quieter.
    """
    candidates = np.ascontiguousarray(candidates, dtype=np.float64)
    n = candidates.shape[0]
    out = np.zeros(n, dtype=np.float64)
    zones = list(zones_rowcol_weight)
    if n == 0 or not zones:
        return out

    zr = np.array([z[0] for z in zones], dtype=np.float64)
    zc = np.array([z[1] for z in zones], dtype=np.float64)
    zw = np.array([z[2] for z in zones], dtype=np.float64)
    if step_m is None:
        step_m = min(map_grid.cell_size_y, map_grid.cell_size_x)

    _rfi_exposure(candidates, elevation, zr, zc, zw,
                  map_grid.cell_size_y, map_grid.cell_size_x,
                  elevation.shape[0], elevation.shape[1],
                  float(step_m), 1.0 / (2.0 * earth_radius_m), out)
    return out
