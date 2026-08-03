"""
Closed-form physics the terrain scan needs but cannot measure from a DEM.

Everything here is analytic. The scan supplies geometry -- where the tau exits, how far
away, through how much local rock -- and these functions supply the physics that
geometry alone does not determine: how much atmosphere a shower has to develop in, how
much Earth a neutrino crossed to get there, how strongly the shower radiates into a
given direction, and how large its radio footprint is when it arrives.

Constants that are genuinely uncertain are parameters with stated defaults rather than
literals buried in an expression. Where a default encodes a convention or an
approximation, the docstring says so.
"""

import math

# Atmosphere: an exponential isothermal approximation, adequate over the few km of
# relief a site search spans.
SEA_LEVEL_DENSITY_KGM3 = 1.225
DENSITY_SCALE_HEIGHT_M = 8400.0
# Refractivity at sea level, (n - 1). Falls with density, hence with altitude.
SEA_LEVEL_REFRACTIVITY = 2.9e-4
# Depth of shower maximum, order of magnitude for the PeV-EeV range
X_MAX_GCM2 = 700.0

# Earth
EARTH_RADIUS_M = 6.371e6
CRUST_DENSITY_GCM3 = 2.65

# Tau energy loss in rock: dE/dX = -(a + bE), with b the term that matters at these
# energies. Values in the literature span roughly 0.4-1.0e-6 cm^2/g and are themselves
# energy-dependent, so this is the least certain number here.
TAU_ENERGY_LOSS_B_CM2G = 0.5e-6

KGM2_TO_GCM2 = 0.1


# ---------------------------------------------------------------- atmosphere

def air_density_kgm3(altitude_m, sea_level_density=SEA_LEVEL_DENSITY_KGM3,
                     scale_height_m=DENSITY_SCALE_HEIGHT_M):
    """Exponential atmosphere, ``rho0 * exp(-h/H)``."""
    return sea_level_density * math.exp(-altitude_m / scale_height_m)


def slant_grammage_gcm2(start_altitude_m, elevation_deg, distance_m,
                        sea_level_density=SEA_LEVEL_DENSITY_KGM3,
                        scale_height_m=DENSITY_SCALE_HEIGHT_M):
    """
    Atmospheric depth along a slanted path, in g/cm^2.

    A shower develops through grammage, not through metres, and air density falls by a
    third between 2000 m and 4500 m. A site search that compares candidates at
    different altitudes while measuring path length in metres is comparing unlike
    things: 20 km at 4000 m is about 1500 g/cm^2, the same 20 km at sea level about
    2450 g/cm^2.

    Integrating ``rho0 exp(-(z0 + l sin(theta))/H) dl`` along the slant path has a
    closed form, so no numerical integration is needed:

        X = rho0 exp(-z0/H) * H/sin(theta) * (1 - exp(-D sin(theta)/H)) / 1

    with the horizontal-path limit ``rho0 exp(-z0/H) * D / cos(theta)`` as theta -> 0.

    Parameters:
    - start_altitude_m: altitude of the near end of the path.
    - elevation_deg: elevation angle of the path, positive upward.
    - distance_m: ground distance covered.
    """
    if distance_m <= 0:
        return 0.0
    theta = math.radians(elevation_deg)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    slant = distance_m / cos_t if abs(cos_t) > 1e-12 else distance_m
    base = sea_level_density * math.exp(-start_altitude_m / scale_height_m)

    if abs(sin_t) < 1e-9:
        kgm2 = base * slant
    else:
        rise = slant * sin_t
        kgm2 = base * (scale_height_m / sin_t) * (1.0 - math.exp(-rise / scale_height_m))
    return kgm2 * KGM2_TO_GCM2


def shower_maturity(grammage_gcm2, x_max_gcm2=X_MAX_GCM2):
    """
    Path grammage as a fraction of the depth of shower maximum.

    Below 1 the shower is still developing when it arrives. What "above 1" means
    depends on what is being detected, and the two cases are not alike:

    *Radio.* Emission comes from the region around shower maximum and then simply
    propagates; air is effectively transparent at 50-200 MHz. Being far beyond maximum
    costs nothing directly, so the criterion is a **threshold**, not a band. The real
    trade at greater distance is amplitude against footprint area, which belongs to the
    footprint term rather than here.

    *Particles.* The charged-particle content peaks at maximum and dies away after, so
    a particle array such as TAMBO does want to sit near it, and there the criterion
    genuinely is a band.

    This is one more reason criteria have to be per-channel rather than global.
    """
    return grammage_gcm2 / x_max_gcm2 if x_max_gcm2 > 0 else 0.0


# ---------------------------------------------------------------- Earth chord

def earth_chord_m(elevation_deg, radius_m=EARTH_RADIUS_M):
    """
    Chord length through the Earth for a ray arriving from below the horizontal.

    A ray making angle theta below the local horizontal cuts a chord of ``2R sin(theta)``:
    zero along the tangent, a full diameter straight down. Directions at or above the
    horizontal return zero.

    This matters because it dwarfs local topography. At -1 degree the chord is about
    220 km, at -3 degrees about 670 km, against the tens of km of mountain a DEM can
    see. The deepest part of even a 670 km chord lies only ~9 km below the surface, so
    it stays in the crust and a constant density is adequate.
    """
    if elevation_deg >= 0:
        return 0.0
    return 2.0 * radius_m * math.sin(math.radians(-elevation_deg))


def earth_chord_gcm2(elevation_deg, radius_m=EARTH_RADIUS_M,
                     density_gcm3=CRUST_DENSITY_GCM3):
    """Column depth of the Earth chord, in g/cm^2."""
    return earth_chord_m(elevation_deg, radius_m) * 100.0 * density_gcm3


def neutrino_survival(elevation_deg, interaction_length_gcm2,
                      radius_m=EARTH_RADIUS_M, density_gcm3=CRUST_DENSITY_GCM3):
    """
    Fraction of neutrinos surviving the Earth chord to reach the exit region.

    ``exp(-X_chord / X_int)``. The interaction length is a parameter because it depends
    on the cross-section at the energy of interest; around an EeV it is of order
    1e8 g/cm^2, which is the same order as the chord at -3 degrees, so the suppression
    across a +/-3 degree window is substantial rather than marginal.
    """
    if interaction_length_gcm2 <= 0:
        return 1.0
    return math.exp(-earth_chord_gcm2(elevation_deg, radius_m, density_gcm3)
                    / interaction_length_gcm2)


# ---------------------------------------------------------------- tau range

def tau_decay_length_m(energy_pev, mass_gev=1.77686, ctau_m=87.03e-6):
    """Lorentz-boosted decay length, ``(E/m) c*tau``."""
    return (energy_pev * 1.0e6 / mass_gev) * ctau_m


def tau_range_gcm2(energy_pev, beta_cm2g=TAU_ENERGY_LOSS_B_CM2G,
                   density_gcm3=CRUST_DENSITY_GCM3):
    """
    Effective tau range in rock, in g/cm^2.

    Two effects compete: the tau decays after a boosted decay length, and it loses
    energy continuously with ``1/beta`` as the loss length. Combining them
    harmonically, ``X_range = X_decay X_loss / (X_decay + X_loss)``, gives a range that
    grows with energy at low energy and saturates once energy loss dominates.

    This is an approximation with a genuinely uncertain constant -- published beta
    values span roughly 0.4-1.0e-6 cm^2/g and are themselves energy-dependent -- so it
    fixes the *scale* of the useful column depth rather than its precise value.
    """
    x_decay = tau_decay_length_m(energy_pev) * 100.0 * density_gcm3
    x_loss = 1.0 / beta_cm2g
    return x_decay * x_loss / (x_decay + x_loss)


def depth_band_from_energy(energy_min_pev, energy_max_pev, low_factor=0.3,
                           high_factor=3.0, beta_cm2g=TAU_ENERGY_LOSS_B_CM2G,
                           density_gcm3=CRUST_DENSITY_GCM3):
    """
    Column-depth band implied by an energy range.

    The tau must be produced, which needs rock, and must escape, which limits how much,
    so the useful depth sits around the tau range and the criterion is a band. The
    range grows with energy, so the band moves too -- a fixed band cannot be right for
    both ends of a wide energy interval.

    ``low_factor`` and ``high_factor`` set how far either side of the range the band
    extends; the defaults are deliberately generous.

    Returns (low_gcm2, high_gcm2).
    """
    lo = tau_range_gcm2(energy_min_pev, beta_cm2g, density_gcm3) * low_factor
    hi = tau_range_gcm2(energy_max_pev, beta_cm2g, density_gcm3) * high_factor
    return (lo, hi) if lo <= hi else (hi, lo)


# ---------------------------------------------------------------- geomagnetic

# Default field for the Peruvian Andes. Provenance differs between the two, and it
# matters:
#
#   declination -6.9 deg  -- IGRF 2026 at Arequipa (16.4 S, 71.5 W), retrieved from
#                            NOAA's geomagnetic calculator.
#   inclination -14.0 deg -- centered-dipole estimate at the same point (see
#                            centered_dipole_inclination), because the IGRF value was
#                            not available. It is approximate.
#
# Both should be replaced with IGRF values for the actual site. They are defaults so
# that the geomagnetic effect is modelled rather than silently omitted, not because
# one field fits all of Peru.
DEFAULT_GEOMAG_DECLINATION_DEG = -6.9
DEFAULT_GEOMAG_INCLINATION_DEG = -14.0

# North geomagnetic pole, recent epoch, for the dipole approximation
GEOMAGNETIC_POLE_LAT_DEG = 80.7
GEOMAGNETIC_POLE_LON_DEG = -72.7


def geomagnetic_latitude_deg(latitude_deg, longitude_deg,
                             pole_lat_deg=GEOMAGNETIC_POLE_LAT_DEG,
                             pole_lon_deg=GEOMAGNETIC_POLE_LON_DEG):
    """Latitude in the centered-dipole frame."""
    lat = math.radians(latitude_deg)
    plat = math.radians(pole_lat_deg)
    dlon = math.radians(longitude_deg - pole_lon_deg)
    sin_mlat = (math.sin(lat) * math.sin(plat)
                + math.cos(lat) * math.cos(plat) * math.cos(dlon))
    return math.degrees(math.asin(max(-1.0, min(1.0, sin_mlat))))


def centered_dipole_inclination(latitude_deg, longitude_deg, **kw):
    """
    Inclination from a centered dipole, ``tan(I) = 2 tan(magnetic latitude)``.

    Useful when no IGRF lookup is available: it captures the dominant behaviour, that
    inclination passes through zero at the magnetic equator and steepens away from it.
    At Arequipa it gives about -14 degrees, the site being roughly 7 degrees south of
    the magnetic equator.

    Note this approximation is only worth using for *inclination*. The dipole
    declination at Arequipa is about -0.2 degrees against an IGRF value of -6.9, so the
    non-dipole terms dominate there and a dipole declination would be misleading.
    """
    mlat = math.radians(geomagnetic_latitude_deg(latitude_deg, longitude_deg, **kw))
    return math.degrees(math.atan(2.0 * math.tan(mlat)))


def geomagnetic_unit_vector(declination_deg, inclination_deg):
    """
    Unit geomagnetic field in local East-North-Up coordinates.

    Declination is measured from geographic north, positive eastward; inclination is
    positive downward, so it enters the Up component with a minus sign.
    """
    d = math.radians(declination_deg)
    i = math.radians(inclination_deg)
    return (math.cos(i) * math.sin(d),      # east
            math.cos(i) * math.cos(d),      # north
            -math.sin(i))                   # up


def geomagnetic_sin_alpha(azimuth_deg, elevation_deg, field_unit_vector):
    """
    ``sin(alpha)`` between a shower axis and the geomagnetic field.

    Radio emission from an air shower is dominantly geomagnetic, with amplitude
    proportional to ``|v x B|``, so a shower travelling along the field radiates very
    little of it. Peru lies near the magnetic equator, where the field is close to
    horizontal and roughly northward: north-south showers are therefore strongly
    suppressed and east-west ones near maximal. The azimuth of a target matters, not
    merely whether one exists.

    The sign of the axis is irrelevant -- ``|(-v) x B| = |v x B|`` -- so the same value
    applies whether the direction is taken as arrival or propagation.
    """
    phi = math.radians(azimuth_deg)
    theta = math.radians(elevation_deg)
    cos_t = math.cos(theta)
    v = (math.sin(phi) * cos_t, math.cos(phi) * cos_t, math.sin(theta))
    bx, by, bz = field_unit_vector
    dot = v[0] * bx + v[1] * by + v[2] * bz
    return math.sqrt(max(0.0, 1.0 - dot * dot))


# ---------------------------------------------------------------- footprint

def refractivity(altitude_m, sea_level_value=SEA_LEVEL_REFRACTIVITY,
                 scale_height_m=DENSITY_SCALE_HEIGHT_M):
    """``n - 1`` at altitude, falling with density."""
    return sea_level_value * math.exp(-altitude_m / scale_height_m)


def cherenkov_angle_rad(altitude_m, **kw):
    """
    Cherenkov angle in air, ``sqrt(2(n-1))`` for small angles.

    About 1.4 degrees at sea level and 1.1 degrees at 4000 m: the cone narrows with
    altitude because the air is thinner.
    """
    return math.sqrt(2.0 * refractivity(altitude_m, **kw))


def cherenkov_footprint_radius_m(altitude_m, distance_m, **kw):
    """
    Radius of the radio footprint on the ground, ``D * theta_C``.

    The consequence for layout is counter-intuitive: a *higher* site has a *smaller*
    footprint, so it needs a *denser* array for the same trigger efficiency. A 1 km
    grid under-samples a footprint of a few hundred metres either way, which is why
    counted antennas are a cost proxy rather than an effective area.
    """
    return distance_m * cherenkov_angle_rad(altitude_m, **kw)


def footprint_sampling(spacing_m, altitude_m, distance_m, **kw):
    """
    Antennas per footprint diameter: ``2 r / spacing``.

    Below 1 the array does not resolve the footprint and triggering relies on a single
    antenna happening to fall inside it.
    """
    radius = cherenkov_footprint_radius_m(altitude_m, distance_m, **kw)
    return (2.0 * radius / spacing_m) if spacing_m > 0 else 0.0
