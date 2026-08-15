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

import numpy as np

_trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))

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


# Shower maximum deepens with energy: more generations are needed before the average
# particle drops below the critical energy. X_MAX_GCM2 is quoted at this reference.
X_MAX_REFERENCE_ENERGY_PEV = 1.0e3               # 1 EeV
# Elongation rate, g/cm^2 per decade of primary energy. ~55 is the usual hadronic
# value; a purely electromagnetic cascade would be nearer the Heitler 85.
ELONGATION_RATE_GCM2_PER_DECADE = 55.0
# Gaisser-Hillas interaction length, setting how fast the profile rises and falls
SHOWER_PROFILE_LAMBDA_GCM2 = 70.0


def shower_maximum_gcm2(energy_pev, x_max_ref_gcm2=X_MAX_GCM2,
                        reference_energy_pev=X_MAX_REFERENCE_ENERGY_PEV,
                        elongation_rate=ELONGATION_RATE_GCM2_PER_DECADE):
    """
    Depth of shower maximum at a given primary energy.

        X_max(E) = X_max(E_ref) + D * log10(E / E_ref)

    Over TAMBO's 3 PeV to 1 EeV this runs from about 560 to 700 g/cm^2, so the energy
    dependence is real but mild -- the band below is set far more by how much of the
    profile is being accepted than by where its peak sits.
    """
    energy = np.asarray(energy_pev, dtype=np.float64)
    return x_max_ref_gcm2 + elongation_rate * np.log10(energy / reference_energy_pev)


def shower_size_fraction(grammage_gcm2, x_max_gcm2,
                         lambda_gcm2=SHOWER_PROFILE_LAMBDA_GCM2):
    """
    Charged-particle content at depth X, as a fraction of the content at maximum.

    Gaisser-Hillas with the first interaction at X_0 = 0:

        N(X)/N_max = (X / X_max)^(X_max/lambda) * exp((X_max - X) / lambda)

    Zero at and below X = 0. This is what makes a particle array's criterion a band
    rather than a threshold: the content rises steeply, peaks, and then dies.
    """
    x = np.asarray(grammage_gcm2, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(x > 0.0, x / x_max_gcm2, 0.0)
        out = np.where(
            x > 0.0,
            ratio ** (x_max_gcm2 / lambda_gcm2) * np.exp((x_max_gcm2 - x) / lambda_gcm2),
            0.0)
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def grammage_band_from_energy(energy_min_pev, energy_max_pev, fraction=0.1,
                              lambda_gcm2=SHOWER_PROFILE_LAMBDA_GCM2, samples=4000,
                              **x_max_kw):
    """
    Atmospheric-depth band over which a particle array still sees a usable shower.

    A particle detector needs the cascade to have grown but not yet died, so the
    criterion is the depth range where the charged-particle content stays above
    ``fraction`` of its own maximum. The band is taken across the requested energy
    range: the low edge at the lowest energy, whose maximum is shallowest, and the high
    edge at the highest, whose maximum is deepest.

    At TAMBO's 3 PeV to 1 EeV and fraction 0.1 this gives roughly 235 to 1300 g/cm^2.
    That matters for siting: a canyon crossing supplies only what its own width of air
    contains -- about 170 g/cm^2 across 2 km of Colca, and ~390 g/cm^2 across its full
    4.5 km rim to rim -- so the criterion selects the *widest* crossings, which is the
    physically right answer and not one the default (X_max, 4*X_max) band could express.

    Returns (low_gcm2, high_gcm2).
    """
    grid = np.linspace(1.0, 5000.0, samples)

    def edges(energy):
        x_max = float(shower_maximum_gcm2(energy, **x_max_kw))
        n = shower_size_fraction(grid, x_max, lambda_gcm2)
        ok = grid[n >= fraction]
        if ok.size == 0:                            # pragma: no cover - degenerate
            return x_max, x_max
        return float(ok.min()), float(ok.max())

    lo = edges(energy_min_pev)[0]
    hi = edges(energy_max_pev)[1]
    return (lo, hi) if lo <= hi else (hi, lo)


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


# ------------------------------------------------------- muon shielding

# Rock overburden a detector wants along the arrival direction so that atmospheric
# muons from that direction cannot reach it. Ref. [2] Fig. 1 annotates >4 km for TAMBO.
DEFAULT_MUON_SHIELDING_KM = 4.0


def muon_shielding_gcm2(thickness_km, density_gcm3=CRUST_DENSITY_GCM3):
    """
    Column depth corresponding to a rock thickness, for muon rejection.

    A mountain in the arrival direction is a muon filter: an air-shower muon from that
    direction would have to cross the whole thickness, which a few km of rock makes
    impossible. Anything detected from behind that much rock is therefore not a
    cosmic-ray muon, which is what lets a surface array claim neutrino purity.

    Unlike the production-and-escape band, this is a **floor**: more rock is always
    better for background rejection, and only the signal side wants an upper limit.

    4 km of standard rock is about 1.06e6 g/cm^2.
    """
    return thickness_km * 1000.0 * 100.0 * density_gcm3


# ---------------------------------------------------------------- tau range

def tau_decay_length_m(energy_pev, mass_gev=1.77686, ctau_m=87.03e-6):
    """Lorentz-boosted decay length, ``(E/m) c*tau``."""
    return (energy_pev * 1.0e6 / mass_gev) * ctau_m


# Charged-current nu-N cross-section, a power-law fit of the standard
# parameterisations: sigma = A (E/GeV)^n cm^2. Good to tens of per cent over
# 1e8-1e10 GeV and an extrapolation above that, where no data constrain it.
SIGMA_CC_COEFF_CM2 = 6.04e-36
SIGMA_CC_INDEX = 0.358
AVOGADRO = 6.022e23


def cc_cross_section_cm2(energy_pev):
    """Charged-current neutrino-nucleon cross-section."""
    return SIGMA_CC_COEFF_CM2 * (energy_pev * 1.0e6) ** SIGMA_CC_INDEX


def neutrino_interaction_length_gcm2(energy_pev):
    """
    Column depth over which a neutrino interacts once, ``1/(N_A sigma)``.

    Falls from about 3.8e8 g/cm^2 at 100 PeV to 7e7 at 10 EeV. Only charged-current
    attenuation is counted; neutral-current regeneration would soften it slightly.
    """
    return 1.0 / (AVOGADRO * cc_cross_section_cm2(energy_pev))


# Tau energy loss is essentially all photonuclear. Bremsstrahlung and pair production
# both scale as 1/m^2, so relative to the muon they are suppressed by
# (m_mu/m_tau)^2 = 1/283 and contribute only a few per cent -- which is exactly why a
# tau outranges a muon. Photonuclear depends on the lepton mass only through the
# virtual-photon flux, logarithmically, so it survives at close to the muon value.
#
# Estimating from muon coefficients in standard rock (brems 1.6, pair 2.0,
# photonuclear 0.4, all 1e-6 cm^2/g) gives beta_tau ~ (0.4-1.0)e-6 cm^2/g. Cross-check:
# 1/beta is then 3.8-9.4 km of rock, bracketing the ~10 km at which the tau range is
# usually quoted to saturate.
#
# Photonuclear cross-sections grow with energy, so beta grows too. The power law below
# reproduces 0.38e-6 at 100 PeV and 0.95e-6 at 10 EeV. It is an ESTIMATE, not a fit to
# published tables; set the index to zero to recover a constant.
BETA_REFERENCE_CM2G = 0.6e-6
BETA_REFERENCE_ENERGY_PEV = 1.0e3        # 1 EeV
BETA_ENERGY_INDEX = 0.20

# Mean inelasticity of a charged-current interaction at these energies: the tau carries
# away roughly (1 - y) of the neutrino energy, with <y> about 0.2.
CC_INELASTICITY = 0.2


def tau_energy_loss_beta(energy_pev, reference=BETA_REFERENCE_CM2G,
                         reference_energy_pev=BETA_REFERENCE_ENERGY_PEV,
                         index=BETA_ENERGY_INDEX):
    """Energy-loss coefficient beta(E), rising with energy as photonuclear does."""
    if index == 0.0:
        return reference
    return reference * (energy_pev / reference_energy_pev) ** index


def tau_range_gcm2(energy_pev, beta_cm2g=None, density_gcm3=CRUST_DENSITY_GCM3):
    """
    Column depth over which a tau's survival probability falls to 1/e.

    Decay and energy loss couple, because losing energy shortens the boosted decay
    length. With ``E(X) = E0 exp(-beta X)`` the decay probability per unit depth is
    ``exp(beta X)/X_decay(E0)``, and integrating the survival gives

        S(X) = exp( -(X_loss/X_decay) (exp(X/X_loss) - 1) ),   X_loss = 1/beta

    so the 1/e point is

        R = X_loss * ln(1 + X_decay/X_loss)

    Note this **grows logarithmically** at high energy rather than saturating at
    ``1/beta``. An earlier version of this module combined the two lengths
    harmonically, which saturates and underestimates the range by a factor of 2 at
    an EeV and 4 at 10 EeV.
    """
    beta = tau_energy_loss_beta(energy_pev) if beta_cm2g is None else beta_cm2g
    x_decay = tau_decay_length_m(energy_pev) * 100.0 * density_gcm3
    x_loss = 1.0 / beta
    return x_loss * math.log1p(x_decay / x_loss)


def tau_survival(depth_gcm2, energy_pev, beta_cm2g=None, density_gcm3=CRUST_DENSITY_GCM3):
    """
    Probability a tau of the given energy crosses ``depth_gcm2`` of rock without decaying.

    The double-exponential form above, which falls far more sharply than a simple
    exponential once the depth exceeds ``1/beta``: the tau is losing energy, so its
    decay length shrinks as it goes.
    """
    beta = tau_energy_loss_beta(energy_pev) if beta_cm2g is None else beta_cm2g
    x_decay = tau_decay_length_m(energy_pev) * 100.0 * density_gcm3
    x_loss = 1.0 / beta
    d = np.asarray(depth_gcm2, dtype=np.float64)
    # exp overflows for very thick rock, where survival is zero anyway
    return np.exp(-(x_loss / x_decay) * np.expm1(np.clip(d / x_loss, 0.0, 700.0)))


def tau_exit_probability(column_depth_gcm2, energy_pev, beta_cm2g=None,
                         density_gcm3=CRUST_DENSITY_GCM3,
                         inelasticity=CC_INELASTICITY, samples=2000):
    """
    Relative probability that a traversing neutrino yields a tau that escapes.

    The neutrino interacts at depth x with probability ``exp(-x/lambda) dx/lambda``;
    the tau, carrying ``(1 - y)`` of the energy, must then cross the remaining
    ``X - x`` and survive:

        P(X) = integral_0^X (dx/lambda) exp(-x/lambda) S(X - x)

    Evaluated numerically, because ``S`` is a double exponential and the integral has
    no clean closed form. Relative, not absolute: normalisation needs the trigger
    response (see aperture.py).
    """
    e_tau = energy_pev * (1.0 - inelasticity)
    lam = neutrino_interaction_length_gcm2(energy_pev)

    def one(X):
        if X <= 0:
            return 0.0
        x = np.linspace(0.0, float(X), samples)
        s = tau_survival(X - x, e_tau, beta_cm2g, density_gcm3)
        return float(_trapezoid(np.exp(-x / lam) / lam * s, x))

    if np.ndim(column_depth_gcm2) == 0:
        return one(column_depth_gcm2)
    return np.array([one(v) for v in np.asarray(column_depth_gcm2)])


def production_escape_optimum_gcm2(energy_pev, beta_cm2g=None,
                                   density_gcm3=CRUST_DENSITY_GCM3,
                                   inelasticity=CC_INELASTICITY, samples=400):
    """
    Column depth maximising :func:`tau_exit_probability`.

    Found on a log grid, since the corrected exit probability has no closed-form peak.
    Rises with energy and then flattens -- about 12 km of standard rock at 100 PeV,
    20 km at 1 EeV, 23 km at 10 EeV -- as the logarithmic growth of the tau range is
    tempered by beta rising.
    """
    grid = np.logspace(4.0, 9.0, samples)
    p = tau_exit_probability(grid, energy_pev, beta_cm2g, density_gcm3, inelasticity)
    return float(grid[int(np.argmax(p))])


def depth_band_from_energy(energy_min_pev, energy_max_pev, fraction=0.5,
                           beta_cm2g=None, density_gcm3=CRUST_DENSITY_GCM3,
                           inelasticity=CC_INELASTICITY, samples=400):
    """
    Column-depth band where the tau exit probability stays above ``fraction`` of peak.

    Very wide -- roughly 5e5 to 2.6e8 g/cm^2 at half maximum across 100 PeV to 10 EeV,
    some two and a half decades -- which is itself the result: column depth is an
    intrinsically weak discriminant and the criterion should not pretend otherwise.

    Returns (low_gcm2, high_gcm2), the low edge taken at the lowest energy and the high
    edge at the highest, so the band spans the whole requested range.
    """
    grid = np.logspace(4.0, 9.0, samples)

    def edges(energy):
        p = tau_exit_probability(grid, energy, beta_cm2g, density_gcm3, inelasticity)
        ok = grid[p >= fraction * p.max()]
        return float(ok.min()), float(ok.max())

    lo = edges(energy_min_pev)[0]
    hi = edges(energy_max_pev)[1]
    return (lo, hi) if lo <= hi else (hi, lo)


def earth_absorption_cutoff_deg(energy_pev, fraction=0.5, radius_m=EARTH_RADIUS_M,
                                density_gcm3=CRUST_DENSITY_GCM3, **kw):
    """
    Elevation below which the Earth chord itself exceeds the useful column depth.

    The chord is ``2R sin(theta)``, hundreds of km for degrees below the horizontal, so
    steep arrival directions carry far more matter than the optimum wants and the
    neutrino is absorbed before reaching the exit region. Setting the chord equal to the
    upper band edge gives the elevation at which acceptance has fallen to ``fraction``
    of peak.

    The result narrows sharply with energy -- about -4.5 degrees at 100 PeV, -2.1 at
    1 EeV, -1.0 at 10 EeV -- so the *effective* arrival window is not a fixed +/-3
    degrees but an energy-dependent one whose lower edge climbs toward the horizon.

    Returns a negative angle, or None when the cut lies below the horizon entirely.
    """
    X = np.logspace(4, 9, 400)
    P = tau_exit_probability(X, energy_pev, **kw)
    upper = float(X[P >= fraction * P.max()].max())
    chord_m = upper / density_gcm3 / 100.0
    sin_theta = chord_m / (2.0 * radius_m)
    if sin_theta >= 1.0:
        return None
    return -math.degrees(math.asin(sin_theta))


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


def default_field_for_site(latitude_deg, longitude_deg,
                           declination_deg=None, inclination_deg=None):
    """
    Geomagnetic field for an arbitrary site, falling back sensibly.

    Inclination is computed from the site's own coordinates with the dipole model, so
    moving the search to another location gets the right inclination automatically --
    it is the quantity that varies most across Peru, from about -5 degrees near Lima to
    -14 near Arequipa.

    Declination cannot be had that way: the dipole gives about -0.2 degrees at Arequipa
    against an IGRF -6.9, so non-dipole terms dominate. It therefore falls back to the
    Arequipa IGRF value, which is right for the prototype region and approximately
    right for the rest of southern Peru. Supply the IGRF declination for anywhere else.

    Returns (declination_deg, inclination_deg).
    """
    if declination_deg is None:
        declination_deg = DEFAULT_GEOMAG_DECLINATION_DEG
    if inclination_deg is None:
        inclination_deg = centered_dipole_inclination(latitude_deg, longitude_deg)
    return float(declination_deg), float(inclination_deg)


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
