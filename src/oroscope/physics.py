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

from __future__ import annotations

import math

import numpy as np


_trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))

# The public surface. Declared so that `from physics import *` and the generated API
# reference both show what this module offers rather than everything it imported.
__all__ = [
    "air_density_kgm3", "slant_grammage_gcm2", "shower_maturity",
    "shower_maximum_gcm2", "shower_size_fraction", "grammage_band_from_energy",
    "earth_chord_m", "earth_chord_gcm2", "neutrino_survival", "muon_shielding_gcm2",
    "tau_decay_length_m", "cc_cross_section_cm2", "neutrino_interaction_length_gcm2",
    "nc_regeneration_factor", "NC_TO_CC_RATIO", "NC_INELASTICITY",
    "tau_energy_loss_beta", "tau_range_gcm2", "tau_survival", "tau_exit_probability",
    "set_tau_energy_loss", "restore_tau_energy_loss", "tau_energy_loss_settings",
    "production_escape_optimum_gcm2", "depth_band_from_energy",
    "earth_absorption_cutoff_deg", "spectrum_weighted_decay_probability",
    "geomagnetic_latitude_deg",
    "centered_dipole_inclination", "default_field_for_site",
    "set_declination_model", "declination_model", "declination_from_grid",
    "geomagnetic_unit_vector", "geomagnetic_sin_alpha", "refractivity",
    "cherenkov_angle_rad", "cherenkov_footprint_radius_m", "footprint_sampling",
]

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

def air_density_kgm3(altitude_m: float,
                     sea_level_density: float = SEA_LEVEL_DENSITY_KGM3,
                     scale_height_m: float = DENSITY_SCALE_HEIGHT_M) -> float:
    """
    Exponential atmosphere, ``rho0 * exp(-h/H)``.

    An isothermal approximation, adequate over the few kilometres of relief a site
    search spans. It is not a substitute for a real profile at large zenith angles.

    Parameters
    ----------
    altitude_m : float
        Altitude above sea level, in metres.
    sea_level_density : float, optional
        Density at sea level, in kg/m^3.
    scale_height_m : float, optional
        Density scale height, in metres.

    Returns
    -------
    float
        Air density in kg/m^3.

    Examples
    --------
    >>> from oroscope import physics
    >>> round(physics.air_density_kgm3(0.0), 3)
    1.225
    >>> round(physics.air_density_kgm3(4000.0), 3)   # a third thinner at Andean altitude
    0.761
    """
    return sea_level_density * math.exp(-altitude_m / scale_height_m)


def slant_grammage_gcm2(start_altitude_m: float, elevation_deg: float,
                        distance_m: float,
                        sea_level_density: float = SEA_LEVEL_DENSITY_KGM3,
                        scale_height_m: float = DENSITY_SCALE_HEIGHT_M) -> float:
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

    Parameters
    ----------
    start_altitude_m : float
        Altitude of the near end of the path, in metres.
    elevation_deg : float
        Elevation angle of the path, in degrees, positive upward.
    distance_m : float
        Ground distance covered, in metres. Zero or less returns zero.
    sea_level_density : float, optional
        Density at sea level, in kg/m^3.
    scale_height_m : float, optional
        Density scale height, in metres.

    Returns
    -------
    float
        Atmospheric depth along the path, in g/cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> horizontal = physics.slant_grammage_gcm2(4000.0, 0.0, 20000.0)
    >>> sea_level = physics.slant_grammage_gcm2(0.0, 0.0, 20000.0)
    >>> f"{horizontal:.0f} vs {sea_level:.0f} g/cm^2"
    '1522 vs 2450 g/cm^2'
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


def shower_maturity(grammage_gcm2: float | np.ndarray,
                    x_max_gcm2: float = X_MAX_GCM2) -> float | np.ndarray:
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

    Parameters
    ----------
    grammage_gcm2 : float or array_like
        Atmospheric depth traversed, in g/cm^2.
    x_max_gcm2 : float, optional
        Depth of shower maximum, in g/cm^2.

    Returns
    -------
    float or ndarray
        Grammage as a fraction of shower maximum. Below 1 the shower is still growing.

    See Also
    --------
    shower_size_fraction : the particle content itself, rather than this ratio.

    Examples
    --------
    >>> from oroscope import physics
    >>> round(physics.shower_maturity(750.0), 3)
    1.071
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


def shower_maximum_gcm2(energy_pev: float | np.ndarray,
                        x_max_ref_gcm2: float = X_MAX_GCM2,
                        reference_energy_pev: float = X_MAX_REFERENCE_ENERGY_PEV,
                        elongation_rate: float = ELONGATION_RATE_GCM2_PER_DECADE
                        ) -> np.ndarray:
    """
    Depth of shower maximum at a given primary energy.

        X_max(E) = X_max(E_ref) + D * log10(E / E_ref)

    Over TAMBO's 3 PeV to 1 EeV this runs from about 560 to 700 g/cm^2, so the energy
    dependence is real but mild -- the band below is set far more by how much of the
    profile is being accepted than by where its peak sits.

    Parameters
    ----------
    energy_pev : float or array_like
        Primary energy, in PeV.
    x_max_ref_gcm2 : float, optional
        Depth of maximum at the reference energy, in g/cm^2.
    reference_energy_pev : float, optional
        Energy at which ``x_max_ref_gcm2`` is quoted, in PeV.
    elongation_rate : float, optional
        Deepening per decade of energy, in g/cm^2. About 55 for a hadronic cascade;
        a purely electromagnetic one is nearer the Heitler value of 85.

    Returns
    -------
    ndarray
        Depth of shower maximum, in g/cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> [f"{float(physics.shower_maximum_gcm2(e)):.0f}" for e in (3.0, 1000.0)]
    ['561', '700']
    """
    energy = np.asarray(energy_pev, dtype=np.float64)
    return x_max_ref_gcm2 + elongation_rate * np.log10(energy / reference_energy_pev)


def shower_size_fraction(grammage_gcm2: float | np.ndarray, x_max_gcm2: float,
                         lambda_gcm2: float = SHOWER_PROFILE_LAMBDA_GCM2
                         ) -> np.ndarray:
    """
    Charged-particle content at depth X, as a fraction of the content at maximum.

    Gaisser-Hillas with the first interaction at X_0 = 0:

        N(X)/N_max = (X / X_max)^(X_max/lambda) * exp((X_max - X) / lambda)

    Zero at and below X = 0. This is what makes a particle array's criterion a band
    rather than a threshold: the content rises steeply, peaks, and then dies.

    Parameters
    ----------
    grammage_gcm2 : float or array_like
        Depth at which to evaluate the profile, in g/cm^2.
    x_max_gcm2 : float
        Depth of shower maximum, in g/cm^2.
    lambda_gcm2 : float, optional
        Gaisser-Hillas interaction length, in g/cm^2, setting how fast the profile
        rises and falls.

    Returns
    -------
    ndarray
        Particle content as a fraction of the content at maximum, in [0, 1].

    Examples
    --------
    >>> from oroscope import physics
    >>> round(float(physics.shower_size_fraction(700.0, 700.0)), 3)   # at maximum
    1.0
    >>> round(float(physics.shower_size_fraction(172.0, 561.0)), 3)   # a 2 km crossing
    0.02
    """
    x = np.asarray(grammage_gcm2, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(x > 0.0, x / x_max_gcm2, 0.0)
        out = np.where(
            x > 0.0,
            ratio ** (x_max_gcm2 / lambda_gcm2) * np.exp((x_max_gcm2 - x) / lambda_gcm2),
            0.0)
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def grammage_band_from_energy(energy_min_pev: float, energy_max_pev: float,
                              fraction: float = 0.1,
                              lambda_gcm2: float = SHOWER_PROFILE_LAMBDA_GCM2,
                              samples: int = 4000,
                              **x_max_kw: float) -> tuple[float, float]:
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

    Parameters
    ----------
    energy_min_pev, energy_max_pev : float
        Ends of the primary-energy range, in PeV.
    fraction : float, optional
        Fraction of peak particle content that still counts as a usable shower. A
        choice about detector capability rather than a property of the shower, and one
        of the parameters a result is most sensitive to.
    lambda_gcm2 : float, optional
        Gaisser-Hillas interaction length, in g/cm^2.
    samples : int, optional
        Points on the depth grid searched for the band edges.
    **x_max_kw
        Passed through to :func:`shower_maximum_gcm2`.

    Returns
    -------
    tuple of float
        ``(low_gcm2, high_gcm2)``, the low edge taken at the lowest energy and the
        high edge at the highest, so the band spans the whole requested range.

    Examples
    --------
    >>> from oroscope import physics
    >>> lo, hi = physics.grammage_band_from_energy(3.0, 1000.0)
    >>> f"{lo:.0f} - {hi:.0f} g/cm^2"
    '236 - 1287 g/cm^2'
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

def earth_chord_m(elevation_deg: float, radius_m: float = EARTH_RADIUS_M) -> float:
    """
    Chord length through the Earth for a ray arriving from below the horizontal.

    A ray making angle theta below the local horizontal cuts a chord of ``2R sin(theta)``:
    zero along the tangent, a full diameter straight down. Directions at or above the
    horizontal return zero.

    This matters because it dwarfs local topography. At -1 degree the chord is about
    220 km, at -3 degrees about 670 km, against the tens of km of mountain a DEM can
    see. The deepest part of even a 670 km chord lies only ~9 km below the surface, so
    it stays in the crust and a constant density is adequate.

    Parameters
    ----------
    elevation_deg : float
        Arrival elevation angle, in degrees. Zero or above returns zero.
    radius_m : float, optional
        Earth radius, in metres. The true radius: the neutrino is not refracted.

    Returns
    -------
    float
        Chord length through the Earth, in metres.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.earth_chord_m(-3.0) / 1000:.0f} km"
    '667 km'
    """
    if elevation_deg >= 0:
        return 0.0
    return 2.0 * radius_m * math.sin(math.radians(-elevation_deg))


def earth_chord_gcm2(elevation_deg: float, radius_m: float = EARTH_RADIUS_M,
                     density_gcm3: float = CRUST_DENSITY_GCM3) -> float:
    """
    Column depth of the Earth chord, in g/cm^2.

    Parameters
    ----------
    elevation_deg : float
        Arrival elevation angle, in degrees.
    radius_m : float, optional
        Earth radius, in metres.
    density_gcm3 : float, optional
        Crust density, in g/cm^3.

    Returns
    -------
    float
        Column depth of the chord, in g/cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.earth_chord_gcm2(-3.0):.2e}"
    '1.77e+08'
    """
    return earth_chord_m(elevation_deg, radius_m) * 100.0 * density_gcm3


# Neutral-current cross-section as a fraction of charged-current, and the mean
# inelasticity of an NC interaction. Both are slowly varying at these energies; these
# are round numbers, and both are parameters of nc_regeneration_factor() so a better
# value can be supplied without touching this.
NC_TO_CC_RATIO = 0.42
NC_INELASTICITY = 0.25


def nc_regeneration_factor(chord_gcm2: float, cc_length_gcm2: float,
                           spectral_index: float = 2.0,
                           nc_to_cc: float = NC_TO_CC_RATIO,
                           inelasticity: float = NC_INELASTICITY) -> float:
    """
    Leading-order enhancement from neutral-current regeneration. **An approximation.**

    A charged-current interaction removes a neutrino from the beam. A neutral-current
    one does not -- it degrades the energy and the neutrino continues. Counting only CC
    absorption, as ``neutrino_survival`` does by default, therefore understates the flux
    arriving at a given energy: neutrinos that started higher up the spectrum scatter
    *down into* the band and partially refill it.

    The size of that refilling follows from the spectrum. For a power law
    ``Phi(E) ~ E**-gamma``, a neutrino observed at ``E`` after one NC interaction of
    inelasticity ``y`` began at ``E' = E/(1-y)``. The flux there is larger by
    ``(1-y)**gamma``, and the Jacobian ``dE'/dE = 1/(1-y)`` gives one more power, so each
    NC scatter contributes ``(1-y)**(gamma-1)`` relative to the unscattered flux. With
    ``tau_nc = (chord/X_cc) * (sigma_NC/sigma_CC)`` NC interactions expected along the
    chord, the leading term is::

        1 + tau_nc * (1 - y)**(gamma - 1)

    **What this is not.** It is the first term of a series, not a solution of the
    cascade equations, and it ignores the ``nu_tau -> tau -> nu_tau`` chain that makes
    the Earth genuinely translucent to tau neutrinos at high energy. It is a correction
    of the right sign and roughly the right size, and it is capped below (see
    :func:`neutrino_survival`) so it cannot manufacture more flux than arrived. Treat a
    result that depends strongly on it as a result that needs a real transport code.

    Parameters
    ----------
    chord_gcm2 : float
        Column depth along the chord, in g/cm^2.
    cc_length_gcm2 : float
        Charged-current interaction length, in g/cm^2.
    spectral_index : float, optional
        ``gamma`` of the assumed power-law flux. A steeper spectrum regenerates *less*:
        the neutrinos scattering down into the band come from ``E/(1-y)``, above it, and
        a steeper spectrum has less flux there.
    nc_to_cc : float, optional
        ``sigma_NC / sigma_CC``.
    inelasticity : float, optional
        Mean fraction of energy lost in one NC interaction.

    Returns
    -------
    float
        Multiplicative enhancement, at least 1.

    Examples
    --------
    A chord one CC interaction length deep, on an E^-2 spectrum:

    >>> from oroscope import physics
    >>> round(physics.nc_regeneration_factor(1.0e8, 1.0e8), 3)
    1.315

    A steeper spectrum regenerates less, since less flux sits above the band:

    >>> round(physics.nc_regeneration_factor(1.0e8, 1.0e8, spectral_index=2.7), 3)
    1.258

    A short chord regenerates nothing:

    >>> round(physics.nc_regeneration_factor(0.0, 1.0e8), 3)
    1.0
    """
    if cc_length_gcm2 <= 0 or chord_gcm2 <= 0:
        return 1.0
    tau_nc = (chord_gcm2 / cc_length_gcm2) * nc_to_cc
    return 1.0 + tau_nc * (1.0 - inelasticity) ** (spectral_index - 1.0)


def neutrino_survival(elevation_deg: float, interaction_length_gcm2: float,
                      radius_m: float = EARTH_RADIUS_M,
                      density_gcm3: float = CRUST_DENSITY_GCM3,
                      nc_regeneration: bool = False,
                      spectral_index: float = 2.0) -> float:
    """
    Fraction of neutrinos surviving the Earth chord to reach the exit region.

    ``exp(-X_chord / X_int)``. The interaction length is a parameter because it depends
    on the cross-section at the energy of interest; around an EeV it is of order
    1e8 g/cm^2, which is the same order as the chord at -3 degrees, so the suppression
    across a +/-3 degree window is substantial rather than marginal.

    With ``nc_regeneration=True`` the result is multiplied by
    :func:`nc_regeneration_factor`, and clipped at 1 so that regeneration can offset
    absorption but never invent flux. **Off by default**, because it is an
    approximation and every published number here was computed without it.

    Parameters
    ----------
    elevation_deg : float
        Arrival elevation angle, in degrees.
    interaction_length_gcm2 : float
        Neutrino interaction length at the energy of interest, in g/cm^2. Zero or
        less disables the attenuation and returns 1.
    radius_m : float, optional
        Earth radius, in metres.
    density_gcm3 : float, optional
        Crust density, in g/cm^3.
    nc_regeneration : bool, optional
        Apply the leading-order neutral-current regeneration correction.
    spectral_index : float, optional
        Power-law index used by that correction.

    Returns
    -------
    float
        Surviving fraction, in [0, 1].

    Examples
    --------
    >>> from oroscope import physics
    >>> lam = physics.neutrino_interaction_length_gcm2(1000.0)
    >>> round(physics.neutrino_survival(-1.0, lam), 3)
    0.7

    Regeneration lifts it, because absorption alone was overstating the suppression:

    >>> round(physics.neutrino_survival(-1.0, lam, nc_regeneration=True), 3)
    0.778
    """
    if interaction_length_gcm2 <= 0:
        return 1.0
    chord = earth_chord_gcm2(elevation_deg, radius_m, density_gcm3)
    survival = math.exp(-chord / interaction_length_gcm2)
    if nc_regeneration:
        survival *= nc_regeneration_factor(chord, interaction_length_gcm2,
                                           spectral_index=spectral_index)
    return min(survival, 1.0)


# ------------------------------------------------------- muon shielding

# Rock overburden a detector wants along the arrival direction so that atmospheric
# muons from that direction cannot reach it. Ref. [2] Fig. 1 annotates >4 km for TAMBO.
DEFAULT_MUON_SHIELDING_KM = 4.0


def muon_shielding_gcm2(thickness_km: float,
                        density_gcm3: float = CRUST_DENSITY_GCM3) -> float:
    """
    Column depth corresponding to a rock thickness, for muon rejection.

    A mountain in the arrival direction is a muon filter: an air-shower muon from that
    direction would have to cross the whole thickness, which a few km of rock makes
    impossible. Anything detected from behind that much rock is therefore not a
    cosmic-ray muon, which is what lets a surface array claim neutrino purity.

    Unlike the production-and-escape band, this is a **floor**: more rock is always
    better for background rejection, and only the signal side wants an upper limit.

    4 km of standard rock is about 1.06e6 g/cm^2.

    Parameters
    ----------
    thickness_km : float
        Rock thickness along the arrival direction, in km.
    density_gcm3 : float, optional
        Rock density, in g/cm^3.

    Returns
    -------
    float
        Column depth, in g/cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.muon_shielding_gcm2(4.0):.2e}"
    '1.06e+06'
    """
    return thickness_km * 1000.0 * 100.0 * density_gcm3


# ---------------------------------------------------------------- tau range

def tau_decay_length_m(energy_pev: float, mass_gev: float = 1.77686,
                       ctau_m: float = 87.03e-6) -> float:
    """
    Lorentz-boosted tau decay length, ``(E/m) c*tau``.

    The quantity that decides whether a tau decays inside a canyon crossing or flies
    through it: 147 m at 3 PeV against 49 km at 1 EeV, over a gap of a few km.

    Parameters
    ----------
    energy_pev : float
        Tau energy, in PeV.
    mass_gev : float, optional
        Tau mass, in GeV.
    ctau_m : float, optional
        Proper decay length ``c * tau``, in metres.

    Returns
    -------
    float
        Decay length in the laboratory frame, in metres.

    Examples
    --------
    >>> from oroscope import physics
    >>> [f"{physics.tau_decay_length_m(e):.0f}" for e in (3.0, 1000.0)]
    ['147', '48980']
    """
    return (energy_pev * 1.0e6 / mass_gev) * ctau_m


# Charged-current nu-N cross-section, a power-law fit of the standard
# parameterisations: sigma = A (E/GeV)^n cm^2. Good to tens of per cent over
# 1e8-1e10 GeV and an extrapolation above that, where no data constrain it.
SIGMA_CC_COEFF_CM2 = 6.04e-36
SIGMA_CC_INDEX = 0.358
AVOGADRO = 6.022e23


def cc_cross_section_cm2(energy_pev: float | np.ndarray) -> float | np.ndarray:
    """
    Charged-current neutrino-nucleon cross-section.

    A power-law fit to the standard parameterisations, good to tens of per cent over
    1e8-1e10 GeV and an extrapolation above that, where no data constrain it.

    Parameters
    ----------
    energy_pev : float or array_like
        Neutrino energy, in PeV.

    Returns
    -------
    float or ndarray
        Cross-section, in cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.cc_cross_section_cm2(1000.0):.2e}"
    '1.01e-32'
    """
    return SIGMA_CC_COEFF_CM2 * (energy_pev * 1.0e6) ** SIGMA_CC_INDEX


def neutrino_interaction_length_gcm2(energy_pev: float | np.ndarray
                                     ) -> float | np.ndarray:
    """
    Column depth over which a neutrino interacts once, ``1/(N_A sigma)``.

    Falls from about 3.8e8 g/cm^2 at 100 PeV to 7e7 at 10 EeV. Only charged-current
    attenuation is counted; neutral-current regeneration would soften it slightly.

    Parameters
    ----------
    energy_pev : float or array_like
        Neutrino energy, in PeV.

    Returns
    -------
    float or ndarray
        Interaction length, in g/cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.neutrino_interaction_length_gcm2(100.0):.2e}"
    '3.76e+08'
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

# The live values, which set_tau_energy_loss() replaces. The BETA_* constants above
# stay as the shipped estimate, so restore_tau_energy_loss() has something to restore
# to and a reader can always see what the default was.
_TAU_BETA = {"reference": BETA_REFERENCE_CM2G,
             "reference_energy_pev": BETA_REFERENCE_ENERGY_PEV,
             "index": BETA_ENERGY_INDEX}


def tau_energy_loss_settings() -> dict:
    """
    The beta parameters currently in force.

    Returns
    -------
    dict
        ``{"reference", "reference_energy_pev", "index"}``, in cm^2/g and PeV.

    Examples
    --------
    >>> from oroscope import physics
    >>> physics.tau_energy_loss_settings()["reference"]
    6e-07
    """
    return dict(_TAU_BETA)


def set_tau_energy_loss(reference: float | None = None,
                        reference_energy_pev: float | None = None,
                        index: float | None = None) -> dict:
    """
    Adopts a different beta, in one place, for every function that uses it.

    beta is the least certain number in this module -- an estimate from mass scaling,
    in the range (0.4-1.0)e-6 cm^2/g -- and the collaboration's own value should replace
    it when there is one. Until now that meant editing the source, which is not
    something a user of an installed package can reasonably do, and which leaves no
    record of what was used.

    **This does not change any site-search result.** beta enters tau *range* and
    *survival* -- production and escape through rock -- and the search does not model
    those: it uses the decay length ``L = (E/m_tau) c*tau``, which is kinematics and
    carries no beta. So this affects :func:`tau_range_gcm2`, :func:`tau_survival` and
    :func:`tau_exit_probability`, which the notebooks and the physics page use, and
    nothing in a search. See ROADMAP 6.38.

    Parameters
    ----------
    reference : float, optional
        beta at ``reference_energy_pev``, in cm^2/g. Left unchanged when None.
    reference_energy_pev : float, optional
        Energy at which ``reference`` applies, in PeV.
    index : float, optional
        Power-law index of the energy dependence. Zero gives a constant beta.

    Returns
    -------
    dict
        The settings now in force, as :func:`tau_energy_loss_settings` returns.

    Examples
    --------
    A collaboration value, adopted once:

    >>> from oroscope import physics
    >>> _ = physics.set_tau_energy_loss(reference=0.8e-6, index=0.0)
    >>> f"{physics.tau_energy_loss_beta(100.0):.2e}"
    '8.00e-07'

    A constant beta means the same value at every energy:

    >>> f"{physics.tau_energy_loss_beta(10000.0):.2e}"
    '8.00e-07'

    And back to the shipped estimate, which does rise with energy:

    >>> _ = physics.restore_tau_energy_loss()
    >>> f"{physics.tau_energy_loss_beta(100.0):.2e}"
    '3.79e-07'
    """
    if reference is not None:
        if reference <= 0.0:
            raise ValueError("beta reference must be positive")
        _TAU_BETA["reference"] = float(reference)
    if reference_energy_pev is not None:
        if reference_energy_pev <= 0.0:
            raise ValueError("beta reference energy must be positive")
        _TAU_BETA["reference_energy_pev"] = float(reference_energy_pev)
    if index is not None:
        _TAU_BETA["index"] = float(index)
    return tau_energy_loss_settings()


def restore_tau_energy_loss() -> dict:
    """
    Restores the shipped estimate, undoing :func:`set_tau_energy_loss`.

    Returns
    -------
    dict
        The settings now in force.

    Examples
    --------
    >>> from oroscope import physics
    >>> physics.restore_tau_energy_loss() == {
    ...     "reference": physics.BETA_REFERENCE_CM2G,
    ...     "reference_energy_pev": physics.BETA_REFERENCE_ENERGY_PEV,
    ...     "index": physics.BETA_ENERGY_INDEX}
    True
    """
    _TAU_BETA.update(reference=BETA_REFERENCE_CM2G,
                     reference_energy_pev=BETA_REFERENCE_ENERGY_PEV,
                     index=BETA_ENERGY_INDEX)
    return tau_energy_loss_settings()

# Mean inelasticity of a charged-current interaction at these energies: the tau carries
# away roughly (1 - y) of the neutrino energy, with <y> about 0.2.
CC_INELASTICITY = 0.2


def tau_energy_loss_beta(energy_pev: float, reference: float | None = None,
                         reference_energy_pev: float | None = None,
                         index: float | None = None) -> float:
    """
    Energy-loss coefficient beta(E), rising with energy as photonuclear does.

    An **estimate**, not a fit to published tables: see the module comments for how it
    was arrived at. It is the least certain number in this module, and it moves the
    production-and-escape optimum in proportion.

    Parameters
    ----------
    energy_pev : float
        Tau energy, in PeV.
    reference : float, optional
        Value of beta at ``reference_energy_pev``, in cm^2/g. Defaults to whatever
        :func:`set_tau_energy_loss` last established, and to
        :data:`BETA_REFERENCE_CM2G` if it has not been called.
    reference_energy_pev : float, optional
        Energy at which ``reference`` applies, in PeV.
    index : float, optional
        Power-law index of the energy dependence. Zero recovers a constant beta.

    Returns
    -------
    float
        Energy-loss coefficient, in cm^2/g.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.tau_energy_loss_beta(100.0):.2e}"
    '3.79e-07'

    An argument overrides the module setting for that one call:

    >>> f"{physics.tau_energy_loss_beta(100.0, reference=1.0e-6, index=0.0):.2e}"
    '1.00e-06'
    """
    live = _TAU_BETA
    reference = live["reference"] if reference is None else reference
    reference_energy_pev = (live["reference_energy_pev"]
                            if reference_energy_pev is None else reference_energy_pev)
    index = live["index"] if index is None else index
    if index == 0.0:
        return reference
    return reference * (energy_pev / reference_energy_pev) ** index


def tau_range_gcm2(energy_pev: float, beta_cm2g: float | None = None,
                   density_gcm3: float = CRUST_DENSITY_GCM3) -> float:
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

    Parameters
    ----------
    energy_pev : float
        Tau energy on entering the rock, in PeV.
    beta_cm2g : float, optional
        Energy-loss coefficient, in cm^2/g. Defaults to :func:`tau_energy_loss_beta`
        at this energy.
    density_gcm3 : float, optional
        Rock density, in g/cm^3.

    Returns
    -------
    float
        Column depth at which survival falls to 1/e, in g/cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.tau_range_gcm2(1000.0) / physics.CRUST_DENSITY_GCM3 / 1e5:.1f} km"
    '13.7 km'
    """
    beta = tau_energy_loss_beta(energy_pev) if beta_cm2g is None else beta_cm2g
    x_decay = tau_decay_length_m(energy_pev) * 100.0 * density_gcm3
    x_loss = 1.0 / beta
    return x_loss * math.log1p(x_decay / x_loss)


def tau_survival(depth_gcm2: float | np.ndarray, energy_pev: float,
                 beta_cm2g: float | None = None,
                 density_gcm3: float = CRUST_DENSITY_GCM3) -> np.ndarray:
    """
    Probability a tau of the given energy crosses ``depth_gcm2`` of rock without decaying.

    The double-exponential form above, which falls far more sharply than a simple
    exponential once the depth exceeds ``1/beta``: the tau is losing energy, so its
    decay length shrinks as it goes.

    Parameters
    ----------
    depth_gcm2 : float or array_like
        Rock traversed, in g/cm^2.
    energy_pev : float
        Tau energy on entering the rock, in PeV.
    beta_cm2g : float, optional
        Energy-loss coefficient, in cm^2/g.
    density_gcm3 : float, optional
        Rock density, in g/cm^3.

    Returns
    -------
    ndarray
        Survival probability, in [0, 1].

    Examples
    --------
    >>> from oroscope import physics
    >>> round(float(physics.tau_survival(0.0, 1000.0)), 3)
    1.0
    """
    beta = tau_energy_loss_beta(energy_pev) if beta_cm2g is None else beta_cm2g
    x_decay = tau_decay_length_m(energy_pev) * 100.0 * density_gcm3
    x_loss = 1.0 / beta
    d = np.asarray(depth_gcm2, dtype=np.float64)
    # exp overflows for very thick rock, where survival is zero anyway
    return np.exp(-(x_loss / x_decay) * np.expm1(np.clip(d / x_loss, 0.0, 700.0)))


def tau_exit_probability(column_depth_gcm2: float | np.ndarray, energy_pev: float,
                         beta_cm2g: float | None = None,
                         density_gcm3: float = CRUST_DENSITY_GCM3,
                         inelasticity: float = CC_INELASTICITY,
                         samples: int = 2000) -> float | np.ndarray:
    """
    Relative probability that a traversing neutrino yields a tau that escapes.

    The neutrino interacts at depth x with probability ``exp(-x/lambda) dx/lambda``;
    the tau, carrying ``(1 - y)`` of the energy, must then cross the remaining
    ``X - x`` and survive:

        P(X) = integral_0^X (dx/lambda) exp(-x/lambda) S(X - x)

    Evaluated numerically, because ``S`` is a double exponential and the integral has
    no clean closed form. Relative, not absolute: normalisation needs the trigger
    response (see :mod:`aperture`).

    Parameters
    ----------
    column_depth_gcm2 : float or array_like
        Rock along the arrival direction, in g/cm^2.
    energy_pev : float
        Neutrino energy, in PeV.
    beta_cm2g : float, optional
        Tau energy-loss coefficient, in cm^2/g.
    density_gcm3 : float, optional
        Rock density, in g/cm^3.
    inelasticity : float, optional
        Mean inelasticity of the charged-current interaction; the tau carries away
        ``1 - y`` of the neutrino energy.
    samples : int, optional
        Points used for the numerical integration over interaction depth, spaced
        logarithmically in the *remaining* depth ``X - x``.

    Returns
    -------
    float or ndarray
        Relative exit probability, matching the shape of ``column_depth_gcm2``.

    Notes
    -----
    **The grid is in ``u = X - x``, not in ``x``, and that is the whole of the
    accuracy.** Only interactions within a few tau ranges of the far surface contribute
    anything -- everything produced deeper is absorbed -- so the integrand is a spike
    against the far end of a range that can be five decades wide. Sampled uniformly in
    ``x``, as this was, the spacing outruns the spike and the trapezoid rule reports the
    area of something it never resolved. Measured at 3 PeV and ``X`` = 10^9 g/cm^2:

    ====================  ==============
    grid                  P(X)
    ====================  ==============
    uniform, 2000 pts     8.884e-05
    uniform, 20,000       1.328e-05
    uniform, 200,000      1.103e-05
    uniform, 2,000,000    1.100e-05
    **log in u, 2000**    **1.1004e-05**
    ====================  ==============

    So the old default was **8x** the converged value, and the substitution reaches that
    value with a thousandth of the points. Worse than the magnitude, it inverted the
    trend: uniform, the exit probability *rose* with depth, giving a spurious maximum at
    the edge of the grid; converged, it falls monotonically, as it must, since more rock
    can only absorb more.

    The error was confined to low energy against deep rock -- at 100 PeV and above the
    worst case over the same grid was 3%, and below 10^7 g/cm^2 it was exact everywhere
    -- which is why nothing else showed it. But
    :func:`depth_band_from_energy` takes its low edge at the *lowest* energy asked for,
    and TAMBO's configured range starts at 3 PeV, so the band inherited the whole of it:
    (1.18e8, 2.89e8) over 3 PeV - 1 EeV where the converged answer is (2.17e4, 1.15e8).
    The published 1 EeV optimum of 5.7e6 g/cm^2 lies inside the second and 20x below the
    first. No published number moved: every config leaves ``depth_band_gcm2`` null, so
    no search has ever called this (roadmap 6.44, 6.56).

    Examples
    --------
    >>> from oroscope import physics
    >>> p = physics.tau_exit_probability(1.0e6, 1000.0)
    >>> 0.0 <= p <= 1.0
    True

    Converged where it used to be eight times too large:

    >>> round(physics.tau_exit_probability(1.0e9, 3.0) * 1.0e5, 3)
    1.1

    and falling with depth, as more rock must:

    >>> p = physics.tau_exit_probability([1.0e8, 3.0e8, 1.0e9, 3.0e9], 3.0)
    >>> bool((p[1:] < p[:-1]).all())
    True
    """
    e_tau = energy_pev * (1.0 - inelasticity)
    lam = neutrino_interaction_length_gcm2(energy_pev)

    def one(X):
        if X <= 0:
            return 0.0
        X = float(X)
        # Substituted u = X - x: u is the rock the tau still has to cross, and S(u) is
        # what kills the integrand. Sampled logarithmically in u so the points land
        # where the weight is, which is within a few tau ranges of the far surface.
        # See the notes: uniform in x, this needed a thousand times as many points.
        u = np.concatenate(([0.0], np.geomspace(max(X * 1.0e-12, 1.0e-9), X, samples)))
        s = tau_survival(u, e_tau, beta_cm2g, density_gcm3)
        return float(_trapezoid(np.exp(-(X - u) / lam) / lam * s, u))

    if np.ndim(column_depth_gcm2) == 0:
        return one(column_depth_gcm2)
    return np.array([one(v) for v in np.asarray(column_depth_gcm2)])


def production_escape_optimum_gcm2(energy_pev: float, beta_cm2g: float | None = None,
                                   density_gcm3: float = CRUST_DENSITY_GCM3,
                                   inelasticity: float = CC_INELASTICITY,
                                   samples: int = 400) -> float:
    """
    Column depth maximising :func:`tau_exit_probability`.

    Found on a log grid, since the corrected exit probability has no closed-form peak.
    Rises with energy and then flattens -- about 12 km of standard rock at 100 PeV,
    20 km at 1 EeV, 23 km at 10 EeV -- as the logarithmic growth of the tau range is
    tempered by beta rising.

    Parameters
    ----------
    energy_pev : float
        Neutrino energy, in PeV.
    beta_cm2g : float, optional
        Tau energy-loss coefficient, in cm^2/g.
    density_gcm3 : float, optional
        Rock density, in g/cm^3.
    inelasticity : float, optional
        Mean charged-current inelasticity.
    samples : int, optional
        Points on the log-spaced depth grid searched for the peak.

    Returns
    -------
    float
        Column depth maximising the exit probability, in g/cm^2.

    Examples
    --------
    >>> from oroscope import physics
    >>> km = physics.production_escape_optimum_gcm2(100.0) / physics.CRUST_DENSITY_GCM3 / 1e5
    >>> 8.0 < km < 16.0
    True
    """
    grid = np.logspace(4.0, 9.0, samples)
    p = tau_exit_probability(grid, energy_pev, beta_cm2g, density_gcm3, inelasticity)
    return float(grid[int(np.argmax(p))])


def depth_band_from_energy(energy_min_pev: float, energy_max_pev: float,
                           fraction: float = 0.5, beta_cm2g: float | None = None,
                           density_gcm3: float = CRUST_DENSITY_GCM3,
                           inelasticity: float = CC_INELASTICITY,
                           samples: int = 400) -> tuple[float, float]:
    """
    Column-depth band where the tau exit probability stays above ``fraction`` of peak.

    Very wide -- roughly 5e5 to 2.6e8 g/cm^2 at half maximum across 100 PeV to 10 EeV,
    some two and a half decades -- which is itself the result: column depth is an
    intrinsically weak discriminant and the criterion should not pretend otherwise.

    Parameters
    ----------
    energy_min_pev, energy_max_pev : float
        Ends of the neutrino energy range, in PeV.
    fraction : float, optional
        Fraction of peak exit probability defining the band edges.
    beta_cm2g : float, optional
        Tau energy-loss coefficient, in cm^2/g.
    density_gcm3 : float, optional
        Rock density, in g/cm^3.
    inelasticity : float, optional
        Mean charged-current inelasticity.
    samples : int, optional
        Points on the log-spaced depth grid.

    Returns
    -------
    tuple of float
        ``(low_gcm2, high_gcm2)``, the low edge taken at the lowest energy and the high
        edge at the highest, so the band spans the whole requested range.

    Examples
    --------
    >>> from oroscope import physics
    >>> lo, hi = physics.depth_band_from_energy(100.0, 10000.0)
    >>> lo < hi
    True
    """
    grid = np.logspace(4.0, 9.0, samples)

    def edges(energy):
        p = tau_exit_probability(grid, energy, beta_cm2g, density_gcm3, inelasticity)
        ok = grid[p >= fraction * p.max()]
        return float(ok.min()), float(ok.max())

    lo = edges(energy_min_pev)[0]
    hi = edges(energy_max_pev)[1]
    return (lo, hi) if lo <= hi else (hi, lo)


def earth_absorption_cutoff_deg(energy_pev: float, fraction: float = 0.5,
                                radius_m: float = EARTH_RADIUS_M,
                                density_gcm3: float = CRUST_DENSITY_GCM3,
                                **kw: float) -> float | None:
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

    Parameters
    ----------
    energy_pev : float
        Neutrino energy, in PeV.
    fraction : float, optional
        Fraction of peak acceptance defining the cut.
    radius_m : float, optional
        Earth radius, in metres.
    density_gcm3 : float, optional
        Crust density, in g/cm^3.
    **kw
        Passed through to :func:`tau_exit_probability`.

    Returns
    -------
    float or None
        A negative elevation angle in degrees, or ``None`` when the cut lies below the
        horizon entirely.

    Examples
    --------
    >>> from oroscope import physics
    >>> cut = physics.earth_absorption_cutoff_deg(1000.0)
    >>> -6.0 < cut < 0.0
    True
    """
    X = np.logspace(4, 9, 400)
    P = tau_exit_probability(X, energy_pev, **kw)
    upper = float(X[P >= fraction * P.max()].max())
    chord_m = upper / density_gcm3 / 100.0
    sin_theta = chord_m / (2.0 * radius_m)
    if sin_theta >= 1.0:
        return None
    return -math.degrees(math.asin(sin_theta))


# Tau flux spectral index. A power law dN/dE ~ E^-gamma; 2 is the canonical
# cosmogenic/astrophysical slope, and the value is a knob because it is an input to the
# search rather than a property of the terrain.
DEFAULT_SPECTRAL_INDEX = 2.0


DECAY_WEIGHTINGS = ("flux", "acceptance", "flux_times_acceptance")


def spectrum_weighted_decay_probability(distance_m, energy_min_pev, energy_max_pev,
                                        spectral_index=DEFAULT_SPECTRAL_INDEX,
                                        shower_development_m=0.0, samples=96,
                                        index_samples=129,
                                        weight_by="flux", response=None):
    r"""
    Probability the tau decays inside the usable gap, folded over a power-law spectrum.

    A tau of energy :math:`E` decays within a usable path :math:`u` with probability
    :math:`1 - e^{-u/L(E)}`, and :math:`L = (E/m_\tau)c\tau` runs over three decades
    across a single experiment's reach. Evaluating that at one representative energy is
    therefore not an approximation but a choice of answer: measured on a real canyon
    search, the reported capacity ran from 10878 detector positions at 3 PeV to zero at
    100 PeV. Weighting by the flux instead gives a number that is a property of the
    terrain and the spectrum rather than of the energy someone picked.

    .. math::

        P(u) = \frac{\int E^{-\gamma}\left(1 - e^{-u/L(E)}\right)\,{\rm d}E}
                    {\int E^{-\gamma}\,{\rm d}E}

    integrated on a log-spaced grid, since the range spans decades.

    **What weights the average is selectable.** An event rate is
    :math:`\int \Phi(E) A(E) P(E)\,{\rm d}E`, and weighting by the flux alone is only
    one of three defensible choices:

    ``"flux"``
        :math:`w(E) = E^{-\gamma}`. The default, and what every published number here
        was computed with. Says: *of the neutrinos that arrive, what fraction decays
        usefully?*
    ``"acceptance"``
        :math:`w(E) = A(E)`. Says: *over the energies this detector actually responds
        to, what fraction decays usefully?* -- a property of the instrument and the
        terrain, with no assumed spectrum. Useful precisely because :math:`\gamma` is an
        assumption and this removes it.
    ``"flux_times_acceptance"``
        :math:`w(E) = E^{-\gamma} A(E)`, the event-rate integrand itself, and the right
        one when both the spectrum and a response table are trusted.

    ``A(E)`` is not something this module can supply; pass a callable, most usefully one
    recovered from a published curve by :func:`aperture.infer_response`, which divides
    that curve by the geometric model and leaves everything else.

    Note a steep spectrum weights low energies heavily, where the tau decays readily --
    so a soft spectrum drives :math:`P` toward 1 and the term stops discriminating. That
    is the physics rather than a defect, but it means the spectral index deserves the
    same scrutiny as any other assumption, and it is the reason ``"acceptance"`` exists.

    Parameters
    ----------
    distance_m : float or array_like
        Distance from the exit point to the detector, in metres.
    energy_min_pev, energy_max_pev : float
        Ends of the tau energy range, in PeV.
    spectral_index : float, optional
        :math:`\gamma` in :math:`{\rm d}N/{\rm d}E \propto E^{-\gamma}`.
    shower_development_m : float, optional
        Path the shower needs after the decay, in metres. Subtracted from the gap, so
        a target closer than this yields nothing usable.
    samples : int, optional
        Points on the log-spaced energy grid.
    index_samples : int, optional
        Points across the spectral-index range, when one is given. Generous by default
        because the index integral is folded into a weight vector computed once, so a
        fine grid costs nothing per candidate. Ignored for a single index.
    weight_by : {'flux', 'acceptance', 'flux_times_acceptance'}, optional
        What weights the average. See the discussion above. ``'acceptance'`` ignores
        ``spectral_index`` entirely, which is the point of it.
    response : callable, optional
        ``A(E)`` in PeV, required by the two acceptance weightings. A response that is
        zero everywhere on the grid leaves nothing to average and raises.

    Returns
    -------
    ndarray
        Flux-weighted decay probability, in [0, 1], matching the shape of
        ``distance_m``.

    Raises
    ------
    ValueError
        If the energy range is inverted, or ``spectral_index`` is neither one value
        nor a pair.

    See Also
    --------
    arrival_scan.decay_probability : the single-energy form, for one baseline window.

    Examples
    --------
    >>> from oroscope import physics
    >>> p = physics.spectrum_weighted_decay_probability(3000.0, 3.0, 1000.0)
    >>> round(float(p), 3)
    0.954

    A harder spectrum puts more weight at high energy, where the tau outruns the gap:

    >>> soft = physics.spectrum_weighted_decay_probability(3000.0, 3.0, 1000.0, 2.7)
    >>> hard = physics.spectrum_weighted_decay_probability(3000.0, 3.0, 1000.0, 1.5)
    >>> bool(hard < soft)
    True

    Or marginalise over the index rather than choosing one, which lands between the
    extremes it spans:

    >>> spread = physics.spectrum_weighted_decay_probability(
    ...     3000.0, 3.0, 1000.0, (1.5, 2.7))
    >>> bool(hard < spread < soft)
    True
    """
    if energy_max_pev <= energy_min_pev:
        raise ValueError("energy_max_pev must exceed energy_min_pev")
    if weight_by not in DECAY_WEIGHTINGS:
        raise ValueError(f"weight_by must be one of {DECAY_WEIGHTINGS}, got {weight_by!r}")
    if weight_by != "flux" and response is None:
        raise ValueError(f"weight_by={weight_by!r} needs a response A(E); pass one, "
                         f"e.g. from aperture.infer_response()")

    # Weighting by acceptance alone is the flux weighting with no spectrum at all, so
    # gamma drops out rather than being quietly retained.
    if weight_by == "acceptance":
        spectral_index = 0.0

    gammas = np.atleast_1d(np.asarray(spectral_index, dtype=np.float64))
    if gammas.size == 2:
        lo, hi = float(gammas.min()), float(gammas.max())
        # Marginalised uniformly over the range: no basis for preferring a value
        # inside it, and a flat prior is the honest way to say so.
        gammas = np.linspace(lo, hi, int(index_samples)) if hi > lo else np.array([lo])
    elif gammas.size != 1:
        raise ValueError("spectral_index must be one value or a (low, high) pair")

    d = np.asarray(distance_m, dtype=np.float64)
    usable = np.clip(d - float(shower_development_m), 0.0, None)

    # Log-spaced, and integrated in ln E: with x = ln E, E^-gamma dE = E^(1-gamma) dx,
    # which keeps the quadrature well conditioned over three decades.
    ln_e = np.linspace(np.log(energy_min_pev), np.log(energy_max_pev), int(samples))
    energies = np.exp(ln_e)
    lengths = tau_decay_length_m(energies)

    # One effective weight, computed once and independent of the candidates.
    #
    # The decay factor depends on the energy but not on the index, so the index
    # integral can be moved inside:
    #
    #     (1/dg) * INT dg [ INT dE decay(E) w(E,g) / Z(g) ]
    #         = INT dE decay(E) * [ (1/dg) INT dg w(E,g)/Z(g) ]
    #
    # so marginalising costs a weight vector rather than a re-integration per index.
    # Done the obvious way it was 45 times the work for a 45-point index grid, which on
    # a real search was 25 s against 1 s; this way a fine grid is free.
    weights = energies[None, :] ** (1.0 - gammas[:, None])
    if response is not None and weight_by != "flux":
        # A(E) multiplies the weight before normalisation, so the average is taken over
        # the energies the detector actually responds to. Outside a tabulated range
        # TabulatedResponse returns zero, which correctly removes those energies rather
        # than extrapolating over them.
        acceptance = np.asarray(response(energies), dtype=np.float64)
        if not np.any(acceptance > 0):
            raise ValueError("the response is zero across the whole energy range, so "
                             "there is nothing to average; check its tabulated range "
                             "against energy_min_pev/energy_max_pev")
        weights = weights * acceptance[None, :]
    normalised = weights / _trapezoid(weights, ln_e, axis=-1)[:, None]
    if gammas.size == 1:
        effective = normalised[0]
    else:
        effective = _trapezoid(normalised, gammas, axis=0) / (gammas[-1] - gammas[0])

    # Chunked over candidates rather than broadcast in one go. A real search carries
    # hundreds of thousands of them, and the full (candidates x energies) outer product
    # would be hundreds of megabytes for no gain; this bounds it to a few tens.
    flat = np.atleast_1d(usable).reshape(-1)
    out = np.empty(flat.size, dtype=np.float64)
    chunk = 65536
    for i in range(0, flat.size, chunk):
        block = flat[i:i + chunk, None]
        decay = -np.expm1(-block / lengths)
        out[i:i + chunk] = _trapezoid(decay * effective, ln_e, axis=-1)

    result = np.clip(out, 0.0, 1.0)
    return result.reshape(usable.shape) if usable.ndim else result[0]


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


def geomagnetic_latitude_deg(latitude_deg: float, longitude_deg: float,
                             pole_lat_deg=GEOMAGNETIC_POLE_LAT_DEG,
                             pole_lon_deg=GEOMAGNETIC_POLE_LON_DEG):
    """
    Latitude in the centered-dipole frame.

    Parameters
    ----------
    latitude_deg, longitude_deg : float
        Geographic coordinates of the site, in degrees.
    pole_lat_deg, pole_lon_deg : float, optional
        Geographic coordinates of the north geomagnetic pole, in degrees.

    Returns
    -------
    float
        Magnetic latitude, in degrees.

    Examples
    --------
    >>> from oroscope import physics
    >>> round(physics.geomagnetic_latitude_deg(-16.4, -71.5), 1)
    -7.1
    """
    lat = math.radians(latitude_deg)
    plat = math.radians(pole_lat_deg)
    dlon = math.radians(longitude_deg - pole_lon_deg)
    sin_mlat = (math.sin(lat) * math.sin(plat)
                + math.cos(lat) * math.cos(plat) * math.cos(dlon))
    return math.degrees(math.asin(max(-1.0, min(1.0, sin_mlat))))


def centered_dipole_inclination(latitude_deg: float, longitude_deg: float,
                                **kw: float) -> float:
    """
    Inclination from a centered dipole, ``tan(I) = 2 tan(magnetic latitude)``.

    Useful when no IGRF lookup is available: it captures the dominant behaviour, that
    inclination passes through zero at the magnetic equator and steepens away from it.
    At Arequipa it gives about -14 degrees, the site being roughly 7 degrees south of
    the magnetic equator.

    Note this approximation is only worth using for *inclination*. The dipole
    declination at Arequipa is about -0.2 degrees against an IGRF value of -6.9, so the
    non-dipole terms dominate there and a dipole declination would be misleading.

    Parameters
    ----------
    latitude_deg, longitude_deg : float
        Geographic coordinates of the site, in degrees.
    **kw
        Passed through to :func:`geomagnetic_latitude_deg`.

    Returns
    -------
    float
        Inclination in degrees, positive downward.

    Examples
    --------
    >>> from oroscope import physics
    >>> round(physics.centered_dipole_inclination(-16.4, -71.5), 1)
    -14.0
    """
    mlat = math.radians(geomagnetic_latitude_deg(latitude_deg, longitude_deg, **kw))
    return math.degrees(math.atan(2.0 * math.tan(mlat)))


_DECLINATION_MODEL = None


def set_declination_model(model):
    """
    Supplies a declination model, so declination can follow the site as inclination does.

    Inclination is computed from the site's coordinates with a centred dipole and is
    good enough. Declination is not: the dipole gives about -0.2 degrees at Arequipa
    against an IGRF -6.9, because the non-dipole terms dominate it. So declination has
    always fallen back to a single constant -- Arequipa's IGRF value -- wherever the DEM
    happened to be, which is right for the prototype region and wrong anywhere else.

    **No IGRF model is shipped, deliberately.** IGRF is a spherical-harmonic expansion
    with a couple of hundred coefficients per epoch, and coefficients typed from memory
    would produce declinations that look entirely plausible and are wrong -- the exact
    failure this project keeps finding. What is provided instead is the socket: pass any
    callable ``model(latitude_deg, longitude_deg) -> declination_deg``, from
    ``ppigrf``/``pyIGRF``, from a NOAA grid via :func:`declination_from_grid`, or from
    the collaboration's own numbers.

    Parameters
    ----------
    model : callable or None
        ``model(latitude_deg, longitude_deg)`` returning degrees east of north.
        ``None`` restores the constant fallback.

    Returns
    -------
    callable or None
        The model now in force.

    Examples
    --------
    >>> from oroscope import physics
    >>> _ = physics.set_declination_model(lambda lat, lon: -6.0 + 0.1 * (lon + 71.5))
    >>> dec, inc = physics.default_field_for_site(-16.4, -71.5)
    >>> f"{dec:.2f}"
    '-6.00'

    It follows the site, which the constant never did:

    >>> dec_east, _ = physics.default_field_for_site(-16.4, -66.5)
    >>> f"{dec_east:.2f}"
    '-5.50'
    >>> _ = physics.set_declination_model(None)
    >>> physics.default_field_for_site(-16.4, -66.5)[0]
    -6.9
    """
    global _DECLINATION_MODEL
    _DECLINATION_MODEL = model
    return _DECLINATION_MODEL


def declination_model() -> object:
    """
    The declination model in force, or ``None`` when the constant fallback is in use.

    Examples
    --------
    >>> from oroscope import physics
    >>> physics.declination_model() is None
    True
    """
    return _DECLINATION_MODEL


def declination_from_grid(latitudes, longitudes, declinations):
    """
    Builds a declination model by bilinear interpolation over a supplied grid.

    The practical way to get real IGRF values in without a spherical-harmonic
    implementation: export a declination grid covering the DEM (NOAA's geomagnetic
    calculator will produce one), and hand the three arrays here.

    Parameters
    ----------
    latitudes : array_like
        Grid latitudes, ascending, in degrees.
    longitudes : array_like
        Grid longitudes, ascending, in degrees.
    declinations : array_like
        ``(len(latitudes), len(longitudes))`` of declination in degrees east of north.

    Returns
    -------
    callable
        ``model(latitude_deg, longitude_deg) -> float``, suitable for
        :func:`set_declination_model`. Queries outside the grid are clamped to its edge,
        which is the right behaviour for a DEM that pokes slightly past its corners and
        a poor one for a grid that does not cover the region at all -- so cover it.

    Examples
    --------
    >>> import numpy as np
    >>> from oroscope import physics
    >>> lats, lons = np.array([-18.0, -14.0]), np.array([-74.0, -70.0])
    >>> dec = np.array([[-5.0, -7.0], [-6.0, -8.0]])
    >>> model = physics.declination_from_grid(lats, lons, dec)
    >>> f"{model(-18.0, -74.0):.2f}"
    '-5.00'

    Halfway in both directions is the average of the four corners:

    >>> f"{model(-16.0, -72.0):.2f}"
    '-6.50'
    """
    lats = np.asarray(latitudes, dtype=float)
    lons = np.asarray(longitudes, dtype=float)
    grid = np.asarray(declinations, dtype=float)
    if grid.shape != (lats.size, lons.size):
        raise ValueError(f"declinations must be {(lats.size, lons.size)}, got {grid.shape}")

    def model(latitude_deg, longitude_deg):
        lat = min(max(float(latitude_deg), lats[0]), lats[-1])
        lon = min(max(float(longitude_deg), lons[0]), lons[-1])
        i = int(np.clip(np.searchsorted(lats, lat) - 1, 0, lats.size - 2))
        j = int(np.clip(np.searchsorted(lons, lon) - 1, 0, lons.size - 2))
        dy = (lat - lats[i]) / (lats[i + 1] - lats[i]) if lats[i + 1] != lats[i] else 0.0
        dx = (lon - lons[j]) / (lons[j + 1] - lons[j]) if lons[j + 1] != lons[j] else 0.0
        return float(grid[i, j] * (1 - dy) * (1 - dx) + grid[i + 1, j] * dy * (1 - dx)
                     + grid[i, j + 1] * (1 - dy) * dx + grid[i + 1, j + 1] * dy * dx)

    return model


def default_field_for_site(latitude_deg: float, longitude_deg: float,
                           declination_deg: float | None = None,
                           inclination_deg: float | None = None
                           ) -> tuple[float, float]:
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

    Parameters
    ----------
    latitude_deg, longitude_deg : float
        Geographic coordinates of the site, in degrees.
    declination_deg : float, optional
        Declination in degrees east of north. Supply the IGRF value for the site; the
        fallback is Arequipa's.
    inclination_deg : float, optional
        Inclination in degrees, positive downward. Derived from the site's coordinates
        when omitted.

    Returns
    -------
    tuple of float
        ``(declination_deg, inclination_deg)``.

    Examples
    --------
    >>> from oroscope import physics
    >>> dec, inc = physics.default_field_for_site(-16.4, -71.5)
    >>> f"{dec:.1f} {inc:.1f}"
    '-6.9 -14.0'
    """
    if declination_deg is None and _DECLINATION_MODEL is not None:
        declination_deg = _DECLINATION_MODEL(latitude_deg, longitude_deg)
    if declination_deg is None:
        declination_deg = DEFAULT_GEOMAG_DECLINATION_DEG
    if inclination_deg is None:
        inclination_deg = centered_dipole_inclination(latitude_deg, longitude_deg)
    return float(declination_deg), float(inclination_deg)


def geomagnetic_unit_vector(declination_deg: float,
                            inclination_deg: float) -> tuple[float, float, float]:
    """
    Unit geomagnetic field in local East-North-Up coordinates.

    Declination is measured from geographic north, positive eastward; inclination is
    positive downward, so it enters the Up component with a minus sign.

    Parameters
    ----------
    declination_deg : float
        Declination, in degrees east of north.
    inclination_deg : float
        Inclination, in degrees, positive downward.

    Returns
    -------
    tuple of float
        ``(east, north, up)`` components of the unit field vector.

    Examples
    --------
    >>> from oroscope import physics
    >>> e, n, u = physics.geomagnetic_unit_vector(0.0, 0.0)
    >>> f"{e:.1f} {n:.1f} {u:.1f}"                  # due north, horizontal
    '0.0 1.0 -0.0'
    """
    d = math.radians(declination_deg)
    i = math.radians(inclination_deg)
    return (math.cos(i) * math.sin(d),      # east
            math.cos(i) * math.cos(d),      # north
            -math.sin(i))                   # up


def geomagnetic_sin_alpha(azimuth_deg: float, elevation_deg: float,
                          field_unit_vector: tuple[float, float, float]) -> float:
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

    Parameters
    ----------
    azimuth_deg : float
        Shower azimuth, in degrees clockwise from north.
    elevation_deg : float
        Shower elevation angle, in degrees.
    field_unit_vector : tuple of float
        Field direction as ``(east, north, up)``, from
        :func:`geomagnetic_unit_vector`.

    Returns
    -------
    float
        ``sin(alpha)``, in [0, 1]. Zero for a shower along the field.

    Examples
    --------
    >>> from oroscope import physics
    >>> B = physics.geomagnetic_unit_vector(0.0, 0.0)      # horizontal, northward
    >>> round(physics.geomagnetic_sin_alpha(0.0, 0.0, B), 3)     # along the field
    0.0
    >>> round(physics.geomagnetic_sin_alpha(90.0, 0.0, B), 3)    # across it
    1.0
    """
    phi = math.radians(azimuth_deg)
    theta = math.radians(elevation_deg)
    cos_t = math.cos(theta)
    v = (math.sin(phi) * cos_t, math.cos(phi) * cos_t, math.sin(theta))
    bx, by, bz = field_unit_vector
    dot = v[0] * bx + v[1] * by + v[2] * bz
    return math.sqrt(max(0.0, 1.0 - dot * dot))


# ---------------------------------------------------------------- footprint

def refractivity(altitude_m: float, sea_level_value: float = SEA_LEVEL_REFRACTIVITY,
                 scale_height_m: float = DENSITY_SCALE_HEIGHT_M) -> float:
    """
    ``n - 1`` at altitude, falling with density.

    Parameters
    ----------
    altitude_m : float
        Altitude above sea level, in metres.
    sea_level_value : float, optional
        Refractivity at sea level.
    scale_height_m : float, optional
        Density scale height, in metres.

    Returns
    -------
    float
        Refractivity ``n - 1``.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.refractivity(4000.0):.2e}"
    '1.80e-04'
    """
    return sea_level_value * math.exp(-altitude_m / scale_height_m)


def cherenkov_angle_rad(altitude_m: float, **kw: float) -> float:
    """
    Cherenkov angle in air, ``sqrt(2(n-1))`` for small angles.

    About 1.4 degrees at sea level and 1.1 degrees at 4000 m: the cone narrows with
    altitude because the air is thinner.

    Parameters
    ----------
    altitude_m : float
        Altitude of the emission point, in metres.
    **kw
        Passed through to :func:`refractivity`.

    Returns
    -------
    float
        Cherenkov angle, in radians.

    Examples
    --------
    >>> import math
    >>> from oroscope import physics
    >>> f"{math.degrees(physics.cherenkov_angle_rad(4000.0)):.2f} deg"
    '1.09 deg'
    """
    return math.sqrt(2.0 * refractivity(altitude_m, **kw))


def cherenkov_footprint_radius_m(altitude_m: float, distance_m: float,
                                 **kw: float) -> float:
    """
    Radius of the radio footprint on the ground, ``D * theta_C``.

    The consequence for layout is counter-intuitive: a *higher* site has a *smaller*
    footprint, so it needs a *denser* array for the same trigger efficiency. A 1 km
    grid under-samples a footprint of a few hundred metres either way, which is why
    counted antennas are a cost proxy rather than an effective area.

    Parameters
    ----------
    altitude_m : float
        Altitude of the emission point, in metres.
    distance_m : float
        Distance from emission point to the ground, in metres.
    **kw
        Passed through to :func:`refractivity`.

    Returns
    -------
    float
        Footprint radius, in metres.

    Examples
    --------
    >>> from oroscope import physics
    >>> f"{physics.cherenkov_footprint_radius_m(4000.0, 20000.0):.0f} m"
    '380 m'
    """
    return distance_m * cherenkov_angle_rad(altitude_m, **kw)


def footprint_sampling(spacing_m: float, altitude_m: float, distance_m: float,
                       **kw: float) -> float:
    """
    Antennas per footprint diameter: ``2 r / spacing``.

    Below 1 the array does not resolve the footprint and triggering relies on a single
    antenna happening to fall inside it.

    Parameters
    ----------
    spacing_m : float
        Detector spacing, in metres. Zero or less returns 0.
    altitude_m : float
        Altitude of the emission point, in metres.
    distance_m : float
        Distance from emission point to the ground, in metres.
    **kw
        Passed through to :func:`refractivity`.

    Returns
    -------
    float
        Antennas spanning one footprint diameter.

    Examples
    --------
    >>> from oroscope import physics
    >>> round(physics.footprint_sampling(1000.0, 4000.0, 20000.0), 2)
    0.76
    """
    radius = cherenkov_footprint_radius_m(altitude_m, distance_m, **kw)
    return (2.0 * radius / spacing_m) if spacing_m > 0 else 0.0
