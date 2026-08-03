"""
Aperture estimation, and validation against what is independently known.

The tool's own geometry gives a site's usable area and the solid angle over which it
sees rock. Turning that into an aperture in m^2 sr needs a detection response, and the
part of that response which depends on tau production and escape is exactly what no
available table supplies (see the roadmap, section 4.10). This module therefore
separates cleanly into:

  * the **geometric aperture**, computed here and fully determined by terrain,
  * the **analytic decay factor**, also computed here and free of free parameters,
  * a **pluggable response**, defaulting to unity and replaceable by a table.

Under that split the absolute normalisation is unknown but the *shape* in energy and
the *ranking* between sites are not. Both are testable, and are tested.

On validating against published apertures. Ref. [1] Fig. 25 and ref. [2] Fig. 3 are
integral quantities over a whole array, all geometries and one site, so they cannot be
applied per pixel. They can anchor the normalisation once supplied as data. Reading
numbers off a published figure by eye is not a measurement, so this module provides
the machinery to compare against a supplied curve rather than a transcription of one.
What *is* validated here are the physical invariants the estimate must satisfy
regardless of normalisation.
"""

import numpy as np

import arrival_scan


def unit_response(energy_pev):
    """Default detection response: energy-independent and equal to one."""
    return np.ones_like(np.asarray(energy_pev, dtype=np.float64))


class TabulatedResponse:
    """
    Detection response interpolated from a supplied table.

    Log-log interpolation, since both response and energy span decades. Supply a
    two-column CSV of energy in PeV against relative response, or the arrays directly.
    """

    def __init__(self, energy_pev, response):
        energy = np.asarray(energy_pev, dtype=np.float64)
        resp = np.asarray(response, dtype=np.float64)
        if energy.ndim != 1 or energy.shape != resp.shape:
            raise ValueError("energy and response must be matching 1-D arrays")
        order = np.argsort(energy)
        self.energy = energy[order]
        self.response = resp[order]

    @classmethod
    def from_csv(cls, path):
        table = np.loadtxt(path, delimiter=",", comments="#")
        return cls(table[:, 0], table[:, 1])

    def __call__(self, energy_pev):
        e = np.asarray(energy_pev, dtype=np.float64)
        safe = np.clip(self.response, 1e-300, None)
        logr = np.interp(np.log10(np.clip(e, 1e-300, None)),
                         np.log10(self.energy), np.log10(safe))
        out = np.power(10.0, logr)
        # Outside the tabulated range the response is unknown, not extrapolated
        return np.where((e < self.energy[0]) | (e > self.energy[-1]), 0.0, out)


def geometric_aperture_m2sr(area_km2, solid_angle_sr):
    """
    Area times accepted solid angle, in m^2 sr.

    The purely geometric part of an aperture: no physics beyond the terrain.
    """
    return np.asarray(area_km2, dtype=np.float64) * 1.0e6 * np.asarray(solid_angle_sr,
                                                                       dtype=np.float64)


def aperture_vs_energy(area_km2, solid_angle_sr, min_dist_m, max_dist_m,
                       energies_pev, response=None):
    """
    Aperture as a function of energy for one site.

        A(E) = area * Omega * P_decay(E) * response(E)

    ``P_decay`` is the probability the tau decays inside the accepted baseline window,
    ``exp(-d_min/L) - exp(-d_max/L)`` with ``L = (E/m_tau) c*tau``. It is exact and
    carries the whole geometric energy dependence: short baselines favour low energies,
    long baselines high ones.

    Returns:
    - ndarray: aperture in m^2 sr, one entry per energy.
    """
    response = response or unit_response
    energies = np.asarray(energies_pev, dtype=np.float64)
    geom = geometric_aperture_m2sr(area_km2, solid_angle_sr)
    decay = np.array([arrival_scan.decay_probability(min_dist_m, max_dist_m, float(e))
                      for e in np.atleast_1d(energies)])
    return geom * decay * np.asarray(response(energies), dtype=np.float64)


def peak_energy_pev(min_dist_m, max_dist_m, energies_pev=None):
    """
    Energy at which the decay factor peaks for a given baseline window.

    A site's geometry therefore predicts which energies it is best suited to, which is
    a check that can be made without any normalisation.
    """
    if energies_pev is None:
        energies_pev = np.logspace(-2, 5, 2000)
    energies = np.asarray(energies_pev, dtype=np.float64)
    decay = np.array([arrival_scan.decay_probability(min_dist_m, max_dist_m, float(e))
                      for e in energies])
    return float(energies[int(np.argmax(decay))])


def summarize_sites(site_details, min_dist_m, max_dist_m, energies_pev, response=None):
    """
    Aperture-versus-energy for every site that carries scan observables, plus the total.

    Uses each site's median accepted solid angle, so a site is credited with the
    acceptance typical of its pixels rather than of its best one.

    Returns:
    - dict with the energy grid, per-site curves, and their sum.
    """
    energies = np.asarray(energies_pev, dtype=np.float64)
    per_site = {}
    total = np.zeros_like(energies)
    for site in site_details:
        obs = site.get("arrival_scan")
        if not obs:
            continue
        curve = aperture_vs_energy(site["area_km2"], obs.get("solid_angle_sr_p50", 0.0),
                                   min_dist_m, max_dist_m, energies, response)
        per_site[int(site["site_id"])] = [float(v) for v in curve]
        total = total + curve
    if not per_site:
        return {}
    return {
        "energies_pev": [float(e) for e in energies],
        "per_site_m2sr": per_site,
        "total_m2sr": [float(v) for v in total],
        "peak_energy_pev": peak_energy_pev(min_dist_m, max_dist_m),
        "note": ("Relative: the geometric aperture and the analytic tau decay factor. "
                 "Absolute normalisation requires a detection response table."),
    }
