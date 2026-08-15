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

from __future__ import annotations

import numpy as np

from oroscope import arrival_scan

__all__ = ["unit_response", "TabulatedResponse", "geometric_aperture_m2sr",
           "aperture_vs_energy", "peak_energy_pev", "infer_response",
           "load_curve_csv", "summarize_sites"]


def unit_response(energy_pev: float | np.ndarray) -> np.ndarray:
    """
    Default detection response: energy-independent and equal to one.

    A placeholder that makes the normalisation explicit rather than hidden. Replace it
    with a :class:`TabulatedResponse` once a real acceptance table is available.

    Parameters
    ----------
    energy_pev : float or array_like
        Energies at which to evaluate the response, in PeV.

    Returns
    -------
    ndarray
        Ones, matching the shape of the input.

    Examples
    --------
    >>> from oroscope import aperture
    >>> aperture.unit_response([1.0, 10.0]).tolist()
    [1.0, 1.0]
    """
    return np.ones_like(np.asarray(energy_pev, dtype=np.float64))


class TabulatedResponse:
    """
    Detection response interpolated from a supplied table.

    Log-log interpolation, since both response and energy span decades. Supply a
    two-column CSV of energy in PeV against relative response, or the arrays directly.

    Outside the tabulated range the response is returned as zero rather than
    extrapolated: beyond the table the response is unknown, and a silent extrapolation
    over decades of energy would be an invention.

    Parameters
    ----------
    energy_pev : array_like
        Tabulated energies, in PeV. Sorted internally, so any order will do.
    response : array_like
        Relative response at each energy, same shape.

    Raises
    ------
    ValueError
        If the two are not matching one-dimensional arrays.

    Examples
    --------
    >>> from oroscope import aperture
    >>> r = aperture.TabulatedResponse([1.0, 10.0, 100.0], [0.1, 1.0, 0.5])
    >>> round(float(r(10.0)), 3)
    1.0
    >>> float(r(1000.0))                 # outside the table: unknown, not extrapolated
    0.0
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
    def from_csv(cls, path: str) -> "TabulatedResponse":
        """
        Builds a response from a two-column CSV of energy in PeV against response.

        Parameters
        ----------
        path : str
            Path to the CSV. ``#`` introduces a comment.

        Returns
        -------
        TabulatedResponse
        """
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


def geometric_aperture_m2sr(area_km2: float | np.ndarray,
                            solid_angle_sr: float | np.ndarray) -> np.ndarray:
    """
    Area times accepted solid angle, in m^2 sr.

    The purely geometric part of an aperture: no physics beyond the terrain.

    Parameters
    ----------
    area_km2 : float or array_like
        Usable area, in km^2.
    solid_angle_sr : float or array_like
        Accepted solid angle, in steradians.

    Returns
    -------
    ndarray
        Geometric aperture, in m^2 sr.

    Examples
    --------
    >>> from oroscope import aperture
    >>> f"{float(aperture.geometric_aperture_m2sr(100.0, 0.5)):.2e}"
    '5.00e+07'
    """
    return np.asarray(area_km2, dtype=np.float64) * 1.0e6 * np.asarray(solid_angle_sr,
                                                                       dtype=np.float64)


def aperture_vs_energy(area_km2: float, solid_angle_sr: float, min_dist_m: float,
                       max_dist_m: float,
                       energies_pev, response=None):
    """
    Aperture as a function of energy for one site.

        A(E) = area * Omega * P_decay(E) * response(E)

    ``P_decay`` is the probability the tau decays inside the accepted baseline window,
    ``exp(-d_min/L) - exp(-d_max/L)`` with ``L = (E/m_tau) c*tau``. It is exact and
    carries the whole geometric energy dependence: short baselines favour low energies,
    long baselines high ones.

    Parameters
    ----------
    area_km2 : float
        Usable area of the site, in km^2.
    solid_angle_sr : float
        Accepted solid angle, in steradians.
    min_dist_m, max_dist_m : float
        Ends of the accepted decay-baseline window, in metres.
    energies_pev : array_like
        Energies at which to evaluate, in PeV.
    response : callable, optional
        Detection response as a function of energy. Defaults to
        :func:`unit_response`, which leaves the normalisation explicit.

    Returns
    -------
    ndarray
        Aperture in m^2 sr, one entry per energy.

    Examples
    --------
    >>> import numpy as np
    >>> from oroscope import aperture
    >>> a = aperture.aperture_vs_energy(100.0, 0.5, 1.0e4, 8.0e4, [1.0, 100.0])
    >>> bool(a[1] > a[0])          # a long baseline favours higher energies
    True
    """
    response = response or unit_response
    energies = np.asarray(energies_pev, dtype=np.float64)
    geom = geometric_aperture_m2sr(area_km2, solid_angle_sr)
    decay = np.array([arrival_scan.decay_probability(min_dist_m, max_dist_m, float(e))
                      for e in np.atleast_1d(energies)])
    return geom * decay * np.asarray(response(energies), dtype=np.float64)


def peak_energy_pev(min_dist_m: float, max_dist_m: float,
                    energies_pev: np.ndarray | None = None) -> float:
    """
    Energy at which the decay factor peaks for a given baseline window.

    A site's geometry therefore predicts which energies it is best suited to, which is
    a check that can be made without any normalisation.

    Parameters
    ----------
    min_dist_m, max_dist_m : float
        Ends of the accepted decay-baseline window, in metres.
    energies_pev : ndarray, optional
        Energies to search, in PeV. Defaults to a wide log-spaced grid.

    Returns
    -------
    float
        Energy at which the decay factor peaks, in PeV.

    Examples
    --------
    >>> from oroscope import aperture
    >>> near = aperture.peak_energy_pev(1.0e3, 5.0e3)
    >>> far = aperture.peak_energy_pev(1.0e4, 8.0e4)
    >>> bool(far > near)           # longer baselines peak at higher energy
    True
    """
    if energies_pev is None:
        energies_pev = np.logspace(-2, 5, 2000)
    energies = np.asarray(energies_pev, dtype=np.float64)
    decay = np.array([arrival_scan.decay_probability(min_dist_m, max_dist_m, float(e))
                      for e in energies])
    return float(energies[int(np.argmax(decay))])


def infer_response(published_energy_pev: np.ndarray, published_value: np.ndarray,
                   area_km2: float, solid_angle_sr: float,
                   min_dist_m: float, max_dist_m: float,
                   min_model_fraction: float = 1e-3
                   ) -> tuple[np.ndarray, np.ndarray]:
    """
    Response function implied by a published curve, given our geometric model.

    This is the useful thing to do with an integral aperture when no differential
    acceptance table exists. Our model supplies the two factors terrain and kinematics
    determine -- geometric aperture and the analytic tau decay probability -- so
    dividing a published curve by them leaves everything else:

        response(E) = published(E) / (area * Omega * P_decay(E))

    What remains is the neutrino interaction and tau exit probability, the trigger
    efficiency, and any normalisation the published configuration carries. Its *shape*
    in energy is then usable as a weight for a site of the same experiment, which is
    exactly the piece section 4.10 had to leave out.

    The caveat is the same one that applies to the published curves themselves: they
    are integral over one array and one site, so the inferred response inherits that
    site's geometry. It is a better weight than a flat response, not a substitute for
    a differential table.

    Where the decay probability is negligible the division is ill-conditioned and the
    ratio explodes -- at 0.35 PeV over a canyon baseline it is of order 1e-8 -- so
    energies whose model value falls below ``min_model_fraction`` of its own peak are
    excluded rather than allowed to dominate the normalisation.

    Parameters
    ----------
    published_energy_pev : array_like
        Energies of the published curve, in PeV.
    published_value : array_like
        Published aperture at each energy, in m^2 sr.
    area_km2 : float
        Usable area of the configuration the curve describes, in km^2.
    solid_angle_sr : float
        Accepted solid angle of that configuration, in steradians.
    min_dist_m, max_dist_m : float
        Ends of its accepted decay-baseline window, in metres.
    min_model_fraction : float, optional
        Energies whose model aperture falls below this fraction of its own peak are
        excluded, since the division there is ill-conditioned.

    Returns
    -------
    energies_pev : ndarray
        The well-conditioned subset of the input energies.
    response : ndarray
        Inferred response over that range, normalised to 1 at its maximum.
    """
    energies = np.asarray(published_energy_pev, dtype=np.float64)
    published = np.asarray(published_value, dtype=np.float64)
    model = aperture_vs_energy(area_km2, solid_angle_sr, min_dist_m, max_dist_m, energies)

    usable = model > min_model_fraction * model.max()
    energies, published, model = energies[usable], published[usable], model[usable]
    if energies.size == 0:
        return energies, model

    with np.errstate(divide="ignore", invalid="ignore"):
        response = np.where(model > 0, published / model, 0.0)
    peak = response.max()
    return energies, (response / peak if peak > 0 else response)


def load_curve_csv(path: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Loads a two-column digitized curve, returning (energy_pev, value).

    The files under ``data/`` store energy in GeV, as the published axes do, and this
    converts to PeV on the way in.

    Parameters
    ----------
    path : str
        Path to a two-column CSV of energy in GeV against value.

    Returns
    -------
    energy_pev : ndarray
        Energies, converted to PeV.
    value : ndarray
        The second column, unchanged.
    """
    table = np.loadtxt(path, delimiter=",", comments="#")
    return table[:, 0] / 1.0e6, table[:, 1]


def summarize_sites(site_details: list[dict], min_dist_m: float, max_dist_m: float,
                    energies_pev: np.ndarray,
                    response: np.ndarray | None = None) -> dict:
    """
    Aperture-versus-energy for every site that carries scan observables, plus the total.

    Uses each site's median accepted solid angle, so a site is credited with the
    acceptance typical of its pixels rather than of its best one.

    Parameters
    ----------
    site_details : list of dict
        Site records from the pipeline. Sites without an ``arrival_scan`` entry are
        skipped.
    min_dist_m, max_dist_m : float
        Ends of the accepted decay-baseline window, in metres.
    energies_pev : array_like
        Energy grid, in PeV.
    response : callable, optional
        Detection response. Defaults to :func:`unit_response`.

    Returns
    -------
    dict
        ``energies_pev``, a ``sites`` mapping of site id to aperture curve, and
        ``total``, their sum.
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
