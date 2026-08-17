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

import math

import numpy as np

from oroscope import arrival_scan

_SIN_60 = math.sin(math.radians(60.0))

__all__ = ["unit_response", "TabulatedResponse", "geometric_aperture_m2sr",
           "aperture_vs_energy", "peak_energy_pev", "infer_response",
           "load_curve_csv", "summarize_sites", "PUBLISHED_ARRAYS",
           "array_scale_factor", "scale_published_curve",
           "absolute_from_published"]

# The array each published curve in ``data/`` was simulated for. A published aperture
# is a property of *that* array at *that* site, and neither factor is separable from
# the other by any operation on the curve. What can be corrected is the array size --
# see :func:`array_scale_factor`, and the assumptions page for what cannot.
#
# ``units`` is the detector count the simulation used; ``spacing_km`` the lattice
# spacing it used them at. Both are needed: a count alone does not say how much ground
# was instrumented, and it is ground that the aperture scales with.
PUBLISHED_ARRAYS = {
    "tambo_aperture_fig3": {
        "units": 5000, "spacing_km": 0.15, "grid_type": "hex",
        "curve_units": "m^2 sr",
        "site": "Colca Canyon",
        "source": "TAMBO Collaboration, Nature Astronomy (2026), Fig. 3",
    },
    "grand_effective_area_fig25": {
        "units": 10000, "spacing_km": 1.0, "grid_type": "hex",
        "curve_units": "cm^2",
        "site": "HotSpot1 (a prototypical site, not a surveyed one)",
        "source": "GRAND Collaboration, arXiv:1810.09994, Fig. 25",
    },
}


def array_scale_factor(target_units, published, target_spacing_km=None,
                       target_grid_type=None):
    """
    How much larger this array is than the one a published curve was simulated for.

    A published aperture or effective area belongs to a specific array at a specific
    site. Oroscope changes both, and only one of them can be corrected for by
    arithmetic. This is that one: the array **size**.

    Aperture scales with **instrumented ground**, not with detector count as such, so
    the factor is

        (N_target * s_target^2) / (N_published * s_published^2)

    which collapses to the ratio of counts when the two spacings agree. Both terms are
    carried because they must be: doubling the count at fixed spacing doubles the
    ground and roughly doubles the aperture, while doubling it at fixed ground only
    makes the array denser, and a denser array past the point where it already samples
    the Cherenkov footprint adds very little. Scaling a densified array by its count
    alone would inflate the answer by exactly the factor by which it was densified.

    Parameters
    ----------
    target_units : int
        Detectors oroscope fits on this ground.
    published : dict or str
        An entry of :data:`PUBLISHED_ARRAYS`, or its key.
    target_spacing_km : float, optional
        Lattice spacing oroscope used. Defaults to the published spacing, which makes
        the factor a plain ratio of counts.
    target_grid_type : {'hex', 'square'}, optional
        Lattice oroscope used. Ground per detector is ``spacing^2`` times sin60 for a
        triangular lattice and 1 for a square one, so the factor cancels only when both
        lattices match. Defaults to the published array's.

    Returns
    -------
    float
        Multiplier to apply to the published curve.

    Raises
    ------
    ValueError
        If the published entry is unknown or either count is not positive.

    Notes
    -----
    **This corrects the array and not the site.** The published simulation carries its
    own terrain -- Colca's walls for TAMBO, a prototypical site for GRAND -- with its
    own distribution of column depth, arrival elevation and target distance, and no
    operation on an integral curve can remove them. A scaled curve therefore reads:
    *"what this many detectors would have achieved on the ground the simulation
    assumed"*, not *"what they will achieve on this ground"*. See
    :doc:`assumptions`.

    Examples
    --------
    >>> from oroscope import aperture
    >>> round(aperture.array_scale_factor(10000, "tambo_aperture_fig3"), 3)
    2.0

    Same detector count, spread twice as far apart, is four times the ground:

    >>> round(aperture.array_scale_factor(5000, "tambo_aperture_fig3",
    ...                                   target_spacing_km=0.30), 3)
    4.0
    """
    if isinstance(published, str):
        if published not in PUBLISHED_ARRAYS:
            raise ValueError(
                f"unknown published array {published!r}; "
                f"expected one of {sorted(PUBLISHED_ARRAYS)}")
        published = PUBLISHED_ARRAYS[published]
    # Ground per detector is spacing^2 times a lattice factor -- sin60 for a triangular
    # lattice, 1 for a square one. It cancels only when both lattices match, and they
    # need not: a square-gridded run against the hex-simulated TAMBO array covers
    # 1/sin60 = 1.1547x the ground per detector, so ignoring it under-scales by 15.5%.
    lattice = {"hex": _SIN_60, "square": 1.0}
    f_pub = lattice.get(str(published.get("grid_type", "hex")).lower(), _SIN_60)
    f_target = lattice.get(str(target_grid_type or
                               published.get("grid_type", "hex")).lower(), _SIN_60)
    n_pub = float(published["units"])
    s_pub = float(published["spacing_km"])
    n_target = float(target_units)
    # `None` means "the published spacing"; 0 does not, and must not be swallowed by a
    # truthiness test into meaning it. A zero spacing is a nonsensical array, and
    # quietly reading it as "same as published" would return a factor of 1 for it.
    s_target = float(s_pub if target_spacing_km is None else target_spacing_km)
    if n_pub <= 0 or n_target <= 0:
        raise ValueError("detector counts must be positive")
    if s_pub <= 0 or s_target <= 0:
        raise ValueError("spacings must be positive")
    return (n_target * s_target ** 2 * f_target) / (n_pub * s_pub ** 2 * f_pub)


def scale_published_curve(values, target_units, published,
                          target_spacing_km=None, target_grid_type=None):
    """
    A published curve rescaled to the array oroscope actually found room for.

    Parameters
    ----------
    values : array_like
        The published curve, in whatever units it came in -- m^2 sr for an aperture,
        cm^2 for an effective area. The scaling is dimensionless, so the units survive.
    target_units : int
        Detectors oroscope fits on this ground.
    published : dict or str
        As :func:`array_scale_factor`.
    target_spacing_km : float, optional
        As :func:`array_scale_factor`.
    target_grid_type : {'hex', 'square'}, optional
        As :func:`array_scale_factor`.

    Returns
    -------
    ndarray
        The curve, scaled.

    Examples
    --------
    Half the published detector count is half the aperture, at the same spacing:

    >>> from oroscope import aperture
    >>> published = [10.0, 100.0, 1000.0]
    >>> aperture.scale_published_curve(published, 2500, "tambo_aperture_fig3")
    array([  5.,  50., 500.])
    """
    factor = array_scale_factor(target_units, published, target_spacing_km,
                                target_grid_type)
    return np.asarray(values, dtype=np.float64) * factor


def absolute_from_published(results, curve_path, published,
                            target_spacing_km=None, target_grid_type=None):
    """
    A run's absolute aperture, by scaling a published curve to the array it found room for.

    This is post-processing, exactly as roadmap §4.10 step 4 intends: the search stores
    what terrain determines, and folding a published curve against it needs no re-run.

    Parameters
    ----------
    results : dict
        A results dictionary, or one loaded from a run's JSON.
    curve_path : str
        Two-column digitized curve, as in ``data/``.
    published : dict or str
        The array that curve was simulated for -- an entry of :data:`PUBLISHED_ARRAYS`
        or its key.
    target_spacing_km : float, optional
        Spacing this run used. Read from the run's own ``spacing_km`` when omitted,
        which is what makes the answer honest: reading it from anywhere else invites
        the density error :func:`array_scale_factor` exists to prevent.
    target_grid_type : {'hex', 'square'}, optional
        Lattice this run used. Read from the run's own ``grid_type`` when omitted.

    Returns
    -------
    dict
        ``{"energies_pev", "aperture", "units", "scale_factor", "detectors",
        "published", "caveat"}``. ``aperture`` carries the published curve's own units.

    Notes
    -----
    **This corrects the array, not the site**, and the returned ``caveat`` says so in
    the artefact itself rather than only here. The published simulation carries its own
    terrain -- its column depths, its target distances, its trigger geometry -- and no
    operation on an integral curve separates those from it. Read the result as *what
    this many detectors would have achieved on the ground that simulation assumed*.

    Examples
    --------
    Paths are resolved by the caller, so this example builds one from the package's own
    location rather than from the working directory -- the test job runs from ``tests/``.

    **This resolves only in a source checkout.** ``data/`` sits beside ``src/`` in the
    repository and is not shipped in the wheel, so from an installed
    ``site-packages/oroscope/aperture.py`` the three ``dirname`` calls land outside
    site-packages entirely. That is not a limitation of the function -- it takes whatever
    path it is given -- but the example below cannot be pasted into a fresh install
    unchanged. Point it at your own copy of the curve, or at a clone.

    >>> import os
    >>> from oroscope import aperture
    >>> repo = os.path.dirname(os.path.dirname(os.path.dirname(aperture.__file__)))
    >>> curve = os.path.join(repo, "data", "tambo_aperture_fig3.csv")
    >>> results = {"results": {"total_capacity": 10000},
    ...            "parameters": {"antenna_spacing_km": 0.15}}
    >>> out = aperture.absolute_from_published(results, curve, "tambo_aperture_fig3")
    >>> out["units"], round(out["scale_factor"], 3), out["detectors"]
    ('m^2 sr', 2.0, 10000)
    """
    if isinstance(published, str):
        published = PUBLISHED_ARRAYS[published]
    params = results.get("parameters", {}) or {}
    # `spacing_km` is what a results file records; `antenna_spacing_km` is the config
    # spelling and never appears there, so reading it left target_spacing_km None and
    # silently fell back to the *published* spacing -- a 44x under-report on a real
    # 1 km GRAND run, which is the density error this argument exists to prevent.
    raw = (results.get("results", {}) or {}).get("total_capacity", 0)
    if not isinstance(raw, (int, float)):
        raise ValueError(
            f"this run reports total_capacity={raw!r} rather than a number, so there "
            f"is no detector count to scale by. Capacity is only counted in "
            f"search_mode 'distributed'.")
    detectors = int(raw)
    if target_spacing_km is None:
        target_spacing_km = params.get("spacing_km",
                                       params.get("antenna_spacing_km"))
    if target_grid_type is None:
        target_grid_type = params.get("grid_type")
    energies, values = load_curve_csv(curve_path)
    factor = array_scale_factor(detectors, published, target_spacing_km,
                                target_grid_type)
    # From the registry, which knows what each curve is, rather than sniffed out of the
    # path: a file moved into a directory named "aperture" relabelled cm^2 as m^2 sr.
    units = published.get("curve_units", "unknown")
    return {
        "energies_pev": [float(e) for e in energies],
        "aperture": [float(v) for v in np.asarray(values) * factor],
        "units": units,
        "scale_factor": float(factor),
        "detectors": detectors,
        "published": dict(published),
        "caveat": ("Scaled for array size only. The published simulation's own site -- "
                   f"{published.get('site', 'unknown')} -- is folded into this curve and "
                   "cannot be divided out: its column depths, target distances and "
                   "trigger geometry all remain. Read as what this many detectors would "
                   "have achieved on the ground that simulation assumed, not on this "
                   "ground."),
    }


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
