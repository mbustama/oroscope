"""
Scoring: turning geometric observables into comparable site quality.

Every criterion returns a value in [0, 1] with a documented shape, so criteria can be
combined and so a site's weakness can be attributed to a named component rather than
disappearing into a single opaque number.

Why scores rather than cuts. Measurement on real terrain showed the binary geometric
test carries almost no discriminating power: in the Andes the whole of a +/-3 degree
arrival window sits below the local horizon, so nearly every direction strikes rock
and a hit/no-hit criterion selects most of the map. What separates sites is *how much*
rock, at *what* distance, over *how much* solid angle -- all continuous quantities.

What these scores are and are not. They rank sites against each other for one
experiment and energy band. They are not apertures. Any factor that is
energy-dependent but site-independent cancels in a ranking, which is precisely why
useful relative scoring is possible without the differential acceptance table that
absolute apertures would need (see the roadmap, section 4.10).
"""

from __future__ import annotations

import numpy as np

import physics

# Composition rules. 'product' is unforgiving -- one bad component sinks the site --
# while 'mean' lets a strong component compensate. 'min' reports the weakest link.
COMPOSITION_MODES = ("product", "mean", "min")

__all__ = ["band_score", "saturating_score", "ramp_score", "compose",
           "score_candidates", "summarize_scores", "COMPOSITION_MODES",
           "DEFAULT_SCORE_CONFIG"]


def band_score(x: float | np.ndarray, lo: float, hi: float,
               soft_lo: float | None = None,
               soft_hi: float | None = None) -> np.ndarray:
    """
    A plateau of 1 between ``lo`` and ``hi``, falling linearly to 0 outside it.

    The shape criteria of this kind want: acceptable over a range, degrading either
    side, rather than a cliff at an arbitrary threshold. ``soft_lo``/``soft_hi`` set
    the width of the falling flanks and default to a quarter of the band width, so a
    value must be well outside the band before it scores zero.

    Column depth is the motivating case: the tau must be produced, which needs rock,
    and must escape, which limits how much -- so its score is a band with an optimum,
    not a floor.

    Parameters
    ----------
    x : float or array_like
        Value or values to score.
    lo, hi : float
        Edges of the plateau. Swapped automatically if given the wrong way round.
    soft_lo, soft_hi : float, optional
        Widths of the falling flanks below ``lo`` and above ``hi``. Default to a
        quarter of the band width each. Zero gives a hard cliff.

    Returns
    -------
    ndarray
        Scores in [0, 1].

    Examples
    --------
    >>> import scoring
    >>> float(scoring.band_score(5.0, 0.0, 10.0))          # inside the plateau
    1.0
    >>> float(scoring.band_score(-5.0, 0.0, 10.0))         # well below it
    0.0
    """
    x = np.asarray(x, dtype=np.float64)
    if hi < lo:
        lo, hi = hi, lo
    span = max(hi - lo, 1e-12)
    soft_lo = span * 0.25 if soft_lo is None else float(soft_lo)
    soft_hi = span * 0.25 if soft_hi is None else float(soft_hi)

    score = np.ones_like(x)
    if soft_lo > 0:
        below = x < lo
        score = np.where(below, np.clip(1.0 - (lo - x) / soft_lo, 0.0, 1.0), score)
    else:
        score = np.where(x < lo, 0.0, score)
    if soft_hi > 0:
        above = x > hi
        score = np.where(above, np.clip(1.0 - (x - hi) / soft_hi, 0.0, 1.0), score)
    else:
        score = np.where(x > hi, 0.0, score)
    return np.clip(score, 0.0, 1.0)


def saturating_score(x: float | np.ndarray, half_value: float) -> np.ndarray:
    """
    ``x / (x + half_value)``: rises from 0, reaching 0.5 at ``half_value``.

    For quantities where more is better with diminishing returns and no natural
    maximum -- accepted solid angle being the case in hand.

    Parameters
    ----------
    x : float or array_like
        Value or values to score. Negative values are clipped to zero.
    half_value : float
        Value scoring 0.5. Choosing it far below the range actually observed saturates
        the term so that it stops discriminating, which is a real failure mode: the
        0.05 sr default suits GRAND and flattens completely against a canyon's
        0.2-1.5 sr.

    Returns
    -------
    ndarray
        Scores in [0, 1).

    Examples
    --------
    >>> import scoring
    >>> float(scoring.saturating_score(0.05, 0.05))
    0.5
    >>> round(float(scoring.saturating_score(0.7, 0.05)), 3)   # saturated
    0.933
    """
    x = np.clip(np.asarray(x, dtype=np.float64), 0.0, None)
    half_value = max(float(half_value), 1e-12)
    return x / (x + half_value)


def ramp_score(x: float | np.ndarray, zero_at: float,
               one_at: float) -> np.ndarray:
    """
    Linear ramp from 0 at ``zero_at`` to 1 at ``one_at``; either order.

    Parameters
    ----------
    x : float or array_like
        Value or values to score.
    zero_at : float
        Value scoring 0.
    one_at : float
        Value scoring 1. May be below ``zero_at``, giving a falling ramp.

    Returns
    -------
    ndarray
        Scores in [0, 1].

    Examples
    --------
    >>> import scoring
    >>> [float(scoring.ramp_score(v, 0.0, 10.0)) for v in (-1.0, 5.0, 99.0)]
    [0.0, 0.5, 1.0]
    """
    x = np.asarray(x, dtype=np.float64)
    if one_at == zero_at:
        return np.where(x >= one_at, 1.0, 0.0)
    return np.clip((x - zero_at) / (one_at - zero_at), 0.0, 1.0)


def compose(components: dict[str, np.ndarray], mode: str = "product",
            weights: dict[str, float] | None = None) -> np.ndarray:
    """
    Combines named component scores into one value in [0, 1].

    Parameters
    ----------
    components : dict
        Component name to array of scores in [0, 1].
    mode : {'product', 'mean', 'min'}, optional
        How to combine them. ``product`` is unforgiving -- one bad component sinks the
        site -- ``mean`` lets a strong component compensate, and ``min`` reports the
        weakest link.
    weights : dict, optional
        Per-component weights, used by ``product`` as exponents and by ``mean`` as
        linear weights. Ignored by ``min``. A weight of 0 excludes a component.

    Returns
    -------
    ndarray
        The composed score, in [0, 1].

    Raises
    ------
    ValueError
        If no components are given, the mode is unknown, a weight is negative, or the
        weights for ``mean`` do not sum to a positive value.

    Notes
    -----
    A product of several components each in [0, 1] concentrates near zero, so a
    threshold on the result sits on a cliff. Measured on a real search, a score cut of
    0.0, 0.35 and 0.5 gave 45928, 2056 and 0 detector positions. Prefer ranking sites
    over thresholding a product.

    Examples
    --------
    >>> import numpy as np, scoring
    >>> parts = {"a": np.array([0.5]), "b": np.array([0.5])}
    >>> float(scoring.compose(parts, "product")[0])
    0.25
    >>> float(scoring.compose(parts, "mean")[0])
    0.5
    """
    if not components:
        raise ValueError("compose() needs at least one component")
    if mode not in COMPOSITION_MODES:
        raise ValueError(f"unknown composition mode {mode!r}; expected one of {COMPOSITION_MODES}")

    names = list(components)
    arrays = [np.clip(np.asarray(components[n], dtype=np.float64), 0.0, 1.0) for n in names]
    w = {n: 1.0 for n in names}
    if weights:
        w.update({n: float(v) for n, v in weights.items() if n in w})

    if mode == "min":
        return np.min(np.stack(arrays), axis=0)
    if mode == "mean":
        total = sum(w[n] for n in names)
        if total <= 0:
            raise ValueError("weights for 'mean' composition must sum to a positive value")
        return sum(a * w[n] for a, n in zip(arrays, names)) / total
    # product: weights act as exponents, so a weight of 0 excludes a component.
    # A component of exactly 0 must sink the total -- physical impossibilities score 0 --
    # so nothing is clipped away from zero here.
    out = np.ones_like(arrays[0])
    for a, n in zip(arrays, names):
        if w[n] < 0:
            raise ValueError(f"negative weight for component {n!r}")
        if w[n] == 0:
            continue
        out = out * (a if w[n] == 1.0 else np.power(a, w[n]))
    return np.clip(out, 0.0, 1.0)


# Defaults are deliberately wide, since the physically motivated column-depth band for
# a given energy range is an open question. Wide defaults rank without pretending to
# encode physics the tool has not been given.
DEFAULT_SCORE_CONFIG = {
    "depth_band_gcm2": (1.0e5, 1.0e7),
    # Shower maturity. For radio the criterion is a threshold: emission comes from
    # around shower maximum and then propagates through transparent air, so being well
    # beyond maximum costs nothing here -- the distance trade is amplitude against
    # footprint area, which the footprint term carries. For a particle array the
    # content dies after maximum, so there it genuinely is a band.
    "grammage_mode": "radio",                  # 'radio' threshold, or 'particle' band
    "grammage_maturity_gcm2": physics.X_MAX_GCM2,
    "grammage_band_gcm2": (physics.X_MAX_GCM2, 4.0 * physics.X_MAX_GCM2),
    # Neutrino interaction length for the Earth-chord attenuation term. None leaves the
    # term out: it is strongly energy-dependent and guessing it would be worse than
    # reporting the chord and letting it be supplied.
    "nu_interaction_length_gcm2": None,
    # Antenna spacing, for comparing the array to the radio footprint
    "spacing_m": None,
    # Rock overburden required along the arrival direction to reject atmospheric muons.
    # A floor, not a band: more rock is always better for background rejection.
    "muon_shielding_km": None,
    "distance_band_m": None,           # defaults to the configured decay-baseline window
    # Tau decay. None leaves the term out, because the probability is strongly
    # energy-dependent and a single number cannot stand in for a spectrum. Supplying an
    # energy asks: given this site's baseline, how often does the tau actually decay in
    # the gap with room left for a shower?
    "decay_energy_pev": None,
    "shower_development_m": 3000.0,
    "solid_angle_half_sr": 0.05,
    "clearance_full_at": 1.0,          # clearance ratio scoring 1 (Fresnel radii)
    "composition": "product",
    "weights": None,
}


def score_candidates(observables: dict[str, np.ndarray],
                     config: dict | None = None,
                     distance_window_m: tuple[float, float] | None = None
                     ) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """
    Scores candidates from their scan observables.

    Components:

    - ``depth``     band on column depth: enough rock to interact, not so much that
                    the tau cannot escape;
    - ``distance``  band on the exit-point distance, defaulting to the configured
                    decay-baseline window;
    - ``solid_angle`` saturating in accepted solid angle, since more acceptance is
                    better with diminishing returns;
    - ``decay``     probability the tau decays in the gap with room left for a shower,
                    present only when ``decay_energy_pev`` is supplied;
    - ``clearance`` present only when a Fresnel frequency was configured.

    Parameters
    ----------
    observables : dict
        Per-candidate arrays from :func:`arrival_scan.scan`. ``cells``,
        ``solid_angle_sr``, ``mean_distance_m`` and ``max_depth_gcm2`` are required;
        the rest enable their components when present.
    config : dict, optional
        Overrides for :data:`DEFAULT_SCORE_CONFIG`. ``None`` values are ignored, so a
        sparse configuration leaves the remaining defaults alone.
    distance_window_m : tuple of float, optional
        Fallback distance band, used when the configuration does not set one. Normally
        the decay-baseline window the scan was run with.

    Returns
    -------
    total : ndarray
        The composed score per candidate. Zero wherever no direction was accepted,
        whatever the components say.
    components : dict
        The individual component scores, by name, so a site's weakness can be
        attributed rather than disappearing into one number.

    Examples
    --------
    >>> import numpy as np, scoring
    >>> obs = {"cells": np.array([4, 0]),
    ...        "solid_angle_sr": np.array([0.4, 0.4]),
    ...        "mean_distance_m": np.array([2.0e4, 2.0e4]),
    ...        "max_depth_gcm2": np.array([1.0e6, 1.0e6])}
    >>> total, parts = scoring.score_candidates(obs)
    >>> float(total[1])                      # no accepted direction scores zero
    0.0
    >>> sorted(parts)
    ['depth', 'solid_angle']
    """
    cfg = dict(DEFAULT_SCORE_CONFIG)
    if config:
        # None means "leave the default alone", so callers can pass a sparse config
        cfg.update({k: v for k, v in config.items() if v is not None})

    depth = np.asarray(observables["max_depth_gcm2"], dtype=np.float64)
    dist = np.asarray(observables["mean_distance_m"], dtype=np.float64)
    omega = np.asarray(observables["solid_angle_sr"], dtype=np.float64)

    dist_band = cfg.get("distance_band_m") or distance_window_m
    components = {
        "depth": band_score(depth, *cfg["depth_band_gcm2"]),
        "solid_angle": saturating_score(omega, cfg["solid_angle_half_sr"]),
    }
    if dist_band:
        components["distance"] = band_score(dist, *dist_band)

    # Muon shielding: a hard floor on column depth. A direction with less rock behind
    # it than this cannot claim neutrino purity, so it scores zero outright.
    shielding_km = cfg.get("muon_shielding_km")
    if shielding_km:
        required = physics.muon_shielding_gcm2(shielding_km)
        components["muon_shielding"] = (depth >= required).astype(np.float64)

    # Geomagnetic: the mean sin(alpha) over accepted directions. A site whose targets
    # all lie along the field radiates little geomagnetic signal however good its
    # terrain is, which no purely geometric measure can see.
    geomag = observables.get("geomag_solid_angle_sr")
    if geomag is not None:
        raw = np.asarray(geomag, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(omega > 0, raw / np.clip(omega, 1e-30, None), 0.0)
        if np.any(ratio < 0.999):        # all-ones means weighting was not applied
            components["geomagnetic"] = np.clip(ratio, 0.0, 1.0)

    # Shower maturity: grammage, not metres. Altitude enters here and nowhere else.
    grammage = observables.get("path_grammage_gcm2")
    if grammage is not None and np.any(np.asarray(grammage) > 0):
        if cfg["grammage_mode"] == "particle":
            components["shower"] = band_score(grammage, *cfg["grammage_band_gcm2"])
        else:
            components["shower"] = ramp_score(grammage, 0.0, cfg["grammage_maturity_gcm2"])

    # Tau decay in the gap. The tau leaves the far surface at `dist` and travels toward
    # the detector; it is only useful if it decays with enough path left for the shower
    # to develop. So the usable stretch is (dist - shower_development), and
    #
    #     P = 1 - exp(-usable / L),   L = the boosted decay length
    #
    # For GRAND this is largely implicit already, since the distance window is derived
    # from L. For a canyon it is not: TAMBO's window comes from the terrain, and at
    # 1 EeV the decay length is ~49 km against a ~3 km crossing, so only a few per cent
    # of taus decay in time. That suppression is invisible to every other term here.
    decay_energy = cfg.get("decay_energy_pev")
    if decay_energy:
        length_m = physics.tau_decay_length_m(decay_energy)
        if length_m > 0:
            usable = np.clip(dist - cfg.get("shower_development_m", 0.0), 0.0, None)
            components["decay"] = 1.0 - np.exp(-usable / length_m)

    # Earth-chord attenuation, only when an interaction length is supplied
    chord = observables.get("earth_chord_gcm2")
    x_int = cfg.get("nu_interaction_length_gcm2")
    if chord is not None and x_int:
        components["nu_survival"] = np.exp(-np.asarray(chord, dtype=np.float64) / x_int)

    # Footprint sampling: a higher site has a narrower Cherenkov cone and so a smaller
    # footprint, and needs a denser array for the same trigger efficiency
    altitude = observables.get("altitude_m")
    spacing = cfg.get("spacing_m")
    if altitude is not None and spacing:
        alt = np.asarray(altitude, dtype=np.float64)
        theta_c = np.sqrt(2.0 * physics.SEA_LEVEL_REFRACTIVITY
                          * np.exp(-alt / physics.DENSITY_SCALE_HEIGHT_M))
        sampling = 2.0 * dist * theta_c / spacing
        components["footprint"] = ramp_score(sampling, 0.0, 1.0)

    clearance = observables.get("best_clearance_ratio")
    if clearance is not None and np.any(np.asarray(clearance) > 0):
        components["clearance"] = ramp_score(clearance, 0.0, cfg["clearance_full_at"])

    # A candidate with no accepted direction scores zero regardless of composition
    accepted = np.asarray(observables["cells"]) > 0
    total = compose(components, cfg["composition"], cfg.get("weights"))
    total = np.where(accepted, total, 0.0)
    return total, components


def summarize_scores(scores: np.ndarray,
                     components: dict[str, np.ndarray]) -> dict:
    """
    Distribution summary of a score set, for reporting and per-site records.

    Storing a distribution rather than a single number is deliberate: it lets a site's
    quality be re-examined later without re-running the terrain analysis.

    Parameters
    ----------
    scores : ndarray
        Composed scores, one per candidate.
    components : dict
        The per-component scores behind them, from :func:`score_candidates`.

    Returns
    -------
    dict
        Count, mean, median and 90th percentile of the total, and the mean of each
        named component so a weak site can be attributed rather than merely ranked.
    """
    if len(scores) == 0:
        return {}
    out = {
        "score_mean": float(np.mean(scores)),
        "score_p50": float(np.median(scores)),
        "score_p90": float(np.percentile(scores, 90)),
        "score_max": float(np.max(scores)),
    }
    for name, values in components.items():
        out[f"{name}_score_p50"] = float(np.median(values))
    return out
