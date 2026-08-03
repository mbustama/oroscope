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

import numpy as np

# Composition rules. 'product' is unforgiving -- one bad component sinks the site --
# while 'mean' lets a strong component compensate. 'min' reports the weakest link.
COMPOSITION_MODES = ("product", "mean", "min")


def band_score(x, lo, hi, soft_lo=None, soft_hi=None):
    """
    A plateau of 1 between ``lo`` and ``hi``, falling linearly to 0 outside it.

    The shape criteria of this kind want: acceptable over a range, degrading either
    side, rather than a cliff at an arbitrary threshold. ``soft_lo``/``soft_hi`` set
    the width of the falling flanks and default to a quarter of the band width, so a
    value must be well outside the band before it scores zero.

    Column depth is the motivating case: the tau must be produced, which needs rock,
    and must escape, which limits how much -- so its score is a band with an optimum,
    not a floor.
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


def saturating_score(x, half_value):
    """
    ``x / (x + half_value)``: rises from 0, reaching 0.5 at ``half_value``.

    For quantities where more is better with diminishing returns and no natural
    maximum -- accepted solid angle being the case in hand.
    """
    x = np.clip(np.asarray(x, dtype=np.float64), 0.0, None)
    half_value = max(float(half_value), 1e-12)
    return x / (x + half_value)


def ramp_score(x, zero_at, one_at):
    """Linear ramp from 0 at ``zero_at`` to 1 at ``one_at``; either order."""
    x = np.asarray(x, dtype=np.float64)
    if one_at == zero_at:
        return np.where(x >= one_at, 1.0, 0.0)
    return np.clip((x - zero_at) / (one_at - zero_at), 0.0, 1.0)


def compose(components, mode="product", weights=None):
    """
    Combines named component scores into one value in [0, 1].

    Parameters:
    - components (dict): name -> array of scores in [0, 1].
    - mode (str): one of COMPOSITION_MODES.
    - weights (dict): optional per-component weights, used by 'product' as exponents
      and by 'mean' as linear weights. Ignored by 'min'.

    Returns:
    - ndarray: the composed score.
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
    "distance_band_m": None,           # defaults to the configured decay-baseline window
    "solid_angle_half_sr": 0.05,
    "clearance_full_at": 1.0,          # clearance ratio scoring 1 (Fresnel radii)
    "composition": "product",
    "weights": None,
}


def score_candidates(observables, config=None, distance_window_m=None):
    """
    Scores candidates from their scan observables.

    Components:
    - ``depth``     band on column depth: enough rock to interact, not so much that
                    the tau cannot escape;
    - ``distance``  band on the exit-point distance, defaulting to the configured
                    decay-baseline window;
    - ``solid_angle`` saturating in accepted solid angle, since more acceptance is
                    better with diminishing returns;
    - ``clearance`` present only when a Fresnel frequency was configured.

    Returns:
    - tuple(ndarray, dict): the composed score and the per-component scores.
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

    clearance = observables.get("best_clearance_ratio")
    if clearance is not None and np.any(np.asarray(clearance) > 0):
        components["clearance"] = ramp_score(clearance, 0.0, cfg["clearance_full_at"])

    # A candidate with no accepted direction scores zero regardless of composition
    accepted = np.asarray(observables["cells"]) > 0
    total = compose(components, cfg["composition"], cfg.get("weights"))
    total = np.where(accepted, total, 0.0)
    return total, components


def summarize_scores(scores, components):
    """Distribution summary of a score set, for reporting and per-site records."""
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
