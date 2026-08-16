"""
Turning a results file into an account of what was found and why.

Everything this module says is already in the results JSON. What is missing there is
the *story*: which constraint did the work, which sites survived and what weakened
them, and which of the numbers on the page are choices rather than measurements. A
reader who assembles that themselves gets it wrong in predictable ways -- most often
by reading the reported area as physics-accepted area, which it is not (see
:doc:`assumptions`, and the 2.29x measured at Colca).

The entry point is :func:`explain_results`, which takes the results dictionary and
returns a string. It runs nothing, opens nothing and needs no DEM, so the pipeline,
the library and a test can all use the same words, and the summary can be regenerated
from an old run's JSON long after the run.

Deliberately plain text: these summaries are meant to be handed to other people, and
ANSI colour does not survive a paste into an email.
"""

from __future__ import annotations

__all__ = ["explain_results", "explain_combination", "binding_constraint",
           "weakest_component", "site_strengths", "constraint_overlap",
           "closing_inflation", "selected_sites", "COMPONENT_MEANING",
           "STAGE_KNOBS", "AREA_INFLATION_AT_COLCA"]

# What each named score component means, and which measured quantity earns it.
#
# The components were named so that a weak site could be attributed. The same naming
# answers the more useful question -- why is this site *good*? -- but only with this
# table, because a component's name says what it is and not what a high score implies.
# Each entry is (label, observable field, what a high score means, unit formatter).
COMPONENT_MEANING = {
    "solid_angle": (
        "accepted sky", "solid_angle_sr",
        "a wide spread of arrival directions reaches usable terrain. This is the "
        "quantity an aperture is proportional to, so it is the closest thing to a "
        "single measure of how good a site is",
        lambda v: f"{v:.2f} sr"),
    "depth": (
        "column depth", "max_depth_gcm2",
        "enough rock behind the exit point to produce a tau, and not so much that it "
        "cannot escape -- a band, not a floor",
        lambda v: f"{v:,.0f} g/cm²"),
    "distance": (
        "exit distance", "mean_distance_m",
        "the terrain a ray strikes sits inside the decay-baseline window, so a tau "
        "has room to decay and its shower room to develop",
        lambda v: f"{v:,.0f} m"),
    "shower": (
        "shower development", "path_grammage_gcm2",
        "the air between the exit point and the detector is deep enough for the "
        "shower to have developed",
        lambda v: f"{v:,.0f} g/cm²"),
    "decay": (
        "tau decay", None,
        "the tau is likely to decay within the gap, folded over the assumed spectrum "
        "rather than evaluated at one energy",
        None),
    "geomagnetic": (
        "geomagnetic angle", None,
        "the accepted directions lie across the geomagnetic field rather than along "
        "it, so the shower radiates. Measured, an east-facing target is worth 3.7x a "
        "north-facing one",
        None),
    "footprint": (
        "footprint sampling", "altitude_m",
        "the array spacing samples the Cherenkov footprint well at this altitude",
        lambda v: f"{v:,.0f} m"),
    "clearance": (
        "Fresnel clearance", "best_clearance_ratio",
        "the radio path clears intervening terrain by enough Fresnel radii",
        lambda v: f"{v:.2f} F1"),
    "muon_shielding": (
        "muon shielding", "max_depth_gcm2",
        "enough rock overburden along the arrival direction to reject atmospheric "
        "muons",
        lambda v: f"{v:,.0f} g/cm²"),
    "nu_survival": (
        "Earth transmission", "earth_chord_gcm2",
        "the neutrino's chord through the Earth is short enough that it is not "
        "absorbed before reaching the exit point",
        lambda v: f"{v:,.3g} g/cm²"),
}

# Bands two experiments must genuinely share before they can stand on the same pixel.
#
# Only properties of the *ground itself* belong here. A pixel has one slope, one
# altitude and one aspect, so both experiments must accept those same values. The
# distance window and the arrival-elevation window are asked of the *view* from that
# pixel, and two experiments looking out from the same hillside at different ranges and
# different elevations are in no conflict whatever -- GRAND scanning 10-40 km within
# +/-3 degrees and TAMBO scanning 2-5 km within +/-20 degrees can both be satisfied
# from one patch of ground. Treating those as shared constraints produced the confident
# and wrong conclusion that two experiments sharing 50 km2 "cannot share ground at all".
_SHARED_BANDS = (
    ("deployable slope", "min_slope_deg", "max_slope_deg", "°"),
    ("altitude", "min_altitude", "max_altitude", " m"),
    ("aspect", "min_aspect_deg", "max_aspect_deg", "°"),
)

# Asked of the view, not of the ground. Reported for contrast, never as a conflict.
_VIEWING_BANDS = (
    ("target distance", "min_dist_km", "max_dist_km", " km"),
    ("arrival elevation", "elev_min_deg", "elev_max_deg", "°"),
)

# Measured with a stride-1 control run at Colca: closing a mask with a 1 km element
# more than doubles the area it reports. Quoted rather than recomputed because it is a
# property of the terrain and the element, not of any one run -- but it is the right
# order of magnitude to warn with whenever closing is enabled at all.
AREA_INFLATION_AT_COLCA = 2.29

# Funnel stages that must not be read as constraints.
#
# Striding is a deliberate subsample, not a filter: it removes four candidates in five
# by construction and the acceptance is unchanged (measured, and unbiased to 0.05%).
# Closing *adds* pixels. And the last stage counts a different thing -- pixels inside
# selected sites, estimated back up from the downsampled map.
_SUBSAMPLE_PREFIX = "kept by stride"
_NOT_A_CONSTRAINT = ("after gap closing", "pixels in selected sites")

# What to reach for when a stage turns out to be the binding one. The funnel's labels
# carry their own thresholds ("slope 20.0-60.0 deg"), so the match is by prefix.
STAGE_KNOBS = (
    ("finite elevation", "the DEM's own nodata voids -- not a parameter"),
    ("slope", "min_slope_deg / max_slope_deg"),
    ("altitude", "min_altitude / max_altitude"),
    ("aspect", "min_aspect_deg / max_aspect_deg"),
    ("road", "max_road_dist_km"),
    ("outside RFI", "rfi_zones"),
    (_SUBSAMPLE_PREFIX, "candidate_stride"),
    ("directions accepted", "the arrival window (elev_min_deg/elev_max_deg), the "
                            "distance window (min_dist_km/max_dist_km), "
                            "min_column_depth_gcm2 and min_target_slope_deg"),
    ("score >=", "min_score -- or switch to --score_percentile, which is rank-based"),
    ("score in top", "score_percentile"),
    ("after gap closing", "gap_close_km"),
    ("after pruning", "min_width_km"),
    ("pixels in selected sites", "min_sub_array_size and target_antennas"),
)


def _knob_for(stage):
    """Names the parameter behind a funnel stage, or None when nothing matches."""
    low = stage.lower()
    for prefix, knob in STAGE_KNOBS:
        if low.startswith(prefix.lower()):
            return knob
    return None


def binding_constraint(funnel):
    """
    Finds the funnel stage that removed the largest share of what reached it.

    This is the single most useful sentence a summary can offer, and it matters most
    when a search returns little or nothing: the stage where the survivor count
    collapses *is* the constraint responsible, and every other explanation is a guess.

    Two stages are excluded by construction. ``kept by stride`` is a deliberate
    subsample whose acceptance is unbiased, so calling it a constraint would name the
    same answer on nearly every run; ``after gap closing`` adds pixels rather than
    removing them.

    Parameters
    ----------
    funnel : dict
        Ordered stage-name to survivor-count mapping, as written to
        ``results["funnel"]``.

    Returns
    -------
    dict or None
        ``{"stage", "survivors", "before", "kept_fraction", "knob", "fatal"}``, or
        ``None`` when the funnel has fewer than two stages. ``fatal`` is True when the
        stage left nothing at all.

    Examples
    --------
    >>> from oroscope import explain
    >>> f = {"DEM pixels": 1000, "slope 3-25 deg": 800, "directions accepted": 40}
    >>> b = explain.binding_constraint(f)
    >>> b["stage"], round(b["kept_fraction"], 3), b["fatal"]
    ('directions accepted', 0.05, False)
    >>> explain.binding_constraint({"DEM pixels": 1000}) is None
    True
    """
    stages = [(name, int(count)) for name, count in funnel.items()]
    if len(stages) < 2:
        return None

    worst = None
    for (_, before), (name, after) in zip(stages, stages[1:]):
        low = name.lower()
        if low.startswith(_SUBSAMPLE_PREFIX) or low.startswith(_NOT_A_CONSTRAINT):
            continue
        if before <= 0:
            continue
        kept = after / before
        # A stage that leaves nothing wins outright, however gentle a stage upstream
        # of it looked: everything after it is zero for a reason that is not its own.
        fatal = after == 0
        candidate = {"stage": name, "survivors": after, "before": before,
                     "kept_fraction": kept, "knob": _knob_for(name), "fatal": fatal}
        if worst is None or fatal and not worst["fatal"] or (
                fatal == worst["fatal"] and kept < worst["kept_fraction"]):
            worst = candidate
        if fatal:
            break
    return worst


def selected_sites(results):
    """
    The sites actually in the result, separated from the ones that merely qualified.

    ``results["results"]["sites"]`` lists every site that cleared the area and capacity
    thresholds, which with ``stop_at_target`` is more than were selected: selection
    walks the capacity-sorted list until the target is met and stops. ``total_sites``,
    ``total_capacity`` and the exported mask all cover the selection only, so summing
    the list over-reports area and site count against every other number in the file.

    Parameters
    ----------
    results : dict
        A results dictionary.

    Returns
    -------
    selected : list of dict
        Sites in the result.
    rejected : list of dict
        Sites that qualified but were not selected. Usually empty.

    Notes
    -----
    Prefers each record's ``selected`` flag. Files written before that flag existed
    fall back to the first ``total_sites`` entries, which is exact: the list is sorted
    by capacity and selection takes it in order.

    Examples
    --------
    >>> from oroscope import explain
    >>> r = {"results": {"total_sites": 1, "sites": [
    ...     {"site_id": 2, "capacity_exact": 252, "selected": True},
    ...     {"site_id": 1, "capacity_exact": 36, "selected": False}]}}
    >>> chosen, rest = explain.selected_sites(r)
    >>> [s["site_id"] for s in chosen], [s["site_id"] for s in rest]
    ([2], [1])
    """
    block = results.get("results", {}) or {}
    sites = block.get("sites", []) or []
    if any("selected" in s for s in sites):
        chosen = [s for s in sites if s.get("selected")]
        return chosen, [s for s in sites if not s.get("selected")]

    total = block.get("total_sites")
    if isinstance(total, int) and 0 <= total < len(sites):
        return sites[:total], sites[total:]
    return list(sites), []


def closing_inflation(funnel, candidate_stride=1):
    """
    How much morphological closing grew this run's mask, measured from its own funnel.

    The 2.29x quoted from Colca is a property of that terrain and a 1 km element, not
    a constant. This run has the number in it: the stage before closing counts the
    accepted candidates, closing counts the pixels after, and the only correction
    needed between them is the stride -- which samples one candidate in
    ``candidate_stride`` and was measured unbiased.

    A ratio below 1 is not a bug and is worth reading carefully: it means the closing
    element was too small to bridge the gaps striding left, so the mask under-reports
    the accepted set rather than inflating it. That happens when the element is a few
    pixels across, as it is at a 100 m detector spacing.

    Parameters
    ----------
    funnel : dict
        Ordered stage-name to survivor-count mapping.
    candidate_stride : int, optional
        The run's ``candidate_stride``, to scale the accepted count back to full
        resolution.

    Returns
    -------
    float or None
        Closed pixels divided by estimated accepted pixels, or ``None`` when the
        funnel does not record both.

    Examples
    --------
    >>> from oroscope import explain
    >>> f = {"directions accepted": 100, "after gap closing": 450}
    >>> round(explain.closing_inflation(f, candidate_stride=5), 2)
    0.9
    """
    stages = list(funnel.items())
    for i, (name, count) in enumerate(stages):
        if not name.lower().startswith("after gap closing") or i == 0:
            continue
        before = stages[i - 1][1] * max(1, int(candidate_stride or 1))
        return count / before if before else None
    return None


def site_strengths(arrival_scan, statistic="p50", threshold=0.75):
    """
    Why a site is good: the criteria it satisfies well, and the measurement behind each.

    The mirror of :func:`weakest_component`, and the more useful half when a site has
    been *selected*. "Site 3555 scored 0.55" says nothing a reader can act on; "it sees
    1.08 sr of usable sky across a 3.1 km gap with 780,000 g/cm² of rock behind it, and
    every criterion but the accepted solid angle is satisfied outright" says what the
    ground is actually like.

    Parameters
    ----------
    arrival_scan : dict
        A site's ``arrival_scan`` record.
    statistic : str, optional
        Which per-site statistic to read, ``"mean"``, ``"p50"`` or ``"p90"``.
    threshold : float, optional
        Score at or above which a component counts as satisfied. 0.75 rather than 1.0
        because a band score falls off smoothly either side of its plateau, so
        insisting on exactly 1 would report nothing on most real sites.

    Returns
    -------
    list of dict
        One entry per satisfied component, strongest first, each with ``name``,
        ``label``, ``score``, ``means`` and -- where the record carries the observable
        behind it -- ``evidence``. Empty when the record has no components.

    Examples
    --------
    >>> from oroscope import explain
    >>> rec = {"score_solid_angle_p50": 0.9, "score_depth_p50": 1.0,
    ...        "solid_angle_sr_p50": 1.08, "max_depth_gcm2_p50": 784440.0}
    >>> [s["name"] for s in explain.site_strengths(rec)]
    ['depth', 'solid_angle']
    >>> explain.site_strengths(rec)[1]["evidence"]
    '1.08 sr'
    """
    suffix = "_" + statistic
    found = []
    for key, value in arrival_scan.items():
        if not key.startswith("score_") or not key.endswith(suffix):
            continue
        name = key[len("score_"):-len(suffix)]
        if not name or float(value) < threshold:
            continue
        label, field, means, fmt = COMPONENT_MEANING.get(
            name, (name.replace("_", " "), None, "", None))
        entry = {"name": name, "label": label, "score": float(value), "means": means}
        if field and fmt:
            measured = arrival_scan.get(f"{field}{suffix}")
            if measured is not None:
                entry["evidence"] = fmt(float(measured))
        found.append(entry)
    return sorted(found, key=lambda e: e["score"], reverse=True)


def constraint_overlap(params_a, params_b, bands=None):
    """
    Where two experiments' screening bands agree, and by how little.

    This is what decides whether two experiments can share ground, and it is decided
    before any arrival geometry is considered: a pixel has one slope, and both
    experiments must accept it. Measured at Colca, that is the whole story --
    GRAND's 3-25 degree deployable band against a canyon's ~40 degree walls leaves a
    20-25 degree sliver, and the joint area follows from that rather than from
    anything about neutrinos.

    Parameters
    ----------
    params_a, params_b : dict
        The two runs' recorded ``parameters`` blocks.
    bands : sequence, optional
        Which bands to compare, as ``(label, low key, high key, unit)``. Defaults to
        the properties of the ground itself -- slope, altitude, aspect -- which are the
        only ones both experiments must agree on. Pass :data:`_VIEWING_BANDS` to
        compare what each asks of the view instead, which need not agree at all.

    Returns
    -------
    list of dict
        One entry per band both runs recorded: ``label``, ``a``, ``b``, ``overlap``
        (a ``(low, high)`` pair or ``None``), ``width``, and ``share_of_narrower`` --
        the fraction of the tighter of the two bands that the overlap covers, which is
        the number that says how much room there is to share.

    Examples
    --------
    >>> from oroscope import explain
    >>> a = {"min_slope_deg": 3.0, "max_slope_deg": 25.0}
    >>> b = {"min_slope_deg": 20.0, "max_slope_deg": 60.0}
    >>> band = explain.constraint_overlap(a, b)[0]
    >>> band["overlap"], round(band["share_of_narrower"], 3)
    ((20.0, 25.0), 0.227)
    """
    out = []
    for label, lo_key, hi_key, unit in (bands or _SHARED_BANDS):
        a_lo, a_hi = _get(params_a, lo_key), _get(params_a, hi_key)
        b_lo, b_hi = _get(params_b, lo_key), _get(params_b, hi_key)
        if None in (a_lo, a_hi, b_lo, b_hi):
            continue                       # unset on one side is not a constraint
        a_parts, width_a = _band_intervals(a_lo, a_hi, label)
        b_parts, width_b = _band_intervals(b_lo, b_hi, label)
        pieces = [(max(x0, y0), min(x1, y1))
                  for x0, x1 in a_parts for y0, y1 in b_parts
                  if min(x1, y1) > max(x0, y0)]
        total = sum(hi - lo for lo, hi in pieces)
        narrower = min(width_a, width_b)
        out.append({
            "label": label, "unit": unit,
            "a": (a_lo, a_hi), "b": (b_lo, b_hi),
            # `overlap` is the widest piece and `width` the total across all of them.
            # For a band that wraps these differ, and the renderer prints them in one
            # sentence -- "shared 350-360 deg (100%)" for a 10 deg piece of a 20 deg
            # overlap. `pieces` carries the whole answer so a consumer can say so.
            "overlap": max(pieces, key=lambda p: p[1] - p[0]) if pieces else None,
            "pieces": pieces,
            "width": total,
            "share_of_narrower": (total / narrower) if pieces and narrower else 0.0,
        })
    return out


# Bands measured on a compass wrap; the others do not. Aspect is the only one here,
# and the screen already handles it -- `min_aspect_deg > max_aspect_deg` means an arc
# through north, and `get_candidates_chunked` reads it that way. This did not, so a
# north-facing window of 350-10 degrees was compared as though it ran *backwards* from
# 350 to 10, and reported no overlap with 0-90 when they plainly share 0-10.
_CIRCULAR_BANDS = ("aspect",)
_FULL_CIRCLE_DEG = 360.0


def _band_intervals(lo, hi, label):
    """
    A band as one or two ordinary intervals, plus its true width.

    A wrapping compass band becomes two pieces, so an intersection can be taken with
    ordinary arithmetic and the pieces summed.
    """
    if label in _CIRCULAR_BANDS and lo > hi:
        return [(lo, _FULL_CIRCLE_DEG), (0.0, hi)], (_FULL_CIRCLE_DEG - lo) + hi
    return [(lo, hi)], hi - lo


def weakest_component(arrival_scan, statistic="p50"):
    """
    Names the score component that held a site back, and its value.

    The score is a product of components each in [0, 1] and each named, so a low total
    can be attributed rather than merely reported. The lowest median component is the
    one to look at first: under a product composition it bounds the total from above.

    Parameters
    ----------
    arrival_scan : dict
        A site's ``arrival_scan`` record.
    statistic : str, optional
        Which per-site statistic to compare, ``"mean"``, ``"p50"`` or ``"p90"``.

    Returns
    -------
    tuple or None
        ``(name, value)`` for the lowest-scoring component, or ``None`` when the
        record carries no components -- which is the case for runs made before they
        were stored.

    Examples
    --------
    >>> from oroscope import explain
    >>> rec = {"score_p50": 0.2, "score_decay_p50": 0.9, "score_shower_p50": 0.22}
    >>> explain.weakest_component(rec)
    ('shower', 0.22)
    >>> explain.weakest_component({"score_p50": 0.2}) is None
    True
    """
    suffix = "_" + statistic
    found = {}
    for key, value in arrival_scan.items():
        if not key.startswith("score_") or not key.endswith(suffix):
            continue
        name = key[len("score_"):-len(suffix)]
        if name:                       # "score_p50" itself is the total, not a component
            found[name] = float(value)
    if not found:
        return None
    name = min(found, key=found.get)
    return name, found[name]


# ==========================================
#            SMALL FORMATTING HELPERS
# ==========================================

_WIDTH = 78


def _banner(text):
    """The one heavy rule at the top, so the sections under it can be lighter."""
    return ["=" * _WIDTH, f" {text}", "=" * _WIDTH]


def _heading(text):
    return [text, "-" * len(text)]


def _field(label, value, pad=16):
    return f"  {label.ljust(pad)} {value}"


def _num(value, places=0):
    """Thousands-separated, or a placeholder when the value is missing."""
    if value is None:
        return "not recorded"
    if places:
        return f"{value:,.{places}f}"
    return f"{value:,}"


def _wrap(text, indent="  ", width=_WIDTH):
    """Greedy wrap. textwrap would do, but this keeps the indent explicit."""
    words, lines, current = text.split(), [], indent
    for word in words:
        if len(current) + len(word) + 1 > width and current.strip():
            lines.append(current.rstrip())
            current = indent + word + " "
        else:
            current += word + " "
    if current.strip():
        lines.append(current.rstrip())
    return lines


def _get(params, *names, default=None):
    """First of ``names`` present in ``params``, searching the nested scan block too."""
    scan = params.get("scan") or {}
    for name in names:
        if name in params and params[name] is not None:
            return params[name]
        if name in scan and scan[name] is not None:
            return scan[name]
    return default


# ==========================================
#                 THE SECTIONS
# ==========================================


def _section_run(results, provenance):
    params = results.get("parameters", {}) or {}
    out = _heading("THE RUN")

    origin = params.get("origin") or [None, None]
    out.append(_field("DEM", params.get("dem", "not recorded")))
    if origin[0] is not None:
        out.append(_field("Origin", f"{origin[0]:.6f}, {origin[1]:.6f}"
                                    f"  ({params.get('origin_source', 'source unrecorded')})"))
    deg = params.get("cell_size_deg")
    if deg:
        out.append(_field("Resolution", f"{deg:.8f}°/px"
                                        f"  =  {params.get('cell_size_y_m', 0):.1f} m N-S"
                                        f" x {params.get('cell_size_x_m', 0):.1f} m E-W"))
    out.append(_field("Layout", f"{params.get('search_mode', 'unknown')} search,"
                                f" {params.get('grid_type', '?')} grid,"
                                f" {_get(params, 'spacing_km', 'antenna_spacing_km')} km spacing"))
    if results.get("timestamp"):
        # Sub-second precision on a wall clock is noise in a summary meant to be read.
        stamp = str(results["timestamp"]).split(".")[0].replace("T", " ")
        out.append(_field("Finished", stamp))

    prov = provenance or {}
    git = prov.get("git") or {}
    if git.get("commit"):
        state = "dirty tree" if git.get("dirty") else "clean tree"
        out.append(_field("Code", f"commit {git['commit'][:7]} on"
                                  f" {git.get('branch', '?')} ({state})"))
    dem = prov.get("dem") or {}
    if dem.get("sha256"):
        out.append(_field("DEM checksum", f"sha256 {dem['sha256'][:16]}…"))
    if prov.get("command"):
        out.append(_field("Command", prov["command"]))
    return out


def _section_headline(results):
    res = results.get("results", {}) or {}
    # The selection, not everything that qualified: total_capacity and the exported
    # mask cover the selection, so summing the full list disagrees with both.
    sites, rejected = selected_sites(results)
    params = results.get("parameters", {}) or {}
    total_area = sum(float(s.get("area_km2", 0.0)) for s in sites)
    capacity = res.get("total_capacity")
    target = _get(params, "target", "target_antennas")

    out = _heading("THE HEADLINE")
    if not sites:
        out += _wrap("No site met all the constraints. That is a result, not a "
                     "failure: read the funnel below, which names the stage where the "
                     "candidates ran out.")
        return out

    line = (f"{len(sites)} site{'s' if len(sites) != 1 else ''} covering "
            f"{total_area:,.1f} km²{' between them' if len(sites) != 1 else ''}")
    if isinstance(capacity, (int, float)):
        line += f", {_num(int(capacity))} detectors"
        if target:
            line += f" against a target of {_num(int(target))}"
    out += _wrap(line + ".")

    best = max(sites, key=lambda s: s.get("capacity_exact", 0))
    where = ""
    if best.get("center_lat") is not None:
        where = (f", centred {best['center_lat']:.4f}, {best['center_lon']:.4f}"
                 f" — paste that into a map")
    out.append("")
    out += _wrap(f"Largest by capacity: site {best.get('site_id')}, "
                 f"{float(best.get('area_km2', 0)):,.2f} km², "
                 f"{_num(int(best.get('capacity_exact', 0)))} detectors, "
                 f"facing {best.get('facing_direction', '?')}{where}.")

    if rejected:
        spare = sum(float(s.get("area_km2", 0.0)) for s in rejected)
        out.append("")
        out += _wrap(f"A further {len(rejected)} site"
                     f"{'s' if len(rejected) != 1 else ''} cleared the thresholds but "
                     f"{'were' if len(rejected) != 1 else 'was'} not selected, holding "
                     f"{spare:,.1f} km² more. Selection stopped at the target "
                     f"(stop_at_target), so those are the next best ground rather than "
                     f"ground that failed. They are in the results file, flagged "
                     f"selected: false, and they are not in the exported mask, the "
                     f"totals above, or this area.")
    if isinstance(capacity, (int, float)) and target and capacity < target:
        out.append("")
        out += _wrap(f"The target of {_num(int(target))} was not reached. The funnel "
                     f"below says why; the shortfall is {_num(int(target - capacity))} "
                     f"detectors, {100.0 * (1 - capacity / target):.0f}% of the target.")
    return out


def _section_funnel(results):
    funnel = results.get("funnel", {}) or {}
    out = _heading("WHERE THE CANDIDATES WENT")
    if not funnel:
        out += _wrap("No funnel was recorded for this run.")
        return out

    stages = list(funnel.items())
    total = max(stages[0][1], 1)
    width = max(len(name) for name, _ in stages)
    out.append(f"  {'stage'.ljust(width)} | {'pixels':>14} | {'of DEM':>8} | {'of prev':>8}")
    out.append("  " + "-" * (width + 39))
    prev = None
    for name, count in stages:
        of_dem = f"{100.0 * count / total:7.3f}%"
        of_prev = "        -" if prev is None else (
            f"{100.0 * count / prev:7.3f}%" if prev else "      n/a")
        out.append(f"  {name.ljust(width)} | {count:>14,} | {of_dem:>8} | {of_prev:>8}")
        prev = count

    out.append("")
    binding = binding_constraint(funnel)
    if binding is None:
        return out

    kept = 100.0 * binding["kept_fraction"]
    if binding["fatal"]:
        out += _wrap(f"The search died at “{binding['stage']}”: "
                     f"{_num(binding['before'])} pixels reached it and none survived. "
                     f"Nothing downstream of that stage could have helped.")
    else:
        out += _wrap(f"The binding constraint is “{binding['stage']}”. It kept "
                     f"{kept:.1f}% of the {_num(binding['before'])} pixels that "
                     f"reached it, a larger cut than any other stage made — so it, "
                     f"not the stages around it, is what set the size of this result.")
    if binding["knob"]:
        out.append("")
        out += _wrap(f"To move it, change: {binding['knob']}.")

    stride_stage = next((s for s in stages
                         if s[0].lower().startswith(_SUBSAMPLE_PREFIX)), None)
    if stride_stage:
        out.append("")
        out += _wrap(f"“{stride_stage[0]}” is not a constraint. Striding subsamples "
                     f"the candidates deliberately, and acceptance was measured "
                     f"identical at strides 1 and 5 on a GRAND-scale search. How well "
                     f"the gaps it leaves are closed again depends on the closing "
                     f"element, which the next section measures for this run.")
    return out


def _section_regions(results):
    regions = results.get("regions", {}) or {}
    if not regions:
        return []
    out = _heading("FROM PIXELS TO SITES")
    out += _wrap(f"{_num(regions.get('labelled_regions'))} connected regions were "
                 f"labelled. {_num(regions.get('passed_area_threshold'))} were large "
                 f"enough (at least {_num(regions.get('required_pixels_per_region'))} "
                 f"pixels), {_num(regions.get('passed_capacity_threshold'))} then held "
                 f"enough detectors (at least "
                 f"{_num(regions.get('capacity_threshold_antennas'))}), and "
                 f"{_num(regions.get('selected'))} were selected.")

    labelled = regions.get("labelled_regions") or 0
    passed = regions.get("passed_area_threshold") or 0
    if labelled and passed < 0.01 * labelled:
        out.append("")
        out += _wrap("Almost every region was too small. That is the signature of an "
                     "accepted set that is scattered rather than contiguous — worth "
                     "checking min_width_km and gap_close_km before reading the area.")
    return out


def _section_sites(results, max_rows=25):
    sites, _ = selected_sites(results)
    if not sites:
        return []
    out = _heading("THE SITES")
    header = (f"  {'id':>6} | {'area km²':>9} | {'capacity':>9} | {'facing':>6}"
              f" | {'score':>6} | weakest component")
    out.append(header)
    out.append("  " + "-" * (len(header) - 2))
    for site in sites[:max_rows]:
        scan = site.get("arrival_scan") or {}
        score = scan.get("score_p50")
        weak = weakest_component(scan)
        weak_text = f"{weak[0]} ({weak[1]:.2f})" if weak else "—"
        out.append(f"  {str(site.get('site_id', '?')):>6}"
                   f" | {float(site.get('area_km2', 0)):>9.2f}"
                   f" | {int(site.get('capacity_exact', 0)):>9,}"
                   f" | {str(site.get('facing_direction', '?')):>6}"
                   f" | {(f'{score:.3f}' if score is not None else '—'):>6}"
                   f" | {weak_text}")
    if len(sites) > max_rows:
        out.append(f"  … and {len(sites) - max_rows} more, in the results JSON.")

    # One line of attribution across the whole selection, which is usually the more
    # actionable statement: if every site is held back by the same component, that
    # component is the result.
    weakest = [weakest_component(s.get("arrival_scan") or {}) for s in sites]
    named = [w[0] for w in weakest if w]
    if named:
        commonest = max(set(named), key=named.count)
        share = 100.0 * named.count(commonest) / len(named)
        out.append("")
        out += _wrap(f"“{commonest}” is the weakest component at "
                     f"{named.count(commonest)} of {len(named)} sites ({share:.0f}%). "
                     f"Under a product composition the weakest component bounds the "
                     f"total from above, so that is where a better site would have to "
                     f"come from.")
    elif sites:
        out.append("")
        out += _wrap("This run recorded no per-component scores, so its sites cannot "
                     "be attributed. Re-running on current code stores them.")
    return out


def _section_why_good(results, max_sites=6):
    """Per site: what the ground is actually like, and which criteria it satisfies."""
    sites, _ = selected_sites(results)
    if not sites:
        return []
    if not any(site_strengths(s.get("arrival_scan") or {}) for s in sites):
        return []

    out = _heading("WHY THESE SITES QUALIFY")
    out += _wrap("What the terrain at each one actually offers. The score is a product "
                 "of named criteria, so a site can be described rather than merely "
                 "ranked — these are the criteria it satisfies, with the measurement "
                 "that earned each.")

    for site in sites[:max_sites]:
        scan = site.get("arrival_scan") or {}
        strong = site_strengths(scan)
        if not strong:
            continue
        out.append("")
        header = (f"  Site {site.get('site_id')} — "
                  f"{float(site.get('area_km2', 0)):,.2f} km², "
                  f"{int(site.get('capacity_exact', 0)):,} detectors, "
                  f"facing {site.get('facing_direction', '?')}")
        if site.get("center_lat") is not None:
            header += f", centred {site['center_lat']:.4f}, {site['center_lon']:.4f}"
        out.append(header)

        # The measured geometry, before any scoring: the reader's own check on it.
        facts = []
        for field, fmt in (("solid_angle_sr", lambda v: f"{v:.2f} sr of accepted sky"),
                           ("mean_distance_m", lambda v: f"targets at {v:,.0f} m"),
                           ("target_slope_deg", lambda v: f"striking {v:.0f}° terrain"),
                           ("max_depth_gcm2", lambda v: f"{v:,.0f} g/cm² of rock behind"),
                           ("altitude_m", lambda v: f"at {v:,.0f} m altitude")):
            value = scan.get(f"{field}_p50")
            if value is not None:
                facts.append(fmt(float(value)))
        if facts:
            out += _wrap("Measured: " + ", ".join(facts) + ".", indent="      ")

        out.append(f"      Satisfies {len(strong)} of "
                   f"{len([k for k in scan if k.startswith('score_') and k.endswith('_p50') and k != 'score_p50'])}"
                   f" criteria:")
        for entry in strong:
            evidence = f" ({entry['evidence']})" if "evidence" in entry else ""
            out.append(f"        • {entry['label']}{evidence} — {entry['score']:.2f}")

        weak = weakest_component(scan)
        if weak and weak[1] < 0.999:
            label = COMPONENT_MEANING.get(weak[0], (weak[0].replace("_", " "),))[0]
            out += _wrap(f"Held back by {label} at {weak[1]:.2f}, which under a product "
                         f"composition is what bounds the total.", indent="      ")

    if len(sites) > max_sites:
        out.append("")
        out += _wrap(f"The remaining {len(sites) - max_sites} sites are in the results "
                     f"file, each with the same per-criterion record.")

    # What the selection has in common, which is a statement about the terrain rather
    # than about any one site.
    shared = None
    for site in sites:
        names = {e["name"] for e in site_strengths(site.get("arrival_scan") or {})}
        shared = names if shared is None else (shared & names)
    if shared:
        labels = sorted(COMPONENT_MEANING.get(n, (n.replace("_", " "),))[0]
                        for n in shared)
        out.append("")
        out += _wrap(f"Every selected site satisfies: {', '.join(labels)}. That is a "
                     f"statement about this terrain, not about any one site — and if a "
                     f"criterion is satisfied everywhere, it is not discriminating "
                     f"between sites and its threshold is worth checking.")
    return out


def _section_reading(results):
    params = results.get("parameters", {}) or {}
    out = _heading("HOW TO READ THESE NUMBERS")

    gap = params.get("gap_close_km")
    spacing = _get(params, "spacing_km", "antenna_spacing_km")
    closing = gap if gap is not None else spacing
    if closing:
        out += _wrap(f"Reported area is not physics-accepted area. The mask is closed "
                     f"morphologically with a {closing} km element before areas are "
                     f"measured. Closing is not wrong — a site has to be a deployable "
                     f"region rather than a scatter of pixels — but it moves the "
                     f"number, and at Colca with a 1 km element it inflated the "
                     f"reported area by {AREA_INFLATION_AT_COLCA:.2f}× against a "
                     f"stride-1 control. Set gap_close_km to 0 to see the accepted "
                     f"area directly.")

        ratio = closing_inflation(results.get("funnel", {}) or {},
                                  params.get("candidate_stride", 1))
        if ratio is not None:
            out.append("")
            if ratio >= 1.05:
                out += _wrap(f"This run's own funnel puts that factor at "
                             f"{ratio:.2f}×: closing turned the stride-corrected "
                             f"accepted set into a mask that much larger. Divide the "
                             f"areas above by it for what the physics accepted.")
            elif ratio <= 0.95:
                out += _wrap(f"This run's own funnel puts that factor at "
                             f"{ratio:.2f}× — below 1, so here closing did not "
                             f"inflate. The element is only a few pixels across at "
                             f"this spacing, too small to bridge the gaps "
                             f"candidate_stride left, so the area above understates "
                             f"the accepted set rather than overstating it. A "
                             f"stride-1 control run is the way to settle it.")
            else:
                out += _wrap(f"This run's own funnel puts that factor at "
                             f"{ratio:.2f}×, so closing barely moved this result.")
        out.append("")

    ds = params.get("downsample_factor") or 1
    if ds > 1:
        out += _wrap(f"Area and capacity are measured on different grids. Area comes "
                     f"from the map downsampled by {ds}, capacity from the "
                     f"full-resolution mask, so a feature a few pixels wide loses area "
                     f"while keeping its detectors — around 30% for a canyon strip. "
                     f"Use downsample_factor 1 for thin features.")
        out.append("")

    out += _wrap("The scores rank sites against each other for one experiment and one "
                 "energy band. They are not apertures, and no flux or detector "
                 "response has been folded in.")
    out.append("")
    out += _wrap("The layout is anchored, not fitted: detectors are placed from each "
                 "site's bounding-box corner rather than optimised, so capacity is an "
                 "estimate for an arbitrarily placed array.")
    return out


def _section_aperture(results):
    """The energy the geometry favours, which no other line of the summary carries."""
    aperture = results.get("aperture", {}) or {}
    energies = aperture.get("energies_pev")
    total = aperture.get("total_m2sr")
    if not energies or not total:
        return []

    out = _heading("WHAT ENERGY THIS GEOMETRY FAVOURS")
    peak = aperture.get("peak_energy_pev")
    best = max(range(len(total)), key=lambda i: total[i])
    out += _wrap(f"Folding each site's accepted solid angle and area against the tau "
                 f"decay length gives a geometric aperture that peaks near "
                 f"{energies[best]:,.0f} PeV"
                 + (f" (analytically {peak:,.0f} PeV)" if peak else "")
                 + f", at {total[best]:,.3g} m² sr.")
    out.append("")

    # A coarse profile: enough to see the shape without a plotting library.
    scale = max(total) or 1.0
    for i in range(0, len(energies), max(1, len(energies) // 8)):
        bar = "#" * int(round(40 * total[i] / scale))
        out.append(f"  {energies[i]:>10,.0f} PeV  {total[i]:>10.3g}  {bar}")

    out.append("")
    out += _wrap("Relative, not absolute: this is geometry and an analytic decay "
                 "factor, with no flux, no cross-section and no detector response. "
                 "The detector acceptance A(E) is the one thing no available table "
                 "supplies, and it is why these are apertures in shape only.")
    return out


def _section_next(results):
    """Concrete things to do next, chosen from what this run actually did."""
    params = results.get("parameters", {}) or {}
    funnel = results.get("funnel", {}) or {}
    sites, _ = selected_sites(results)
    suggestions = []

    binding = binding_constraint(funnel)
    if binding and binding.get("knob"):
        first = binding["knob"].split(" ")[0].split("/")[0].strip(",")
        suggestions.append(
            (f"Test the constraint that bound this run — “{binding['stage']}”",
             f"oroscope-sensitivity <config> --sweep {first} <values>"))

    if params.get("min_score") and params.get("score_percentile") is None:
        suggestions.append(
            ("Replace the absolute score cut with a rank, which is scale-free",
             f"oroscope <config> --score_percentile 10   "
             f"(instead of --min_score {params['min_score']:g})"))

    gap = params.get("gap_close_km")
    if gap is None or gap:
        suggestions.append(
            ("See the area the physics accepted, without morphological closing",
             "oroscope <config> --gap_close_km 0"))

    if (params.get("candidate_stride") or 1) > 1:
        suggestions.append(
            ("Confirm the stride is unbiased at this closing element",
             f"oroscope <config> --candidate_stride 1   "
             f"(currently {params['candidate_stride']})"))

    if not sites:
        suggestions.append(
            ("Loosen the binding constraint before anything else — nothing survived",
             "read the funnel above, then change the parameter it names"))

    if not suggestions:
        return []

    out = _heading("WHAT TO TRY NEXT")
    out += _wrap("Chosen from what this run did, not a generic list.")
    for why, how in suggestions:
        out.append("")
        out += _wrap(why, indent="  ")
        out.append(f"      {how}")
    return out


def _section_assumptions(results):
    params = results.get("parameters", {}) or {}
    out = _heading("WHICH OF THESE ARE ASSUMPTIONS")
    out += _wrap("Choices, not measurements. Check them before quoting a result.")
    out.append("")

    items = []

    percentile = params.get("score_percentile")
    min_score = params.get("min_score")
    if percentile is not None:
        items.append(("score_percentile", f"{percentile:g}",
                      "Rank-based, so scale-free. The preferred form: it does not "
                      "move when the composition or the number of components changes."))
    elif min_score:
        items.append(("min_score", f"{min_score:g}",
                      "A cut on a product of components, whose distribution piles up "
                      "near zero, so any threshold in the middle sits on a cliff. "
                      "Measured on one search: 0.0, 0.35 and 0.5 gave 45928, 2056 and "
                      "zero detectors. This is the dominant assumption — prefer "
                      "--score_percentile."))

    index = params.get("decay_spectral_index")
    if index is not None:
        pretty = (f"{index[0]:g}–{index[1]:g}, marginalised"
                  if isinstance(index, (list, tuple)) and len(index) == 2
                  else f"{index if not isinstance(index, (list, tuple)) else index[0]:g}, pinned")
        items.append(("decay_spectral_index", pretty,
                      "The flux slope the decay term is folded against. Capacity "
                      "varies by 1.46× across a plausible range."))
    elif params.get("decay_energy_pev"):
        items.append(("decay_energy_pev", f"{params['decay_energy_pev']:g} PeV",
                      "A single energy for a term that runs over three decades: the "
                      "decay length goes from 147 m at 3 PeV to 49 km at 1 EeV. This "
                      "chooses the answer rather than approximating it — give "
                      "decay_energy_min_pev/max_pev instead and fold over a spectrum."))

    slope = _get(params, "min_target_slope_deg")
    if slope is not None:
        items.append(("min_target_slope_deg", f"{slope:g}°",
                      "The floor separating a canyon wall from a hillside. Strongly "
                      "selective: 0°, 25° and 35° gave 7442, 2056 and zero detectors."))

    frac = params.get("grammage_band_fraction")
    if frac is not None:
        items.append(("grammage_band_fraction", f"{frac:g}",
                      "How far down the shower profile still counts as usable. A "
                      "claim about detector capability, not about the shower."))

    decl = params.get("geomag_declination_deg")
    if decl is not None and _get(params, "use_geomagnetic", default=True):
        items.append(("geomag_declination_deg", f"{decl:g}°",
                      "Constant across the DEM unless a declination model was supplied "
                      "— the centred dipole that gives inclination is unreliable for "
                      "declination (-0.2° against a measured -6.9° at Arequipa), so it "
                      "falls back to Arequipa's IGRF value. Right for southern Peru, "
                      "wrong elsewhere. physics.set_declination_model() takes any "
                      "callable, and declination_from_grid() builds one from a NOAA "
                      "export."))

    if not _get(params, "max_range_km"):
        items.append(("max_range_km", "unset",
                      "Column depth accumulates over the whole profile walk, which "
                      "then stops at max_dist_km — so the reported depth is a property "
                      "of where the walk stopped, not of the target's thickness. "
                      "Measured on TAMBO at Colca: walking to 20 km instead of the "
                      "5 km distance window raised the reported depth 6.4× and changed "
                      "the selection not at all. Read this depth as a lower bound. "
                      "Do not simply set it large, though — at 60 km the same run kept "
                      "only 6.0% of directions against 17.5%, so the walk length is a "
                      "parameter to check rather than to maximise."))


    for name, value, why in items:
        out.append(f"  • {name} = {value}")
        out += _wrap(why, indent="      ")
    out.append("")
    out += _wrap("Not modelled at all: tau production and escape through rock (so β, "
                 "the energy-loss constant, does not enter these numbers — the decay "
                 "length used here is kinematics, E/m·cτ, and carries no β), "
                 "neutral-current regeneration (so Earth-chord suppression is "
                 "overstated), shower simulation, detector response and trigger, and "
                 "any geology beyond one standard rock density. Nothing here has been "
                 "checked against an external simulation.")
    return out


def explain_combination(report, runs=None):
    """
    Explains an overlay of two or more searches: who can share ground with whom, and why.

    The combined report gives a joint area and a Jaccard index. Neither says *why* the
    number is what it is, and the reason is usually not about neutrinos at all: a pixel
    has one slope, and every experiment deployed on it must accept that slope. Measured
    at Colca, the entire co-location result follows from GRAND's 3-25 degree deployable
    band against a canyon's ~40 degree walls.

    So where the runs' own parameters are available, this compares their screening
    bands and names the one that limits the sharing.

    Parameters
    ----------
    report : dict
        A ``combined_report.json``, as :mod:`combine_experiments` writes it.
    runs : dict, optional
        Label to results dictionary, for the runs being combined. Supplying them adds
        the constraint comparison; without them only the areas are explained.

    Returns
    -------
    str
        The summary, as plain text.

    Examples
    --------
    >>> from oroscope import explain
    >>> text = explain.explain_combination({
    ...     "runs": [{"label": "A", "area_km2": 100.0, "pixels": 10,
    ...               "area_in_joint_km2": 20.0, "fraction_of_own_area_in_joint": 0.2,
    ...               "reported_sites": 1, "reported_capacity": 500}],
    ...     "joint": {"area_km2": 20.0}, "union": {"area_km2": 100.0},
    ...     "joint_requires": ["A"], "pairwise_overlap": {}})
    >>> "WHERE THESE EXPERIMENTS CAN SHARE GROUND" in text
    True
    """
    if not isinstance(report, dict):
        raise TypeError("explain_combination takes a combined report dictionary, "
                        f"not {type(report).__name__}")
    runs = runs or {}
    entries = report.get("runs", []) or []
    lines = _banner("WHERE THESE EXPERIMENTS CAN SHARE GROUND, AND WHY")

    # --- what each brings
    lines += ["", *_heading("WHAT EACH BRINGS")]
    width = max([len(e.get("label", "?")) for e in entries] + [10])
    lines.append(f"  {'experiment'.ljust(width)} | {'area km²':>12} | {'sites':>7}"
                 f" | {'detectors':>10} | {'in joint':>9}")
    lines.append("  " + "-" * (width + 48))
    for entry in entries:
        sites = entry.get("reported_sites")
        capacity = entry.get("reported_capacity")
        lines.append(f"  {str(entry.get('label', '?')).ljust(width)}"
                     f" | {float(entry.get('area_km2', 0)):>12,.1f}"
                     f" | {(f'{sites:,}' if sites is not None else '—'):>7}"
                     f" | {(f'{capacity:,}' if capacity is not None else '—'):>10}"
                     f" | {100 * float(entry.get('fraction_of_own_area_in_joint', 0)):>8.1f}%")

    joint = float((report.get("joint") or {}).get("area_km2", 0.0))
    union = float((report.get("union") or {}).get("area_km2", 0.0))
    required = report.get("joint_requires") or [e.get("label") for e in entries]
    lines.append("")
    lines += _wrap(f"Ground satisfying {' and '.join(str(r) for r in required)} at once: "
                   f"{joint:,.1f} km². Ground useful to any of them: {union:,.1f} km². "
                   f"The first is the number that matters for one site, one road and "
                   f"one power feed serving two experiments; the second is what the "
                   f"programme as a whole could use.")

    # --- why the joint is the size it is
    labels = [e.get("label") for e in entries]
    if len(labels) == 2 and all(lab in runs for lab in labels):
        a, b = labels
        bands = constraint_overlap((runs[a].get("parameters") or {}),
                                   (runs[b].get("parameters") or {}))
        if bands:
            lines += ["", *_heading("WHAT DECIDES IT")]
            lines += _wrap("A pixel has one slope, one altitude and one aspect, and "
                           "both experiments must accept those same values to stand on "
                           "it. So co-location is settled by whichever of those bands "
                           "they share least of — before any arrival geometry is "
                           "considered.")
            lines.append("")
            lines.append(f"  {'band'.ljust(18)} | {a:>16} | {b:>16} | shared")
            lines.append("  " + "-" * 68)
            for band in bands:
                unit = band["unit"]
                a_s = f"{band['a'][0]:g}–{band['a'][1]:g}{unit}"
                b_s = f"{band['b'][0]:g}–{band['b'][1]:g}{unit}"
                if band["overlap"]:
                    shared = (f"{band['overlap'][0]:g}–{band['overlap'][1]:g}{unit}"
                              f"  ({100 * band['share_of_narrower']:.0f}%)")
                else:
                    shared = "none — disjoint"
                lines.append(f"  {band['label'].ljust(18)} | {a_s:>16} | {b_s:>16}"
                             f" | {shared}")

            tightest = min(bands, key=lambda x: x["share_of_narrower"])
            lines.append("")
            if not tightest["overlap"]:
                lines += _wrap(f"Their {tightest['label']} bands are disjoint, so no "
                               f"pixel can satisfy both. Any joint area reported above "
                               f"is then an artefact worth investigating.")
            else:
                lines += _wrap(
                    f"“{tightest['label']}” is what limits the sharing: they overlap "
                    f"only over {tightest['overlap'][0]:g}–{tightest['overlap'][1]:g}"
                    f"{tightest['unit']}, which is "
                    f"{100 * tightest['share_of_narrower']:.0f}% of the narrower of the "
                    f"two bands. Co-location is decided there, not by anything about "
                    f"the physics of either experiment.")

            # The viewing criteria, explicitly *not* a conflict. Worth stating, because
            # they look like conflicting requirements written side by side and are not.
            views = constraint_overlap((runs[a].get("parameters") or {}),
                                       (runs[b].get("parameters") or {}),
                                       bands=_VIEWING_BANDS)
            if views:
                lines.append("")
                lines += _wrap("What they ask of the *view* differs, and that is no "
                               "obstacle: two experiments can look out from the same "
                               "hillside at different ranges and different elevations "
                               "without conflict.")
                lines.append("")
                for band in views:
                    unit = band["unit"]
                    lines.append(
                        f"  {band['label'].ljust(18)} | "
                        f"{a}: {band['a'][0]:g}–{band['a'][1]:g}{unit}, "
                        f"{b}: {band['b'][0]:g}–{band['b'][1]:g}{unit}")

    # --- the overlaps themselves
    pairwise = report.get("pairwise_overlap") or {}
    if pairwise:
        lines += ["", *_heading("PAIRWISE")]
        for pair, stats in pairwise.items():
            lines += _wrap(f"{pair}: {float(stats.get('area_km2', 0)):,.1f} km² shared, "
                           f"Jaccard {float(stats.get('jaccard', 0)):.4f}.")
            for key, value in stats.items():
                if key.startswith("fraction_of_"):
                    lines.append(f"      {100 * float(value):5.1f}% of "
                                 f"{key[len('fraction_of_'):]}'s own ground")

    lines += ["", *_heading("HOW TO READ THIS")]
    lines += _wrap("Every caveat on the individual runs applies here and compounds: "
                   "these are the same masks, so morphological closing has already "
                   "moved each area, and the joint of two inflated masks is inflated "
                   "twice over. The overlay is exact — the masks are pixel-aligned and "
                   "the alignment is checked, not assumed — but what it overlays is "
                   "only as good as each run.")
    lines.append("")
    lines += _wrap("And co-location is a question about ground, not about detectors. "
                   "Two experiments sharing a hillside still need their own arrays, "
                   "their own spacing and their own trigger; what they share is the "
                   "site, the access and the power.")
    return "\n".join(lines) + "\n"


def explain_results(results, provenance=None):
    """
    Writes a human-readable account of one search: what was found, and why.

    Everything it reports is already in the results dictionary. The value added is the
    reading: which stage of the funnel actually set the size of the answer, what held
    the surviving sites back, and which of the numbers are assumptions. Those are the
    three things a reader gets wrong without help, and a run that is going to be handed
    to someone else needs all three on the same page as the result.

    Pure: it opens no files, runs nothing, and needs no DEM. An old run's JSON can be
    explained months later, and the pipeline, the library and the tests all get the
    same words.

    Parameters
    ----------
    results : dict
        A results dictionary as returned by
        :func:`site_searcher.find_grand_regions_interactive` and written to the run's
        results JSON. Missing sections are tolerated: each is reported as absent rather
        than raising, so a partial or older file still explains.
    provenance : dict, optional
        The matching ``provenance.json`` contents -- git commit, DEM checksum, package
        versions. Adds the reproducibility block when given.

    Returns
    -------
    str
        The summary, as plain text with no ANSI colour, ready to print or to save
        beside the results.

    Examples
    --------
    >>> from oroscope import explain
    >>> text = explain.explain_results({
    ...     "funnel": {"DEM pixels": 1000, "slope 3-25 deg": 900,
    ...                "directions accepted": 12},
    ...     "results": {"total_sites": 0, "total_capacity": 0, "sites": []}})
    >>> "WHERE THE CANDIDATES WENT" in text
    True
    >>> "directions accepted" in text
    True
    """
    if not isinstance(results, dict):
        raise TypeError("explain_results takes the results dictionary, "
                        f"not {type(results).__name__}")

    blocks = [
        _banner("WHAT THIS SEARCH FOUND, AND WHY"),
        _section_run(results, provenance),
        _section_headline(results),
        _section_funnel(results),
        _section_regions(results),
        _section_sites(results),
        _section_why_good(results),
        _section_aperture(results),
        _section_reading(results),
        _section_assumptions(results),
        _section_next(results),
    ]
    lines = []
    for block in blocks:
        if not block:
            continue
        if lines:
            lines.append("")
        lines.extend(block)
    return "\n".join(lines) + "\n"
