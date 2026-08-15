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

__all__ = ["explain_results", "binding_constraint", "weakest_component",
           "closing_inflation", "STAGE_KNOBS", "AREA_INFLATION_AT_COLCA"]

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
    >>> import explain
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
    >>> import explain
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
    >>> import explain
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
    sites = res.get("sites", []) or []
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

    line = (f"{len(sites)} site{'s' if len(sites) != 1 else ''}, "
            f"{total_area:,.1f} km² of mapped area summed over them")
    if isinstance(capacity, (int, float)):
        line += f", {_num(int(capacity))} detectors"
        if target:
            line += f" against a target of {_num(int(target))}"
    out += _wrap(line + ".")

    best = max(sites, key=lambda s: s.get("capacity_exact", 0))
    out.append("")
    out += _wrap(f"Largest by capacity: site {best.get('site_id')}, "
                 f"{float(best.get('area_km2', 0)):,.2f} km², "
                 f"{_num(int(best.get('capacity_exact', 0)))} detectors, "
                 f"facing {best.get('facing_direction', '?')}.")
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
    sites = (results.get("results", {}) or {}).get("sites", []) or []
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
                      "Does NOT follow the site — it falls back to Arequipa's IGRF "
                      "value wherever the DEM is. Inclination does follow, through a "
                      "centred dipole. Supply the IGRF declination per site."))

    if not _get(params, "max_range_km"):
        items.append(("max_range_km", "unset",
                      "Column depth accumulates over the whole profile walk, which "
                      "then stops at max_dist_km — so the reported depth is a property "
                      "of where the walk stopped, not of the target's thickness."))

    items.append(("β, the tau energy-loss constant", "0.6e-6 cm²/g",
                  "Estimated from mass scaling, in the range (0.4–1.0)e-6. Moves the "
                  "production-and-escape optimum in proportion. Not yet pinned to a "
                  "collaboration value."))

    for name, value, why in items:
        out.append(f"  • {name} = {value}")
        out += _wrap(why, indent="      ")
    out.append("")
    out += _wrap("Not modelled at all: neutral-current regeneration (so Earth-chord "
                 "suppression is overstated), shower simulation, detector response and "
                 "trigger, and any geology beyond one standard rock density. Nothing "
                 "here has been checked against an external simulation.")
    return out


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
    >>> import explain
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
        _section_reading(results),
        _section_assumptions(results),
    ]
    lines = []
    for block in blocks:
        if not block:
            continue
        if lines:
            lines.append("")
        lines.extend(block)
    return "\n".join(lines) + "\n"
