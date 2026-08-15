"""
Oroscope: terrain site search for particle-astrophysics observatories.

Greek *oros*, mountain, and *skopein*, to look at. Oroscope searches digital elevation
models for ground that could host an observatory, and answers one structural question —
*from this patch of ground, is there a target surface at the right range, in the right
direction, at the right relative orientation, with the right matter behind it?* That one
question is what lets a single engine serve experiments that look nothing alike: GRAND
and TAMBO are configurations, not code paths.

Run a search::

    import oroscope

    results = oroscope.find_grand_regions_interactive(
        dem_path="input/dem/colca.tif", run_output_dir="output/colca",
        min_slope_deg=3.0, max_slope_deg=25.0,
        min_dist_km=10.0, max_dist_km=40.0)

    print(results["results"]["total_capacity"])
    print(results["explanation"])

Everything the ``oroscope`` command line can do, this library can do — the console
scripts are argument parsing and file placement over these same functions. The names
below are the whole public surface; the submodules they come from stay importable when
you want a narrower namespace::

    from oroscope import physics, explain
    physics.tau_decay_length_m(100.0)

Where to start:

``find_grand_regions_interactive``
    The pipeline. Screens terrain, scans arrival directions, scores, cleans up, labels
    sites, packs detectors, writes everything, and **returns its results**.
``load_config``, ``default_config``, ``generate_config``
    Configuration as data. Prefer starting from ``default_config()`` and overriding:
    the pipeline signature's own defaults differ from the template's for five
    parameters.
``explain_results``
    What a run found and why, as text. A pure function of the results dictionary, so an
    old run can be explained from its JSON with no DEM and nothing re-run.
``preflight_memory``
    What a search will cost, and a cap so one that outgrows the machine fails with
    ``MemoryError`` rather than inviting the OOM killer to choose a victim.
"""

from __future__ import annotations

__version__ = "0.1.0"

# Submodules, so `from oroscope import physics` and `oroscope.physics.x` both work.
# figures, fetch_dem and generate_env are deliberately not imported here: the first
# pulls in a plotting stack for pictures the library itself never needs, and the other
# two are one-shot setup tools. All three remain importable by name.
from oroscope import (  # noqa: F401
    aperture,
    arrival_scan,
    combine_experiments,
    crop_dem,
    explain,
    physics,
    scoring,
    sensitivity,
    site_searcher,
)

# --- the pipeline, and everything needed to drive it
from oroscope.site_searcher import (  # noqa: F401
    find_grand_regions_interactive,
    main,
    # configuration
    load_config, generate_config, default_config, CONFIG_PRESETS,
    validate_parameters, parse_score_weights, explicitly_passed,
    # memory
    preflight_memory, estimate_peak_memory_gb, apply_memory_cap, available_memory_gb,
    # geometry and the DEM
    resolve_grid_geometry, read_dem_geometry, read_dem_origin, resolve_origin,
    build_elevation_cache, load_dem_and_init_buffers, MapGrid,
    # terrain and screening
    terrain_gradients, terrain_derivatives, slope_band_gradient_sq,
    slope_baseline_pixels, get_candidates_chunked,
    # the scan, morphology and capacity
    run_arrival_scan, summarize_observables_by_site, clean_shape_artifacts,
    apply_morphology_pingpong, separable_closing, separable_opening,
    analyze_sites_and_capacity, count_grid_capacity,
    # outputs and accounting
    create_world_file, generate_kml_file, generate_visualizations_and_outputs,
    collect_provenance, emit_explanation, Funnel,
    find_results_json, RESULTS_PREFIX, LEGACY_RESULTS_PREFIX,
    is_point_in_poly, apply_poly_mask_numba,
)

# --- reading a run
from oroscope.explain import (  # noqa: F401
    explain_results, explain_combination, binding_constraint, weakest_component,
    site_strengths, constraint_overlap, closing_inflation, selected_sites,
    COMPONENT_MEANING, STAGE_KNOBS, AREA_INFLATION_AT_COLCA,
)

# --- the closed-form physics, usable with no terrain at all
from oroscope.physics import (  # noqa: F401
    air_density_kgm3, slant_grammage_gcm2, shower_maturity, shower_maximum_gcm2,
    shower_size_fraction, grammage_band_from_energy, earth_chord_m, earth_chord_gcm2,
    neutrino_survival, muon_shielding_gcm2, tau_decay_length_m, cc_cross_section_cm2,
    neutrino_interaction_length_gcm2, tau_energy_loss_beta, tau_range_gcm2,
    tau_survival, tau_exit_probability, production_escape_optimum_gcm2,
    depth_band_from_energy, earth_absorption_cutoff_deg,
    spectrum_weighted_decay_probability, geomagnetic_latitude_deg,
    centered_dipole_inclination, default_field_for_site, geomagnetic_unit_vector,
    geomagnetic_sin_alpha, refractivity, cherenkov_angle_rad,
    cherenkov_footprint_radius_m, footprint_sampling,
)

# --- the scan kernel
from oroscope.arrival_scan import (  # noqa: F401
    scan, rfi_exposure, earth_radius_for_k, azimuth_fan, balanced_order,
    energy_pev_for_decay_length, decay_probability, distance_window_from_energy,
    STANDARD_ROCK_DENSITY, TRUE_EARTH_RADIUS_M, RADIO_EARTH_RADIUS_M,
)

# --- scoring, and the aperture it feeds
from oroscope.scoring import (  # noqa: F401
    band_score, saturating_score, ramp_score, compose, score_candidates,
    summarize_scores, COMPOSITION_MODES, DEFAULT_SCORE_CONFIG,
)
from oroscope.aperture import (  # noqa: F401
    unit_response, TabulatedResponse, geometric_aperture_m2sr, aperture_vs_energy,
    peak_energy_pev, infer_response, load_curve_csv, summarize_sites,
)

# --- the tools, as functions
from oroscope.crop_dem import crop, read_geo  # noqa: F401
from oroscope.combine_experiments import (  # noqa: F401
    load_run, check_alignment, read_world_file, pixel_area_km2, capacity_of,
)
from oroscope.sensitivity import run_once, summarise  # noqa: F401

__all__ = [
    "__version__",
    # submodules
    "aperture", "arrival_scan", "combine_experiments", "crop_dem", "explain",
    "physics", "scoring", "sensitivity", "site_searcher",
    # the pipeline
    "find_grand_regions_interactive", "main",
    "load_config", "generate_config", "default_config", "CONFIG_PRESETS",
    "validate_parameters", "parse_score_weights", "explicitly_passed",
    "preflight_memory", "estimate_peak_memory_gb", "apply_memory_cap",
    "available_memory_gb",
    "resolve_grid_geometry", "read_dem_geometry", "read_dem_origin", "resolve_origin",
    "build_elevation_cache", "load_dem_and_init_buffers", "MapGrid",
    "terrain_gradients", "terrain_derivatives", "slope_band_gradient_sq",
    "slope_baseline_pixels", "get_candidates_chunked",
    "run_arrival_scan", "summarize_observables_by_site", "clean_shape_artifacts",
    "apply_morphology_pingpong", "separable_closing", "separable_opening",
    "analyze_sites_and_capacity", "count_grid_capacity",
    "create_world_file", "generate_kml_file", "generate_visualizations_and_outputs",
    "collect_provenance", "emit_explanation", "Funnel",
    "find_results_json", "RESULTS_PREFIX", "LEGACY_RESULTS_PREFIX",
    "is_point_in_poly", "apply_poly_mask_numba",
    # reading a run
    "explain_results", "explain_combination", "binding_constraint",
    "weakest_component", "site_strengths", "constraint_overlap", "closing_inflation",
    "selected_sites", "COMPONENT_MEANING", "STAGE_KNOBS", "AREA_INFLATION_AT_COLCA",
    # physics
    "air_density_kgm3", "slant_grammage_gcm2", "shower_maturity",
    "shower_maximum_gcm2", "shower_size_fraction", "grammage_band_from_energy",
    "earth_chord_m", "earth_chord_gcm2", "neutrino_survival", "muon_shielding_gcm2",
    "tau_decay_length_m", "cc_cross_section_cm2", "neutrino_interaction_length_gcm2",
    "tau_energy_loss_beta", "tau_range_gcm2", "tau_survival", "tau_exit_probability",
    "production_escape_optimum_gcm2", "depth_band_from_energy",
    "earth_absorption_cutoff_deg", "spectrum_weighted_decay_probability",
    "geomagnetic_latitude_deg", "centered_dipole_inclination", "default_field_for_site",
    "geomagnetic_unit_vector", "geomagnetic_sin_alpha", "refractivity",
    "cherenkov_angle_rad", "cherenkov_footprint_radius_m", "footprint_sampling",
    # the scan
    "scan", "rfi_exposure", "earth_radius_for_k", "azimuth_fan", "balanced_order",
    "energy_pev_for_decay_length", "decay_probability", "distance_window_from_energy",
    "STANDARD_ROCK_DENSITY", "TRUE_EARTH_RADIUS_M", "RADIO_EARTH_RADIUS_M",
    # scoring and aperture
    "band_score", "saturating_score", "ramp_score", "compose", "score_candidates",
    "summarize_scores", "COMPOSITION_MODES", "DEFAULT_SCORE_CONFIG",
    "unit_response", "TabulatedResponse", "geometric_aperture_m2sr",
    "aperture_vs_energy", "peak_energy_pev", "infer_response", "load_curve_csv",
    "summarize_sites",
    # the tools
    "crop", "read_geo",
    "load_run", "check_alignment", "read_world_file", "pixel_area_km2", "capacity_of",
    "run_once", "summarise",
]
