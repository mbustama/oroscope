The command line
================

Five console scripts, installed by ``pip install oroscope``, and the **complete option
reference** for the search itself at the bottom of this page.

The command line is a convenience, not the primary interface. Everything here is a thin
wrapper over a function, and :doc:`quickstart` shows the same work done from Python,
which is how most of it gets used.

.. contents::
   :local:
   :depth: 2


Everything here can be done from code
-------------------------------------

**There is no CLI-only behaviour.** Every command below is a thin wrapper over a
function, and the wrapper's job is argument parsing and file placement, nothing else.
That is a deliberate property, not a coincidence: a parameter sweep, a notebook or a
service driving the pipeline in a loop should not have to shell out, and anything the
command line can reach that the library cannot is a bug.

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - Command line
     - Code
   * - ``oroscope --config_path run.json``
     - ``oroscope.find_grand_regions_interactive(**params)``
   * - ``oroscope --generate_config t.json --config_preset arequipa``
     - ``oroscope.generate_config("t.json", "arequipa")``
   * - reading that file back
     - ``oroscope.load_config("t.json")``
   * - the defaults it contains
     - ``oroscope.default_config("arequipa")``
   * - ``--max_memory_gb 6``
     - ``max_memory_gb=6.0``, or ``oroscope.preflight_memory(dem, ...)`` on its own
   * - ``--no_explain``
     - ``explain=False``
   * - the summary a run prints
     - ``oroscope.explain_results(results)``, on any results dictionary
   * - ``oroscope-crop src dst --north … --east …``
     - ``oroscope.crop(src, dst, north, south, west, east)``
   * - ``oroscope-combine a b --labels A B``
     - ``oroscope.load_run(...)`` + ``oroscope.explain_combination(report, runs)``
   * - ``oroscope-sensitivity run.json --sweep min_score 0 0.35``
     - ``oroscope.run_once(config, out_dir)`` and ``oroscope.summarise(results)``
   * - ``oroscope-fetch-dem --open_topography_api_key KEY``
     - ``oroscope.fetch_dem.download_dem(region, bounds, key, out_dir)``

The pipeline **returns its results dictionary**, so a caller never has to find and
re-read the JSON it has just written:

.. code-block:: python

   import oroscope

   config = oroscope.load_config("config/tambo_colca_config.json")
   results = oroscope.find_grand_regions_interactive(
       dem_path="input/dem/colca.tif", run_output_dir="output/scan",
       **{k: v for k, v in config.items()
          if not k.startswith("_")
          and k not in ("dem_path", "print_info",
                        "output_directory_base_with_given_json")})

   print(results["results"]["total_capacity"])
   print(results["explanation"])

.. note::

   **A parameter's default is the same wherever you read it** — off the function
   signature, off ``oroscope --help``, or out of ``oroscope.default_config()``. They
   disagreed on ten parameters once, so omitting one meant different things depending
   on which door you came in by; a test now pins all three together.


``oroscope`` — run a search
---------------------------

.. code-block:: shell

   oroscope --config_path config/grand_colca_config.json

An explicitly typed option beats the configuration file, and says so when it does:

.. code-block:: shell

   oroscope --config_path config/grand_colca_config.json --min_slope_deg 5

Generate a template naming every key, then edit it:

.. code-block:: shell

   oroscope --generate_config arequipa.json --config_preset arequipa
   oroscope --config_path arequipa.json

Resume a run that died after the expensive part:

.. code-block:: shell

   oroscope --config_path run.json --resume --resume_dir output/run

Every run prints a plain-language summary of itself and saves it as
``explanation.txt``; ``--no_explain`` suppresses that.


``oroscope-crop`` — cut a window out of a DEM
---------------------------------------------

.. code-block:: shell

   oroscope-crop input/dem/arequipa_SRTMGL1.tif input/dem/colca.tif \
       --north -15.30 --south -15.85 --west -72.40 --east -71.55

The crop carries **its own** north-west corner, so it stands alone as a georeferenced
file. The window is the smallest pixel-aligned box *containing* what you asked for: the
start is floored and the stop ceiled, so you never get less than you requested.


``oroscope-combine`` — overlay two or more searches
----------------------------------------------------

.. code-block:: shell

   oroscope-combine output/grand_colca_config output/tambo_colca_config \
       --labels GRAND TAMBO --out output/combined_colca

Reports terrain viable for each, for **both** (co-location) and for **either**, with a
membership raster and an overview map — and an account of *why* the joint area is the
size it is, saved as ``combination_explanation.txt``. Co-location is usually decided by
slope: a pixel has one slope and both experiments must accept it.

The inputs must be pixel-aligned — same shape, same pixel size, same corner. That is
checked and **refused** rather than resampled, because two runs on differently-cropped
DEMs would silently compare the wrong ground.


``oroscope-sensitivity`` — how firm is the answer?
---------------------------------------------------

.. code-block:: shell

   oroscope-sensitivity config/tambo_colca_config.json \
       --sweep min_score 0.0 0.2 0.35 0.5 \
       --sweep min_target_slope_deg 0 15 25 35

Varies one parameter at a time about a baseline and tabulates what moves. Each point
runs in its own subprocess, so memory is reclaimed between them and one failed point
reports a failed row instead of ending the sweep.

Read this before quoting a capacity. Several criteria sit near cliffs.


``oroscope-fetch-dem`` — download the elevation models
--------------------------------------------------------

.. code-block:: shell

   oroscope-fetch-dem --open_topography_api_key YOUR_KEY

Fetches the bundled regions — Lima and Arequipa — into ``input/dem/`` and writes a
ready-to-run configuration for each. The key is free from
`OpenTopography <https://portal.opentopography.org/myopentopo>`_.

Run it from ``src/``: it writes to ``../input/dem/`` and ``../config/``, both relative
to the working directory. Unlike the search, this one has not been made
config-relative — there is no configuration file to be relative *to* — so it still
needs a directory one level below the repository root.


``tools/run_arequipa_full.py`` — the full-DEM run
---------------------------------------------------

Not a console script — it lives in the repository, because it is about *this* project's
full-DEM run rather than about searching in general. It runs GRAND, then TAMBO, then
the combination, over the whole Arequipa DEM, and stores the small artefacts for
:doc:`notebook 8 <notebooks>` to read.

.. code-block:: shell

   python tools/run_arequipa_full.py --dry-run

**What ``--dry-run`` does.** It reports what the real run would cost and then stops,
without starting a search, writing a file or touching the store:

.. code-block:: text

   DEM:       input/dem/arequipa_SRTMGL1.tif
   estimate:  5.08 GiB at downsample_factor 4
   available: 6.4 GiB
   would run: grand, tambo, then combine
   expected:  ~25 min for grand, ~1 min for tambo
   store:     results/arequipa_full

Each line answers a question worth answering before committing an hour:

``DEM``
   which file will be searched, and whether it is present at all. A missing DEM is
   reported here rather than after the first search has started.

``estimate`` and ``available``
   the pre-flight memory estimate against what the system reports free. This is the
   number that decides ``downsample_factor`` and ``candidate_stride``: the same DEM
   needs 7.2 GiB at 1 and 5.1 GiB at 4. Downsampling scales the labelling arrays as
   its inverse square but leaves the candidates untouched -- they are taken on the
   native grid -- and at this scale the candidates dominate, so striding is the
   stronger lever. The estimate deliberately excludes the memory-mapped DEM, which is
   file-backed and evictable.

``would run``
   which searches, honouring ``--only``. ``--only grand`` runs one and skips the
   combination.

``expected``
   the wall time to expect per search, so a run is not started ten minutes before you
   need the machine.

``store``
   where the artefacts will land — the results JSON, provenance and explanation, a few
   hundred kilobytes that notebook 8 reads.

**In a dry run no memory cap is applied**, since nothing is allocated. The real run
caps the process's address space so a search that outgrows the machine fails with
``MemoryError`` naming itself rather than letting the kernel's OOM killer choose a
victim.

Then, for real:

.. code-block:: shell

   python tools/run_arequipa_full.py

or one experiment at a time:

.. code-block:: shell

   python tools/run_arequipa_full.py --only grand

**Run it when a configuration changes, not otherwise.** The store carries a manifest
naming the configurations and the time, so a stale one is detectable rather than merely
suspected — which matters, because the whole premise of storing rather than recomputing
is that nobody looks again.


Where things are resolved from
------------------------------

Four sources, first match wins:

#. **An option you actually typed.** This beats everything, and announces itself when it
   overrides a configuration file.
#. **The configuration file** given with ``--config_path``.
#. **``config/fallbacks.json``**, if present. Every value taken from here is announced,
   because a fallback is the least visible input the tool has.
#. **The built-in default.**

Only the DEM is genuinely required: ``origin_lat``/``origin_lon`` are read from the
file's own GeoTIFF tiepoint, and a supplied origin disagreeing with the file by more
than ~100 m is reported rather than silently honoured.

.. note::

   **Paths in a configuration are relative to the configuration file**, not to the
   working directory, so a search runs the same from anywhere. ``"dem_path":
   "../input/dem/colca.tif"`` in ``config/`` means ``input/dem/colca.tif`` in the
   repository whether you are standing in the root, in ``src/`` or elsewhere; the
   outputs follow the same rule. Absolute paths are left alone, and a path that
   resolves only against the working directory still works, with a warning — the old
   behaviour, kept so this does not break a setup that relied on it.

   This replaces the long-standing requirement to ``cd src`` first.


Every option, in full
---------------------

The search accepts 82 options. Every one is also a parameter of
``oroscope.find_grand_regions_interactive`` under the same name, except the five
negative-form flags — ``--require_sky``, ``--nearest_sampling``, ``--no_geomagnetic``,
``--include_near_field``, ``--no_print_info``, ``--no_explain`` — whose positive forms
(``require_terrain``, ``bilinear_sampling``, ``use_geomagnetic``,
``exclude_near_field``, ``print_info``, ``explain``) are what the function takes.

Generated from the parser, so it cannot drift from the code: a test asserts that every
option below exists and that every option that exists appears below.

.. list-table:: Every option (82 of them)
   :header-rows: 1
   :widths: 26 14 14 46

   * - Option
     - Type
     - Default
     - What it does
   * - ``--dem_path``
     - string
     - ``—``
     - Path to the Digital Elevation Model (.tif) file.
   * - ``--origin_lat``
     - float
     - ``—``
     - Reference origin latitude (e.g., -10.228).
   * - ``--origin_lon``
     - float
     - ``—``
     - Reference origin longitude (e.g., -78.076).
   * - ``--target_antennas``
     - int
     - ``10000``
     - Total target capacity for the array (default: 10000).
   * - ``--min_width_km``
     - float
     - ``2.0``
     - Minimum acceptable width of the array site in km (default: 2.0).
   * - ``--min_altitude``
     - float
     - ``—``
     - Minimum allowable altitude in meters (optional).
   * - ``--max_altitude``
     - float
     - ``—``
     - Maximum allowable altitude in meters (optional).
   * - ``--antenna_spacing_km``
     - float
     - ``1.0``
     - Distance between antennas in km (default: 1.0).
   * - ``--min_dist_km``
     - float
     - ``10.0``
     - Minimum required distance to target mountain in km (default: 10.0).
   * - ``--max_dist_km``
     - float
     - ``80.0``
     - Maximum required distance to target mountain in km (default: 80.0).
   * - ``--grid_type``
     - square / hex
     - ``hex``
     - Antenna layout grid type (default: 'hex').
   * - ``--min_slope_deg``
     - float
     - ``3.0``
     - Minimum terrain steepness in degrees (default: 3.0).
   * - ``--max_slope_deg``
     - float
     - ``25.0``
     - Maximum terrain steepness in degrees (default: 25.0).
   * - ``--downsample_factor``
     - int
     - ``4``
     - Internal capacity mask downsampling factor for processing speed (default: 4).
   * - ``--cell_size_deg``
     - float
     - ``—``
     - Map resolution in degrees per pixel. Defaults to reading the DEM's GeoTIFF tags.
   * - ``--slope_baseline_m``
     - float
     - ``—``
     - Ground distance in metres over which slope is measured. Default: the DEM's native resolution.
   * - ``--candidate_stride``
     - int
     - ``5``
     - Keep every Nth candidate pixel before ray tracing (default: 5). Use 1 for no thinning.
   * - ``--tile_size``
     - int
     - ``2048``
     - Size of the square memory chunk for RAM management (default: 2048).
   * - ``--num_cores``
     - int
     - ``-1``
     - Number of CPU cores to use. Set to -1 to use all available cores (default: -1).
   * - ``--energy_min_pev``
     - float
     - ``—``
     - Lower tau energy in PeV. With --energy_max_pev, derives the decay-baseline distance window.
   * - ``--energy_max_pev``
     - float
     - ``—``
     - Upper tau energy in PeV.
   * - ``--n_azimuths``
     - int
     - ``9``
     - Azimuths scanned per candidate in scan mode (default: 9).
   * - ``--azimuth_half_width_deg``
     - float
     - ``60.0``
     - Half-width of the azimuth fan about the aspect. Use -1 for a full 360 sweep (default: 60).
   * - ``--elev_min_deg``
     - float
     - ``-3.0``
     - Lower edge of the accepted arrival elevation window (default: -3).
   * - ``--elev_max_deg``
     - float
     - ``3.0``
     - Upper edge of the accepted arrival elevation window (default: +3).
   * - ``--n_elev_bins``
     - int
     - ``12``
     - Elevation bins across the window (default: 12). Nearly free: cost scales with azimuths.
   * - ``--min_column_depth_gcm2``
     - float
     - ``0.0``
     - Column depth a direction must have to count, in g/cm2 (default: 0).
   * - ``--require_sky``
     - flag
     - ``off``
     - Invert the test: accept directions that reach clear sky, for cosmic-ray style channels.
   * - ``--fresnel_frequency_mhz``
     - float
     - ``—``
     - Radio band for the Fresnel clearance measurement, e.g. 50. Omitted skips the second pass.
   * - ``--antenna_height_m``
     - float
     - ``2.0``
     - Antenna height above ground, for the Fresnel measurement (default: 2).
   * - ``--include_near_field``
     - flag
     - ``off (feature on)``
     - Measure Fresnel clearance from the antenna outward instead of skipping the near field. Included for study: the result is then dominated by ground beside the antenna rather than by intervening terrain.
   * - ``--fresnel_near_field_m``
     - float
     - ``500.0``
     - Skip this much of the path when measuring Fresnel clearance (default: 500). Below ~500 m the measure is dominated by ground beside the antenna rather than by intervening terrain.
   * - ``--nearest_sampling``
     - flag
     - ``off (feature on)``
     - Sample terrain profiles at pixel centres instead of interpolating. Faster, but treats terrain as blocky, which over-estimates how much it blocks a ray.
   * - ``--muon_shielding_km``
     - float
     - ``—``
     - Rock overburden required along the arrival direction to reject atmospheric muons, in km (TAMBO quotes >4). A floor on column depth, not a band.
   * - ``--geomag_declination_deg``
     - float
     - ``—``
     - Geomagnetic declination, degrees east of north. Defaults to the Arequipa IGRF 2026 value (-6.9); supply the IGRF value for other regions.
   * - ``--geomag_inclination_deg``
     - float
     - ``—``
     - Geomagnetic inclination, degrees, positive downward. Defaults to a centered-dipole estimate at the DEM's own centre, so it follows the site automatically.
   * - ``--no_geomagnetic``
     - flag
     - ``off (feature on)``
     - Ignore the geomagnetic angle and weight all directions equally.
   * - ``--grammage_mode``
     - radio / particle
     - ``radio``
     - How atmospheric depth is scored. 'radio' is a maturity threshold, since emission comes from shower maximum and then propagates through transparent air. 'particle' is a band, since particle content dies after maximum (default: radio).
   * - ``--grammage_band_gcm2``
     - float ×2
     - ``—``
     - Atmospheric depth band scoring 1 in 'particle' mode, in g/cm2. Defaults to (X_max, 4*X_max) = (700, 2800), which suits a long path to a distant target. A short crossing gives far less: Colca supplies about 170 g/cm2, so a detector there sees a shower that is still developing and this band must be lowered or nothing scores.
   * - ``--grammage_maturity_gcm2``
     - float
     - ``—``
     - Atmospheric depth at which the 'radio' maturity ramp reaches 1, in g/cm2 (default: X_max = 700).
   * - ``--decay_energy_pev``
     - float
     - ``—``
     - Tau energy, in PeV, at which to score the probability that it decays in the gap with room left for a shower. Left out by default because the probability is strongly energy-dependent and one number cannot stand in for a spectrum. Matters most across a canyon: at 1 EeV the decay length is ~49 km against a ~3 km crossing.
   * - ``--max_range_km``
     - float
     - ``—``
     - How far to walk each profile, in km. Defaults to max_dist_km. Worth setting larger for a short-range search: column depth accumulates over the whole walk, so tying the two makes the reported depth a property of where the walk stopped rather than of the target's thickness.
   * - ``--score_percentile``
     - float
     - ``—``
     - Keep this percentage of viable candidates, ranked by score, instead of cutting at an absolute --min_score. Preferred: the default score is a product whose distribution piles up near zero, so an absolute threshold sits on a cliff, while a percentile is scale-free.
   * - ``--stop_at_target``
     - flag
     - ``off``
     - In distributed mode, stop selecting sites once target_antennas is reached. Sites are ranked by capacity, so this reports the best sites for the array actually wanted rather than every patch of qualifying ground.
   * - ``--max_memory_gb``
     - float
     - ``—``
     - Ceiling on this process's address space, in GiB. Defaults to 80%% of what the system reports available, so a search that outgrows the machine fails with MemoryError instead of inviting the OOM killer to choose a victim. 0 disables the cap.
   * - ``--decay_energy_min_pev``
     - float
     - ``—``
     - Lower end of the tau energy range for the decay term. With --decay_energy_max_pev this folds the decay probability over a power-law spectrum, which is the defensible form: the probability runs over three decades across one experiment's reach, so a single energy chooses the answer rather than approximating it.
   * - ``--decay_energy_max_pev``
     - float
     - ``—``
     - Upper end of that range, in PeV.
   * - ``--decay_spectral_index``
     - float ×1–2
     - ``—``
     - Spectral index gamma in dN/dE ~ E^-gamma for the folded decay term (default: 2.0). Give one value to pin the spectrum, or two to marginalise uniformly over that range when the index is not known -- a flat prior says so rather than pretending to a value. A softer spectrum weights low energies, where the tau decays readily, so it drives the term toward 1.
   * - ``--shower_development_m``
     - float
     - ``3000.0``
     - Path the shower needs after the tau decays, in metres (default: 3000). Used both by the decay term and as the far endpoint of the Fresnel clearance measurement.
   * - ``--gap_close_km``
     - float
     - ``—``
     - Size of the morphological closing element that fills gaps between accepted pixels, in km. Defaults to antenna_spacing_km, which couples two unrelated things. Closing more than doubles the reported area on real terrain (measured 2.29x at Colca), so this is worth setting deliberately; 0 disables it.
   * - ``--min_target_slope_deg``
     - float
     - ``—``
     - Require the terrain a ray strikes to be at least this steep, measured along the arrival azimuth. Unset by default, which asks only that rock is present -- true almost everywhere in the Andes. TAMBO's tau exits a canyon *wall*, so this is what separates a canyon from a hillside.
   * - ``--max_target_slope_deg``
     - float
     - ``—``
     - Upper bound on the struck terrain's slope along the arrival azimuth. Unset by default. Note a ceiling does not empty the result: a flat valley floor passes any ceiling, so this removes walls rather than everything.
   * - ``--grammage_band_fraction``
     - float
     - ``—``
     - When the shower band is derived from an energy range, the fraction of peak particle content that still counts as a usable shower (default: 0.1). Lower admits younger and older showers, so it widens the band and accepts narrower canyons.
   * - ``--shower_elongation_rate_gcm2``
     - float
     - ``—``
     - How much deeper shower maximum sits per decade of primary energy, in g/cm2 (default: 55, the usual hadronic value; a purely electromagnetic cascade is nearer 85).
   * - ``--shower_lambda_gcm2``
     - float
     - ``—``
     - Gaisser-Hillas interaction length setting how fast the shower profile rises and falls, in g/cm2 (default: 70).
   * - ``--solid_angle_half_sr``
     - float
     - ``—``
     - Accepted solid angle scoring 0.5, in steradians (default: 0.05). This is a GRAND-scale value: an experiment looking across a canyon sees far more sky, and leaving it at 0.05 saturates the term so it stops discriminating.
   * - ``--distance_band_m``
     - float ×2
     - ``—``
     - Exit-point distance band scoring 1, in metres. Defaults to the configured decay-baseline window.
   * - ``--clearance_full_at``
     - float
     - ``—``
     - Fresnel clearance ratio, in first-Fresnel radii, that scores 1 (default: 1.0).
   * - ``--score_weights``
     - string
     - ``—``
     - Per-component weights for --score_composition weighted, as name=value pairs, e.g. 'shower=2,solid_angle=1,depth=0.5'. Components not named default to weight 1.
   * - ``--nu_interaction_length_gcm2``
     - float
     - ``—``
     - Neutrino interaction length for the Earth-chord attenuation term, g/cm2 (order 1e8 near an EeV). Omitted reports the chord without weighting by it.
   * - ``--refraction_k``
     - float
     - ``—``
     - Refraction k-factor for the RADIO path only (default: 4/3). Particle trajectories always use the true Earth radius, since neutrinos and taus are not refracted.
   * - ``--depth_band_gcm2``
     - float ×2
     - ``—``
     - Column depth band scoring 1, in g/cm2. The tau must be produced and must escape, so this is a band, not a floor.
   * - ``--score_composition``
     - product / mean / min
     - ``product``
     - How component scores combine (default: product).
   * - ``--min_score``
     - float
     - ``0.0``
     - Discard candidates scoring below this (default: 0, keep all).
   * - ``--rfi_zones``
     - string
     - ``none``
     - Can be preset ('lima', 'arequipa') or a valid JSON string outlining custom exclusion zones.
   * - ``--road_map_path``
     - string
     - ``—``
     - Path to a raster mapping distance-to-roads (optional).
   * - ``--max_road_dist_km``
     - float
     - ``20.0``
     - Maximum distance allowed from a road in km (default: 20.0).
   * - ``--search_mode``
     - single / distributed
     - ``distributed``
     - 'single' finds one monolithic site, 'distributed' allows sub-arrays.
   * - ``--min_sub_array_size``
     - int
     - ``500``
     - Minimum required capacity for a sub-array to be considered valid (default: 500).
   * - ``--min_aspect_deg``
     - float
     - ``—``
     - Minimum bound for site facing direction in degrees (0-360).
   * - ``--max_aspect_deg``
     - float
     - ``—``
     - Maximum bound for site facing direction in degrees (0-360).
   * - ``--region_name``
     - string
     - ``—``
     - Cosmetic region name to print on the final visualization chart.
   * - ``--generate_kml``
     - flag
     - ``off``
     - Include this flag to generate a Google Earth KML file of the findings.
   * - ``--no_print_info``
     - flag
     - ``off (feature on)``
     - Include this flag to skip printing the detailed explanatory text.
   * - ``--no_explain``
     - flag
     - ``off (feature on)``
     - Skip the plain-language summary of the run. It is printed by default, and saved as explanation.txt beside the results: what was found, which constraint set the size of the answer, what held the surviving sites back, and which of the numbers are assumptions rather than measurements. A results file can be re-explained at any time with explain.explain_results().
   * - ``--config_path``
     - string
     - ``—``
     - Path to external JSON configuration file.
   * - ``--output_directory_base_with_given_json``
     - string
     - ``../output/``
     - Base directory for outputs when a JSON config is supplied (default: ../output/).
   * - ``--output_image_format``
     - string
     - ``png``
     - Format of the saved map visual, e.g., png, pdf, svg (default: png).
   * - ``--resume``
     - flag
     - ``off``
     - Include this flag to resume a previous run from the ray-tracing checkpoint.
   * - ``--resume_dir``
     - string
     - ``—``
     - Path to an output folder from a previously failed run to resume from the ray-tracing checkpoint.
   * - ``--generate_config``
     - string
     - ``—``
     - Supply a filepath to generate a default JSON config template and exit.
   * - ``--config_preset``
     - default / lima / arequipa
     - ``default``
     - Optional presets to inject when using --generate_config.
