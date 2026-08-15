The command line
================

Five console scripts, installed by ``pip install oroscope``. This page is task-shaped —
*how do I do X* — rather than exhaustive; the complete list of the search's 83 options,
grouped by what they do, is in the `README
<https://github.com/mbustama/oroscope#4-parameter-configuration-hierarchy>`_.

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
     - ``site_searcher.find_grand_regions_interactive(**params)``
   * - ``oroscope --generate_config t.json --config_preset arequipa``
     - ``site_searcher.generate_config("t.json", "arequipa")``
   * - reading that file back
     - ``site_searcher.load_config("t.json")``
   * - the defaults it contains
     - ``site_searcher.default_config("arequipa")``
   * - ``--max_memory_gb 6``
     - ``max_memory_gb=6.0``, or ``site_searcher.preflight_memory(dem, ...)`` on its own
   * - ``--no_explain``
     - ``explain=False``
   * - the summary a run prints
     - ``explain.explain_results(results)``, on any results dictionary
   * - ``oroscope-crop src dst --north … --east …``
     - ``crop_dem.crop(src, dst, north, south, west, east)``
   * - ``oroscope-combine a b --labels A B``
     - ``combine_experiments.load_run(...)`` + ``explain.explain_combination(report, runs)``
   * - ``oroscope-sensitivity run.json --sweep min_score 0 0.35``
     - ``sensitivity.run_once(config, out_dir)`` and ``sensitivity.summarise(results)``
   * - ``oroscope-fetch-dem --open_topography_api_key KEY``
     - ``fetch_dem.download_dem(region, bounds, key, out_dir)``

The pipeline **returns its results dictionary**, so a caller never has to find and
re-read the JSON it has just written:

.. code-block:: python

   import site_searcher as ss

   config = ss.load_config("config/tambo_colca_config.json")
   results = ss.find_grand_regions_interactive(
       dem_path="input/dem/colca.tif", run_output_dir="output/scan",
       **{k: v for k, v in config.items()
          if not k.startswith("_")
          and k not in ("dem_path", "print_info",
                        "output_directory_base_with_given_json")})

   print(results["results"]["total_capacity"])
   print(results["explanation"])

.. warning::

   The function's defaults are **not** the configuration template's defaults. Five
   parameters differ — ``search_mode``, ``grid_type``, ``target_antennas``,
   ``min_dist_km`` and ``min_sub_array_size`` — so omitting one means different things
   depending on which entry point you used. Start from ``ss.default_config()`` and
   override, rather than relying on the signature's defaults.


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

Run it from ``src/``: both paths are relative to it.


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
   estimate:  2.32 GiB at downsample_factor 4
   available: 5.4 GiB
   would run: grand, tambo, then combine
   expected:  ~25-30 minutes each
   store:     results/arequipa_full

Each line answers a question worth answering before committing an hour:

``DEM``
   which file will be searched, and whether it is present at all. A missing DEM is
   reported here rather than after the first search has started.

``estimate`` and ``available``
   the pre-flight memory estimate against what the system reports free. This is the
   number that decides ``downsample_factor``: the same DEM needs 4.5 GiB at 1 and
   2.3 GiB at 4, and the labelling arrays scale as its inverse square. The estimate
   deliberately excludes the memory-mapped DEM, which is file-backed and evictable.

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

   The search resolves several paths relative to the working directory, so the bundled
   configurations expect to be run from ``src/``. That is a known wart, not a design:
   the console scripts are the first half of removing it, and making those paths
   relative to the configuration file is the second.
