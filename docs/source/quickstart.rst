Quickstart
==========

.. jupyter-execute::
   :hide-code:
   :hide-output:

   # jupyter-sphinx runs these blocks in a plain kernel, which renders a Figure as its
   # text repr rather than as a PNG unless the inline backend is switched on. Without
   # this every diagram below silently became "<Figure size 1020x415 with 2 Axes>".
   %matplotlib inline


A search in three commands, then what each knob does.

Cut a region
------------

.. code-block:: shell

   oroscope-crop input/dem/arequipa_SRTMGL1.tif input/dem/colca.tif \
       --north -15.30 --south -15.85 --west -72.40 --east -71.55

prints the corner coordinates to paste into a configuration file.

Run a search
------------

.. code-block:: shell

   oroscope --config_path config/grand_colca_config.json

Every parameter can also be given on the command line, and **an explicitly typed
option beats the configuration file**, announcing itself when it does:

.. code-block:: shell

   oroscope --config_path config/grand_colca_config.json --min_slope_deg 5

Read what it found
------------------

Every run ends with a plain-language summary of itself, and saves it as
``explanation.txt`` beside the results. It is on by default; ``--no_explain``
suppresses it.

It says four things the results file contains but does not spell out:

* **What was found** — sites, area, capacity against the target.
* **Which constraint bound.** The funnel stage that removed the largest share of what
  reached it, and the parameter behind it. When a search returns little or nothing,
  this is the answer.
* **What held the sites back.** The score is a product of named components, so the
  lowest one is reported per site: under a product it bounds the total from above.
* **Which numbers are assumptions**, with the measured sensitivity of each. Including
  the one every reader gets wrong unaided: reported area is not physics-accepted
  area, because morphological closing inflated it 2.29× at Colca.

Any results file can be re-explained later, with no DEM and nothing re-run:

.. code-block:: python

   import json, explain

   with open("output/tambo_colca_config/oroscope_results_colca.json") as f:
       results = json.load(f)
   print(explain.explain_results(results))

Combine two experiments
-----------------------

.. code-block:: shell

   oroscope --config_path config/tambo_colca_config.json
   oroscope-combine output/grand_colca_config output/tambo_colca_config \
       --labels GRAND TAMBO --out output/combined

reports terrain viable for each, for both (co-location) and for either, with a
membership raster and an overview map.

Using it as a library
---------------------

The physics layer takes arrays and returns values, with no side effects, so it is
usable on its own:

.. jupyter-execute::

   import physics

   print(f"tau decay length at 100 PeV: {physics.tau_decay_length_m(100.0):,.0f} m")
   print(f"X_max at 100 PeV:            {float(physics.shower_maximum_gcm2(100.0)):.0f} g/cm^2")

   lo, hi = physics.grammage_band_from_energy(3.0, 1000.0)
   print(f"particle shower band:        {lo:.0f} - {hi:.0f} g/cm^2")

The scan itself takes a terrain array and a list of candidate pixels:

.. jupyter-execute::

   import numpy as np, arrival_scan, site_searcher as ss

   # A ridge to the east of a flat plain, so a westward-facing candidate sees it
   n = 300
   cell_y, cell_x = 30.7, 29.8
   cols = np.arange(n)[None, :].repeat(n, 0)
   z = (2000.0 + 1200.0 * np.exp(-((cols - 220) / 18.0) ** 2)).astype(np.float32)

   grid = ss.resolve_grid_geometry("none.tif", -15.6, cell_size_deg=1 / 3600)
   candidates = np.array([[150.0, 60.0, 90.0]])        # row, col, aspect east

   out = arrival_scan.scan(candidates, z, grid, n_azimuths=1, half_width_deg=0.0,
                           elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
                           min_dist_km=1.0, max_dist_km=10.0, max_range_m=10000.0)

   print(f"accepted directions: {int(out['cells'][0])}")
   print(f"mean exit distance:  {out['mean_distance_m'][0]:,.0f} m")
   print(f"horizon:             {out['horizon_deg'][0]:.2f} deg")

And so is the pipeline. Anything the command line does, the library does: it reads
and writes configuration files, estimates and caps its own memory, explains itself,
and hands back its results rather than making the caller find the file it just wrote.

.. code-block:: python

   import site_searcher as ss

   config = ss.load_config("config/tambo_colca_config.json")
   config["min_score"] = 0.2

   results = ss.find_grand_regions_interactive(
       dem_path="input/dem/colca.tif", run_output_dir="output/scan",
       max_memory_gb=6.0, explain=False,
       **{k: v for k, v in config.items()
          if k not in ("dem_path", "print_info",
                       "output_directory_base_with_given_json")})

   print(results["results"]["total_capacity"])
   print(results["explanation"])

``ss.generate_config(path, preset)`` writes a template naming every knob, and
``ss.preflight_memory(dem_path, downsample_factor=...)`` gives the estimate and the
address-space cap on their own — worth calling before a loop, since it is a sweep
that once reached 6.9 GB and was killed by the kernel.

Choosing parameters
-------------------

The criteria are per experiment. The two bundled configurations differ in these ways,
and the differences are the whole of what makes one GRAND and the other TAMBO:

.. list-table::
   :header-rows: 1
   :widths: 30 22 22 26

   * - Parameter
     - GRAND
     - TAMBO
     - Why
   * - ``min_slope_deg`` / ``max_slope_deg``
     - 3 / 25
     - 20 / 60
     - Deployable ground against a canyon wall
   * - ``min_target_slope_deg``
     - unset
     - 25
     - TAMBO needs the *far* wall to be a wall
   * - ``min_dist_km`` / ``max_dist_km``
     - 10 / 40
     - 2 / 5
     - To the horizon, or across a canyon
   * - ``elev_min_deg`` / ``elev_max_deg``
     - −3 / +3
     - −20 / +20
     - Earth-skimming, against a wall that subtends tens of degrees
   * - ``antenna_spacing_km``
     - 1.0
     - 0.1
     - Radio antennas against particle detectors
   * - ``min_width_km``
     - 2.0
     - 0.0
     - A compact array against a strip along a wall
   * - ``fresnel_frequency_mhz``
     - 50
     - null
     - Radio propagation, or none
   * - ``use_geomagnetic``
     - true
     - false
     - Radio emission goes as :math:`|v \times B|`; particles do not care
   * - ``grammage_mode``
     - ``radio``
     - ``particle``
     - A maturity threshold against a band

When a search returns nothing
-----------------------------

Read the run's own summary, which names the stage where the count collapsed and the
parameter behind it. It is printed at the end of every run and saved as
``explanation.txt``. Behind it is the funnel, in the results JSON: the survivor count
after each filter, where the constraint responsible is the line where the count
collapses.

Then check how firm the answer is:

.. code-block:: shell

   oroscope-sensitivity config/tambo_colca_config.json \
       --sweep min_score 0.0 0.2 0.35 0.5 \
       --sweep decay_energy_pev 3 55 1000

Several criteria sit near cliffs, and a result that moves by an order of magnitude
across a plausible range of an assumption is a result about the assumption. See
:doc:`assumptions`.
