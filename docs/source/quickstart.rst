Quickstart
==========

.. jupyter-execute::
   :hide-code:
   :hide-output:

   # jupyter-sphinx runs these blocks in a plain kernel, which renders a Figure as its
   # text repr rather than as a PNG unless the inline backend is switched on. Without
   # this every diagram below silently became "<Figure size 1020x415 with 2 Axes>".
   %matplotlib inline


A search, and then what each knob does. **One import**: everything below is on
``oroscope``, and the submodules stay available when a narrower namespace reads better.

.. code-block:: shell

   pip install oroscope

Cut a region
------------

.. code-block:: python

   import oroscope

   info = oroscope.crop("input/dem/arequipa_SRTMGL1.tif", "input/dem/colca.tif",
                        north=-15.30, south=-15.85, west=-72.40, east=-71.55)
   print(info["origin_lat"], info["origin_lon"], info["rows"], info["cols"])

The crop carries its own north-west corner, so it stands alone as a georeferenced file
and the search reads its geometry back out of it.

Run a search
------------

.. code-block:: python

   config = oroscope.load_config("config/grand_colca_config.json")

   results = oroscope.find_grand_regions_interactive(
       dem_path="input/dem/colca.tif", run_output_dir="output/colca",
       **{k: v for k, v in config.items()
          if not k.startswith("_")
          and k not in ("dem_path", "print_info",
                        "output_directory_base_with_given_json")})

It **returns its results**, so nothing has to find and re-read the file it just wrote:

.. code-block:: python

   print(results["results"]["total_sites"])
   print(results["results"]["total_capacity"])

   for site in oroscope.selected_sites(results)[0]:
       print(site["site_id"], site["area_km2"], site["center_lat"], site["center_lon"])

Starting from ``oroscope.default_config()`` and overriding is the clearer habit — it
puts every knob in front of you — but it is no longer a correctness matter: the
function signature, ``oroscope --help`` and the template state the same default for
every parameter, and a test keeps all three in step.

The same search from a shell is ``oroscope --config_path config/grand_colca_config.json``
— see :doc:`cli` for every option, and for the other four console scripts.

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
  area, because morphological closing inflated it 2.35× at Colca.

Any results file can be re-explained later, with no DEM and nothing re-run:

.. code-block:: python

   import json, explain

   with open("output/tambo_colca_config/oroscope_results_colca.json") as f:
       results = json.load(f)
   print(explain.explain_results(results))

Combine two experiments
-----------------------

.. code-block:: python

   grand = oroscope.load_run("output/grand_colca_config")
   tambo = oroscope.load_run("output/tambo_colca_config")
   oroscope.check_alignment([grand, tambo])   # refuses to overlay the wrong ground

   print(oroscope.explain_combination(report, {"GRAND": grand_results,
                                               "TAMBO": tambo_results}))

reports terrain viable for each, for both (co-location) and for either — and says
*which screening band decides that*, which is usually slope: a pixel has one slope and
both experiments must accept it. ``oroscope-combine`` does the same and writes a
membership raster and an overview map.

The physics on its own
----------------------

The physics layer takes arrays and returns values, with no side effects and no terrain:

.. jupyter-execute::

   from oroscope import physics

   print(f"tau decay length at 100 PeV: {physics.tau_decay_length_m(100.0):,.0f} m")
   print(f"X_max at 100 PeV:            {float(physics.shower_maximum_gcm2(100.0)):.0f} g/cm^2")

   lo, hi = physics.grammage_band_from_energy(3.0, 1000.0)
   print(f"particle shower band:        {lo:.0f} - {hi:.0f} g/cm^2")

The scan itself takes a terrain array and a list of candidate pixels:

.. jupyter-execute::

   import numpy as np
   from oroscope import arrival_scan, site_searcher as ss

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

Before a long run
-----------------

``oroscope.preflight_memory`` gives the estimate and the address-space cap on their
own — worth calling before a loop, since it was a sweep that once reached 6.9 GB and
was killed by the kernel:

.. code-block:: python

   report = oroscope.preflight_memory("input/dem/arequipa_SRTMGL1.tif",
                                      downsample_factor=4, candidate_stride=5)
   print(report["estimate_gb"], report["available_gb"])

Passing ``max_memory_gb`` to a search caps its address space, so one that outgrows the
machine fails with ``MemoryError`` naming itself rather than inviting the kernel's OOM
killer to choose a victim.

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
