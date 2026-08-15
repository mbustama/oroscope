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

Read the funnel, which every run prints and stores in its results JSON. It gives the
survivor count after each filter, so the constraint responsible is the line where the
count collapses.

Then check how firm the answer is:

.. code-block:: shell

   oroscope-sensitivity config/tambo_colca_config.json \
       --sweep min_score 0.0 0.2 0.35 0.5 \
       --sweep decay_energy_pev 3 55 1000

Several criteria sit near cliffs, and a result that moves by an order of magnitude
across a plausible range of an assumption is a result about the assumption. See
:doc:`assumptions`.
