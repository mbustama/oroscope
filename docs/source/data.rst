Getting the data
================

The tool needs one input and takes two optional ones. This page says where each comes
from, at what resolution, and what the limits are.

.. contents::
   :local:
   :depth: 2


Elevation models
----------------

**A digital elevation model is the only required input.** Everything else is context.

Where to get one
^^^^^^^^^^^^^^^^

`OpenTopography <https://portal.opentopography.org/>`_ serves the global datasets through
an API, and ``oroscope-fetch-dem`` wraps it:

.. code-block:: shell

   export OPENTOPOGRAPHY_API_KEY=...          # see below
   cd src
   oroscope-fetch-dem --region arequipa

**Getting a key.** Free, and it takes a minute.

1. Register at `portal.opentopography.org/myopentopo
   <https://portal.opentopography.org/myopentopo>`_ and sign in.
2. Open **myOpenTopo Authorizations and API Key** from the account menu.
3. Copy the key.

Pass it as ``--open_topography_api_key``, or set ``OPENTOPOGRAPHY_API_KEY`` in the
environment — which keeps it out of your shell history and out of any file that might be
committed by accident.

Which resolution
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 16 12 12 60

   * - Dataset
     - Spacing
     - Coverage
     -
   * - ``SRTMGL1``
     - 1 arc-sec, ~30 m
     - 60° N – 56° S
     - **The default choice.** Void-filled, so glaciated peaks are present rather than
       holes. Every department run in this project uses it.
   * - ``SRTMGL3``
     - 3 arc-sec, ~90 m
     - 60° N – 56° S
     - For areas too large for 30 m. A survey instrument: it resolves plateaux and
       ranges, **not canyon walls**.
   * - ``COP30`` / ``COP90``
     - 1 / 3 arc-sec
     - global
     - Copernicus. Newer than SRTM and often better in steep terrain; not used here only
       because switching would break comparability with existing runs.
   * - ``AW3D30``
     - 1 arc-sec
     - global
     - ALOS. Was used for Lima and was replaced by SRTMGL1, because a dataset difference
       between regions is indistinguishable from a difference in the ground.

.. warning::

   **Requests are capped by area, per dataset:** 450,000 km² for every 30 m dataset and
   4,050,000 km² for the 90 m ones. Peru's bounding box is about 2.86 million km² — inside
   the 90 m limit, six times over the 30 m one. That, and memory, is why the national
   survey is 3 arc-seconds.

**Choose the resolution for the feature you are looking for, not for the area.** TAMBO
selects on canyon walls: Colca's floor is ~1 km wide, so at 90 m a canyon is eleven pixels
across and the wall the array stands on is a handful. A 90 m grid has averaged away the
thing being screened for. GRAND's plateaux survive it.

Bundled regions
^^^^^^^^^^^^^^^

``--region`` fetches one of four boxes, defined in
:data:`oroscope.fetch_dem.REGIONS`. The three departments share a dataset deliberately,
so that runs over them are comparable.

.. list-table::
   :header-rows: 1
   :widths: 14 14 12 60

   * - Region
     - Dataset
     - Size
     -
   * - ``arequipa``
     - SRTMGL1
     - 129 Mpx
     - The high plateau. Most published numbers here come from it or a crop of it.
   * - ``ancash``
     - SRTMGL1
     - 69 Mpx
     - The Cordillera Blanca. Bounds from OpenStreetMap's administrative boundary.
   * - ``lima``
     - SRTMGL1
     - 105 Mpx
     - Coastal desert rising to the Andean flank.
   * - ``peru``
     - SRTMGL3
     - 339 Mpx
     - The whole country, at the only resolution that fits.

.. note::

   **Regions are bounding boxes and departments are not rectangles.** Boxes overlap —
   ``ancash`` and ``lima`` share 9,198 km² — and ground found in one may lie
   administratively in another. The largest joint patch of the Ancash run reverse-geocodes
   to Cajatambo, **Lima**. File results by box; read them by geography.

Anywhere else
^^^^^^^^^^^^^

For a region with no entry, download the tiles from the OpenTopography portal, merge them
into one GeoTIFF if the area spans several, and cut the window you want:

.. code-block:: shell

   oroscope-crop big.tif window.tif --north -8.80 --south -9.90 --west -78.00 --east -77.20

The crop carries **its own** north-west corner, so it stands alone as a georeferenced
file. The window is the smallest pixel-aligned box *containing* what you asked for.

Cropping is not only for convenience. A crop small enough to run at
``downsample_factor 1`` and ``candidate_stride 1`` is free of both sampling biases, and is
the only way to get a TAMBO area that is not a lower bound.

What the tool needs from the file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A GeoTIFF carrying ``ModelPixelScaleTag`` and ``ModelTiepointTag``. Both are read
automatically, so ``origin_lat``, ``origin_lon`` and ``cell_size_deg`` can all be left
``null`` in a configuration — which is the recommended use. A supplied origin that
disagrees with the tiepoint is reported rather than silently preferred.


Roads and settlements
---------------------

Optional context for the maps, from `OpenStreetMap <https://www.openstreetmap.org>`_ via
the Overpass API. **No search number moves because of them** — they are drawn, not
applied.

.. code-block:: shell

   cd src
   python -m oroscope.fetch_roads --dem ../input/dem/ancash_SRTMGL1.tif --places

or, if the console script is installed, ``oroscope-fetch-roads`` with the same arguments.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Option
     - What it does
   * - ``--dem PATH``
     - Take the bounding box from a DEM, which is the usual way.
   * - ``--bbox S N W E``
     - Give the box explicitly instead.
   * - ``--places``
     - Also fetch populated places — cities, towns and villages, with population where
       OSM has it.
   * - ``--places_only``
     - Skip the roads.
   * - ``--classes``
     - Which road classes to keep. Default: motorway, trunk, primary, secondary,
       tertiary.
   * - ``--step_deg``
     - Tile size for the Overpass queries. A large box is split; lower it if requests
       time out.
   * - ``--out PATH``
     - Where to write. Places go to ``<stem>_places.geojson`` beside it.

Then pass them to a run or to a combination:

.. code-block:: shell

   oroscope --config_path run.json \
       --roads_geojson input/roads/ancash_SRTMGL1.geojson \
       --settlements input/roads/ancash_SRTMGL1_places.geojson

.. note::

   Overpass is a shared public service and rate-limits. A department-sized box takes
   several minutes and occasionally stalls; the places query is the slower half. **Fetch
   the context before starting a search**, because a run resolves its map inputs once at
   the beginning and will draw without them if they arrive late.

Roads can also be a *criterion* rather than context, through ``road_map_path`` and
``max_road_dist_km`` — but that path wants an aligned distance-to-road raster rather than
a GeoJSON, and none of the runs here uses it.


Radio-frequency interference zones
-----------------------------------

``rfi_zones`` excludes ground near settlements and industry. It accepts a preset name
(``"arequipa"``, ``"lima"``), ``"none"``, a JSON string, or an explicit list of
``('circle', lat, lon, radius_km, name)`` and ``('poly', [(lat, lon), ...], name)``
entries.

The presets are **hand-curated lists of specific places**, not a transferable rule. There
is no preset for a region that has not had one written, and applying another region's list
excludes nothing while looking as though it did. Where a run has no list of its own, these
docs say so and quantify the difference: Arequipa's five circles cover ~3,500 km² of a
~120,000 km² box, about 2.9%.

A future improvement is to build zones from the OSM places above, which are already
downloaded and carry population. That is not implemented.
