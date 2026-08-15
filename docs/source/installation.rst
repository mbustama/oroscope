Installation
============

Requirements
------------

Python 3.9 or newer. Every dependency is installed automatically; none is optional in
practice, which is why they are all base dependencies rather than extras.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Package
     - What needs it
   * - ``numpy``
     - Everything.
   * - ``numba``
     - The arrival scan is a compiled parallel kernel. Without it the profile walk
       falls back to interpreted Python and a real DEM becomes unusable.
   * - ``scipy``
     - Morphology and connected-component labelling.
   * - ``tifffile``
     - Reads and writes the GeoTIFF DEMs and mask rasters.
   * - ``imagecodecs``
     - Not imported by name anywhere, and required anyway: the DEMs are LZW-compressed
       and ``tifffile`` delegates that codec to it. Omitting it produces a decode
       failure that is hard to trace back to a missing package.
   * - ``matplotlib``
     - The overview maps; the KML contours are extracted from its paths.
   * - ``tqdm``
     - Progress bars over the tiled passes.

From PyPI
---------

.. code-block:: shell

   pip install oroscope

From a clone
------------

Preferred if you want the notebooks, the example configurations, the benchmark harness
or the test suite, none of which ship in the wheel:

.. code-block:: shell

   git clone https://github.com/mbustama/oroscope.git
   cd oroscope
   pip install -e .

Optional extras
---------------

.. code-block:: shell

   pip install -e '.[test]'        # pytest and coverage
   pip install -e '.[docs]'        # everything `cd docs && make html` needs
   pip install -e '.[notebooks]'   # enough to execute notebooks/ unaided
   pip install -e '.[diagnostics]' # psutil, for RAM reporting during a run

Getting a DEM
-------------

A search needs a digital elevation model in GeoTIFF form, at roughly 30 m resolution.
`OpenTopography <https://opentopography.org/>`_ serves SRTM and ALOS AW3D30, both of
which work.

.. code-block:: shell

   oroscope-fetch-dem --open_topography_api_key YOUR_KEY

downloads the tiles for the bundled regions and writes matching configuration files.
For anywhere else, download the tiles yourself, merge them if the region spans several,
and cut the window you want:

.. code-block:: shell

   oroscope-crop big.tif colca.tif --north -15.30 --south -15.85 \
       --west -72.40 --east -71.55

which prints the ``origin_lat`` and ``origin_lon`` to paste into a configuration.

.. note::

   Two searches can only be combined if they cover **identical** ground at the same
   ``downsample_factor``. ``oroscope-combine`` checks the world files and refuses
   rather than resampling, because two runs on differently-cropped DEMs would
   silently compare the wrong terrain.

Verifying the install
---------------------

.. code-block:: shell

   python -c "import arrival_scan, physics, scoring; print('ok')"
   cd tests && python -m unittest discover

The suite is standard-library ``unittest`` and needs nothing beyond the package. Tests
that need a real DEM skip automatically when ``input/dem/`` is absent.
