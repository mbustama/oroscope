"""
Cutting a window out of a DEM.

A crop that takes the wrong pixels, or writes the right pixels with the wrong corner,
does not fail — it produces a perfectly ordinary GeoTIFF describing the wrong ground,
and every search run on it is quietly about somewhere else. That is the same failure
mode the origin check exists for, one step earlier in the chain.

So these assert against arithmetic: a synthetic DEM whose elevation encodes its own
pixel coordinates, so the identity of every cropped pixel is checkable.
"""

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
import tifffile as tiff

from _support import quiet  # noqa: F401  (also sets up sys.path)

from oroscope import crop_dem

# A DEM whose north-west corner is (-15.0, -72.0), one arc-second pixels, 400 x 600.
LAT0, LON0 = -15.0, -72.0
CELL = 1 / 3600
ROWS, COLS = 400, 600


def write_dem(path, rows=ROWS, cols=COLS, cell=CELL, lat0=LAT0, lon0=LON0):
    """Elevation encodes the pixel's own indices: z = row * 1000 + col."""
    r = np.arange(rows)[:, None].repeat(cols, 1)
    c = np.arange(cols)[None, :].repeat(rows, 0)
    z = (r * 1000 + c).astype(np.float32)
    tiff.imwrite(path, z, extratags=[
        (33550, "d", 3, (cell, cell, 0.0)),
        (33922, "d", 6, (0.0, 0.0, 0.0, lon0, lat0, 0.0)),
    ])
    return path


class CropCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="oroscope_crop_")
        self.src = write_dem(os.path.join(self.tmp, "src.tif"))
        self.dst = os.path.join(self.tmp, "out.tif")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


class TestReadGeo(CropCase):
    def test_it_reads_the_pixel_size_and_the_corner(self):
        cell_x, cell_y, lon0, lat0, rows, cols = crop_dem.read_geo(self.src)
        self.assertAlmostEqual(cell_x, CELL, places=12)
        self.assertAlmostEqual(cell_y, CELL, places=12)
        self.assertAlmostEqual(lon0, LON0, places=9)
        self.assertAlmostEqual(lat0, LAT0, places=9)
        self.assertEqual((rows, cols), (ROWS, COLS))


class TestTheWindowTakesTheRightPixels(CropCase):
    """
    Rows run north to south and columns west to east, so a window 100 rows below the
    top and 50 columns in starts at elevation 100*1000 + 50.

    The window is deliberately requested at half-pixel offsets. ``crop`` floors the
    start and ceils the stop, so the result is the smallest pixel-aligned box that
    *contains* the request — asking on an exact pixel boundary leaves it to floating
    point whether the neighbouring pixel comes too, which is a poor thing to assert on.
    """

    # north/south/west/east at .5 of a pixel: floor and ceil are then unambiguous.
    NORTH, SOUTH = LAT0 - 100.5 * CELL, LAT0 - 150.5 * CELL
    WEST, EAST = LON0 + 50.5 * CELL, LON0 + 90.5 * CELL

    def crop_it(self):
        return crop_dem.crop(self.src, self.dst,
                             self.NORTH, self.SOUTH, self.WEST, self.EAST)

    def test_the_first_pixel_is_the_one_the_window_asked_for(self):
        info = self.crop_it()
        window = tiff.imread(self.dst)
        self.assertEqual(float(window[0, 0]), 100 * 1000 + 50)
        self.assertEqual(info["rows"], 51)          # rows 100..150 inclusive
        self.assertEqual(info["cols"], 41)          # cols 50..90 inclusive

    def test_the_crop_contains_the_requested_window(self):
        """The contract of flooring the start and ceiling the stop."""
        info = self.crop_it()
        self.assertGreaterEqual(info["origin_lat"], self.NORTH,
                                "the crop's top edge must be at or north of the request")
        self.assertLessEqual(info["south"], self.SOUTH)
        self.assertLessEqual(info["origin_lon"], self.WEST)
        self.assertGreaterEqual(info["east"], self.EAST)

    def test_it_is_no_more_than_one_pixel_larger_than_asked_on_each_side(self):
        info = self.crop_it()
        self.assertLess(info["origin_lat"] - self.NORTH, CELL)
        self.assertLess(self.SOUTH - info["south"], CELL)

    def test_the_crop_carries_its_own_corner_not_the_source_corner(self):
        """Standing alone as a georeferenced file is the whole point."""
        info = self.crop_it()
        self.assertAlmostEqual(info["origin_lat"], LAT0 - 100 * CELL, places=9)
        self.assertAlmostEqual(info["origin_lon"], LON0 + 50 * CELL, places=9)

        _, cell_y, lon0, lat0, _, _ = crop_dem.read_geo(self.dst)
        self.assertAlmostEqual(lat0, info["origin_lat"], places=9)
        self.assertAlmostEqual(lon0, info["origin_lon"], places=9)
        self.assertAlmostEqual(cell_y, CELL, places=12)

    def test_the_reported_south_and_east_edges_match_the_pixels_written(self):
        info = self.crop_it()
        self.assertAlmostEqual(info["south"],
                               info["origin_lat"] - info["rows"] * CELL, places=9)
        self.assertAlmostEqual(info["east"],
                               info["origin_lon"] + info["cols"] * CELL, places=9)

    def test_the_elevation_range_is_the_window_not_the_whole_dem(self):
        info = self.crop_it()
        self.assertEqual(info["z_min"], 100 * 1000 + 50)
        self.assertEqual(info["z_max"], 150 * 1000 + 90)

    def test_a_window_larger_than_the_dem_is_clipped_to_it(self):
        info = crop_dem.crop(self.src, self.dst, LAT0 + 1.0, LAT0 - 1.0,
                             LON0 - 1.0, LON0 + 1.0)
        self.assertEqual((info["rows"], info["cols"]), (ROWS, COLS))
        self.assertAlmostEqual(info["origin_lat"], LAT0, places=9)
        self.assertAlmostEqual(info["origin_lon"], LON0, places=9)

    def test_a_window_outside_the_dem_is_refused_rather_than_returning_nothing(self):
        with self.assertRaises(SystemExit) as caught:
            crop_dem.crop(self.src, self.dst, LAT0 + 2.0, LAT0 + 1.0,
                          LON0 + 1.0, LON0 + 2.0)
        self.assertIn("does not overlap", str(caught.exception))

    def test_the_crop_is_readable_by_the_searcher(self):
        """The point of writing the tags: the next stage reads the corner back."""
        from oroscope import site_searcher as ss
        info = self.crop_it()
        lat, lon = ss.read_dem_origin(self.dst)
        self.assertAlmostEqual(lat, info["origin_lat"], places=9)
        self.assertAlmostEqual(lon, info["origin_lon"], places=9)


class TestCropCli(CropCase):
    def run_main(self, *args):
        argv = sys.argv
        sys.argv = ["crop_dem.py", *args]
        try:
            with quiet():
                crop_dem.main()
        finally:
            sys.argv = argv

    def test_it_crops_from_the_command_line(self):
        self.run_main(self.src, self.dst,
                      "--north", repr(LAT0 - 10.5 * CELL),
                      "--south", repr(LAT0 - 60.5 * CELL),
                      "--west", repr(LON0 + 5.5 * CELL),
                      "--east", repr(LON0 + 45.5 * CELL))
        self.assertTrue(os.path.exists(self.dst))
        # rows 10..60 and cols 5..45 inclusive, by the floor/ceil contract
        self.assertEqual(tiff.imread(self.dst).shape, (51, 41))


if __name__ == "__main__":
    unittest.main()
