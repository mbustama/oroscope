"""
Map geometry: resolution detection and the two metric pixel sizes.

A geographic DEM is square in degrees but not in metres, and every physical
conversion in the pipeline depends on getting that right.
"""

import math
import os
import shutil
import tempfile
import unittest

import numpy as np

from _support import ss
import synthetic


class TestResolutionDetection(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="sitesearch_geom_")
        self.z = np.zeros((64, 64), dtype=np.float32)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, cell_deg):
        path = os.path.join(self.tmp, name)
        return synthetic.write_geotiff(path, self.z, -16.0, -72.0, cell_size_deg=cell_deg)

    def test_reads_pixel_size_from_geotiff_tags(self):
        path = self._write("arcsec.tif", synthetic.CELL_DEG)
        grid = ss.resolve_grid_geometry(path, -16.0)
        self.assertAlmostEqual(grid.cell_size_deg, synthetic.CELL_DEG, places=12)
        self.assertIn("auto-detected", grid.source)

    def test_reads_a_non_arcsecond_resolution(self):
        """Three arc-seconds (SRTMGL3) must not be silently treated as one."""
        path = self._write("3arcsec.tif", 3.0 / 3600.0)
        grid = ss.resolve_grid_geometry(path, -16.0)
        self.assertAlmostEqual(grid.cell_size_deg, 3.0 / 3600.0, places=12)
        self.assertAlmostEqual(grid.cell_size_y, 3 * synthetic.CELL_DEG * 110.6 * 1000.0, places=6)

    def test_explicit_override_wins_over_the_file(self):
        path = self._write("override.tif", synthetic.CELL_DEG)
        grid = ss.resolve_grid_geometry(path, -16.0, cell_size_deg=0.001)
        self.assertAlmostEqual(grid.cell_size_deg, 0.001, places=12)
        self.assertEqual(grid.source, "user-specified")

    def test_falls_back_when_the_file_cannot_be_read(self):
        grid = ss.resolve_grid_geometry(os.path.join(self.tmp, "absent.tif"), -16.0)
        self.assertAlmostEqual(grid.cell_size_deg, ss.DEFAULT_CELL_SIZE_DEG, places=12)
        self.assertIn("assumed", grid.source)


class TestPixelAnisotropy(unittest.TestCase):
    def test_east_west_pixels_are_shorter_than_north_south_in_the_south(self):
        grid = ss.resolve_grid_geometry("nonexistent.tif", -16.0, cell_size_deg=synthetic.CELL_DEG)
        self.assertLess(grid.cell_size_x, grid.cell_size_y)
        self.assertAlmostEqual(grid.cell_size_x / grid.cell_size_y,
                               (111.32 / 110.6) * math.cos(math.radians(grid.center_lat)),
                               places=9)

    def test_pixels_are_nearly_square_at_the_equator(self):
        grid = ss.resolve_grid_geometry("nonexistent.tif", 0.0, cell_size_deg=synthetic.CELL_DEG)
        self.assertAlmostEqual(grid.cell_size_x / grid.cell_size_y, 111.32 / 110.6, places=6)

    def test_longitude_scale_is_taken_at_the_dem_centre_not_its_top_edge(self):
        """
        Evaluating cos(latitude) at the top edge biases the whole map; the centre
        spreads the residual error evenly. A tall DEM makes the difference visible.
        """
        tmp = tempfile.mkdtemp(prefix="sitesearch_centre_")
        try:
            z = np.zeros((10000, 16), dtype=np.float32)   # ~2.8 degrees tall
            path = synthetic.write_geotiff(os.path.join(tmp, "tall.tif"), z, -14.5, -72.0)
            grid = ss.resolve_grid_geometry(path, -14.5)
            self.assertLess(grid.center_lat, -15.8)
            self.assertGreater(grid.center_lat, -16.1)
            at_edge = synthetic.CELL_DEG * 111.32 * math.cos(math.radians(-14.5)) * 1000.0
            self.assertLess(grid.cell_size_x, at_edge)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestNonSquarePixelWarning(unittest.TestCase):
    def test_reports_the_latitude_scale_for_non_square_pixels(self):
        tmp = tempfile.mkdtemp(prefix="sitesearch_nonsq_")
        try:
            import tifffile as tiff
            path = os.path.join(tmp, "nonsquare.tif")
            tiff.imwrite(path, np.zeros((32, 32), dtype=np.float32),
                         extratags=[(33550, "d", 3, (0.0005, 0.00025, 0.0)),
                                    (33922, "d", 6, (0.0, 0.0, 0.0, -72.0, -16.0, 0.0))])
            import contextlib, io
            with contextlib.redirect_stdout(io.StringIO()):
                deg, rows = ss.read_dem_geometry(path)
            self.assertAlmostEqual(deg, 0.00025, places=12)
            self.assertEqual(rows, 32)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
