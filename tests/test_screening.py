"""
Topographic screening: slope/aspect recovery, exclusion zones, candidate thinning
and the funnel accounting.
"""

import contextlib
import io
import unittest

import numpy as np

from _support import ss
import synthetic


def screen(elevation, grid, **kwargs):
    """Runs the screening stage with console output suppressed."""
    with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
        return ss.get_candidates_chunked(elevation, grid, kwargs.pop("rfi_zones", None),
                                         kwargs.pop("origin_lat", -16.0),
                                         kwargs.pop("origin_lon", -72.0), **kwargs)


def make_grid(latitude=-16.0):
    return ss.resolve_grid_geometry("nonexistent.tif", latitude, cell_size_deg=synthetic.CELL_DEG)


class TestSlopeAndAspectScreening(unittest.TestCase):
    def setUp(self):
        self.grid = make_grid()

    def test_plane_inside_the_slope_band_is_kept_entirely(self):
        z = synthetic.planar(200, 10.0, 90.0, self.grid.cell_size_y, self.grid.cell_size_x)
        cands = screen(z, self.grid, tile_size=200, candidate_stride=1,
                       min_slope_deg=3.0, max_slope_deg=25.0)
        # Every interior pixel qualifies; edges use one-sided differences
        self.assertGreater(len(cands), 200 * 200 * 0.97)

    def test_plane_outside_the_slope_band_is_rejected_entirely(self):
        for slope_deg in (1.0, 40.0):
            z = synthetic.planar(200, slope_deg, 90.0, self.grid.cell_size_y, self.grid.cell_size_x)
            cands = screen(z, self.grid, tile_size=200, candidate_stride=1,
                           min_slope_deg=3.0, max_slope_deg=25.0)
            self.assertEqual(len(cands), 0, msg=f"slope {slope_deg} should be excluded")

    def test_reported_aspect_matches_the_plane(self):
        for aspect_deg in (0.0, 90.0, 180.0, 270.0):
            z = synthetic.planar(200, 10.0, aspect_deg, self.grid.cell_size_y, self.grid.cell_size_x)
            cands = screen(z, self.grid, tile_size=200, candidate_stride=1)
            got = cands[:, 2]
            offset = (got - aspect_deg + 180) % 360 - 180
            self.assertLess(float(np.abs(offset).max()), 0.01, msg=f"aspect {aspect_deg}")

    def test_aspect_window_selects_only_matching_terrain(self):
        z = synthetic.planar(200, 10.0, 270.0, self.grid.cell_size_y, self.grid.cell_size_x)
        kept = screen(z, self.grid, tile_size=200, candidate_stride=1,
                      min_aspect_deg=250.0, max_aspect_deg=290.0)
        dropped = screen(z, self.grid, tile_size=200, candidate_stride=1,
                         min_aspect_deg=0.0, max_aspect_deg=45.0)
        self.assertGreater(len(kept), 0)
        self.assertEqual(len(dropped), 0)

    def test_altitude_bounds_bind(self):
        z = synthetic.planar(200, 10.0, 90.0, self.grid.cell_size_y, self.grid.cell_size_x, base=2000.0)
        below = screen(z, self.grid, tile_size=200, candidate_stride=1, max_alt=float(z.min()) - 1.0)
        above = screen(z, self.grid, tile_size=200, candidate_stride=1, min_alt=float(z.max()) + 1.0)
        self.assertEqual(len(below), 0)
        self.assertEqual(len(above), 0)


class TestCandidateStride(unittest.TestCase):
    def test_stride_thins_candidates_proportionally(self):
        grid = make_grid()
        z = synthetic.planar(200, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        full = len(screen(z, grid, tile_size=200, candidate_stride=1))
        for stride in (2, 5, 10):
            thinned = len(screen(z, grid, tile_size=200, candidate_stride=stride))
            self.assertAlmostEqual(thinned / full, 1.0 / stride, delta=0.02,
                                   msg=f"stride {stride}")


class TestExclusionZones(unittest.TestCase):
    """An RFI zone must be a circle on the ground, not in pixel space."""

    def setUp(self):
        self.grid = make_grid()
        self.n = 1200
        self.origin_lat, self.origin_lon = -16.0, -72.0
        self.radius_km = 8.0
        self.z = synthetic.planar(self.n, 10.0, 90.0, self.grid.cell_size_y, self.grid.cell_size_x)
        self.zone_lat = self.origin_lat - (self.n / 2) * synthetic.CELL_DEG
        self.zone_lon = self.origin_lon + (self.n / 2) * synthetic.CELL_DEG

    def _mask(self):
        zones = [("circle", self.zone_lat, self.zone_lon, self.radius_km, "Test")]
        cands = screen(self.z, self.grid, rfi_zones=zones, origin_lat=self.origin_lat,
                       origin_lon=self.origin_lon, tile_size=600, candidate_stride=1)
        mask = np.zeros((self.n, self.n), dtype=bool)
        mask[cands[:, 0].astype(int), cands[:, 1].astype(int)] = True
        return mask

    def test_excluded_region_has_the_requested_radius_on_both_axes(self):
        mask = self._mask()
        mid = self.n // 2
        ew_km = np.count_nonzero(~mask[mid, :]) * self.grid.cell_size_x / 2000.0
        ns_km = np.count_nonzero(~mask[:, mid]) * self.grid.cell_size_y / 2000.0
        self.assertAlmostEqual(ew_km, self.radius_km, delta=0.05)
        self.assertAlmostEqual(ns_km, self.radius_km, delta=0.05)

    def test_zone_is_an_ellipse_in_pixel_space(self):
        """Equal ground radii span more columns than rows, by exactly cell_y/cell_x."""
        mask = self._mask()
        mid = self.n // 2
        cols = np.count_nonzero(~mask[mid, :])
        rows = np.count_nonzero(~mask[:, mid])
        self.assertAlmostEqual(cols / rows, self.grid.cell_size_y / self.grid.cell_size_x,
                               delta=0.01)


class TestFunnelAccounting(unittest.TestCase):
    def test_funnel_records_a_monotonic_screening_cascade(self):
        grid = make_grid()
        z = synthetic.planar(200, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        funnel = ss.Funnel()
        screen(z, grid, tile_size=200, candidate_stride=1, funnel=funnel)
        counts = funnel.as_dict()
        self.assertEqual(counts["DEM pixels"], 200 * 200)
        self.assertEqual(counts["finite elevation"], 200 * 200)
        self.assertLessEqual(counts["slope 3.0-25.0 deg"], counts["finite elevation"])
        self.assertGreater(counts["slope 3.0-25.0 deg"], 0)

    def test_funnel_stride_count_matches_returned_candidates(self):
        grid = make_grid()
        z = synthetic.planar(200, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        funnel = ss.Funnel()
        cands = screen(z, grid, tile_size=64, candidate_stride=5, funnel=funnel)
        self.assertEqual(funnel.get("kept by stride 5"), len(cands))

    def test_funnel_renders_without_stages(self):
        self.assertEqual(ss.Funnel().render(), "")


if __name__ == "__main__":
    unittest.main()
