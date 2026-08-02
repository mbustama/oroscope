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


class TestTilingInvariance(unittest.TestCase):
    """
    Screening must not depend on how the map is cut into tiles.

    Derivatives are computed on a haloed block and cropped, so a pixel's slope is the
    same whether or not a tile boundary happens to fall next to it. Without the halo,
    np.gradient falls back to one-sided differences at every tile edge.
    """

    def candidates(self, z, grid, tile_size, **kwargs):
        cands = screen(z, grid, tile_size=tile_size, candidate_stride=1, **kwargs)
        return set(map(tuple, cands[:, :2].astype(int)))

    def test_tiled_screening_matches_untiled(self):
        grid = make_grid(-15.6)
        n = 600
        z = synthetic.ridge_and_slope(n, grid.cell_size_x)
        untiled = self.candidates(z, grid, n)
        for tile_size in (64, 128, 256):
            self.assertEqual(self.candidates(z, grid, tile_size), untiled,
                             msg=f"tile_size {tile_size} disagrees with the untiled result")

    def test_tiled_screening_matches_untiled_with_a_slope_baseline(self):
        grid = make_grid(-15.6)
        n = 600
        z = synthetic.ridge_and_slope(n, grid.cell_size_x)
        untiled = self.candidates(z, grid, n, slope_baseline_m=500.0)
        for tile_size in (128, 256):
            self.assertEqual(self.candidates(z, grid, tile_size, slope_baseline_m=500.0), untiled,
                             msg=f"tile_size {tile_size} disagrees at a 500 m baseline")


class TestSlopeBaseline(unittest.TestCase):
    """Slope is scale-dependent; the baseline makes the scale explicit."""

    def test_baseline_converts_to_a_per_axis_pixel_window(self):
        grid = make_grid(-16.0)
        ny, nx = ss.slope_baseline_pixels(grid, 1000.0)
        self.assertEqual(ny, round(1000.0 / grid.cell_size_y))
        self.assertEqual(nx, round(1000.0 / grid.cell_size_x))
        self.assertNotEqual(ny, nx, "anisotropic pixels give different windows per axis")

    def test_no_baseline_means_native_resolution(self):
        grid = make_grid(-16.0)
        self.assertEqual(ss.slope_baseline_pixels(grid, None), (0, 0))
        self.assertEqual(ss.slope_baseline_pixels(grid, 0), (0, 0))

    def test_a_plane_keeps_its_slope_at_every_baseline(self):
        """Smoothing must not bias slope on terrain that has no curvature."""
        grid = make_grid(-16.0)
        z = synthetic.planar(300, 12.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        for baseline in (None, 200.0, 1000.0):
            smooth = ss.slope_baseline_pixels(grid, baseline)
            slope, aspect = ss.terrain_derivatives(z, grid.cell_size_y, grid.cell_size_x, *smooth)
            interior = slope[50:-50, 50:-50]
            self.assertAlmostEqual(float(interior.mean()), 12.0, places=3,
                                   msg=f"baseline {baseline}")

    def test_longer_baselines_smooth_rough_terrain(self):
        """On real-shaped terrain a wider window lowers the measured slope."""
        grid = make_grid(-16.0)
        rng = np.random.default_rng(0)
        z = synthetic.planar(400, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        z = z + rng.normal(0, 15.0, z.shape).astype(np.float32)
        medians = []
        for baseline in (None, 250.0, 1000.0):
            smooth = ss.slope_baseline_pixels(grid, baseline)
            slope, _ = ss.terrain_derivatives(z, grid.cell_size_y, grid.cell_size_x, *smooth)
            medians.append(float(np.median(slope[20:-20, 20:-20])))
        self.assertEqual(medians, sorted(medians, reverse=True),
                         msg=f"expected decreasing slope with baseline, got {medians}")


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
