"""
Verifies the synthetic fixtures themselves.

A fixture that does not have the geometry it claims cannot validate anything built on
it, so these run before the tests that depend on them.
"""

import unittest

import numpy as np

import _support  # noqa: F401  (path setup)
import synthetic


class TestPlanarFixture(unittest.TestCase):
    def test_plane_has_exactly_the_requested_slope_and_aspect(self):
        cell_y, cell_x = synthetic.cell_sizes(-16.0)
        for slope_deg in (3.0, 12.5, 24.0, 40.0):
            for aspect_deg in (0.0, 45.0, 90.0, 180.0, 270.0, 315.0):
                z = synthetic.planar(64, slope_deg, aspect_deg, cell_y, cell_x)
                dy, dx = np.gradient(z.astype(np.float64), cell_y, cell_x)
                slope = np.degrees(np.arctan(np.sqrt(dx ** 2 + dy ** 2)))
                aspect = np.degrees(np.arctan2(-dx, dy)) % 360

                # Interior only: np.gradient uses one-sided differences at the edges.
                # Tolerance is 1e-4 deg: the fixture is stored float32 like a real DEM,
                # which limits agreement to ~1e-6 deg, far below anything physical.
                self.assertAlmostEqual(float(slope[1:-1, 1:-1].mean()), slope_deg, places=4,
                                       msg=f"slope {slope_deg} aspect {aspect_deg}")
                got = float(aspect[1:-1, 1:-1].mean())
                self.assertAlmostEqual((got - aspect_deg + 180) % 360 - 180, 0.0, places=4,
                                       msg=f"slope {slope_deg} aspect {aspect_deg}")


class TestCanyonFixture(unittest.TestCase):
    def setUp(self):
        _, self.cell_x = synthetic.cell_sizes(-15.6)
        self.floor_w = 1000.0
        self.depth = 1500.0
        self.wall = 35.0
        self.n = 400
        self.z = synthetic.canyon(self.n, self.cell_x, self.floor_w, self.depth, self.wall)

    def test_depth_matches_request(self):
        self.assertAlmostEqual(float(self.z.max() - self.z.min()), self.depth, places=3)

    def test_rim_to_rim_separation_matches_closed_form(self):
        expected = synthetic.canyon_rim_separation_m(self.floor_w, self.depth, self.wall)
        profile = self.z[0]
        rim = self.z.max()
        below = np.where(profile < rim - 1e-3)[0]
        measured = (below[-1] - below[0] + 1) * self.cell_x
        # Within one pixel of the analytic value
        self.assertLess(abs(measured - expected), self.cell_x * 1.5,
                        msg=f"measured {measured:.1f} m vs expected {expected:.1f} m")

    def test_walls_have_the_requested_slope(self):
        dy, dx = np.gradient(self.z.astype(np.float64), self.cell_x, self.cell_x)
        slope = np.degrees(np.arctan(np.sqrt(dx ** 2 + dy ** 2)))
        # Sample well inside the west wall, away from floor and rim breaks
        col = int((self.n / 2) - (self.floor_w / 2 + self.depth / 2 / np.tan(np.radians(self.wall))) / self.cell_x)
        self.assertAlmostEqual(float(slope[10, col]), self.wall, places=3)

    def test_colca_like_parameters_reproduce_published_separation(self):
        """~1.5 km deep with 35 degree walls gives roughly Colca's 4.5 km rim separation."""
        sep = synthetic.canyon_rim_separation_m(1000.0, 1500.0, 35.0)
        self.assertGreater(sep, 4000.0)
        self.assertLess(sep, 5500.0)


class TestPeakFixture(unittest.TestCase):
    def test_block_is_where_it_was_asked_for(self):
        z = synthetic.flat_with_peak(200, 100, 150, height=4000.0, half_width=2)
        self.assertEqual(float(z[100, 150]), 4000.0)
        self.assertEqual(float(z[100, 150 - 3]), 0.0)
        self.assertEqual(int(np.count_nonzero(z)), 25)


if __name__ == "__main__":
    unittest.main()
