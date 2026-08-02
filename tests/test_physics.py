"""
Ray-casting kernel: does the ray go where the aspect says, and do the bounds bind?

The A/B against an isotropic pixel scale pins the bug fixed in commit a9843f9, where
stepping in raw pixels skewed every ray away from its stated aspect.
"""

import unittest

import numpy as np

from _support import ss
import synthetic


class RayCastingBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cell_y, cls.cell_x = synthetic.cell_sizes(-16.0)
        cls.n = 1400
        cls.r0, cls.c0 = 700, 300
        cls.dist_m = 20000.0
        # Targets are several pixels across, as real ridges are. A single-pixel target
        # makes the test sensitive to the kernel's truncation convention rather than
        # to its direction, which produces false failures.
        cls.half = 3

    def cast(self, elevation, aspect_deg, min_km=5.0, max_km=40.0,
             cell_y=None, cell_x=None, fresnel=200.0):
        cand = np.array([[self.r0, self.c0, aspect_deg]], dtype=np.float64)
        hits_r, _ = ss.check_physics_chunk(
            cand, elevation,
            self.cell_y if cell_y is None else cell_y,
            self.cell_x if cell_x is None else cell_x,
            self.n, self.n, fresnel, min_km, max_km,
        )
        return len(hits_r) == 1


class TestRayDirection(RayCastingBase):
    def setUp(self):
        self.col_east = self.c0 + int(self.dist_m / self.cell_x)
        self.row_north = self.r0 - int(self.dist_m / self.cell_y)
        self.east_target = synthetic.flat_with_peak(self.n, self.r0, self.col_east,
                                                    height=5000.0, half_width=self.half)
        self.north_target = synthetic.flat_with_peak(self.n, self.row_north, self.c0,
                                                     height=5000.0, half_width=self.half)

    def test_finds_target_due_east_when_facing_east(self):
        self.assertTrue(self.cast(self.east_target, 90.0))

    def test_ignores_target_due_east_when_facing_north(self):
        self.assertFalse(self.cast(self.east_target, 0.0))

    def test_finds_target_due_north_when_facing_north(self):
        self.assertTrue(self.cast(self.north_target, 0.0))

    def test_ignores_target_due_north_when_facing_east(self):
        self.assertFalse(self.cast(self.north_target, 90.0))

    def test_isotropic_pixel_scale_misses_the_eastward_target(self):
        """
        Regression for a9843f9. Using the north-south scale on both axes points the
        ray at a column ~22 px short of the target, so it must miss what the
        corrected kernel finds.
        """
        wrong = self.cast(self.east_target, 90.0, cell_x=self.cell_y)
        self.assertTrue(self.cast(self.east_target, 90.0))
        self.assertFalse(wrong)


class TestRayBounds(RayCastingBase):
    def setUp(self):
        col_east = self.c0 + int(self.dist_m / self.cell_x)
        self.target = synthetic.flat_with_peak(self.n, self.r0, col_east,
                                               height=5000.0, half_width=self.half)

    def test_target_beyond_max_distance_is_rejected(self):
        self.assertFalse(self.cast(self.target, 90.0, min_km=5.0, max_km=15.0))

    def test_target_inside_min_distance_is_rejected(self):
        self.assertFalse(self.cast(self.target, 90.0, min_km=25.0, max_km=40.0))

    def test_target_below_height_threshold_is_rejected(self):
        """Must clear detector + 1 km interaction depth + Fresnel buffer."""
        low = self.target.copy()
        low[low > 0] = 900.0
        self.assertFalse(self.cast(low, 90.0))

    def test_raising_the_fresnel_buffer_can_reject_a_marginal_target(self):
        marginal = self.target.copy()
        marginal[marginal > 0] = 1300.0        # clears 1000 + 200, but not 1000 + 500
        self.assertTrue(self.cast(marginal, 90.0, fresnel=200.0))
        self.assertFalse(self.cast(marginal, 90.0, fresnel=500.0))

    def test_earth_curvature_lowers_distant_targets(self):
        """
        A target just above threshold at short range fails at long range once the
        d^2/2R drop is applied. At 20 km the drop is ~24 m; at 80 km, ~376 m.
        """
        near_col = self.c0 + int(10000.0 / self.cell_x)
        far_col = self.c0 + int(70000.0 / self.cell_x)
        height = 1000.0 + 200.0 + 150.0        # clears threshold by 150 m before curvature
        near = synthetic.flat_with_peak(self.n, self.r0, near_col, height, self.half)
        far = synthetic.flat_with_peak(self.n, self.r0, far_col, height, self.half)
        self.assertTrue(self.cast(near, 90.0, min_km=5.0, max_km=80.0))
        self.assertFalse(self.cast(far, 90.0, min_km=60.0, max_km=80.0))


class TestNoDataHandling(RayCastingBase):
    def test_nan_samples_do_not_count_as_targets(self):
        elevation = np.full((self.n, self.n), 0.0, dtype=np.float32)
        elevation[:, self.c0 + 100:] = np.nan
        self.assertFalse(self.cast(elevation, 90.0))


if __name__ == "__main__":
    unittest.main()
