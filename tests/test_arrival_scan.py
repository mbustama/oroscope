"""
Arrival-direction scanning against terrain whose answer is known in closed form.

The fixtures are chosen so the expected horizon angle, exit distance and column depth
can be worked out with trigonometry rather than read off a previous run.
"""

import math
import unittest

import numpy as np

import _support  # noqa: F401  (path setup)
import arrival_scan as scan_mod
import synthetic
from _support import ss


def grid_at(latitude=-15.6):
    return ss.resolve_grid_geometry("nonexistent.tif", latitude, cell_size_deg=synthetic.CELL_DEG)


class TestGroundUnderfoot(unittest.TestCase):
    """
    A detector standing on the ground has every downward direction blocked by the
    ground at its own feet.

    Over flat terrain a ray tilted below horizontal goes underground within a pixel or
    two and stays there, so its "exit point" is metres away and its column depth is the
    whole traced path. This is correct and it is why the decay-baseline window matters:
    only a site whose local terrain falls away can use the lower half of an acceptance
    window. Without a ``min_dist_km`` these near-field hits swamp everything.
    """

    def setUp(self):
        self.grid = grid_at()
        self.elevation = np.zeros((400, 400), dtype=np.float32)
        self.cands = np.array([[200.0, 200.0, 90.0]])

    def scan(self, **kw):
        params = dict(n_azimuths=1, half_width_deg=0.0,
                      elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
                      max_range_m=5000.0, min_dist_km=0.0, max_dist_km=5.0)
        params.update(kw)
        return scan_mod.scan(self.cands, self.elevation, self.grid, **params)

    def test_downward_directions_strike_the_ground_immediately(self):
        out = self.scan()
        self.assertGreater(int(out["cells"][0]), 0)
        self.assertLess(float(out["mean_distance_m"][0]), 200.0)

    def test_a_decay_baseline_window_removes_them(self):
        self.assertEqual(int(self.scan(min_dist_km=1.0)["cells"][0]), 0)

    def test_upward_directions_are_unobstructed(self):
        """Only the sub-horizontal half is blocked; above the horizon is clear."""
        out = self.scan(elev_min_deg=0.5, elev_max_deg=3.0, n_elev_bins=5)
        self.assertEqual(int(out["cells"][0]), 0)


class TestOpenSky(unittest.TestCase):
    """Flat terrain: nothing to strike, so a neutrino channel sees nothing."""

    def setUp(self):
        self.grid = grid_at()
        self.n = 400
        self.elevation = np.zeros((self.n, self.n), dtype=np.float32)
        self.cands = np.array([[200.0, 200.0, 90.0]])

    def test_no_terrain_means_no_accepted_directions(self):
        out = scan_mod.scan(self.cands, self.elevation, self.grid,
                            max_range_m=5000.0, min_dist_km=1.0, max_dist_km=5.0)
        self.assertEqual(int(out["cells"][0]), 0)
        self.assertEqual(float(out["max_depth_gcm2"][0]), 0.0)

    def test_flat_ground_horizon_sits_just_below_horizontal(self):
        """Curvature alone pulls a flat surface below the horizontal."""
        out = scan_mod.scan(self.cands, self.elevation, self.grid,
                            max_range_m=5000.0, min_dist_km=1.0, max_dist_km=5.0)
        self.assertLess(float(out["horizon_deg"][0]), 0.0)
        self.assertGreater(float(out["horizon_deg"][0]), -0.1)

    def test_cosmic_ray_mode_accepts_the_whole_window_over_flat_ground(self):
        out = scan_mod.scan(self.cands, self.elevation, self.grid,
                            elev_min_deg=0.0, elev_max_deg=30.0, n_elev_bins=10,
                            n_azimuths=8, half_width_deg=None, use_aspect=False,
                            require_terrain=False, max_range_m=5000.0)
        self.assertEqual(int(out["cells"][0]), 8 * 10)
        self.assertGreater(float(out["solid_angle_sr"][0]), 0.0)


class TestHorizonAngle(unittest.TestCase):
    """A wall at a known distance and height subtends a computable angle."""

    def setUp(self):
        self.grid = grid_at()
        self.n = 1200
        self.r0, self.c0 = 600, 200
        self.wall_dist_m = 10000.0
        self.wall_height = 800.0
        self.elevation = np.zeros((self.n, self.n), dtype=np.float32)
        wall_col = self.c0 + int(self.wall_dist_m / self.grid.cell_size_x)
        self.elevation[:, wall_col:wall_col + 40] = self.wall_height
        self.cands = np.array([[float(self.r0), float(self.c0), 90.0]])

    def test_horizon_matches_the_analytic_angle(self):
        out = scan_mod.scan(self.cands, self.elevation, self.grid,
                            n_azimuths=1, half_width_deg=0.0,
                            max_range_m=20000.0, max_dist_km=20.0,
                            elev_min_deg=-3.0, elev_max_deg=10.0, n_elev_bins=26)
        drop = self.wall_dist_m ** 2 / (2 * scan_mod.DEFAULT_EARTH_RADIUS_M)
        expected = math.degrees(math.atan((self.wall_height - drop) / self.wall_dist_m))
        self.assertAlmostEqual(float(out["horizon_deg"][0]), expected, delta=0.05)

    def test_directions_above_the_horizon_are_not_accepted(self):
        """Restricting the window to angles above the wall leaves nothing to hit."""
        out = scan_mod.scan(self.cands, self.elevation, self.grid,
                            n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=6.0, elev_max_deg=10.0, n_elev_bins=8,
                            max_range_m=20000.0, max_dist_km=20.0)
        self.assertEqual(int(out["cells"][0]), 0)

    def test_directions_below_the_horizon_strike_the_wall(self):
        out = scan_mod.scan(self.cands, self.elevation, self.grid,
                            n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
                            max_range_m=20000.0, min_dist_km=1.0, max_dist_km=20.0)
        self.assertGreater(int(out["cells"][0]), 0)
        self.assertAlmostEqual(float(out["mean_distance_m"][0]), self.wall_dist_m,
                               delta=2 * self.grid.cell_size_x)


class TestColumnDepth(unittest.TestCase):
    """A block of known thickness gives a column depth that can be checked by hand."""

    def setUp(self):
        self.grid = grid_at()
        self.n = 1400
        self.r0, self.c0 = 700, 100
        self.elevation = np.zeros((self.n, self.n), dtype=np.float32)
        self.block_start_m = 8000.0
        self.block_len_m = 4000.0
        start = self.c0 + int(self.block_start_m / self.grid.cell_size_x)
        end = self.c0 + int((self.block_start_m + self.block_len_m) / self.grid.cell_size_x)
        self.elevation[:, start:end] = 3000.0        # tall enough to fill the window
        self.cands = np.array([[float(self.r0), float(self.c0), 90.0]])

    def scan(self, **kw):
        params = dict(n_azimuths=1, half_width_deg=0.0,
                      elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                      max_range_m=30000.0, min_dist_km=1.0, max_dist_km=30.0)
        params.update(kw)
        return scan_mod.scan(self.cands, self.elevation, self.grid, **params)

    def test_column_depth_matches_the_block_thickness(self):
        out = self.scan()
        expected_m = self.block_len_m
        expected = expected_m * scan_mod.STANDARD_ROCK_DENSITY * scan_mod.KGM2_TO_GCM2
        self.assertAlmostEqual(float(out["max_depth_gcm2"][0]), expected,
                               delta=0.02 * expected)

    def test_depth_scales_with_rock_density(self):
        light = self.scan(rock_density=1000.0)["max_depth_gcm2"][0]
        heavy = self.scan(rock_density=2000.0)["max_depth_gcm2"][0]
        self.assertAlmostEqual(heavy / light, 2.0, delta=0.01)

    def test_a_depth_threshold_rejects_a_thin_target(self):
        thick = float(self.scan()["max_depth_gcm2"][0])
        accepted = self.scan(min_depth_gcm2=thick * 0.5)
        rejected = self.scan(min_depth_gcm2=thick * 2.0)
        self.assertGreater(int(accepted["cells"][0]), 0)
        self.assertEqual(int(rejected["cells"][0]), 0)

    def test_exit_point_is_the_near_face_not_the_far_one(self):
        out = self.scan()
        self.assertAlmostEqual(float(out["mean_distance_m"][0]), self.block_start_m,
                               delta=3 * self.grid.cell_size_x)

    def test_two_ridges_accumulate_both_chords(self):
        """A ray crossing two ridges traverses the rock of both."""
        one = float(self.scan()["max_depth_gcm2"][0])
        second_start = self.c0 + int(16000.0 / self.grid.cell_size_x)
        second_end = self.c0 + int(18000.0 / self.grid.cell_size_x)
        self.elevation[:, second_start:second_end] = 3000.0
        two = float(self.scan()["max_depth_gcm2"][0])
        extra = 2000.0 * scan_mod.STANDARD_ROCK_DENSITY * scan_mod.KGM2_TO_GCM2
        self.assertAlmostEqual(two - one, extra, delta=0.05 * extra)


class TestDistanceWindow(unittest.TestCase):
    def setUp(self):
        self.grid = grid_at()
        self.n = 1400
        self.r0, self.c0 = 700, 100
        self.elevation = np.zeros((self.n, self.n), dtype=np.float32)
        start = self.c0 + int(15000.0 / self.grid.cell_size_x)
        self.elevation[:, start:start + 200] = 3000.0
        self.cands = np.array([[float(self.r0), float(self.c0), 90.0]])

    def scan(self, **kw):
        params = dict(n_azimuths=1, half_width_deg=0.0,
                      elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                      max_range_m=40000.0)
        params.update(kw)
        return scan_mod.scan(self.cands, self.elevation, self.grid, **params)

    def test_target_inside_the_window_is_accepted(self):
        self.assertGreater(int(self.scan(min_dist_km=10.0, max_dist_km=20.0)["cells"][0]), 0)

    def test_target_beyond_the_window_is_rejected(self):
        self.assertEqual(int(self.scan(min_dist_km=1.0, max_dist_km=10.0)["cells"][0]), 0)

    def test_target_nearer_than_the_window_is_rejected(self):
        self.assertEqual(int(self.scan(min_dist_km=20.0, max_dist_km=40.0)["cells"][0]), 0)


class TestAzimuthFan(unittest.TestCase):
    def test_fan_is_symmetric_about_the_aspect(self):
        offsets = scan_mod.azimuth_fan(9, 60.0)
        self.assertEqual(len(offsets), 9)
        self.assertAlmostEqual(offsets[0], -60.0)
        self.assertAlmostEqual(offsets[-1], 60.0)
        self.assertAlmostEqual(float(offsets.sum()), 0.0, places=9)

    def test_full_sweep_covers_the_compass_without_repeating(self):
        offsets = scan_mod.azimuth_fan(8, None)
        self.assertEqual(len(offsets), 8)
        self.assertAlmostEqual(offsets[0], 0.0)
        self.assertAlmostEqual(offsets[-1], 315.0)

    def test_a_wider_fan_finds_targets_a_narrow_one_misses(self):
        """The measured cost of the old single-ray screen, on a controlled fixture."""
        grid = grid_at()
        n = 1200
        elevation = np.zeros((n, n), dtype=np.float32)
        r0, c0 = 600, 300
        # Ridge placed 40 degrees off the candidate's aspect
        bearing = math.radians(40.0)
        dist = 12000.0
        tr = r0 - int(dist * math.cos(bearing) / grid.cell_size_y)
        tc = c0 + int(dist * math.sin(bearing) / grid.cell_size_x)
        elevation[tr - 60:tr + 60, tc - 60:tc + 60] = 3000.0
        cands = np.array([[float(r0), float(c0), 0.0]])     # aspect due north

        common = dict(elev_min_deg=-3.0, elev_max_deg=3.0, n_elev_bins=12,
                      max_range_m=25000.0, min_dist_km=5.0, max_dist_km=20.0)
        narrow = scan_mod.scan(cands, elevation, grid, n_azimuths=1,
                               half_width_deg=0.0, **common)
        wide = scan_mod.scan(cands, elevation, grid, n_azimuths=9,
                             half_width_deg=60.0, **common)
        self.assertEqual(int(narrow["cells"][0]), 0)
        self.assertGreater(int(wide["cells"][0]), 0)


class TestCanyonGeometry(unittest.TestCase):
    """TAMBO-style: a candidate on one wall looking across at the other."""

    def setUp(self):
        self.grid = grid_at()
        self.n = 600
        _, cell_x = synthetic.cell_sizes(-15.6)
        self.elevation = synthetic.colca_like(self.n, cell_x)
        # Sit on the west rim looking east across the canyon
        self.rim_col = int(((self.n * cell_x) / 2.0
                            - synthetic.canyon_rim_separation_m(**{
                                k: synthetic.COLCA[k] for k in
                                ("floor_width_m", "depth_m", "wall_slope_deg")}) / 2.0)
                           / cell_x)
        self.cands = np.array([[300.0, float(self.rim_col), 90.0]])

    def test_opposite_wall_is_found_at_the_published_separation(self):
        out = scan_mod.scan(self.cands, self.elevation, self.grid,
                            n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-30.0, elev_max_deg=0.0, n_elev_bins=30,
                            max_range_m=12000.0, min_dist_km=0.5, max_dist_km=10.0)
        self.assertGreater(int(out["cells"][0]), 0)
        # Colca's rim-to-rim separation is 4.5 km; the far wall is reached within that
        self.assertLess(float(out["mean_distance_m"][0]), 5000.0)
        self.assertGreater(float(out["mean_distance_m"][0]), 500.0)

    def test_looking_away_from_the_canyon_finds_nothing(self):
        away = np.array([[300.0, float(self.rim_col), 270.0]])
        out = scan_mod.scan(away, self.elevation, self.grid,
                            n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-30.0, elev_max_deg=0.0, n_elev_bins=30,
                            max_range_m=12000.0, min_dist_km=0.5, max_dist_km=10.0)
        self.assertEqual(int(out["cells"][0]), 0)


class TestNoDataAndEdges(unittest.TestCase):
    def test_nan_candidate_elevation_yields_nothing(self):
        grid = grid_at()
        elevation = np.zeros((200, 200), dtype=np.float32)
        elevation[100, 100] = np.nan
        out = scan_mod.scan(np.array([[100.0, 100.0, 90.0]]), elevation, grid,
                            max_range_m=2000.0, min_dist_km=1.0, max_dist_km=2.0)
        self.assertEqual(int(out["cells"][0]), 0)

    def test_scan_walks_off_the_map_without_error(self):
        grid = grid_at()
        elevation = np.zeros((200, 200), dtype=np.float32)
        out = scan_mod.scan(np.array([[5.0, 5.0, 225.0]]), elevation, grid,
                            max_range_m=50000.0, min_dist_km=1.0, max_dist_km=50.0)
        self.assertEqual(int(out["cells"][0]), 0)

    def test_empty_candidate_set_is_handled(self):
        grid = grid_at()
        out = scan_mod.scan(np.zeros((0, 3)), np.zeros((10, 10), dtype=np.float32), grid)
        self.assertEqual(len(out["cells"]), 0)


if __name__ == "__main__":
    unittest.main()
