"""
Arrival-direction scanning against terrain whose answer is known in closed form.

The fixtures are chosen so the expected horizon angle, exit distance and column depth
can be worked out with trigonometry rather than read off a previous run.
"""

import math
import unittest

import numpy as np

import _support  # noqa: F401  (path setup)
from oroscope import arrival_scan as scan_mod
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


class TestTauPhysicsHelpers(unittest.TestCase):
    """Closed-form quantities that need no simulation input."""

    def test_decay_length_reproduces_the_published_tambo_range(self):
        """Ref. [2] Fig. 1 annotates a 50 m - 5 km range for 1-100 PeV."""
        self.assertAlmostEqual(scan_mod.tau_decay_length_m(1.0), 49.0, delta=1.0)
        self.assertAlmostEqual(scan_mod.tau_decay_length_m(100.0), 4900.0, delta=50.0)

    def test_decay_length_is_linear_in_energy(self):
        self.assertAlmostEqual(scan_mod.tau_decay_length_m(20.0)
                               / scan_mod.tau_decay_length_m(2.0), 10.0, places=9)

    def test_inherited_grand_window_corresponds_to_sub_eev_energies(self):
        """The hardcoded 10-80 km default silently encoded an energy assumption."""
        lo = scan_mod.energy_pev_for_decay_length(10000.0)
        hi = scan_mod.energy_pev_for_decay_length(80000.0)
        self.assertAlmostEqual(lo, 204.0, delta=5.0)     # ~0.2 EeV
        self.assertAlmostEqual(hi, 1633.0, delta=20.0)   # ~1.6 EeV

    def test_energy_and_decay_length_round_trip(self):
        for energy in (1.0, 37.0, 1000.0):
            length = scan_mod.tau_decay_length_m(energy)
            self.assertAlmostEqual(scan_mod.energy_pev_for_decay_length(length),
                                   energy, places=6)

    def test_decay_probability_is_bounded_and_peaks_at_the_decay_length(self):
        length = scan_mod.tau_decay_length_m(100.0)
        p_wide = scan_mod.decay_probability(0.0, 1e9, 100.0)
        self.assertAlmostEqual(p_wide, 1.0, places=6)
        p_one_l = scan_mod.decay_probability(0.0, length, 100.0)
        self.assertAlmostEqual(p_one_l, 1 - math.exp(-1), places=6)
        self.assertEqual(scan_mod.decay_probability(1000.0, 1000.0, 100.0), 0.0)

    def test_energy_window_sets_a_few_km_baseline_for_tambo(self):
        lo, hi = scan_mod.distance_window_from_energy(1.0, 100.0)
        self.assertLess(lo, 100.0)
        self.assertGreater(hi, 4000.0)
        self.assertLess(hi, 10000.0)


class TestOrientationIsSubsumedByDepth(unittest.TestCase):
    """
    Column depth already encodes what a "target must face me" test was proxying for.

    A face sloping away from the observer is never struck at all, and among faces that
    are struck, a shallower one presents more rock for the same summit height. So the
    depth measurement subsumes the orientation criterion (roadmap 4.6) rather than
    needing it alongside.

    Note the direction of the effect: a *gentler* face gives greater column depth,
    because a near-horizontal ray runs further underground before reaching the summit.
    That is one reason the depth criterion wants an optimum band and not a floor --
    the tau has to escape as well as be produced.
    """

    def setUp(self):
        self.grid = grid_at()
        self.n = 1600
        self.r0, self.c0 = 800, 100
        self.summit = 2000.0
        self.start_m = 6000.0

    def ramp(self, face_slope_deg, descending=False):
        """
        A ridge of fixed summit height whose near face has the requested slope.

        The far side drops away, so the rock a ray traverses is set by the face's
        horizontal run rather than by how far the scan happens to reach.
        """
        z = np.zeros((self.n, self.n), dtype=np.float32)
        start = self.c0 + int(self.start_m / self.grid.cell_size_x)
        run_px = max(1, int(self.summit / math.tan(math.radians(face_slope_deg))
                            / self.grid.cell_size_x))
        cols = np.arange(self.n)
        if descending:
            profile = -np.clip(cols - start, 0, None) * self.grid.cell_size_x \
                      * math.tan(math.radians(face_slope_deg))
        else:
            profile = np.clip((cols - start) / run_px, 0.0, 1.0) * self.summit
            profile[cols > start + run_px] = 0.0
        z[:, :] = profile[None, :]
        return z

    def scan(self, elevation, **kw):
        params = dict(n_azimuths=1, half_width_deg=0.0,
                      elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                      max_range_m=45000.0, min_dist_km=1.0, max_dist_km=45.0)
        params.update(kw)
        cands = np.array([[float(self.r0), float(self.c0), 90.0]])
        return scan_mod.scan(cands, elevation, self.grid, **params)

    def test_a_face_sloping_away_is_never_struck(self):
        out = self.scan(self.ramp(20.0, descending=True))
        self.assertEqual(int(out["cells"][0]), 0)

    def test_a_face_sloping_toward_the_observer_is_struck(self):
        out = self.scan(self.ramp(20.0))
        self.assertGreater(int(out["cells"][0]), 0)

    def test_a_gentler_face_presents_more_rock(self):
        shallow = float(self.scan(self.ramp(5.0))["max_depth_gcm2"][0])
        steep = float(self.scan(self.ramp(30.0))["max_depth_gcm2"][0])
        self.assertGreater(shallow, steep)


class TestFresnelClearance(unittest.TestCase):
    """
    Clearance is measured only when a band is given, and against intervening terrain
    rather than against the ground beside the antenna.
    """

    def setUp(self):
        """
        A candidate on a shoulder whose ground falls away, looking at a distant wall.

        Flat ground would be the wrong fixture: a near-horizontal ray skims it for
        kilometres, so the ground itself fills the first Fresnel zone and every path
        scores poorly. Real sites score well precisely because their terrain drops
        away, which is what this reproduces.
        """
        self.grid = grid_at()
        self.n = 1400
        self.r0, self.c0 = 700, 100
        cols = np.arange(self.n)
        profile = np.clip(800.0 - (cols - self.c0) * self.grid.cell_size_x * 0.8, 0.0, 800.0)
        profile[:self.c0] = 800.0
        self.elevation = np.repeat(profile[None, :], self.n, axis=0).astype(np.float32)
        start = self.c0 + int(15000.0 / self.grid.cell_size_x)
        self.elevation[:, start:start + 300] = 3000.0
        self.cands = np.array([[float(self.r0), float(self.c0), 90.0]])

    def scan(self, **kw):
        params = dict(n_azimuths=1, half_width_deg=0.0,
                      elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                      max_range_m=30000.0, min_dist_km=5.0, max_dist_km=30.0)
        params.update(kw)
        return scan_mod.scan(self.cands, self.elevation, self.grid, **params)

    def test_not_measured_without_a_frequency(self):
        out = self.scan()
        self.assertEqual(float(out["best_clearance_ratio"][0]), 0.0)

    def test_measured_when_a_frequency_is_given(self):
        out = self.scan(frequency_mhz=50.0)
        self.assertGreater(float(out["best_clearance_ratio"][0]), 0.0)

    def test_a_clear_path_clears_many_fresnel_radii(self):
        """Flat ground to a distant wall should not obstruct anything."""
        out = self.scan(frequency_mhz=50.0)
        self.assertGreater(float(out["best_clearance_ratio"][0]), 1.0)

    def test_an_intervening_ridge_reduces_the_clearance(self):
        clear = float(self.scan(frequency_mhz=50.0)["best_clearance_ratio"][0])
        # 650 m at 8 km subtends about -1.1 deg from the 800 m shoulder, just below the
        # window, so it obstructs the path without ever becoming a target itself.
        blocked = self.elevation.copy()
        mid = self.c0 + int(8000.0 / self.grid.cell_size_x)
        blocked[:, mid:mid + 10] = 650.0
        out = scan_mod.scan(self.cands, blocked, self.grid,
                            n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                            max_range_m=30000.0, min_dist_km=5.0, max_dist_km=30.0,
                            frequency_mhz=50.0)
        self.assertLess(float(out["best_clearance_ratio"][0]), clear)

    def test_higher_frequency_needs_less_clearance(self):
        """r1 shrinks as sqrt(lambda), so the same terrain clears more radii."""
        low = float(self.scan(frequency_mhz=50.0)["best_clearance_ratio"][0])
        high = float(self.scan(frequency_mhz=200.0)["best_clearance_ratio"][0])
        self.assertAlmostEqual(high / low, 2.0, delta=0.1)

    def test_near_field_exclusion_removes_the_mast_height_sensitivity(self):
        """
        Without it the measure is dominated by ground beside the antenna and swings by
        more than an order of magnitude with mast height, which has nothing to do with
        site quality. On real Andean terrain the spread falls from 28x to about 2x once
        the first 500 m are skipped.
        """
        def spread(near_field):
            vals = [float(self.scan(frequency_mhz=50.0, antenna_height_m=h,
                                    near_field_m=near_field)["best_clearance_ratio"][0])
                    for h in (0.0, 5.0, 20.0)]
            return max(vals) / max(min(vals), 1e-9)

        self.assertLess(spread(500.0), spread(0.0))


class TestTwoRadii(unittest.TestCase):
    """
    Particle trajectories and radio signals do not curve the same way.

    Neutrinos and taus are not refracted, so the geometry that decides where the tau
    exits uses the true Earth radius. The radio signal is refracted by the tropospheric
    density gradient, which the 4/3 convention absorbs -- and that applies to the
    Fresnel clearance of the signal path and to nothing else.
    """

    def test_particle_geometry_defaults_to_the_true_radius(self):
        self.assertAlmostEqual(scan_mod.TRUE_EARTH_RADIUS_M, 6.371e6, places=1)
        self.assertAlmostEqual(scan_mod.DEFAULT_EARTH_RADIUS_M,
                               scan_mod.TRUE_EARTH_RADIUS_M, places=1)

    def test_radio_path_keeps_the_four_thirds_convention(self):
        self.assertAlmostEqual(scan_mod.RADIO_EARTH_RADIUS_M, 8.5e6, places=1)

    def test_the_radio_radius_does_not_affect_where_the_tau_exits(self):
        """Changing the radio radius must leave the particle geometry untouched."""
        grid = grid_at()
        n = 1400
        elevation = np.zeros((n, n), dtype=np.float32)
        col = 200 + int(20000.0 / grid.cell_size_x)
        elevation[:, col:col + 60] = 2000.0
        cands = np.array([[700.0, 200.0, 90.0]])
        common = dict(n_azimuths=1, half_width_deg=0.0, elev_min_deg=-1.0,
                      elev_max_deg=1.0, n_elev_bins=4, max_range_m=30000.0,
                      min_dist_km=5.0, max_dist_km=30.0)
        a = scan_mod.scan(cands, elevation, grid, radio_earth_radius_m=8.5e6, **common)
        b = scan_mod.scan(cands, elevation, grid, radio_earth_radius_m=6.371e6, **common)
        self.assertEqual(float(a["horizon_deg"][0]), float(b["horizon_deg"][0]))
        self.assertEqual(float(a["mean_distance_m"][0]), float(b["mean_distance_m"][0]))


class TestFresnelDefaults(unittest.TestCase):
    def test_antenna_height_defaults_to_two_metres(self):
        import inspect
        sig = inspect.signature(scan_mod.scan)
        self.assertEqual(sig.parameters["antenna_height_m"].default, 2.0)

    def test_near_field_is_excluded_by_default(self):
        import inspect
        sig = inspect.signature(scan_mod.scan)
        self.assertEqual(sig.parameters["near_field_m"].default, 500.0)

    def test_near_field_can_be_included_by_setting_it_to_zero(self):
        grid = grid_at()
        n = 1400
        cols = np.arange(n)
        profile = np.clip(800.0 - (cols - 100) * grid.cell_size_x * 0.8, 0.0, 800.0)
        profile[:100] = 800.0
        elevation = np.repeat(profile[None, :], n, axis=0).astype(np.float32)
        start = 100 + int(15000.0 / grid.cell_size_x)
        elevation[:, start:start + 300] = 3000.0
        cands = np.array([[700.0, 100.0, 90.0]])
        common = dict(n_azimuths=1, half_width_deg=0.0, elev_min_deg=-1.0,
                      elev_max_deg=1.0, n_elev_bins=4, max_range_m=30000.0,
                      min_dist_km=5.0, max_dist_km=30.0, frequency_mhz=50.0)
        excluded = scan_mod.scan(cands, elevation, grid, near_field_m=500.0, **common)
        included = scan_mod.scan(cands, elevation, grid, near_field_m=0.0, **common)
        # Including it can only lower the measure: more of the path is examined
        self.assertLessEqual(float(included["best_clearance_ratio"][0]),
                             float(excluded["best_clearance_ratio"][0]))


class TestRefractionKFactor(unittest.TestCase):
    def test_k_factor_maps_to_the_conventional_radius(self):
        self.assertAlmostEqual(scan_mod.earth_radius_for_k(1.0), 6371000.0, places=3)
        self.assertAlmostEqual(scan_mod.earth_radius_for_k(4.0 / 3.0), 8494666.67, places=1)

    def test_true_geometry_lowers_distant_terrain_more_than_the_radio_convention(self):
        """At 80 km the drop is 502 m at k=1 against 376 m at k=4/3."""
        d = 80000.0
        drop_true = d ** 2 / (2 * scan_mod.earth_radius_for_k(1.0))
        drop_radio = d ** 2 / (2 * scan_mod.earth_radius_for_k(4.0 / 3.0))
        self.assertAlmostEqual(drop_true, 502.0, delta=2.0)
        self.assertAlmostEqual(drop_radio, 376.0, delta=2.0)

    def test_a_smaller_radius_makes_a_marginal_target_disappear(self):
        grid = grid_at()
        n = 2200
        r0, c0 = 1100, 100
        elevation = np.zeros((n, n), dtype=np.float32)
        col = c0 + int(70000.0 / grid.cell_size_x)
        elevation[:, col:col + 60] = 450.0     # clears the k=4/3 drop, not the k=1 drop
        cands = np.array([[float(r0), float(c0), 90.0]])
        common = dict(n_azimuths=1, half_width_deg=0.0, elev_min_deg=-1.0,
                      elev_max_deg=1.0, n_elev_bins=4, max_range_m=80000.0,
                      min_dist_km=5.0, max_dist_km=80.0)
        radio = scan_mod.scan(cands, elevation, grid,
                              earth_radius_m=scan_mod.earth_radius_for_k(4 / 3), **common)
        true_geom = scan_mod.scan(cands, elevation, grid,
                                  earth_radius_m=scan_mod.earth_radius_for_k(1.0), **common)
        self.assertGreater(float(radio["horizon_deg"][0]), float(true_geom["horizon_deg"][0]))


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


class TestBinningBoundary(unittest.TestCase):
    """
    Regression: terrain just below the acceptance window must not be counted.

    The original binned by ``k = int((theta - elev_min)/bin)`` and kept samples with
    ``k >= 0``. C-style int() truncates toward zero, so a sample up to one bin *below*
    elev_min gave k = 0 and was added to the lowest bin, inflating its column depth by
    as much as a third. Working in slope against pre-computed tangent edges removes
    the boundary entirely.
    """

    def test_terrain_below_the_window_is_excluded_from_the_depth(self):
        grid = grid_at()
        n = 1200
        r0, c0 = 600, 100
        # A slab whose top subtends an angle just below the window's lower edge
        elevation = np.zeros((n, n), dtype=np.float32)
        start = c0 + int(6000.0 / grid.cell_size_x)
        drop = -6000.0 * math.tan(math.radians(1.4))     # about -1.4 deg, below -1.0
        elevation[:, start:start + 400] = drop
        cands = np.array([[float(r0), float(c0), 90.0]])
        out = scan_mod.scan(cands, elevation, grid, n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                            max_range_m=20000.0, min_dist_km=1.0, max_dist_km=20.0)
        self.assertEqual(int(out["cells"][0]), 0)
        self.assertEqual(float(out["max_depth_gcm2"][0]), 0.0)

    def test_a_target_above_the_window_credits_the_upward_bins(self):
        """
        A wall subtending 20 degrees is far above a +/-1 degree window, so it lands in
        the overflow bin. The inclusive suffix sum means overflow counts toward every
        bin, so the upward bins are still credited with its rock.

        Only two of the four bins are accepted: the two downward ones strike the ground
        underfoot well inside min_dist and are rejected there, not here.
        """
        grid = grid_at()
        n = 1200
        elevation = np.zeros((n, n), dtype=np.float32)
        col = 100 + int(8000.0 / grid.cell_size_x)
        elevation[:, col:col + 200] = 3000.0
        cands = np.array([[600.0, 100.0, 90.0]])
        out = scan_mod.scan(cands, elevation, grid, n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                            max_range_m=20000.0, min_dist_km=1.0, max_dist_km=20.0)
        self.assertEqual(int(out["cells"][0]), 2)
        self.assertGreater(float(out["max_depth_gcm2"][0]), 0.0)
        self.assertAlmostEqual(float(out["mean_distance_m"][0]), 8000.0,
                               delta=3 * grid.cell_size_x)


class TestProfileSampling(unittest.TestCase):
    """
    Sub-pixel sampling of the terrain profile.

    Nearest-neighbour quantises the profile to pixel centres and treats terrain as
    piecewise constant, which over-estimates how much it blocks a ray: the ray is
    stopped by a whole pixel's worth of the nearby maximum. It also samples through
    ``int()``, which truncates toward zero and so biases the sample point back toward
    the candidate by up to half a pixel, asymmetrically in azimuth.
    """

    def test_bilinear_is_exact_on_a_plane(self):
        """
        Bilinear interpolation reproduces linear functions exactly, so on a planar DEM
        the sampled profile is the true one and the horizon is analytic.
        """
        grid = grid_at()
        n = 800
        slope_deg = 4.0
        z = synthetic.planar(n, slope_deg, 90.0, grid.cell_size_y, grid.cell_size_x,
                             base=5000.0)
        # Look uphill: due west, since aspect 90 means the surface falls to the east
        cands = np.array([[400.0, 600.0, 270.0]])
        out = scan_mod.scan(cands, z, grid, n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=0.0, elev_max_deg=6.0, n_elev_bins=12,
                            max_range_m=8000.0, min_dist_km=0.5, max_dist_km=8.0,
                            bilinear=True)
        # Curvature pulls the far end down slightly, so the horizon sits just under the
        # plane's own slope
        self.assertLess(float(out["horizon_deg"][0]), slope_deg)
        self.assertGreater(float(out["horizon_deg"][0]), slope_deg - 0.2)

    def test_both_modes_agree_on_a_plane(self):
        """With no sub-pixel structure there is nothing for interpolation to recover."""
        grid = grid_at()
        z = synthetic.planar(800, 4.0, 90.0, grid.cell_size_y, grid.cell_size_x, base=5000.0)
        cands = np.array([[400.0, 600.0, 270.0]])
        kw = dict(n_azimuths=1, half_width_deg=0.0, elev_min_deg=0.0, elev_max_deg=6.0,
                  n_elev_bins=12, max_range_m=8000.0, min_dist_km=0.5, max_dist_km=8.0)
        a = scan_mod.scan(cands, z, grid, bilinear=False, **kw)
        b = scan_mod.scan(cands, z, grid, bilinear=True, **kw)
        self.assertAlmostEqual(float(a["horizon_deg"][0]), float(b["horizon_deg"][0]),
                               delta=0.05)

    def test_nearest_sampling_blocks_more_than_bilinear(self):
        """
        A narrow ridge stops a nearest-neighbour ray over a wider band of directions,
        because the ray is blocked by a whole pixel rather than by the interpolated
        surface it actually grazes.
        """
        grid = grid_at()
        n = 1000
        z = np.zeros((n, n), dtype=np.float32)
        col = 200 + int(6000.0 / grid.cell_size_x)
        z[:, col:col + 2] = 400.0                  # a two-pixel ridge
        cands = np.array([[500.0, 200.0, 90.0]])
        kw = dict(n_azimuths=1, half_width_deg=0.0, elev_min_deg=-1.0, elev_max_deg=5.0,
                  n_elev_bins=24, max_range_m=12000.0, min_dist_km=1.0, max_dist_km=12.0)
        near = scan_mod.scan(cands, z, grid, bilinear=False, **kw)
        bil = scan_mod.scan(cands, z, grid, bilinear=True, **kw)
        # Tolerance, because on this fixture the two agree to seven digits and the last
        # of them is not a property of the sampling. Compiled with fastmath the nearest
        # form came out 7e-7 deg the larger; interpreted, under NUMBA_DISABLE_JIT=1, the
        # bilinear form did -- so an exact >= turned a floating-point reassociation into
        # a failing test, and made the only route to coverage of these kernels look like
        # a regression. The claim is that nearest never sees *further*, not that it
        # differs at the eighth digit.
        self.assertGreaterEqual(float(near["horizon_deg"][0]),
                                float(bil["horizon_deg"][0]) - 1.0e-5)

    def test_nodata_neighbours_fall_back_to_nearest(self):
        grid = grid_at()
        z = np.zeros((600, 600), dtype=np.float32)
        z[:, 400:] = np.nan
        cands = np.array([[300.0, 100.0, 90.0]])
        out = scan_mod.scan(cands, z, grid, n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                            max_range_m=10000.0, min_dist_km=1.0, max_dist_km=10.0,
                            bilinear=True)
        self.assertEqual(int(out["cells"][0]), 0)     # no crash, nothing spurious


class TestTargetSlopeCriterion(unittest.TestCase):
    """
    Is the target a wall, or merely ground?

    The scan asks whether rock lies at the right range and bearing, which on real
    terrain is nearly always true somewhere -- 92% of Andean candidates passed a
    canyon-shaped criterion before this existed. TAMBO's geometry needs more than the
    presence of rock: the tau exits a *wall*, and a gentle rise at the same distance is
    not the same target. The slope is measured along the arrival azimuth, so an
    obliquely-viewed wall counts as the tau would actually cross it.

    The canyon fixture has walls of known slope, so the recovered value is checkable
    against arithmetic rather than against a previous run.

    Note the un-filtered mean is deliberately *not* the wall slope: rays aimed lower
    strike the flat canyon floor, whose slope really is zero, so the mean over all
    accepted directions is a mixture. Filtering is what isolates the wall.
    """

    FLOOR_W, DEPTH = 1000.0, 1500.0

    def setUp(self):
        self.n = 400
        self.cell_y, self.cell_x = synthetic.cell_sizes(-15.6)
        self.grid = grid_at()

    def scan_across(self, wall_slope_deg, **kw):
        """One candidate part-way down the west wall, looking east across the canyon."""
        z = synthetic.canyon(self.n, self.cell_x, floor_width_m=self.FLOOR_W,
                             depth_m=self.DEPTH, wall_slope_deg=wall_slope_deg)
        col = int((self.n * self.cell_x / 2.0 - 1200.0) / self.cell_x)
        cands = np.array([[self.n // 2, col, 90.0]], dtype=np.float64)
        params = dict(elev_min_deg=-25.0, elev_max_deg=25.0, n_elev_bins=25,
                      # use_aspect, so the single azimuth is the candidate's own 90 deg
                      # -- due east, across the canyon. With use_aspect off the fan
                      # returns 0 deg and looks north, along a canyon that is uniform
                      # north-south and so has no wall to find.
                      n_azimuths=1, half_width_deg=0.0, use_aspect=True,
                      min_dist_km=0.3, max_dist_km=8.0, max_range_m=8000.0)
        params.update(kw)
        return scan_mod.scan(cands, z, self.grid, **params)

    def test_recovers_the_wall_slope_it_was_built_with(self):
        """Filtered to wall hits, the measured slope is the fixture's own parameter."""
        for wall in (15.0, 25.0, 35.0, 45.0):
            out = self.scan_across(wall, min_target_slope_deg=wall - 5.0)
            self.assertGreater(int(out["cells"][0]), 0, f"{wall} deg wall: nothing accepted")
            got = float(out["target_slope_deg"][0])
            self.assertAlmostEqual(got, wall, delta=0.5,
                                   msg=f"built a {wall} deg wall, measured {got:.1f}")

    def test_unfiltered_mean_sits_between_the_floor_and_the_wall(self):
        wall = 35.0
        out = self.scan_across(wall)
        got = float(out["target_slope_deg"][0])
        self.assertGreater(got, 0.0)
        self.assertLess(got, wall, "rays reaching the flat floor must pull the mean down")

    def test_a_steep_requirement_rejects_a_shallow_wall(self):
        self.assertEqual(int(self.scan_across(15.0, min_target_slope_deg=30.0)["cells"][0]), 0)

    def test_the_same_requirement_accepts_a_steep_wall(self):
        self.assertGreater(int(self.scan_across(45.0, min_target_slope_deg=30.0)["cells"][0]), 0)

    def test_the_criterion_only_ever_removes_directions(self):
        base = self.scan_across(45.0)
        cut = self.scan_across(45.0, min_target_slope_deg=30.0)
        self.assertLessEqual(int(cut["cells"][0]), int(base["cells"][0]))
        self.assertLessEqual(float(cut["solid_angle_sr"][0]),
                             float(base["solid_angle_sr"][0]) + 1e-12)

    def test_unset_by_default_so_grand_is_unaffected(self):
        implicit = self.scan_across(35.0)
        explicit = self.scan_across(35.0, min_target_slope_deg=None, max_target_slope_deg=None)
        self.assertEqual(int(implicit["cells"][0]), int(explicit["cells"][0]))
        self.assertAlmostEqual(float(implicit["target_slope_deg"][0]),
                               float(explicit["target_slope_deg"][0]), places=9)

    def test_an_upper_bound_excludes_the_wall_but_keeps_the_floor(self):
        """
        A ceiling does not empty the result: the canyon floor is flatter than any
        ceiling worth setting, so what it removes is the wall itself.
        """
        base = self.scan_across(45.0)
        capped = self.scan_across(45.0, max_target_slope_deg=30.0)
        self.assertGreater(int(capped["cells"][0]), 0)
        self.assertLess(int(capped["cells"][0]), int(base["cells"][0]))
        self.assertLessEqual(float(capped["target_slope_deg"][0]), 30.0)

    def test_the_flat_floor_is_reported_as_flat(self):
        """A direction that strikes the canyon floor must measure ~0, not the wall."""
        out = self.scan_across(35.0, max_target_slope_deg=5.0)
        self.assertGreater(int(out["cells"][0]), 0, "floor hits should still be accepted")
        self.assertLess(abs(float(out["target_slope_deg"][0])), 6.0)
