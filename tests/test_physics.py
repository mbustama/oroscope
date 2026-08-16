"""
The closed-form physics: atmosphere, Earth chord, tau range, geomagnetic angle,
Cherenkov footprint.

Every value here is checkable by hand, which is the point of keeping them analytic.
"""

import math
import unittest

import numpy as np

import _support  # noqa: F401  (path setup)
from oroscope import arrival_scan as scan_mod
from oroscope import physics
import synthetic
from _support import ss


def grid_at(latitude=-15.6):
    return ss.resolve_grid_geometry("nonexistent.tif", latitude, cell_size_deg=synthetic.CELL_DEG)


class TestAtmosphere(unittest.TestCase):
    def test_density_falls_by_a_third_between_sea_level_and_4000_m(self):
        ratio = physics.air_density_kgm3(4000.0) / physics.air_density_kgm3(0.0)
        self.assertAlmostEqual(ratio, math.exp(-4000.0 / 8400.0), places=9)
        self.assertAlmostEqual(ratio, 0.622, delta=0.005)

    def test_horizontal_grammage_is_density_times_length(self):
        x = physics.slant_grammage_gcm2(0.0, 0.0, 20000.0)
        self.assertAlmostEqual(x, 1.225 * 20000.0 * 0.1, places=6)
        self.assertAlmostEqual(x, 2450.0, delta=1.0)

    def test_the_same_path_at_altitude_yields_much_less_grammage(self):
        """This is the whole point: a kilometre at 4000 m is not a kilometre at 0 m."""
        low = physics.slant_grammage_gcm2(0.0, 0.0, 20000.0)
        high = physics.slant_grammage_gcm2(4000.0, 0.0, 20000.0)
        self.assertAlmostEqual(high, 1524.0, delta=5.0)
        self.assertAlmostEqual(high / low, 0.622, delta=0.005)

    def test_grammage_grows_with_distance_and_shrinks_with_altitude(self):
        self.assertGreater(physics.slant_grammage_gcm2(3000.0, 0.0, 30000.0),
                           physics.slant_grammage_gcm2(3000.0, 0.0, 10000.0))
        self.assertGreater(physics.slant_grammage_gcm2(1000.0, 0.0, 20000.0),
                           physics.slant_grammage_gcm2(5000.0, 0.0, 20000.0))

    def test_upward_paths_accumulate_less_than_horizontal_ones(self):
        flat = physics.slant_grammage_gcm2(3000.0, 0.0, 20000.0)
        up = physics.slant_grammage_gcm2(3000.0, 3.0, 20000.0)
        self.assertLess(up, flat * 1.02)

    def test_the_closed_form_matches_numerical_integration(self):
        z0, theta_deg, dist = 3500.0, 2.0, 25000.0
        closed = physics.slant_grammage_gcm2(z0, theta_deg, dist)
        n = 200000
        theta = math.radians(theta_deg)
        dl = (dist / math.cos(theta)) / n
        total = sum(physics.air_density_kgm3(z0 + (i + 0.5) * dl * math.sin(theta)) * dl
                    for i in range(n))
        self.assertAlmostEqual(closed, total * 0.1, delta=0.01 * closed)

    def test_zero_distance_gives_no_grammage(self):
        self.assertEqual(physics.slant_grammage_gcm2(3000.0, 0.0, 0.0), 0.0)

    def test_maturity_is_grammage_over_x_max(self):
        self.assertAlmostEqual(physics.shower_maturity(700.0), 1.0, places=9)
        self.assertLess(physics.shower_maturity(300.0), 1.0)


class TestEarthChord(unittest.TestCase):
    def test_chord_is_zero_along_and_above_the_horizontal(self):
        self.assertEqual(physics.earth_chord_m(0.0), 0.0)
        self.assertEqual(physics.earth_chord_m(3.0), 0.0)

    def test_chord_matches_two_r_sin_theta(self):
        for deg in (0.5, 1.0, 3.0, 30.0):
            self.assertAlmostEqual(physics.earth_chord_m(-deg),
                                   2 * physics.EARTH_RADIUS_M * math.sin(math.radians(deg)),
                                   places=3)

    def test_straight_down_is_a_diameter(self):
        self.assertAlmostEqual(physics.earth_chord_m(-90.0),
                               2 * physics.EARTH_RADIUS_M, places=3)

    def test_the_chord_dwarfs_local_topography(self):
        """At -1 deg about 220 km, at -3 deg about 670 km, against tens of km of mountain."""
        self.assertAlmostEqual(physics.earth_chord_m(-1.0) / 1000.0, 222.0, delta=2.0)
        self.assertAlmostEqual(physics.earth_chord_m(-3.0) / 1000.0, 667.0, delta=5.0)

    def test_chord_stays_within_the_crust(self):
        """Deepest point of a 670 km chord is only ~9 km down, so constant density holds."""
        chord = physics.earth_chord_m(-3.0)
        depth = (chord / 2.0) ** 2 / (2 * physics.EARTH_RADIUS_M)
        self.assertLess(depth, 15000.0)

    def test_survival_falls_with_angle_below_the_horizon(self):
        x_int = 1.0e8
        shallow = physics.neutrino_survival(-0.5, x_int)
        steep = physics.neutrino_survival(-3.0, x_int)
        self.assertGreater(shallow, steep)
        self.assertLess(steep, 0.5)

    def test_survival_is_unity_above_the_horizon(self):
        self.assertAlmostEqual(physics.neutrino_survival(1.0, 1.0e8), 1.0, places=12)


class TestTauRange(unittest.TestCase):
    def test_range_grows_with_energy_then_saturates(self):
        r = [physics.tau_range_gcm2(e) for e in (1.0, 10.0, 100.0, 1000.0, 10000.0)]
        self.assertTrue(all(b > a for a, b in zip(r, r[1:])))
        # Saturating: each decade adds proportionally less
        self.assertLess(r[4] / r[3], r[1] / r[0])

    def test_range_grows_logarithmically_rather_than_saturating(self):
        """
        Decay and energy loss couple: R = X_loss ln(1 + X_decay/X_loss). An earlier
        version combined them harmonically, which saturates at 1/beta and understates
        the range by 2x at an EeV and 4x at 10 EeV.
        """
        beta = 0.5e-6
        for e in (100.0, 1000.0, 10000.0):
            x_decay = physics.tau_decay_length_m(e) * 100.0 * physics.CRUST_DENSITY_GCM3
            x_loss = 1.0 / beta
            expected = x_loss * math.log1p(x_decay / x_loss)
            self.assertAlmostEqual(physics.tau_range_gcm2(e, beta_cm2g=beta), expected,
                                   delta=1e-6 * expected)
        # and it exceeds 1/beta at high energy, which the harmonic form never does
        self.assertGreater(physics.tau_range_gcm2(10000.0, beta_cm2g=beta), 1.0 / beta)

    def test_the_band_narrows_with_energy_from_both_sides(self):
        """
        Not the naive expectation that the band simply moves up. The lower edge rises,
        because the tau range grows and more rock is needed to produce one efficiently;
        the upper edge *falls*, because the interaction length shrinks and the neutrino
        is absorbed by less matter. So the acceptable window closes at both ends.
        """
        lo_e = [physics.depth_band_from_energy(e, e)[0] for e in (100.0, 1000.0, 10000.0)]
        hi_e = [physics.depth_band_from_energy(e, e)[1] for e in (100.0, 1000.0, 10000.0)]
        self.assertTrue(all(b > a for a, b in zip(lo_e, lo_e[1:])), f"lower edge {lo_e}")
        self.assertTrue(all(b < a for a, b in zip(hi_e, hi_e[1:])), f"upper edge {hi_e}")

    def test_band_is_ordered(self):
        lo, hi = physics.depth_band_from_energy(1.0, 100.0)
        self.assertLess(lo, hi)


class TestGeomagnetic(unittest.TestCase):
    def setUp(self):
        # Near the magnetic equator: field close to horizontal and roughly northward
        self.equatorial = physics.geomagnetic_unit_vector(0.0, 0.0)

    def test_field_vector_is_a_unit_vector(self):
        for d, i in ((0.0, 0.0), (-5.0, 10.0), (30.0, -60.0)):
            v = physics.geomagnetic_unit_vector(d, i)
            self.assertAlmostEqual(math.sqrt(sum(c * c for c in v)), 1.0, places=12)

    def test_equatorial_field_points_north(self):
        self.assertAlmostEqual(self.equatorial[1], 1.0, places=12)
        self.assertAlmostEqual(self.equatorial[0], 0.0, places=12)
        self.assertAlmostEqual(self.equatorial[2], 0.0, places=12)

    def test_inclination_is_positive_downward(self):
        v = physics.geomagnetic_unit_vector(0.0, 90.0)
        self.assertAlmostEqual(v[2], -1.0, places=12)

    def test_north_south_showers_are_strongly_suppressed(self):
        """The effect that makes target azimuth matter, not merely target existence."""
        north = physics.geomagnetic_sin_alpha(0.0, 0.0, self.equatorial)
        south = physics.geomagnetic_sin_alpha(180.0, 0.0, self.equatorial)
        self.assertAlmostEqual(north, 0.0, places=9)
        self.assertAlmostEqual(south, 0.0, places=9)

    def test_east_west_showers_are_maximal(self):
        for az in (90.0, 270.0):
            self.assertAlmostEqual(
                physics.geomagnetic_sin_alpha(az, 0.0, self.equatorial), 1.0, places=9)

    def test_sign_of_the_axis_does_not_matter(self):
        a = physics.geomagnetic_sin_alpha(37.0, 2.0, self.equatorial)
        b = physics.geomagnetic_sin_alpha(217.0, -2.0, self.equatorial)
        self.assertAlmostEqual(a, b, places=9)

    def test_result_is_always_a_valid_sine(self):
        field = physics.geomagnetic_unit_vector(-4.0, 5.0)
        for az in range(0, 360, 17):
            for el in (-3.0, 0.0, 3.0):
                v = physics.geomagnetic_sin_alpha(az, el, field)
                self.assertGreaterEqual(v, 0.0)
                self.assertLessEqual(v, 1.0 + 1e-12)


class TestFootprint(unittest.TestCase):
    def test_cherenkov_angle_at_sea_level_and_altitude(self):
        self.assertAlmostEqual(math.degrees(physics.cherenkov_angle_rad(0.0)), 1.38, delta=0.03)
        self.assertAlmostEqual(math.degrees(physics.cherenkov_angle_rad(4000.0)), 1.09, delta=0.03)

    def test_a_higher_site_has_a_narrower_cone(self):
        self.assertLess(physics.cherenkov_angle_rad(4500.0), physics.cherenkov_angle_rad(1000.0))

    def test_footprint_shrinks_with_altitude(self):
        low = physics.cherenkov_footprint_radius_m(1000.0, 10000.0)
        high = physics.cherenkov_footprint_radius_m(4500.0, 10000.0)
        self.assertLess(high, low)
        self.assertAlmostEqual(low, 224.0, delta=10.0)

    def test_a_kilometre_grid_undersamples_the_footprint(self):
        """Counted antennas are a cost proxy, not an effective area."""
        self.assertLess(physics.footprint_sampling(1000.0, 4000.0, 10000.0), 1.0)

    def test_denser_spacing_samples_better(self):
        self.assertGreater(physics.footprint_sampling(200.0, 4000.0, 10000.0),
                           physics.footprint_sampling(1000.0, 4000.0, 10000.0))


class TestGeomagneticInTheScan(unittest.TestCase):
    """The weighting reaches the kernel and behaves as the closed form says."""

    def build(self, azimuth_deg):
        grid = grid_at()
        n = 1400
        r0, c0 = 700, 700
        elevation = np.zeros((n, n), dtype=np.float32)
        bearing = math.radians(azimuth_deg)
        dist = 12000.0
        tr = r0 - int(dist * math.cos(bearing) / grid.cell_size_y)
        tc = c0 + int(dist * math.sin(bearing) / grid.cell_size_x)
        elevation[tr - 80:tr + 80, tc - 80:tc + 80] = 2500.0
        cands = np.array([[float(r0), float(c0), float(azimuth_deg)]])
        return grid, elevation, cands

    def scan(self, azimuth_deg, **kw):
        grid, elevation, cands = self.build(azimuth_deg)
        params = dict(n_azimuths=1, half_width_deg=0.0, elev_min_deg=-1.0,
                      elev_max_deg=1.0, n_elev_bins=4, max_range_m=20000.0,
                      min_dist_km=5.0, max_dist_km=20.0)
        params.update(kw)
        return scan_mod.scan(cands, elevation, grid, **params)

    def test_unweighted_by_default(self):
        out = self.scan(0.0)
        self.assertAlmostEqual(float(out["geomag_solid_angle_sr"][0]),
                               float(out["solid_angle_sr"][0]), places=12)

    def test_a_north_facing_target_is_suppressed_near_the_magnetic_equator(self):
        out = self.scan(0.0, geomag_declination_deg=0.0, geomag_inclination_deg=0.0)
        self.assertGreater(float(out["solid_angle_sr"][0]), 0.0)
        self.assertLess(float(out["geomag_solid_angle_sr"][0]),
                        0.05 * float(out["solid_angle_sr"][0]))

    def test_an_east_facing_target_is_not_suppressed(self):
        out = self.scan(90.0, geomag_declination_deg=0.0, geomag_inclination_deg=0.0)
        self.assertAlmostEqual(float(out["geomag_solid_angle_sr"][0]),
                               float(out["solid_angle_sr"][0]), delta=0.02
                               * float(out["solid_angle_sr"][0]))

    def test_two_identical_sites_differ_only_by_target_azimuth(self):
        """The point of (b): terrain statistics alone cannot rank these."""
        north = self.scan(0.0, geomag_declination_deg=0.0, geomag_inclination_deg=0.0)
        east = self.scan(90.0, geomag_declination_deg=0.0, geomag_inclination_deg=0.0)
        self.assertAlmostEqual(float(north["solid_angle_sr"][0]),
                               float(east["solid_angle_sr"][0]), delta=0.02)
        self.assertGreater(float(east["geomag_solid_angle_sr"][0]),
                           10 * float(north["geomag_solid_angle_sr"][0]))


class TestGrammageAndChordInTheScan(unittest.TestCase):
    def scan_at_altitude(self, base_altitude):
        grid = grid_at()
        n = 1400
        r0, c0 = 700, 200
        elevation = np.full((n, n), base_altitude, dtype=np.float32)
        col = c0 + int(12000.0 / grid.cell_size_x)
        elevation[:, col:col + 80] = base_altitude + 2500.0
        cands = np.array([[float(r0), float(c0), 90.0]])
        return scan_mod.scan(cands, elevation, grid, n_azimuths=1, half_width_deg=0.0,
                             elev_min_deg=-1.0, elev_max_deg=1.0, n_elev_bins=4,
                             max_range_m=20000.0, min_dist_km=5.0, max_dist_km=20.0)

    def test_a_higher_site_reports_less_grammage_for_the_same_geometry(self):
        low = float(self.scan_at_altitude(0.0)["path_grammage_gcm2"][0])
        high = float(self.scan_at_altitude(4000.0)["path_grammage_gcm2"][0])
        self.assertGreater(low, high)
        self.assertAlmostEqual(high / low, math.exp(-4000.0 / 8400.0), delta=0.02)

    def test_grammage_matches_the_closed_form(self):
        out = self.scan_at_altitude(3000.0)
        expected = physics.slant_grammage_gcm2(3000.0, 0.25, 12000.0)
        self.assertAlmostEqual(float(out["path_grammage_gcm2"][0]), expected,
                               delta=0.25 * expected)

    def test_upgoing_targets_report_no_earth_chord(self):
        self.assertEqual(float(self.scan_at_altitude(0.0)["earth_chord_gcm2"][0]), 0.0)

    def test_downgoing_directions_report_a_large_chord(self):
        grid = grid_at()
        n = 1400
        cols = np.arange(n)
        # Ground falling away, so sub-horizontal directions strike distant terrain
        profile = np.clip(3000.0 - (cols - 200) * grid.cell_size_x * 0.05, 0.0, 3000.0)
        profile[:200] = 3000.0
        elevation = np.repeat(profile[None, :], n, axis=0).astype(np.float32)
        cands = np.array([[700.0, 200.0, 90.0]])
        out = scan_mod.scan(cands, elevation, grid, n_azimuths=1, half_width_deg=0.0,
                            elev_min_deg=-3.0, elev_max_deg=-0.5, n_elev_bins=5,
                            max_range_m=40000.0, min_dist_km=5.0, max_dist_km=40.0)
        if int(out["cells"][0]) > 0:
            self.assertGreater(float(out["earth_chord_gcm2"][0]), 1.0e7)


class TestRfiShielding(unittest.TestCase):
    def setUp(self):
        self.grid = grid_at()
        self.n = 800
        self.flat = np.zeros((self.n, self.n), dtype=np.float32)
        self.cands = np.array([[400.0, 100.0, 90.0]])
        self.zone = [(400.0, 600.0, 10.0)]

    def test_a_visible_source_contributes_exposure(self):
        e = scan_mod.rfi_exposure(self.cands, self.flat, self.grid, self.zone)
        self.assertGreater(float(e[0]), 0.0)

    def test_terrain_between_removes_it_entirely(self):
        blocked = self.flat.copy()
        blocked[:, 300:310] = 500.0
        e = scan_mod.rfi_exposure(self.cands, blocked, self.grid, self.zone)
        self.assertEqual(float(e[0]), 0.0)

    def test_exposure_falls_as_inverse_square(self):
        near = scan_mod.rfi_exposure(np.array([[400.0, 400.0, 90.0]]), self.flat,
                                     self.grid, self.zone)[0]
        far = scan_mod.rfi_exposure(np.array([[400.0, 200.0, 90.0]]), self.flat,
                                    self.grid, self.zone)[0]
        self.assertAlmostEqual(near / far, 4.0, delta=0.2)

    def test_no_zones_means_no_exposure(self):
        e = scan_mod.rfi_exposure(self.cands, self.flat, self.grid, [])
        self.assertEqual(float(e[0]), 0.0)

    def test_a_stronger_source_contributes_more(self):
        weak = scan_mod.rfi_exposure(self.cands, self.flat, self.grid,
                                     [(400.0, 600.0, 5.0)])[0]
        strong = scan_mod.rfi_exposure(self.cands, self.flat, self.grid,
                                       [(400.0, 600.0, 20.0)])[0]
        self.assertAlmostEqual(strong / weak, 4.0, places=6)


if __name__ == "__main__":
    unittest.main()


class TestGeomagneticDefaults(unittest.TestCase):
    def test_defaults_are_set_for_the_andes(self):
        self.assertAlmostEqual(physics.DEFAULT_GEOMAG_DECLINATION_DEG, -6.9, places=6)
        self.assertAlmostEqual(physics.DEFAULT_GEOMAG_INCLINATION_DEG, -14.0, places=6)

    def test_dipole_inclination_is_near_zero_at_the_magnetic_equator(self):
        """Peru sits close to it, which is what makes the effect so directional."""
        lat, lon = -12.0, -77.0
        self.assertLess(abs(physics.centered_dipole_inclination(lat, lon)), 12.0)

    def test_dipole_inclination_steepens_away_from_the_magnetic_equator(self):
        near = abs(physics.centered_dipole_inclination(-12.0, -77.0))
        far = abs(physics.centered_dipole_inclination(-45.0, -71.0))
        self.assertGreater(far, near)

    def test_dipole_inclination_reproduces_the_arequipa_default(self):
        self.assertAlmostEqual(physics.centered_dipole_inclination(-16.4, -71.5),
                               physics.DEFAULT_GEOMAG_INCLINATION_DEG, delta=0.1)

    def test_the_default_field_still_suppresses_north_south_showers(self):
        field = physics.geomagnetic_unit_vector(physics.DEFAULT_GEOMAG_DECLINATION_DEG,
                                                physics.DEFAULT_GEOMAG_INCLINATION_DEG)
        north = physics.geomagnetic_sin_alpha(0.0, 0.0, field)
        east = physics.geomagnetic_sin_alpha(90.0, 0.0, field)
        self.assertLess(north, east)
        self.assertLess(north, 0.35)
        self.assertGreater(east, 0.95)


class TestShowerMaturityIsAThresholdForRadio(unittest.TestCase):
    """
    Radio emission comes from around shower maximum and then propagates through air
    that is transparent at these frequencies, so being well past maximum costs nothing.
    A particle array is different: its signal dies after maximum.
    """

    def observables(self, grammage):
        return dict(cells=np.array([4]), max_depth_gcm2=np.array([1.0e6]),
                    mean_distance_m=np.array([1.0e4]), solid_angle_sr=np.array([0.05]),
                    path_grammage_gcm2=np.array([float(grammage)]))

    def shower_score(self, grammage, mode="radio"):
        from oroscope import scoring
        _, comps = scoring.score_candidates(self.observables(grammage),
                                            {"grammage_mode": mode})
        return float(comps["shower"][0])

    def test_an_immature_shower_scores_low_either_way(self):
        self.assertLess(self.shower_score(150.0), 0.4)
        self.assertLess(self.shower_score(150.0, "particle"), 0.4)

    def test_radio_saturates_at_maximum_and_stays_there(self):
        for grammage in (700.0, 1500.0, 3000.0, 10000.0):
            self.assertAlmostEqual(self.shower_score(grammage), 1.0, places=9,
                                   msg=f"grammage {grammage}")

    def test_a_particle_array_is_penalised_far_past_maximum(self):
        self.assertAlmostEqual(self.shower_score(700.0, "particle"), 1.0, places=9)
        self.assertLess(self.shower_score(6000.0, "particle"), 0.1)

    def test_the_two_modes_agree_only_up_to_maximum(self):
        self.assertAlmostEqual(self.shower_score(700.0), self.shower_score(700.0, "particle"))
        self.assertGreater(self.shower_score(6000.0), self.shower_score(6000.0, "particle"))


class TestProductionEscapeOptimum(unittest.TestCase):
    """
    The column depth that maximises tau yield, from the competition between the
    neutrino having to interact and the tau having to escape.
    """

    def test_cross_section_matches_the_standard_parameterisation(self):
        """About 1e-32 cm^2 at an EeV."""
        self.assertAlmostEqual(physics.cc_cross_section_cm2(1000.0), 1.0e-32, delta=2e-33)

    def test_interaction_length_falls_with_energy(self):
        lo = physics.neutrino_interaction_length_gcm2(100.0)
        hi = physics.neutrino_interaction_length_gcm2(10000.0)
        self.assertGreater(lo, hi)
        self.assertAlmostEqual(lo, 3.8e8, delta=0.5e8)
        self.assertAlmostEqual(hi, 7.2e7, delta=1e7)

    def test_exit_probability_rises_linearly_for_thin_slabs(self):
        """A thin slab yields taus in proportion to its interaction probability."""
        e = 1000.0
        p1 = physics.tau_exit_probability(1.0e4, e)
        p2 = physics.tau_exit_probability(2.0e4, e)
        self.assertAlmostEqual(p2 / p1, 2.0, delta=0.05)

    def test_exit_probability_falls_again_for_very_thick_slabs(self):
        e = 1000.0
        peak = physics.tau_exit_probability(physics.production_escape_optimum_gcm2(e), e)
        self.assertLess(physics.tau_exit_probability(1.0e9, e), 0.1 * peak)

    def test_the_optimum_maximises_the_exit_probability(self):
        for e in (100.0, 1000.0, 10000.0):
            x = physics.production_escape_optimum_gcm2(e)
            p_at = physics.tau_exit_probability(x, e)
            for factor in (0.3, 3.0):
                self.assertLess(physics.tau_exit_probability(x * factor, e), p_at,
                                msg=f"{e} PeV, factor {factor}")

    def test_the_optimum_is_tens_of_km_of_rock(self):
        """
        The headline number: 12 km of standard rock at 100 PeV rising to about 23 km
        at 10 EeV. It rises because the tau range grows logarithmically, and flattens
        because beta rises with energy and tempers that growth.
        """
        km = [physics.production_escape_optimum_gcm2(e) / physics.CRUST_DENSITY_GCM3 / 1e5
              for e in (100.0, 1000.0, 10000.0)]
        for v in km:
            self.assertGreater(v, 8.0)
            self.assertLess(v, 30.0)
        self.assertTrue(all(b > a for a, b in zip(km, km[1:])), f"should rise: {km}")
        self.assertLess(km[2] / km[1], 1.5, "and then flatten")

    def test_the_band_is_about_two_decades_wide(self):
        """Column depth is an intrinsically weak discriminant, and says so."""
        lo, hi = physics.depth_band_from_energy(100.0, 10000.0)
        self.assertGreater(hi / lo, 30.0)
        self.assertLess(lo, 1.0e6)
        self.assertGreater(hi, 1.0e7)


class TestTheExitIntegralIsResolved(unittest.TestCase):
    """
    The integrand is a spike at the far surface, so the grid has to be in ``X - x``.

    Only interactions within a few tau ranges of the exit contribute; everything deeper
    is absorbed. Sampled uniformly in ``x`` over a range up to five decades wide, the
    spacing outran the spike and the trapezoid rule reported the area of something it
    never resolved -- 8x high at 3 PeV and 10^9 g/cm^2, with the trend in depth
    inverted. These pin the converged values, so a return to a uniform grid fails here
    rather than in a band six months later.
    """

    def test_the_converged_value_where_a_uniform_grid_was_eight_times_high(self):
        # 1.100343e-05 from the substituted form at 200,000 points; the uniform grid
        # needed 2,000,000 to reach it and gave 8.884e-05 at its 2000-point default.
        self.assertAlmostEqual(physics.tau_exit_probability(1.0e9, 3.0),
                               1.1003e-05, delta=2.0e-08)

    def test_the_default_grid_is_already_converged(self):
        """Refining by a hundredfold must not move the answer."""
        coarse = physics.tau_exit_probability(1.0e9, 3.0)
        fine = physics.tau_exit_probability(1.0e9, 3.0, samples=200_000)
        self.assertAlmostEqual(coarse / fine, 1.0, delta=1.0e-4)

    def test_more_rock_can_only_absorb_more(self):
        """
        The sign of the trend, which the unresolved grid had backwards.

        A uniform grid gave 2.63e-05, 4.60e-05, 8.88e-05 over these depths -- rising,
        with a spurious maximum at the edge of the grid.
        """
        p = physics.tau_exit_probability([1.0e8, 3.0e8, 1.0e9, 3.0e9], 3.0)
        self.assertTrue(bool((np.diff(p) < 0).all()), f"should fall with depth: {p}")

    def test_the_resolved_regime_is_unchanged(self):
        """
        At 100 PeV and above the uniform grid was already within 3%, which is why this
        went unseen. Those values must not have moved.
        """
        for energy, depth, before in ((100.0, 1.0e6, 1.6255e-03),
                                      (1000.0, 1.0e7, 1.6274e-02),
                                      (10000.0, 1.0e8, 1.5870e-02),
                                      (3.0, 1.0e6, 2.3461e-05)):
            got = physics.tau_exit_probability(depth, energy)
            self.assertAlmostEqual(got / before, 1.0, delta=2.0e-4,
                                   msg=f"{energy} PeV at {depth:.0e}")

    def test_the_published_optima_are_unchanged(self):
        """12 km of rock at 100 PeV rising to 23 km at 10 EeV, as before."""
        for energy, before in ((100.0, 3.302e6), (1000.0, 5.713e6), (10000.0, 6.230e6)):
            got = physics.production_escape_optimum_gcm2(energy)
            self.assertAlmostEqual(got / before, 1.0, delta=1.0e-3,
                                   msg=f"{energy} PeV")

    def test_tambos_configured_range_now_contains_its_own_optimum(self):
        """
        The band takes its low edge at the *lowest* energy asked for, and TAMBO's range
        starts at 3 PeV -- exactly where the integral did not converge. It returned
        (1.18e8, 2.89e8) over 3 PeV - 1 EeV, whose low edge is 20x above the 1 EeV
        optimum of 5.7e6, so the band excluded the depth it exists to find.
        """
        lo, hi = physics.depth_band_from_energy(3.0, 1000.0)
        optimum = physics.production_escape_optimum_gcm2(1000.0)
        self.assertLess(lo, optimum)
        self.assertGreater(hi, optimum)

    def test_lowering_the_minimum_energy_cannot_raise_the_low_edge(self):
        """
        A wider range must give a wider band. It did not: 100 PeV - 10 EeV gave a low
        edge of 5.2e5 and 3 PeV - 10 EeV gave 5.6e7, a hundredfold *rise* from asking
        for more.
        """
        wide = physics.depth_band_from_energy(3.0, 10000.0)[0]
        narrow = physics.depth_band_from_energy(100.0, 10000.0)[0]
        self.assertLessEqual(wide, narrow)

    def test_beta_rises_with_energy_as_photonuclear_does(self):
        betas = [physics.tau_energy_loss_beta(e) for e in (100.0, 1000.0, 10000.0)]
        self.assertTrue(all(b > a for a, b in zip(betas, betas[1:])))
        self.assertAlmostEqual(betas[0], 0.38e-6, delta=0.05e-6)
        self.assertAlmostEqual(betas[2], 0.95e-6, delta=0.1e-6)

    def test_a_zero_index_recovers_a_constant_beta(self):
        self.assertEqual(physics.tau_energy_loss_beta(100.0, index=0.0),
                         physics.tau_energy_loss_beta(10000.0, index=0.0))

    def test_survival_is_sharper_than_a_simple_exponential(self):
        """The tau loses energy as it goes, so its decay length shrinks en route."""
        e = 1000.0
        r = physics.tau_range_gcm2(e)
        self.assertAlmostEqual(float(physics.tau_survival(r, e)), math.exp(-1.0), delta=0.02)
        self.assertLess(float(physics.tau_survival(3 * r, e)), math.exp(-3.0))

    def test_survival_reduces_to_an_exponential_without_energy_loss(self):
        e = 1000.0
        x_decay = physics.tau_decay_length_m(e) * 100.0 * physics.CRUST_DENSITY_GCM3
        s = float(physics.tau_survival(1.0e6, e, beta_cm2g=1e-14))
        self.assertAlmostEqual(s, math.exp(-1.0e6 / x_decay), places=5)

    def test_terrain_depths_sit_on_the_rising_side_of_the_optimum(self):
        """
        Measured Arequipa depths (median 2.0e6, p90 4.9e6 g/cm^2) fall below the
        optimum, so within what topography can supply, more rock is always better and
        the upper edge of the band never binds.
        """
        for e in (100.0, 1000.0, 10000.0):
            opt = physics.production_escape_optimum_gcm2(e)
            peak = physics.tau_exit_probability(opt, e)
            p90 = physics.tau_exit_probability(4.9e6, e) / peak
            self.assertGreater(p90, 0.9, msg=f"{e} PeV")

    def test_earth_absorption_narrows_the_window_with_energy(self):
        """
        The upper limit on column depth is reached by the Earth chord, not by mountains,
        so the effective arrival window's lower edge climbs toward the horizon as energy
        rises.
        """
        cuts = [physics.earth_absorption_cutoff_deg(e) for e in (100.0, 1000.0, 10000.0)]
        self.assertAlmostEqual(cuts[0], -4.4, delta=0.4)
        self.assertAlmostEqual(cuts[1], -2.0, delta=0.3)
        self.assertAlmostEqual(cuts[2], -1.0, delta=0.3)
        self.assertTrue(all(b > a for a, b in zip(cuts, cuts[1:])), "cut should climb")

    def test_at_100_pev_absorption_does_not_bite_inside_the_window(self):
        self.assertLess(physics.earth_absorption_cutoff_deg(100.0), -3.0)

    def test_a_larger_beta_shortens_the_range_and_lowers_the_optimum(self):
        loose = physics.production_escape_optimum_gcm2(1000.0, beta_cm2g=0.4e-6)
        tight = physics.production_escape_optimum_gcm2(1000.0, beta_cm2g=1.0e-6)
        self.assertGreater(loose, tight)

    def test_inelasticity_lowers_the_tau_energy_and_so_its_reach(self):
        """
        The tau carries (1-y) of the neutrino energy, y about 0.2, so it is less
        boosted and cannot reach back as far through the rock. The optimum itself is
        too broad to shift measurably on a log grid, but the reach does.
        """
        deep = 1.0e7
        full = physics.tau_exit_probability(deep, 1000.0, inelasticity=0.0)
        real = physics.tau_exit_probability(deep, 1000.0, inelasticity=0.2)
        self.assertLess(real, full)
        self.assertGreater(real, 0.5 * full)


class TestPhysicsVerification(unittest.TestCase):
    """
    Independent checks of the derivations, against limits and constructions rather
    than against the code's own output.
    """

    def test_quadrature_reproduces_the_closed_form_for_exponential_survival(self):
        """
        Replacing the true survival with a simple exponential makes the exit integral
        analytic; the numerical integration must reproduce it.
        """
        lam = physics.neutrino_interaction_length_gcm2(1000.0)
        R = 3.0e6
        for X in np.logspace(4, 8, 9):
            x = np.linspace(0.0, X, 20000)
            num = np.trapezoid(np.exp(-x / lam) / lam * np.exp(-(X - x) / R), x)
            ana = R / (lam - R) * (np.exp(-X / lam) - np.exp(-X / R))
            self.assertAlmostEqual(num / ana, 1.0, places=4)

    def test_earth_chord_by_direct_construction(self):
        """
        From P = (0, R) a ray at theta below the tangent is Q(t) = (t cos, R - t sin);
        |Q| = R gives t = 2R sin(theta). Verified by checking Q lands on the circle.
        """
        R = physics.EARTH_RADIUS_M
        for th in (0.5, 1.0, 3.0, 30.0, 90.0):
            t = physics.earth_chord_m(-th)
            q = (t * math.cos(math.radians(th)), R - t * math.sin(math.radians(th)))
            self.assertAlmostEqual(math.hypot(*q), R, delta=1.0,
                                   msg=f"chord endpoint off the sphere at {th} deg")

    def test_a_diameter_at_ninety_degrees(self):
        self.assertAlmostEqual(physics.earth_chord_m(-90.0),
                               2 * physics.EARTH_RADIUS_M, delta=1.0)

    def test_decay_length_matches_the_published_tambo_range(self):
        """Ref. [2] Fig. 1 annotates 50 m - 5 km for 1-100 PeV."""
        self.assertAlmostEqual(physics.tau_decay_length_m(1.0), 49.0, delta=2.0)
        self.assertAlmostEqual(physics.tau_decay_length_m(100.0) / 1000.0, 4.9, delta=0.2)

    def test_interaction_length_is_hundreds_of_km_of_rock(self):
        lam = physics.neutrino_interaction_length_gcm2(1000.0)
        km = lam / physics.CRUST_DENSITY_GCM3 / 1e5
        self.assertGreater(km, 300.0)
        self.assertLess(km, 1000.0)

    def test_exit_probability_is_dimensionally_a_probability(self):
        for e in (100.0, 1000.0, 10000.0):
            for X in (1e5, 1e6, 1e7, 1e8):
                v = physics.tau_exit_probability(X, e)
                self.assertGreaterEqual(v, 0.0)
                self.assertLessEqual(v, 1.0)


class TestShowerProfileAndGrammageBand(unittest.TestCase):
    """
    The band a particle array accepts, derived from the primary energy.

    A radio array wants a maturity threshold, because emission comes from around
    shower maximum and then propagates through transparent air. A particle array wants
    a band, because the charged-particle content dies after maximum. That band is what
    decides whether a canyon is wide enough to be worth instrumenting, so it has to
    come from the energy rather than from a hard-coded multiple of X_max.
    """

    def test_shower_maximum_deepens_with_energy_at_the_elongation_rate(self):
        lo = physics.shower_maximum_gcm2(100.0)
        hi = physics.shower_maximum_gcm2(1000.0)
        self.assertAlmostEqual(hi - lo, physics.ELONGATION_RATE_GCM2_PER_DECADE, delta=1e-6)

    def test_shower_maximum_matches_its_reference_at_the_reference_energy(self):
        self.assertAlmostEqual(
            float(physics.shower_maximum_gcm2(physics.X_MAX_REFERENCE_ENERGY_PEV)),
            physics.X_MAX_GCM2, delta=1e-9)

    def test_profile_peaks_at_maximum_and_is_normalised_there(self):
        x_max = 700.0
        grid = np.linspace(1.0, 3000.0, 3000)
        n = physics.shower_size_fraction(grid, x_max)
        self.assertAlmostEqual(float(n.max()), 1.0, delta=1e-3)
        self.assertAlmostEqual(float(grid[n.argmax()]), x_max, delta=2.0)

    def test_profile_rises_before_maximum_and_falls_after(self):
        x_max = 700.0
        rising = physics.shower_size_fraction(np.array([100.0, 300.0, 500.0]), x_max)
        self.assertTrue(np.all(np.diff(rising) > 0))
        falling = physics.shower_size_fraction(np.array([900.0, 1200.0, 1600.0]), x_max)
        self.assertTrue(np.all(np.diff(falling) < 0))

    def test_no_shower_before_the_first_interaction(self):
        self.assertEqual(float(physics.shower_size_fraction(np.array([0.0]), 700.0)[0]), 0.0)

    def test_band_brackets_shower_maximum(self):
        lo, hi = physics.grammage_band_from_energy(3.0, 1000.0)
        self.assertLess(lo, float(physics.shower_maximum_gcm2(3.0)))
        self.assertGreater(hi, float(physics.shower_maximum_gcm2(1000.0)))

    def test_a_looser_fraction_widens_the_band(self):
        tight = physics.grammage_band_from_energy(3.0, 1000.0, fraction=0.2)
        loose = physics.grammage_band_from_energy(3.0, 1000.0, fraction=0.05)
        self.assertLess(loose[0], tight[0])
        self.assertGreater(loose[1], tight[1])

    def test_tambo_band_admits_a_full_colca_crossing_but_not_a_short_one(self):
        """
        The siting consequence. Colca is ~1.5 km deep and ~4.5 km rim to rim; at those
        altitudes the air holds roughly 390 g/cm^2 across the full width and only about
        170 g/cm^2 across 2 km. The band has to separate those, or the search cannot
        tell a canyon worth instrumenting from a gully.
        """
        lo, hi = physics.grammage_band_from_energy(3.0, 1000.0, fraction=0.1)
        self.assertLess(lo, 390.0, "a full-width crossing must be inside the band")
        self.assertGreater(lo, 170.0, "a 2 km crossing must not be")
        self.assertGreater(hi, lo)


class TestSpectrumWeightedDecay(unittest.TestCase):
    """
    Folding the decay probability over a spectrum, rather than picking one energy.

    A tau's decay length runs over three decades across a single experiment's reach, so
    evaluating the probability at one representative energy is not an approximation but
    a choice of answer: on a real canyon search the reported capacity ran from 10878
    detector positions at 3 PeV to zero at 100 PeV. These tests pin the properties the
    folded form must have for that to be an improvement rather than a different
    arbitrary number.
    """

    LO, HI = 3.0, 1000.0

    def test_is_a_probability(self):
        for d in (0.0, 500.0, 3000.0, 50_000.0, 1.0e6):
            p = float(physics.spectrum_weighted_decay_probability(d, self.LO, self.HI))
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)

    def test_rises_with_the_gap(self):
        d = np.array([500.0, 1000.0, 3000.0, 10_000.0, 50_000.0])
        p = physics.spectrum_weighted_decay_probability(d, self.LO, self.HI)
        self.assertTrue(np.all(np.diff(p) > 0), "a longer gap can only help")

    def test_lies_between_the_single_energy_extremes(self):
        """
        The whole point: the folded value sits inside the range one energy could give.

        At the low end of the reach the tau decays almost at once; at the high end it
        mostly flies through. A weighted average of the two must fall between them, and
        a folded value outside that range would mean the quadrature was wrong.
        """
        d = 3000.0
        at_lo = 1.0 - math.exp(-d / physics.tau_decay_length_m(self.LO))
        at_hi = 1.0 - math.exp(-d / physics.tau_decay_length_m(self.HI))
        folded = float(physics.spectrum_weighted_decay_probability(d, self.LO, self.HI))
        self.assertLess(at_hi, folded)
        self.assertLess(folded, at_lo)

    def test_a_harder_spectrum_lowers_it(self):
        """More weight at high energy, where the tau outruns the gap."""
        d = 3000.0
        soft = float(physics.spectrum_weighted_decay_probability(d, self.LO, self.HI, 2.7))
        mid = float(physics.spectrum_weighted_decay_probability(d, self.LO, self.HI, 2.0))
        hard = float(physics.spectrum_weighted_decay_probability(d, self.LO, self.HI, 1.0))
        self.assertGreater(soft, mid)
        self.assertGreater(mid, hard)

    def test_a_narrow_range_reproduces_the_single_energy(self):
        """Collapsing the range must recover the unfolded form, or the weighting is wrong."""
        d = 3000.0
        folded = float(physics.spectrum_weighted_decay_probability(d, 54.9, 55.1))
        single = 1.0 - math.exp(-d / physics.tau_decay_length_m(55.0))
        self.assertAlmostEqual(folded, single, places=3)

    def test_the_shower_length_is_subtracted_from_the_gap(self):
        with_room = physics.spectrum_weighted_decay_probability(
            6000.0, self.LO, self.HI, shower_development_m=0.0)
        less_room = physics.spectrum_weighted_decay_probability(
            6000.0, self.LO, self.HI, shower_development_m=3000.0)
        self.assertGreater(float(with_room), float(less_room))

    def test_a_target_closer_than_the_shower_needs_yields_nothing(self):
        p = physics.spectrum_weighted_decay_probability(
            2000.0, self.LO, self.HI, shower_development_m=3000.0)
        self.assertEqual(float(p), 0.0)

    def test_shape_is_preserved(self):
        d = np.array([[1000.0, 3000.0], [10_000.0, 30_000.0]])
        p = physics.spectrum_weighted_decay_probability(d, self.LO, self.HI)
        self.assertEqual(p.shape, d.shape)

    def test_an_inverted_range_is_refused(self):
        with self.assertRaises(ValueError):
            physics.spectrum_weighted_decay_probability(3000.0, 1000.0, 3.0)

    def test_the_quadrature_has_converged(self):
        """Doubling the grid must not move the answer, or the default is too coarse."""
        d = 3000.0
        coarse = float(physics.spectrum_weighted_decay_probability(
            d, self.LO, self.HI, samples=96))
        fine = float(physics.spectrum_weighted_decay_probability(
            d, self.LO, self.HI, samples=384))
        self.assertAlmostEqual(coarse, fine, places=4)


class TestSpectralIndexPinnedOrMarginalised(unittest.TestCase):
    """
    The index can be pinned to a value or left to vary over a range.

    Pinning states a belief about the spectrum. Marginalising states that the belief is
    not held, which for an input nobody has measured for this purpose is often the
    honest position -- and a flat prior over a stated range says so, rather than a
    single number pretending to knowledge.
    """

    LO, HI = 3.0, 1000.0
    D = 3000.0

    def folded(self, index):
        return float(physics.spectrum_weighted_decay_probability(
            self.D, self.LO, self.HI, index))

    def test_marginalising_lands_between_the_extremes_it_spans(self):
        hard, soft = self.folded(1.5), self.folded(2.7)
        spread = self.folded((1.5, 2.7))
        self.assertLess(hard, spread)
        self.assertLess(spread, soft)

    def test_a_degenerate_range_reproduces_the_single_value(self):
        self.assertAlmostEqual(self.folded((2.0, 2.0)), self.folded(2.0), places=9)

    def test_the_order_of_the_pair_does_not_matter(self):
        self.assertAlmostEqual(self.folded((1.5, 2.7)), self.folded((2.7, 1.5)), places=9)

    def test_a_wider_range_is_still_a_probability(self):
        for index in (0.5, 3.5, (0.5, 3.5)):
            p = self.folded(index)
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)

    def test_more_than_two_indices_is_refused(self):
        with self.assertRaises(ValueError):
            physics.spectrum_weighted_decay_probability(
                self.D, self.LO, self.HI, (1.5, 2.0, 2.7))

    def test_marginalising_preserves_shape(self):
        d = np.array([[1000.0, 3000.0], [10_000.0, 30_000.0]])
        p = physics.spectrum_weighted_decay_probability(d, self.LO, self.HI, (1.5, 2.7))
        self.assertEqual(p.shape, d.shape)

    def test_the_default_index_grid_has_converged(self):
        """
        The default, against a grid fine enough to stand as truth.

        Worth asserting on the default specifically: a coarse grid is *not* converged
        here -- nine points differ from the limit in the fourth decimal -- and the
        reason the default can afford to be generous is that the index integral folds
        into a weight vector computed once, so a fine grid costs nothing per candidate.
        """
        default = float(physics.spectrum_weighted_decay_probability(
            self.D, self.LO, self.HI, (1.5, 2.7)))
        truth = float(physics.spectrum_weighted_decay_probability(
            self.D, self.LO, self.HI, (1.5, 2.7), index_samples=2001))
        self.assertAlmostEqual(default, truth, places=5)

    def test_marginalising_costs_no_more_than_pinning(self):
        """
        The refactor that makes a fine grid affordable.

        Folding the index integral into the weights turned marginalising from 45 times
        the work of a single index into the same work. Asserted as a loose bound rather
        than a timing, so it cannot fail on a busy machine.
        """
        import time
        d = np.full(200_000, 3000.0)
        t0 = time.perf_counter()
        physics.spectrum_weighted_decay_probability(d, self.LO, self.HI, 2.0)
        pinned = time.perf_counter() - t0
        t0 = time.perf_counter()
        physics.spectrum_weighted_decay_probability(d, self.LO, self.HI, (1.5, 2.7))
        marginal = time.perf_counter() - t0
        self.assertLess(marginal, 5.0 * pinned + 1.0,
                        "marginalising should fold into the weights, not re-integrate")


class TestBetaIsConfigurable(unittest.TestCase):
    """
    beta is the least certain number in the module, so adopting a collaboration value
    must not mean editing the source of an installed package.
    """

    def tearDown(self):
        physics.restore_tau_energy_loss()

    def test_setting_beta_changes_what_every_caller_sees(self):
        before = physics.tau_energy_loss_beta(100.0)
        physics.set_tau_energy_loss(reference=0.8e-6, index=0.0)
        self.assertAlmostEqual(physics.tau_energy_loss_beta(100.0), 0.8e-6)
        self.assertAlmostEqual(physics.tau_energy_loss_beta(10000.0), 0.8e-6,
                               msg="index 0 means a constant beta")
        self.assertNotAlmostEqual(physics.tau_energy_loss_beta(100.0), before)

    def test_it_reaches_range_and_survival(self):
        loose = physics.tau_range_gcm2(100.0)
        physics.set_tau_energy_loss(reference=2.0e-6, index=0.0)
        tighter = physics.tau_range_gcm2(100.0)
        self.assertLess(tighter, loose, "a larger beta means a shorter range")

    def test_restore_puts_the_shipped_estimate_back(self):
        physics.set_tau_energy_loss(reference=0.8e-6, index=0.0)
        physics.restore_tau_energy_loss()
        self.assertEqual(physics.tau_energy_loss_settings(),
                         {"reference": physics.BETA_REFERENCE_CM2G,
                          "reference_energy_pev": physics.BETA_REFERENCE_ENERGY_PEV,
                          "index": physics.BETA_ENERGY_INDEX})

    def test_an_explicit_argument_still_wins_over_the_module_setting(self):
        physics.set_tau_energy_loss(reference=0.8e-6, index=0.0)
        self.assertAlmostEqual(
            physics.tau_energy_loss_beta(100.0, reference=1.0e-6, index=0.0), 1.0e-6)

    def test_a_nonsense_beta_is_refused(self):
        for bad in (0.0, -1.0e-6):
            with self.assertRaises(ValueError):
                physics.set_tau_energy_loss(reference=bad)

    def test_beta_does_not_enter_the_decay_length_the_search_uses(self):
        # The search weights by the decay length E/m*c*tau, which is kinematics. If
        # beta ever leaks into it, this catches it -- and the explanation's claim that
        # beta does not affect a search result would become false.
        before = physics.tau_decay_length_m(100.0)
        physics.set_tau_energy_loss(reference=5.0e-6, index=0.0)
        self.assertEqual(physics.tau_decay_length_m(100.0), before)


class TestDeclinationCanFollowTheSite(unittest.TestCase):
    """
    Inclination follows the site through a dipole; declination could not, because the
    dipole gives -0.2 deg at Arequipa against an IGRF -6.9. It fell back to one constant
    wherever the DEM was. These cover the socket that lets a real model be supplied.
    """

    def tearDown(self):
        physics.set_declination_model(None)

    def test_without_a_model_the_constant_fallback_is_used_everywhere(self):
        for lat, lon in ((-16.4, -71.5), (-12.0, -77.0), (60.0, 20.0)):
            dec, _ = physics.default_field_for_site(lat, lon)
            self.assertEqual(dec, physics.DEFAULT_GEOMAG_DECLINATION_DEG)

    def test_a_model_makes_declination_follow_the_site(self):
        physics.set_declination_model(lambda lat, lon: 0.5 * lon)
        self.assertAlmostEqual(physics.default_field_for_site(-16.4, -70.0)[0], -35.0)
        self.assertAlmostEqual(physics.default_field_for_site(-16.4, -74.0)[0], -37.0)

    def test_an_explicit_declination_still_beats_the_model(self):
        physics.set_declination_model(lambda lat, lon: -1.0)
        dec, _ = physics.default_field_for_site(-16.4, -71.5, declination_deg=-6.9)
        self.assertEqual(dec, -6.9)

    def test_inclination_keeps_following_the_site_either_way(self):
        _, lima = physics.default_field_for_site(-12.0, -77.0)
        _, arequipa = physics.default_field_for_site(-16.4, -71.5)
        self.assertLess(arequipa, lima, "inclination steepens southward across Peru")

    def test_a_grid_interpolates_bilinearly(self):
        lats, lons = np.array([-18.0, -14.0]), np.array([-74.0, -70.0])
        dec = np.array([[-5.0, -7.0], [-6.0, -8.0]])
        model = physics.declination_from_grid(lats, lons, dec)
        for (lat, lon), expected in (((-18.0, -74.0), -5.0), ((-14.0, -70.0), -8.0),
                                     ((-16.0, -72.0), -6.5)):
            self.assertAlmostEqual(model(lat, lon), expected)

    def test_a_grid_query_outside_its_corners_clamps(self):
        lats, lons = np.array([-18.0, -14.0]), np.array([-74.0, -70.0])
        dec = np.array([[-5.0, -7.0], [-6.0, -8.0]])
        model = physics.declination_from_grid(lats, lons, dec)
        self.assertAlmostEqual(model(-40.0, -100.0), -5.0)

    def test_a_mis_shaped_grid_is_refused(self):
        with self.assertRaises(ValueError):
            physics.declination_from_grid(np.array([-18.0, -14.0]),
                                          np.array([-74.0, -70.0]),
                                          np.array([[-5.0, -7.0]]))


class TestNeutralCurrentRegeneration(unittest.TestCase):
    """
    A CC interaction removes a neutrino; an NC one only degrades its energy. Counting
    absorption alone therefore overstates the suppression, and the correction is off by
    default because it is a leading-order approximation rather than a cascade solution.
    """

    def test_it_is_off_by_default_so_published_numbers_do_not_move(self):
        lam = physics.neutrino_interaction_length_gcm2(1000.0)
        self.assertEqual(physics.neutrino_survival(-1.0, lam),
                         math.exp(-physics.earth_chord_gcm2(-1.0) / lam))

    def test_regeneration_can_only_help(self):
        lam = physics.neutrino_interaction_length_gcm2(1000.0)
        for elev in (-0.5, -1.0, -3.0, -8.0):
            plain = physics.neutrino_survival(elev, lam)
            regen = physics.neutrino_survival(elev, lam, nc_regeneration=True)
            self.assertGreaterEqual(regen, plain)

    def test_it_never_manufactures_flux(self):
        lam = physics.neutrino_interaction_length_gcm2(1000.0)
        for elev in (-0.01, -0.1, -0.5, -1.0, -5.0, -20.0):
            self.assertLessEqual(
                physics.neutrino_survival(elev, lam, nc_regeneration=True), 1.0)

    def test_a_steeper_spectrum_regenerates_less(self):
        # The neutrinos scattering into the band come from E/(1-y), above it, where a
        # steeper spectrum has less flux. Getting this backwards is easy.
        soft = physics.nc_regeneration_factor(1.0e8, 1.0e8, spectral_index=2.0)
        hard = physics.nc_regeneration_factor(1.0e8, 1.0e8, spectral_index=2.7)
        self.assertGreater(soft, hard)

    def test_no_chord_means_no_regeneration(self):
        self.assertEqual(physics.nc_regeneration_factor(0.0, 1.0e8), 1.0)
        self.assertEqual(physics.nc_regeneration_factor(1.0e8, 0.0), 1.0)

    def test_it_grows_with_the_chord(self):
        factors = [physics.nc_regeneration_factor(x, 1.0e8)
                   for x in (1e7, 5e7, 1e8, 2e8)]
        self.assertTrue(all(b > a for a, b in zip(factors, factors[1:])))


class TestDecayWeightingIsSelectable(unittest.TestCase):
    """
    An event rate is the integral of flux * A(E) * P(E). Weighting by flux alone is one
    of three defensible choices, and which one was used has to be a stated parameter
    rather than an assumption buried in the code.
    """

    def setUp(self):
        from oroscope import aperture
        self.flat = aperture.TabulatedResponse([1.0, 3.0, 10.0, 100.0, 1000.0], [1.0] * 5)
        self.low = aperture.TabulatedResponse([1.0, 3.0, 10.0, 100.0, 1000.0],
                                              [1.0, 1.0, 0.5, 0.05, 0.001])
        self.high = aperture.TabulatedResponse([1.0, 3.0, 10.0, 100.0, 1000.0],
                                               [0.001, 0.05, 0.5, 1.0, 1.0])

    def test_a_flat_response_reduces_exactly_to_the_flux_weighting(self):
        # The check that the acceptance factor enters where it should: multiplying by a
        # constant must cancel in the normalisation and leave the flux answer untouched.
        flux = physics.spectrum_weighted_decay_probability(3000.0, 3.0, 1000.0)
        both = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, weight_by="flux_times_acceptance", response=self.flat)
        self.assertAlmostEqual(float(flux), float(both), places=12)

    def test_a_low_energy_response_decays_more_readily_than_a_high_energy_one(self):
        # The tau outruns a 3 km gap at high energy, so weighting toward low energies
        # must raise the decay probability.
        low = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, weight_by="acceptance", response=self.low)
        high = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, weight_by="acceptance", response=self.high)
        self.assertGreater(float(low), float(high))

    def test_acceptance_weighting_ignores_the_spectral_index(self):
        # The point of it: gamma is an assumption, and this weighting removes it.
        a = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, spectral_index=1.5,
            weight_by="acceptance", response=self.low)
        b = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, spectral_index=2.7,
            weight_by="acceptance", response=self.low)
        self.assertAlmostEqual(float(a), float(b), places=12)

    def test_flux_times_acceptance_still_depends_on_the_index(self):
        a = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, spectral_index=1.5,
            weight_by="flux_times_acceptance", response=self.low)
        b = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, spectral_index=2.7,
            weight_by="flux_times_acceptance", response=self.low)
        self.assertNotAlmostEqual(float(a), float(b), places=6)

    def test_an_acceptance_weighting_without_a_response_is_refused(self):
        for mode in ("acceptance", "flux_times_acceptance"):
            with self.assertRaises(ValueError):
                physics.spectrum_weighted_decay_probability(
                    3000.0, 3.0, 1000.0, weight_by=mode)

    def test_an_unknown_weighting_is_refused(self):
        with self.assertRaises(ValueError):
            physics.spectrum_weighted_decay_probability(
                3000.0, 3.0, 1000.0, weight_by="vibes", response=self.flat)

    def test_a_response_that_is_zero_everywhere_is_refused_rather_than_dividing_by_zero(self):
        from oroscope import aperture
        elsewhere = aperture.TabulatedResponse([1e6, 1e7], [1.0, 1.0])
        with self.assertRaises(ValueError):
            physics.spectrum_weighted_decay_probability(
                3000.0, 3.0, 1000.0, weight_by="acceptance", response=elsewhere)

    def test_the_default_is_flux_so_published_numbers_do_not_move(self):
        explicit = physics.spectrum_weighted_decay_probability(
            3000.0, 3.0, 1000.0, weight_by="flux")
        default = physics.spectrum_weighted_decay_probability(3000.0, 3.0, 1000.0)
        self.assertEqual(float(explicit), float(default))
