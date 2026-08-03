"""
Score shapes, composition, and the aperture estimate's physical invariants.

The absolute normalisation of an aperture cannot be checked without a response table,
but the invariants it must satisfy regardless of normalisation can be, and are.
"""

import math
import os
import unittest

import numpy as np

import _support  # noqa: F401  (path setup)
import aperture
import arrival_scan
import scoring


class TestBandScore(unittest.TestCase):
    def test_plateau_inside_the_band(self):
        s = scoring.band_score([2.0, 5.0, 8.0], 2.0, 8.0)
        np.testing.assert_allclose(s, [1.0, 1.0, 1.0])

    def test_falls_off_outside_and_reaches_zero(self):
        s = scoring.band_score([0.0, 10.0], 2.0, 8.0, soft_lo=1.0, soft_hi=1.0)
        np.testing.assert_allclose(s, [0.0, 0.0])

    def test_flanks_are_linear(self):
        s = scoring.band_score([1.5], 2.0, 8.0, soft_lo=1.0, soft_hi=1.0)
        self.assertAlmostEqual(float(s[0]), 0.5, places=9)

    def test_default_flanks_are_a_quarter_of_the_band(self):
        s = scoring.band_score([2.0 - 0.75], 2.0, 8.0)     # quarter of 6 is 1.5
        self.assertAlmostEqual(float(s[0]), 0.5, places=9)

    def test_reversed_bounds_are_tolerated(self):
        np.testing.assert_allclose(scoring.band_score([5.0], 8.0, 2.0), [1.0])

    def test_hard_edges_when_flanks_are_zero(self):
        s = scoring.band_score([1.99, 2.0, 8.0, 8.01], 2.0, 8.0, soft_lo=0.0, soft_hi=0.0)
        np.testing.assert_allclose(s, [0.0, 1.0, 1.0, 0.0])


class TestOtherShapes(unittest.TestCase):
    def test_saturating_score_reaches_half_at_the_half_value(self):
        self.assertAlmostEqual(float(scoring.saturating_score([0.05], 0.05)[0]), 0.5)
        self.assertAlmostEqual(float(scoring.saturating_score([0.0], 0.05)[0]), 0.0)
        self.assertGreater(float(scoring.saturating_score([1.0], 0.05)[0]), 0.9)

    def test_saturating_score_is_monotonic(self):
        v = scoring.saturating_score([0.0, 0.01, 0.1, 1.0], 0.05)
        self.assertTrue(np.all(np.diff(v) > 0))

    def test_ramp_score_clamps_at_both_ends(self):
        v = scoring.ramp_score([-1.0, 0.0, 0.5, 1.0, 2.0], 0.0, 1.0)
        np.testing.assert_allclose(v, [0.0, 0.0, 0.5, 1.0, 1.0])


class TestComposition(unittest.TestCase):
    def setUp(self):
        self.components = {"a": np.array([1.0, 0.5, 0.0]), "b": np.array([1.0, 1.0, 1.0])}

    def test_product_is_unforgiving(self):
        v = scoring.compose(self.components, "product")
        np.testing.assert_allclose(v, [1.0, 0.5, 0.0], atol=1e-9)

    def test_mean_lets_a_strong_component_compensate(self):
        v = scoring.compose(self.components, "mean")
        np.testing.assert_allclose(v, [1.0, 0.75, 0.5])

    def test_min_reports_the_weakest_link(self):
        v = scoring.compose(self.components, "min")
        np.testing.assert_allclose(v, [1.0, 0.5, 0.0])

    def test_weights_act_as_exponents_for_product(self):
        v = scoring.compose({"a": np.array([0.25])}, "product", {"a": 0.5})
        self.assertAlmostEqual(float(v[0]), 0.5, places=6)

    def test_zero_weight_excludes_a_component_from_the_product(self):
        v = scoring.compose(self.components, "product", {"a": 0.0})
        np.testing.assert_allclose(v, [1.0, 1.0, 1.0], atol=1e-9)

    def test_a_zero_component_sinks_the_product(self):
        """A physical impossibility must score exactly zero, not merely small."""
        v = scoring.compose({"a": np.array([0.0]), "b": np.array([1.0])}, "product")
        self.assertEqual(float(v[0]), 0.0)

    def test_negative_weights_are_rejected(self):
        with self.assertRaises(ValueError):
            scoring.compose(self.components, "product", {"a": -1.0})

    def test_unknown_mode_is_rejected(self):
        with self.assertRaises(ValueError):
            scoring.compose(self.components, "geometric-ish")

    def test_no_components_is_rejected(self):
        with self.assertRaises(ValueError):
            scoring.compose({}, "product")

    def test_all_modes_stay_in_range(self):
        rng = np.random.default_rng(0)
        comps = {k: rng.random(200) for k in "abc"}
        for mode in scoring.COMPOSITION_MODES:
            v = scoring.compose(comps, mode)
            self.assertGreaterEqual(float(v.min()), 0.0)
            self.assertLessEqual(float(v.max()), 1.0)


class TestScoreCandidates(unittest.TestCase):
    def observables(self, **kw):
        base = dict(cells=np.array([4, 4, 0]),
                    max_depth_gcm2=np.array([1e6, 1e6, 1e6]),
                    mean_distance_m=np.array([10000.0, 10000.0, 10000.0]),
                    solid_angle_sr=np.array([0.05, 0.05, 0.05]))
        base.update(kw)
        return base

    def test_candidates_with_no_accepted_direction_score_zero(self):
        total, _ = scoring.score_candidates(self.observables())
        self.assertEqual(float(total[2]), 0.0)
        self.assertGreater(float(total[0]), 0.0)

    def test_depth_outside_the_band_lowers_the_score(self):
        inside, _ = scoring.score_candidates(self.observables())
        outside, _ = scoring.score_candidates(
            self.observables(max_depth_gcm2=np.array([1e9, 1e9, 1e9])))
        self.assertGreater(float(inside[0]), float(outside[0]))

    def test_more_solid_angle_scores_higher(self):
        low, _ = scoring.score_candidates(
            self.observables(solid_angle_sr=np.array([0.005, 0.005, 0.005])))
        high, _ = scoring.score_candidates(
            self.observables(solid_angle_sr=np.array([0.5, 0.5, 0.5])))
        self.assertGreater(float(high[0]), float(low[0]))

    def test_components_are_reported_separately(self):
        _, comps = scoring.score_candidates(self.observables())
        self.assertIn("depth", comps)
        self.assertIn("solid_angle", comps)

    def test_clearance_enters_only_when_measured(self):
        _, without = scoring.score_candidates(self.observables())
        self.assertNotIn("clearance", without)
        _, with_c = scoring.score_candidates(
            self.observables(best_clearance_ratio=np.array([0.5, 0.5, 0.5])))
        self.assertIn("clearance", with_c)

    def test_distance_band_defaults_to_the_configured_window(self):
        _, comps = scoring.score_candidates(self.observables(),
                                            distance_window_m=(5000.0, 25000.0))
        self.assertIn("distance", comps)
        self.assertAlmostEqual(float(comps["distance"][0]), 1.0)

    def test_a_sparse_config_leaves_other_defaults_alone(self):
        """None means 'use the default', so a partial config is safe to pass."""
        total, comps = scoring.score_candidates(
            self.observables(), {"depth_band_gcm2": None, "composition": "mean"})
        self.assertIn("depth", comps)
        self.assertGreater(float(total[0]), 0.0)

    def test_an_explicit_depth_band_overrides_the_default(self):
        narrow, _ = scoring.score_candidates(
            self.observables(), {"depth_band_gcm2": (1.0e9, 2.0e9)})
        self.assertEqual(float(narrow[0]), 0.0)

    def test_scores_stay_in_range(self):
        rng = np.random.default_rng(1)
        obs = dict(cells=rng.integers(0, 5, 500),
                   max_depth_gcm2=rng.random(500) * 1e8,
                   mean_distance_m=rng.random(500) * 5e4,
                   solid_angle_sr=rng.random(500))
        total, _ = scoring.score_candidates(obs)
        self.assertGreaterEqual(float(total.min()), 0.0)
        self.assertLessEqual(float(total.max()), 1.0)


class TestApertureInvariants(unittest.TestCase):
    """
    What must hold whatever the unknown normalisation is.
    """

    ENERGIES = np.logspace(0, 4, 40)

    def test_geometric_aperture_is_area_times_solid_angle(self):
        self.assertAlmostEqual(aperture.geometric_aperture_m2sr(1.0, 1.0), 1.0e6)
        self.assertAlmostEqual(aperture.geometric_aperture_m2sr(2.0, 0.5), 1.0e6)

    def test_aperture_scales_linearly_with_area(self):
        a1 = aperture.aperture_vs_energy(100.0, 0.05, 5e3, 2.5e4, self.ENERGIES)
        a2 = aperture.aperture_vs_energy(200.0, 0.05, 5e3, 2.5e4, self.ENERGIES)
        np.testing.assert_allclose(a2 / np.clip(a1, 1e-300, None), 2.0, rtol=1e-9)

    def test_aperture_scales_linearly_with_solid_angle(self):
        a1 = aperture.aperture_vs_energy(100.0, 0.02, 5e3, 2.5e4, self.ENERGIES)
        a2 = aperture.aperture_vs_energy(100.0, 0.06, 5e3, 2.5e4, self.ENERGIES)
        np.testing.assert_allclose(a2 / np.clip(a1, 1e-300, None), 3.0, rtol=1e-9)

    def test_a_short_baseline_favours_lower_energies(self):
        """Geometry alone predicts the energies a site is suited to."""
        short = aperture.peak_energy_pev(50.0, 5.0e3)        # a canyon
        long = aperture.peak_energy_pev(1.0e4, 8.0e4)        # GRAND's inherited window
        self.assertLess(short, long)

    def test_a_canyon_baseline_peaks_in_tambo_energy_range(self):
        """Colca's 4.5 km separation should favour the PeV band TAMBO targets."""
        peak = aperture.peak_energy_pev(500.0, 4.5e3)
        self.assertGreater(peak, 1.0)
        self.assertLess(peak, 300.0)

    def test_the_inherited_grand_window_peaks_near_an_eev(self):
        peak = aperture.peak_energy_pev(1.0e4, 8.0e4)
        self.assertGreater(peak, 200.0)
        self.assertLess(peak, 5000.0)

    def test_aperture_falls_away_either_side_of_the_peak(self):
        """
        Below the window the tau decays before arriving; far above it seldom decays at
        all. Neither limit is a hard zero -- the high-energy tail falls as 1/E, since
        the decay probability over a fixed window goes as (d_max - d_min)/L.
        """
        peak_e = aperture.peak_energy_pev(1.0e4, 8.0e4)
        peak = float(aperture.aperture_vs_energy(100.0, 0.05, 1.0e4, 8.0e4,
                                                 np.array([peak_e]))[0])
        low, high = aperture.aperture_vs_energy(100.0, 0.05, 1.0e4, 8.0e4,
                                                np.array([1.0e-3, 1.0e9]))
        self.assertAlmostEqual(float(low), 0.0, places=9)
        self.assertLess(float(high) / peak, 1e-4)

    def test_high_energy_tail_falls_inversely_with_energy(self):
        e = np.array([1.0e7, 1.0e8])
        curve = aperture.aperture_vs_energy(100.0, 0.05, 1.0e4, 8.0e4, e)
        self.assertAlmostEqual(curve[0] / curve[1], 10.0, delta=0.1)

    def test_decay_factor_carries_the_whole_energy_dependence(self):
        """With a unit response, the shape is exactly the decay probability."""
        curve = aperture.aperture_vs_energy(100.0, 0.05, 5e3, 2.5e4, self.ENERGIES)
        decay = np.array([arrival_scan.decay_probability(5e3, 2.5e4, float(e))
                          for e in self.ENERGIES])
        geom = aperture.geometric_aperture_m2sr(100.0, 0.05)
        np.testing.assert_allclose(curve, geom * decay, rtol=1e-12)


class TestTabulatedResponse(unittest.TestCase):
    def setUp(self):
        self.table = aperture.TabulatedResponse([1.0, 10.0, 100.0], [1.0, 10.0, 100.0])

    def test_interpolates_log_log(self):
        self.assertAlmostEqual(float(self.table(np.array([3.1622776]))[0]), 3.1622776, places=4)

    def test_reproduces_tabulated_points(self):
        np.testing.assert_allclose(self.table(np.array([1.0, 10.0, 100.0])),
                                   [1.0, 10.0, 100.0], rtol=1e-9)

    def test_returns_zero_outside_the_table_rather_than_extrapolating(self):
        np.testing.assert_allclose(self.table(np.array([0.1, 1000.0])), [0.0, 0.0])

    def test_mismatched_arrays_are_rejected(self):
        with self.assertRaises(ValueError):
            aperture.TabulatedResponse([1.0, 2.0], [1.0])

    def test_a_response_multiplies_the_curve(self):
        energies = np.array([10.0])
        plain = aperture.aperture_vs_energy(100.0, 0.05, 5e3, 2.5e4, energies)
        scaled = aperture.aperture_vs_energy(100.0, 0.05, 5e3, 2.5e4, energies,
                                             response=self.table)
        self.assertAlmostEqual(float(scaled[0]) / float(plain[0]), 10.0, places=6)


class TestSiteSummary(unittest.TestCase):
    def test_totals_sum_the_per_site_curves(self):
        sites = [
            {"site_id": 1, "area_km2": 100.0, "arrival_scan": {"solid_angle_sr_p50": 0.05}},
            {"site_id": 2, "area_km2": 50.0, "arrival_scan": {"solid_angle_sr_p50": 0.05}},
        ]
        out = aperture.summarize_sites(sites, 5e3, 2.5e4, np.logspace(0, 4, 20))
        total = np.array(out["total_m2sr"])
        summed = np.array(out["per_site_m2sr"][1]) + np.array(out["per_site_m2sr"][2])
        np.testing.assert_allclose(total, summed, rtol=1e-12)

    def test_sites_without_observables_are_skipped(self):
        self.assertEqual(aperture.summarize_sites([{"site_id": 1, "area_km2": 1.0}],
                                                  5e3, 2.5e4, [10.0]), {})


if __name__ == "__main__":
    unittest.main()


class TestDigitizedCurves(unittest.TestCase):
    """
    The published curves, hand-digitized from the figures into data/.

    These are transcriptions of vector figures, not the collaborations' tabulated
    values, so the tests check internal consistency and documented properties rather
    than exact numbers.
    """

    DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

    def test_tambo_curve_loads_and_covers_the_published_range(self):
        e, a = aperture.load_curve_csv(os.path.join(self.DATA, "tambo_aperture_fig3.csv"))
        self.assertLess(e.min(), 1.0)          # below a PeV
        self.assertGreater(e.max(), 1000.0)    # above an EeV
        self.assertTrue(np.all(np.diff(e) > 0))

    def test_tambo_aperture_rises_monotonically(self):
        """Non-decreasing: the running maximum used to clean the trace leaves plateaus."""
        _, a = aperture.load_curve_csv(os.path.join(self.DATA, "tambo_aperture_fig3.csv"))
        self.assertTrue(np.all(np.diff(a) >= 0))
        self.assertGreater(a[-1] / a[0], 100.0)

    def test_tambo_aperture_flattens_above_an_eev(self):
        """The paper states the aperture flattens above 1 EeV."""
        e, a = aperture.load_curve_csv(os.path.join(self.DATA, "tambo_aperture_fig3.csv"))
        high = a[e > 1000.0]
        self.assertGreater(len(high), 3)
        self.assertLess(high.max() / high.min(), 2.0)
        self.assertAlmostEqual(high.max(), 6.7e4, delta=2e4)

    def test_tambo_crosses_icecube_scale_near_a_few_pev(self):
        """Above ~3 PeV the published aperture exceeds IceCube's, of order 1e2-1e3."""
        e, a = aperture.load_curve_csv(os.path.join(self.DATA, "tambo_aperture_fig3.csv"))
        at_3pev = 10 ** np.interp(np.log10(3.0), np.log10(e), np.log10(a))
        self.assertGreater(at_3pev, 1.0e2)
        self.assertLess(at_3pev, 1.0e4)

    def test_grand_effective_area_loads_and_rises(self):
        e, a = aperture.load_curve_csv(
            os.path.join(self.DATA, "grand_effective_area_fig25.csv"))
        self.assertGreater(e.min(), 50.0)         # 1e8 GeV = 100 PeV
        self.assertLess(e.min(), 200.0)
        self.assertTrue(np.all(np.diff(a) >= 0))
        self.assertGreater(a[-1] / a[0], 10.0)

    def test_grand_effective_area_is_of_the_published_order(self):
        """Fig. 25 shows ~1e10 cm^2 around an EeV for GRAND10k."""
        e, a = aperture.load_curve_csv(
            os.path.join(self.DATA, "grand_effective_area_fig25.csv"))
        at_1eev = 10 ** np.interp(np.log10(1000.0), np.log10(e), np.log10(a))
        self.assertGreater(at_1eev, 1e9)
        self.assertLess(at_1eev, 1e11)


class TestInferredResponse(unittest.TestCase):
    """
    Dividing a published curve by our geometric and decay factors leaves the response.

    This is what makes an integral published curve useful without a differential
    acceptance table: whatever remains after the terrain and kinematics are divided out
    is the neutrino-interaction, tau-exit and trigger physics, and its energy shape can
    weight a site of the same experiment.
    """

    DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

    def setUp(self):
        e, a = aperture.load_curve_csv(os.path.join(self.DATA, "tambo_aperture_fig3.csv"))
        self.energy, self.response = aperture.infer_response(
            e, a, area_km2=100.0, solid_angle_sr=0.05,
            min_dist_m=300.0, max_dist_m=4500.0)

    def test_ill_conditioned_energies_are_excluded(self):
        """Where the decay probability is ~1e-8 the ratio is meaningless."""
        self.assertGreater(self.energy.min(), 0.5)

    def test_response_is_normalised_to_its_peak(self):
        self.assertAlmostEqual(float(self.response.max()), 1.0, places=9)

    def test_response_rises_with_energy_where_it_is_well_conditioned(self):
        well = self.energy > 60.0
        r = self.response[well]
        self.assertGreater(len(r), 5)
        self.assertTrue(np.all(np.diff(r) > 0), "response should rise above ~60 PeV")

    def test_the_rise_is_roughly_a_power_law(self):
        """Slope near E^1.2, consistent with a rising cross-section and tau range."""
        well = self.energy > 60.0
        slope = np.polyfit(np.log10(self.energy[well]), np.log10(self.response[well]), 1)[0]
        self.assertGreater(slope, 0.7)
        self.assertLess(slope, 2.0)

    def test_a_flat_published_curve_infers_the_inverse_decay_shape(self):
        """Sanity: dividing out a known model returns what was divided."""
        e = np.logspace(1, 4, 30)
        model = aperture.aperture_vs_energy(100.0, 0.05, 300.0, 4500.0, e)
        _, r = aperture.infer_response(e, model, 100.0, 0.05, 300.0, 4500.0)
        np.testing.assert_allclose(r, np.ones_like(r), rtol=1e-9)
