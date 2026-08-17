"""
The run summary, and the parity that lets a library caller get one.

Two things are being pinned here. First, that :mod:`explain` says the *right* thing
about a results dictionary -- names the constraint that actually bound, attributes a
weak site to the component that weakened it, and does not quietly promise more than
the numbers support. Second, that it is reachable: on by default, suppressible, and
composable from an old results file with no DEM and no pipeline.

The summary is prose, so these assert on the claims it makes rather than on its
wording. A test that pins the sentences would fail on every edit and pass on every
lie.
"""

import json
import os
import shutil
import tempfile
import unittest

from _support import quiet, run_pipeline   # noqa: F401  (also sets up sys.path)
import synthetic

from oroscope import explain
from oroscope import site_searcher as ss

ORIGIN_LAT, ORIGIN_LON = -15.6, -72.3


def a_funnel(**stages):
    return dict(stages)


class TestBindingConstraint(unittest.TestCase):
    """
    The single most useful sentence a summary can offer, and the one a reader is most
    likely to get wrong unaided: which stage set the size of the answer.
    """

    def test_the_stage_with_the_largest_cut_is_named(self):
        funnel = {"DEM pixels": 1000, "slope 3-25 deg": 900,
                  "directions accepted": 90}
        self.assertEqual(explain.binding_constraint(funnel)["stage"],
                         "directions accepted")

    def test_a_gentle_late_stage_does_not_outrank_a_severe_early_one(self):
        funnel = {"DEM pixels": 1000, "slope 3-25 deg": 100,
                  "directions accepted": 95}
        self.assertEqual(explain.binding_constraint(funnel)["stage"], "slope 3-25 deg")

    def test_a_stage_that_leaves_nothing_wins_outright(self):
        """
        Everything downstream of an empty stage is zero for a reason that is not its
        own, so the fatal stage is the answer however severe a later ratio looks.
        """
        funnel = {"DEM pixels": 1000, "slope 3-25 deg": 800,
                  "directions accepted": 0, "after pruning (< 2.0 km wide)": 0}
        binding = explain.binding_constraint(funnel)
        self.assertEqual(binding["stage"], "directions accepted")
        self.assertTrue(binding["fatal"])

    def test_striding_is_never_the_binding_constraint(self):
        """
        It removes four candidates in five by construction, so it would otherwise be
        named on nearly every run -- and it is not a constraint: acceptance was
        measured identical at strides 1 and 5.
        """
        funnel = {"DEM pixels": 1000, "kept by stride 5": 200,
                  "directions accepted": 150}
        self.assertEqual(explain.binding_constraint(funnel)["stage"],
                         "directions accepted")

    def test_gap_closing_is_never_the_binding_constraint(self):
        """It adds pixels; a stage that grows the set cannot be what shrank it."""
        funnel = {"DEM pixels": 1000, "directions accepted": 300,
                  "after gap closing": 900}
        self.assertEqual(explain.binding_constraint(funnel)["stage"],
                         "directions accepted")

    def test_a_funnel_too_short_to_have_a_constraint_says_so(self):
        self.assertIsNone(explain.binding_constraint({}))
        self.assertIsNone(explain.binding_constraint({"DEM pixels": 10}))

    def test_the_parameter_behind_the_stage_is_named(self):
        """Naming the stage without naming its knob leaves the reader stuck."""
        funnel = {"DEM pixels": 1000, "score >= 0.35": 4}
        self.assertIn("min_score", explain.binding_constraint(funnel)["knob"])

    def test_an_unrecognised_stage_name_has_no_knob_rather_than_a_wrong_one(self):
        funnel = {"DEM pixels": 1000, "something new": 4}
        self.assertIsNone(explain.binding_constraint(funnel)["knob"])


class TestSelectedSites(unittest.TestCase):
    """
    The results file lists every site that cleared the thresholds. With
    ``stop_at_target`` that is more than were selected, and only the selection is in
    ``total_sites``, ``total_capacity`` and the exported mask — so a summary that
    totals the raw list disagrees with every other number in the file. It did:
    measured on a synthetic run, 2 sites and 243.9 km² reported against a mask holding
    1 site and 215.7.
    """

    def test_the_flag_decides(self):
        results = {"results": {"total_sites": 1, "sites": [
            {"site_id": 2, "selected": True}, {"site_id": 1, "selected": False}]}}
        chosen, rejected = explain.selected_sites(results)
        self.assertEqual([s["site_id"] for s in chosen], [2])
        self.assertEqual([s["site_id"] for s in rejected], [1])

    def test_a_file_without_the_flag_falls_back_to_the_count(self):
        """
        Exact for an older file: the list is capacity-sorted and selection walks it in
        order, so the first ``total_sites`` entries are the selected ones.
        """
        results = {"results": {"total_sites": 2, "sites": [
            {"site_id": 5}, {"site_id": 4}, {"site_id": 3}]}}
        chosen, rejected = explain.selected_sites(results)
        self.assertEqual([s["site_id"] for s in chosen], [5, 4])
        self.assertEqual([s["site_id"] for s in rejected], [3])

    def test_nothing_is_discarded_when_the_file_says_nothing(self):
        results = {"results": {"sites": [{"site_id": 1}, {"site_id": 2}]}}
        chosen, rejected = explain.selected_sites(results)
        self.assertEqual(len(chosen), 2)
        self.assertEqual(rejected, [])

    def test_a_count_larger_than_the_list_does_not_truncate(self):
        results = {"results": {"total_sites": 9, "sites": [{"site_id": 1}]}}
        self.assertEqual(len(explain.selected_sites(results)[0]), 1)

    def test_the_headline_counts_the_selection_not_the_shortlist(self):
        results = {"results": {"total_sites": 1, "total_capacity": 252, "sites": [
            {"site_id": 2, "area_km2": 215.69, "capacity_exact": 252,
             "facing_direction": "W", "selected": True},
            {"site_id": 1, "area_km2": 28.19, "capacity_exact": 36,
             "facing_direction": "E", "selected": False}]}}
        text = explain.explain_results(results)
        self.assertIn("1 site covering 215.7 km²", text)
        self.assertNotIn("243.9", text)
        # The shortlist is still worth knowing about, just not counted in.
        self.assertIn("not selected", text)


class TestSiteStrengths(unittest.TestCase):
    """
    Why a site is *good* — the mirror of attribution, and the more useful half once a
    site has been selected. "Scored 0.55" is not something a reader can act on.
    """

    RECORD = {
        "score_p50": 0.52,
        "score_solid_angle_p50": 0.90, "solid_angle_sr_p50": 1.08,
        "score_depth_p50": 1.00, "max_depth_gcm2_p50": 784440.0,
        "score_distance_p50": 1.00, "mean_distance_m_p50": 3137.0,
        "score_decay_p50": 0.30,
    }

    def test_it_reports_the_satisfied_criteria_strongest_first(self):
        names = [s["name"] for s in explain.site_strengths(self.RECORD)]
        self.assertEqual(names[:2], ["depth", "distance"])
        self.assertIn("solid_angle", names)

    def test_a_weak_component_is_not_a_strength(self):
        self.assertNotIn("decay", [s["name"] for s in explain.site_strengths(self.RECORD)])

    def test_each_strength_carries_the_measurement_that_earned_it(self):
        by_name = {s["name"]: s for s in explain.site_strengths(self.RECORD)}
        self.assertEqual(by_name["solid_angle"]["evidence"], "1.08 sr")
        self.assertIn("784,440", by_name["depth"]["evidence"])

    def test_each_strength_says_what_it_means_physically(self):
        for entry in explain.site_strengths(self.RECORD):
            self.assertTrue(entry["means"], f"{entry['name']} has no explanation")

    def test_the_threshold_can_be_moved(self):
        strict = explain.site_strengths(self.RECORD, threshold=1.0)
        self.assertEqual({s["name"] for s in strict}, {"depth", "distance"})

    def test_a_component_without_a_stored_observable_still_reports(self):
        """``decay`` has no single observable behind it; it must not be dropped."""
        entries = explain.site_strengths({"score_decay_p50": 0.95})
        self.assertEqual([e["name"] for e in entries], ["decay"])
        self.assertNotIn("evidence", entries[0])

    def test_a_record_without_components_yields_nothing(self):
        self.assertEqual(explain.site_strengths({"score_p50": 0.9}), [])


class TestConstraintOverlap(unittest.TestCase):
    """
    What decides whether two experiments can share ground. A pixel has one slope, and
    both must accept it.
    """

    GRAND = {"min_slope_deg": 3.0, "max_slope_deg": 25.0}
    TAMBO = {"min_slope_deg": 20.0, "max_slope_deg": 60.0}

    def test_it_finds_the_shared_interval(self):
        band = explain.constraint_overlap(self.GRAND, self.TAMBO)[0]
        self.assertEqual(band["overlap"], (20.0, 25.0))

    def test_the_share_is_of_the_narrower_band(self):
        """GRAND's band is 22 degrees wide, TAMBO's 40; 5/22 is the tighter squeeze."""
        band = explain.constraint_overlap(self.GRAND, self.TAMBO)[0]
        self.assertAlmostEqual(band["share_of_narrower"], 5.0 / 22.0, places=6)

    def test_disjoint_bands_report_no_overlap(self):
        band = explain.constraint_overlap({"min_slope_deg": 3.0, "max_slope_deg": 10.0},
                                          self.TAMBO)[0]
        self.assertIsNone(band["overlap"])
        self.assertEqual(band["share_of_narrower"], 0.0)

    def test_a_band_unset_on_either_side_is_not_a_constraint(self):
        bands = explain.constraint_overlap({"min_altitude": 2000.0}, {"max_altitude": 5000})
        self.assertEqual([b["label"] for b in bands], [])

    def test_the_viewing_windows_are_not_treated_as_shared_constraints(self):
        """
        Two experiments looking out from the same hillside at different ranges and
        different elevations are in no conflict whatever. Treating those windows as
        shared produced the confident and wrong conclusion that GRAND and TAMBO, which
        demonstrably share 50 km², "cannot share ground at all".
        """
        grand = dict(self.GRAND, min_dist_km=10.0, max_dist_km=40.0)
        tambo = dict(self.TAMBO, min_dist_km=2.0, max_dist_km=5.0)
        labels = [b["label"] for b in explain.constraint_overlap(grand, tambo)]
        self.assertIn("deployable slope", labels)
        self.assertNotIn("target distance", labels)


class TestExplainCombination(unittest.TestCase):
    """The overlay, explained: what each brings and what limits the sharing."""

    def setUp(self):
        self.report = {
            "runs": [
                {"label": "GRAND", "area_km2": 4580.2, "pixels": 5005057,
                 "reported_sites": 1, "reported_capacity": 5317,
                 "area_in_joint_km2": 50.1, "fraction_of_own_area_in_joint": 0.011},
                {"label": "TAMBO", "area_km2": 83.6, "pixels": 91332,
                 "reported_sites": 15, "reported_capacity": 9717,
                 "area_in_joint_km2": 50.1, "fraction_of_own_area_in_joint": 0.599},
            ],
            "joint": {"area_km2": 50.1}, "union": {"area_km2": 4613.7},
            "joint_requires": ["GRAND", "TAMBO"],
            "pairwise_overlap": {"GRAND & TAMBO": {"area_km2": 50.1, "jaccard": 0.0109,
                                                   "fraction_of_GRAND": 0.011,
                                                   "fraction_of_TAMBO": 0.599}},
        }
        self.runs = {
            "GRAND": {"parameters": {"min_slope_deg": 3.0, "max_slope_deg": 25.0,
                                     "min_dist_km": 10.0, "max_dist_km": 40.0}},
            "TAMBO": {"parameters": {"min_slope_deg": 20.0, "max_slope_deg": 60.0,
                                     "min_dist_km": 2.0, "max_dist_km": 5.0}},
        }

    def test_it_reports_what_each_experiment_brings(self):
        text = explain.explain_combination(self.report)
        self.assertIn("GRAND", text)
        self.assertIn("4,580.2", text)
        self.assertIn("9,717", text)

    def test_it_reports_the_joint_and_the_union(self):
        text = explain.explain_combination(self.report)
        self.assertIn("50.1", text)
        self.assertIn("4,613.7", text)

    def test_it_names_the_band_that_limits_the_sharing(self):
        text = explain.explain_combination(self.report, self.runs)
        self.assertIn("deployable slope", text)
        self.assertIn("20–25", text)

    def test_it_does_not_call_the_viewing_windows_a_conflict(self):
        text = explain.explain_combination(self.report, self.runs)
        self.assertNotIn("cannot share ground at all", text)
        self.assertIn("no obstacle", text)

    def test_disjoint_ground_bands_are_reported_as_such(self):
        self.runs["TAMBO"]["parameters"].update(min_slope_deg=40.0, max_slope_deg=60.0)
        text = explain.explain_combination(self.report, self.runs)
        self.assertIn("disjoint", text)

    def test_it_works_without_the_runs(self):
        """The report alone still explains its areas; only the 'why' needs parameters."""
        text = explain.explain_combination(self.report)
        self.assertIn("WHERE THESE EXPERIMENTS CAN SHARE GROUND", text)

    def test_it_warns_that_the_caveats_compound(self):
        text = explain.explain_combination(self.report)
        self.assertIn("compounds", text)

    def test_it_refuses_something_that_is_not_a_report(self):
        with self.assertRaises(TypeError):
            explain.explain_combination("combined_report.json")


class TestClosingInflation(unittest.TestCase):
    """
    The gap between accepted pixels and reported area, measured from the run rather
    than quoted from Colca. On the GRAND Colca config it gives 2.25x against the 2.35x
    a stride-1 control measured, which is the cross-check that makes it worth printing.
    """

    def test_it_corrects_for_the_stride(self):
        funnel = {"directions accepted": 100, "after gap closing": 1000}
        self.assertAlmostEqual(explain.closing_inflation(funnel, 5), 2.0)

    def test_a_small_element_against_a_coarse_stride_reports_below_one(self):
        """
        Not a bug: an element a few pixels across cannot bridge the gaps striding
        leaves, so the mask understates the accepted set instead of inflating it.
        """
        funnel = {"directions accepted": 100, "after gap closing": 300}
        self.assertLess(explain.closing_inflation(funnel, 5), 1.0)

    def test_a_funnel_without_a_closing_stage_declines(self):
        self.assertIsNone(explain.closing_inflation({"directions accepted": 100}, 5))

    def test_it_does_not_divide_by_an_empty_stage(self):
        funnel = {"directions accepted": 0, "after gap closing": 0}
        self.assertIsNone(explain.closing_inflation(funnel, 5))


class TestWeakestComponent(unittest.TestCase):
    """Attribution: which named criterion held a site back."""

    def test_the_lowest_component_is_returned_with_its_value(self):
        record = {"score_p50": 0.1, "score_decay_p50": 0.8,
                  "score_shower_p50": 0.13, "score_distance_p50": 0.9}
        self.assertEqual(explain.weakest_component(record), ("shower", 0.13))

    def test_the_total_is_not_mistaken_for_a_component(self):
        """``score_p50`` is the product, not a criterion; it is usually the lowest."""
        record = {"score_p50": 0.01, "score_decay_p50": 0.8}
        self.assertEqual(explain.weakest_component(record)[0], "decay")

    def test_a_record_without_components_declines(self):
        self.assertIsNone(explain.weakest_component({"score_p50": 0.2}))

    def test_the_statistic_can_be_chosen(self):
        record = {"score_decay_mean": 0.4, "score_shower_mean": 0.9,
                  "score_decay_p50": 0.95, "score_shower_p50": 0.2}
        self.assertEqual(explain.weakest_component(record, "mean")[0], "decay")
        self.assertEqual(explain.weakest_component(record, "p50")[0], "shower")


class TestExplainText(unittest.TestCase):
    """What the composed summary must contain, whatever words it uses."""

    def setUp(self):
        self.results = {
            "timestamp": "2026-08-15T22:43:07.123456",
            "mode": "distributed",
            "parameters": {
                "dem": "colca.tif", "origin": [-15.3, -72.4],
                "origin_source": "auto-detected from the GeoTIFF tiepoint",
                "cell_size_deg": 0.000277, "cell_size_y_m": 30.7,
                "cell_size_x_m": 29.8, "search_mode": "distributed",
                "grid_type": "hex", "spacing_km": 0.1, "target": 10000,
                "downsample_factor": 4, "gap_close_km": 1.0, "min_score": 0.35,
                "scan": {"min_target_slope_deg": 25.0},
            },
            "results": {
                "total_sites": 2, "total_capacity": 900,
                "sites": [
                    {"site_id": 7, "area_km2": 5.0, "capacity_exact": 600,
                     "grid_type": "hex", "facing_direction": "SE",
                     "arrival_scan": {"score_p50": 0.4, "score_decay_p50": 0.9,
                                      "score_shower_p50": 0.45}},
                    {"site_id": 9, "area_km2": 2.0, "capacity_exact": 300,
                     "grid_type": "hex", "facing_direction": "W",
                     "arrival_scan": {"score_p50": 0.3, "score_decay_p50": 0.8,
                                      "score_shower_p50": 0.35}},
                ],
            },
            "funnel": {"DEM pixels": 6000000, "slope 20.0-60.0 deg": 2400000,
                       "kept by stride 5": 480000, "directions accepted": 84000},
            "regions": {"labelled_regions": 900, "passed_area_threshold": 2,
                        "passed_capacity_threshold": 2, "selected": 2,
                        "required_pixels_per_region": 2700,
                        "capacity_threshold_antennas": 250},
        }
        self.text = explain.explain_results(self.results)

    def test_it_reports_what_was_found(self):
        self.assertIn("2 sites", self.text)
        self.assertIn("900", self.text)

    def test_it_names_the_binding_constraint(self):
        self.assertIn("directions accepted", self.text)

    def test_it_warns_that_reported_area_is_not_accepted_area(self):
        """
        The error a reader makes unaided, and the reason the warning defaults to on:
        morphological closing inflated the reported area 2.35x at Colca.

        Asserted against the constant rather than a literal. The literal was "2.29"
        here, which is a copy of a measurement kept in a second place: re-measuring the
        control moves the constant and leaves this test asserting the old value still
        appears, so the one test guarding the warning would have failed for the right
        reason and been "fixed" by pasting the new number in beside the old mistake.
        """
        self.assertIn(f"{explain.AREA_INFLATION_AT_COLCA:.2f}", self.text)

    def test_it_attributes_the_sites_to_a_component(self):
        self.assertIn("shower", self.text)

    def test_it_flags_min_score_as_the_dominant_assumption(self):
        self.assertIn("min_score", self.text)
        self.assertIn("score_percentile", self.text)

    def test_it_says_nothing_has_been_checked_externally(self):
        self.assertIn("external simulation", self.text)

    def test_it_carries_no_ansi_escapes(self):
        """It is meant to be pasted into an email, where colour becomes noise."""
        self.assertNotIn("\x1b[", self.text)

    def test_the_downsampling_caveat_appears_only_when_it_applies(self):
        self.assertIn("downsampled by 4", self.text)
        self.results["parameters"]["downsample_factor"] = 1
        self.assertNotIn("downsampled by 1", explain.explain_results(self.results))

    def test_provenance_is_included_when_supplied(self):
        text = explain.explain_results(self.results, {
            "git": {"commit": "34887d999342c4fa", "branch": "dev", "dirty": False},
            "dem": {"sha256": "29676fa74e05ff83cff690f158a0c1091"},
        })
        self.assertIn("34887d9", text)
        self.assertIn("29676fa74e05ff83", text)


class TestExplainSurvivesThinInput(unittest.TestCase):
    """
    It must never be the thing that fails. A summary is worth less than the search it
    describes, so a missing section is reported as missing rather than raised.
    """

    def test_an_empty_dictionary_still_explains(self):
        self.assertIn("WHAT THIS SEARCH FOUND", explain.explain_results({}))

    def test_a_search_that_found_nothing_says_so_and_still_shows_the_funnel(self):
        text = explain.explain_results({
            "results": {"total_sites": 0, "total_capacity": 0, "sites": []},
            "funnel": {"DEM pixels": 1000, "slope 3-25 deg": 900,
                       "directions accepted": 0}})
        self.assertIn("No site met all the constraints", text)
        self.assertIn("directions accepted", text)

    def test_an_older_file_without_components_says_it_cannot_attribute(self):
        text = explain.explain_results({
            "results": {"total_sites": 1, "total_capacity": 10, "sites": [
                {"site_id": 1, "area_km2": 1.0, "capacity_exact": 10,
                 "facing_direction": "N", "arrival_scan": {"score_p50": 0.5}}]}})
        self.assertIn("no per-component scores", text)

    def test_it_refuses_something_that_is_not_a_results_dictionary(self):
        with self.assertRaises(TypeError):
            explain.explain_results("the results file")


class TestExplainIsWiredIntoThePipeline(unittest.TestCase):
    """
    On by default, and reaching both the console and a file beside the results. The
    owner asked for the default explicitly: a run that has to be asked to explain
    itself is one that mostly does not.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="oroscope_explain_")
        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(700, grid_x)
        cls.dem = synthetic.write_geotiff(os.path.join(cls.tmp, "ridge.tif"), z,
                                          ORIGIN_LAT, ORIGIN_LON)
        cls.out = os.path.join(cls.tmp, "out")
        cls.results = run_pipeline(cls.dem, cls.out, ORIGIN_LAT, ORIGIN_LON)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_the_pipeline_returns_its_results(self):
        """
        It returned None, so every caller re-read the JSON it had just written.
        """
        self.assertIsInstance(self.results, dict)
        self.assertIn("funnel", self.results)
        self.assertIn("results", self.results)

    def test_what_is_returned_matches_what_was_written(self):
        with open(ss.find_results_json(self.out)) as f:
            on_disk = json.load(f)
        self.assertEqual(self.results["results"]["total_sites"],
                         on_disk["results"]["total_sites"])
        self.assertEqual(self.results["funnel"], on_disk["funnel"])

    def test_the_explanation_is_produced_without_being_asked(self):
        self.assertIn("explanation", self.results)
        self.assertIn("WHAT THIS SEARCH FOUND", self.results["explanation"])

    def test_the_explanation_is_saved_beside_the_results(self):
        path = os.path.join(self.out, "explanation.txt")
        self.assertTrue(os.path.exists(path), "explanation.txt should be written")
        with open(path, encoding="utf-8") as f:
            self.assertIn("WHERE THE CANDIDATES WENT", f.read())

    def test_the_sites_carry_their_named_score_components(self):
        """Without these the summary can report a weak site but not attribute it."""
        sites = self.results["results"]["sites"]
        self.assertTrue(sites)
        scan = sites[0]["arrival_scan"]
        components = [k for k in scan if k.startswith("score_") and k.endswith("_p50")
                      and k != "score_p50"]
        self.assertTrue(components, f"no score components in {sorted(scan)}")

    def test_it_can_be_suppressed(self):
        out = os.path.join(self.tmp, "quiet")
        results = run_pipeline(self.dem, out, ORIGIN_LAT, ORIGIN_LON, explain=False)
        self.assertFalse(os.path.exists(os.path.join(out, "explanation.txt")))
        # Still composed, so a caller that wants the text has it without re-running.
        self.assertIn("explanation", results)

    def test_the_summarised_area_is_the_area_of_the_exported_mask(self):
        """
        The invariant that was violated. Whatever the mode, the area the summary adds
        up must be the area of the raster the run actually wrote -- that raster is
        what `oroscope-combine` measures and what anyone opens in a GIS.
        """
        from oroscope import combine_experiments as ce
        for label, kw in (("plain", {}),
                          ("stop_at_target", dict(stop_at_target=True,
                                                  target_antennas=50,
                                                  min_sub_array_size=5)),
                          ("downsampled", dict(downsample_factor=3)),
                          ("single", dict(search_mode="single"))):
            with self.subTest(mode=label):
                out = os.path.join(self.tmp, "area_" + label)
                results = run_pipeline(self.dem, out, ORIGIN_LAT, ORIGIN_LON, **kw)
                run = ce.load_run(out)
                top = run["world"][5]
                centre = top + 0.5 * run["mask"].shape[0] * run["world"][3]
                mask_km2 = run["mask"].sum() * ce.pixel_area_km2(run["world"], centre)
                listed = results["results"]["sites"]
                self.assertTrue(all("selected" in s for s in listed),
                                "every site record must say whether it was selected; "
                                "without it a reader can only guess")
                chosen, _ = explain.selected_sites(results)
                summary_km2 = sum(s["area_km2"] for s in chosen)
                self.assertEqual(len(chosen), results["results"]["total_sites"],
                                 "summary counted a different number of sites")
                self.assertAlmostEqual(summary_km2, mask_km2, delta=0.01 * mask_km2,
                                       msg="summary area disagrees with the mask")

    def test_an_old_results_file_can_be_explained_without_a_dem(self):
        """The point of keeping it a pure function of the results dictionary."""
        with open(ss.find_results_json(self.out)) as f:
            on_disk = json.load(f)
        self.assertIn("WHERE THE CANDIDATES WENT", explain.explain_results(on_disk))


class TestLibraryParity(unittest.TestCase):
    """
    Everything the command line can do, the library must do too. These were the
    measured gaps.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="oroscope_parity_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_a_config_template_can_be_generated_without_the_cli(self):
        path = os.path.join(self.tmp, "nested", "config.json")
        written = ss.generate_config(path, "arequipa")
        self.assertTrue(os.path.exists(path))
        with open(path) as f:
            self.assertEqual(json.load(f), written)

    def test_every_preset_is_accepted_and_an_unknown_one_is_not(self):
        for preset in ss.CONFIG_PRESETS:
            self.assertIn("dem_path", ss.default_config(preset))
        with self.assertRaises(ValueError):
            ss.default_config("patagonia")

    def test_the_presets_differ_where_they_should(self):
        lima = ss.default_config("lima")
        arequipa = ss.default_config("arequipa")
        self.assertNotEqual(lima["origin_lat"], arequipa["origin_lat"])
        self.assertEqual(arequipa["rfi_zones"], "arequipa")

    def test_a_config_can_be_loaded_without_the_cli(self):
        path = os.path.join(self.tmp, "config.json")
        ss.generate_config(path)
        self.assertEqual(ss.load_config(path)["min_slope_deg"], 3.0)

    def test_a_missing_config_loads_as_empty_rather_than_raising(self):
        """Matching how --config_path has always treated one."""
        self.assertEqual(ss.load_config(os.path.join(self.tmp, "absent.json")), {})
        self.assertEqual(ss.load_config(None), {})

    def test_a_malformed_config_is_a_failure_rather_than_a_silent_empty(self):
        path = os.path.join(self.tmp, "broken.json")
        with open(path, "w") as f:
            f.write("{not json")
        with self.assertRaises(json.JSONDecodeError):
            ss.load_config(path)

    def test_the_generated_template_names_every_pipeline_parameter_it_claims_to(self):
        """
        A template with holes falls back silently for whatever it omits, and the
        fallbacks are the least visible input the tool has.
        """
        import inspect
        accepted = set(inspect.signature(
            ss.find_grand_regions_interactive).parameters)
        config = ss.default_config()
        # Keys that are the config file's own business: output placement, console
        # verbosity, and the negative-form spelling main() inverts.
        cli_only = {"print_info", "output_directory_base_with_given_json",
                    "require_sky"}
        unknown = set(config) - accepted - cli_only
        self.assertEqual(unknown, set(), f"template keys the pipeline ignores: {unknown}")

    def test_the_memory_preflight_is_callable_from_the_library(self):
        """
        It ran only inside main(), so the caller most likely to need it -- a sweep,
        which is what reached 6.9 GB and was killed -- did not get it.
        """
        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(200, grid_x)
        dem = synthetic.write_geotiff(os.path.join(self.tmp, "small.tif"), z,
                                      ORIGIN_LAT, ORIGIN_LON)
        with quiet():
            report = ss.preflight_memory(dem, downsample_factor=2, max_memory_gb=0)
        self.assertGreater(report["estimate_gb"], 0.0)
        self.assertFalse(report["capped"], "max_memory_gb=0 must disable the cap")

    def test_the_preflight_says_when_the_cap_cannot_protect_the_machine(self):
        """
        A cap above the available memory is not a cap: the OOM killer arrives before
        RLIMIT_AS fires. Easy to set by accident, because the estimate counts only
        anonymous memory, so on a large DEM the cap must clear it by the size of the
        memory-mapped file -- and raising it for that reason can carry it past what the
        machine has. A 339 Mpx search capped at 13.0 GiB against 8.0 GiB available took
        the machine down; nothing said a word.
        """
        import resource

        keep = resource.getrlimit(resource.RLIMIT_AS)
        try:
            with quiet():
                # Far above any real machine, so the comparison is deterministic, and
                # high enough that applying it constrains nothing.
                report = ss.preflight_memory("nonexistent.tif", max_memory_gb=1.0e6)
            self.assertTrue(report["capped"])
            self.assertTrue(report["cap_exceeds_available"],
                            "a cap of 1e6 GiB must be reported as unprotective")
        finally:
            resource.setrlimit(resource.RLIMIT_AS, keep)

        with quiet():
            report = ss.preflight_memory("nonexistent.tif", max_memory_gb=0)
        self.assertFalse(report["cap_exceeds_available"],
                         "no cap applied means nothing to warn about")

    def test_the_preflight_counts_the_map_as_well_as_the_search(self):
        """
        Three runs in one session finished their searches and then died drawing the
        picture — the JSON and the GeoTIFF already written, only the map lost. A
        pre-flight that models the search alone sizes a cap that cannot survive the run.
        """
        viz = ss.estimate_visualisation_memory_gb(3961, 2881, 1)
        self.assertGreater(viz, 0.5, "a 2.85 Mpx viz raster costs more than half a GiB")
        # It renders at downsample_factor * 2, so the cost falls as its square.
        self.assertLess(ss.estimate_visualisation_memory_gb(3961, 2881, 4), viz)

        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(200, grid_x)
        dem = synthetic.write_geotiff(os.path.join(self.tmp, "viz.tif"), z,
                                      ORIGIN_LAT, ORIGIN_LON)
        with quiet():
            report = ss.preflight_memory(dem, max_memory_gb=0)
        self.assertGreater(report["visualisation_gb"], 0.0)
        self.assertAlmostEqual(report["estimate_gb"],
                               report["search_gb"] + report["visualisation_gb"],
                               places=6, msg="the reported total must be the sum")

    def test_the_combination_is_the_larger_figure_and_is_now_counted(self):
        """
        The overlay renders at the mask's own resolution, not at ``viz_ds``, so it is
        four times the pixels. It was outside the estimate entirely, ran last, and is
        what actually kept failing: measured 1.96 GiB at the Arequipa mask's 2551x3151
        against the 0.48 GiB budgeted for a map. Composited in float32 it is 1.31 GiB,
        and the model is fitted to four measured sizes.
        """
        rows, cols = 10204, 12603
        for ds in (1, 2, 4):
            solo = ss.estimate_visualisation_memory_gb(rows, cols, ds)
            joint = ss.estimate_visualisation_memory_gb(rows, cols, ds, combine=True)
            self.assertGreater(joint, solo, msg=f"downsample_factor {ds}")

        # Against the bench: 8.04 Mpx of raster measured 1.31 GiB peak RSS.
        self.assertAlmostEqual(
            ss.estimate_visualisation_memory_gb(10204, 12603, 4, combine=True),
            1.31, delta=0.06)

    def test_the_preflight_judges_the_largest_stage_not_the_sum(self):
        """
        The search, its map and the combination happen one after another, so what has
        to fit is the biggest of them. The map is the exception -- it is drawn while the
        search's arrays are still live -- so those two are added and the combination is
        compared against that total.
        """
        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(200, grid_x)
        dem = synthetic.write_geotiff(os.path.join(self.tmp, "combine.tif"), z,
                                      ORIGIN_LAT, ORIGIN_LON)
        with quiet():
            without = ss.preflight_memory(dem, max_memory_gb=0)
            with_it = ss.preflight_memory(dem, max_memory_gb=0, combine=True)

        self.assertIsNone(without["combine_gb"], "not asked for, not reported")
        self.assertIsNotNone(with_it["combine_gb"])
        self.assertGreaterEqual(with_it["estimate_gb"], without["estimate_gb"],
                                "counting a stage cannot lower the estimate")
        self.assertEqual(with_it["estimate_gb"],
                         max(with_it["search_gb"] + with_it["visualisation_gb"],
                             with_it["combine_gb"]))

    def test_the_preflight_can_refuse_instead_of_merely_warning(self):
        """
        Warning is what it did while three runs died anyway. ``refuse=True`` is the
        gate, and it fires before anything is allocated.
        """
        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(200, grid_x)
        dem = synthetic.write_geotiff(os.path.join(self.tmp, "refuse.tif"), z,
                                      ORIGIN_LAT, ORIGIN_LON)
        real = ss.available_memory_gb

        ss.available_memory_gb = lambda: 0.001          # anything is too much
        try:
            with self.assertRaises(MemoryError):
                with quiet():
                    ss.preflight_memory(dem, max_memory_gb=0, refuse=True)
            with quiet():                                # the default still only warns
                ss.preflight_memory(dem, max_memory_gb=0)
        finally:
            ss.available_memory_gb = real

    def test_the_preflight_survives_a_dem_it_cannot_measure(self):
        with quiet():
            report = ss.preflight_memory("nonexistent.tif", max_memory_gb=0)
        self.assertIsNone(report["estimate_gb"])

    def test_the_pipeline_creates_its_own_output_directory(self):
        """
        main() created it, so the library raised FileNotFoundError from inside numpy's
        open_memmap, naming a scratch buffer rather than the directory that was missing.
        """
        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(300, grid_x)
        dem = synthetic.write_geotiff(os.path.join(self.tmp, "r.tif"), z,
                                      ORIGIN_LAT, ORIGIN_LON)
        out = os.path.join(self.tmp, "does", "not", "exist")
        with quiet():
            results = ss.find_grand_regions_interactive(
                dem_path=dem, origin_lat=ORIGIN_LAT, origin_lon=ORIGIN_LON,
                run_output_dir=out, downsample_factor=2, tile_size=256,
                num_cores=2, generate_kml=False, explain=False)
        self.assertTrue(os.path.isdir(out))
        self.assertIsInstance(results, dict)

    def test_the_pipeline_accepts_the_ceiling_as_a_parameter(self):
        """The one flag with no library equivalent at all."""
        import inspect
        params = inspect.signature(ss.find_grand_regions_interactive).parameters
        self.assertIn("max_memory_gb", params)
        self.assertIn("explain", params)
        self.assertIs(params["explain"].default, True)


if __name__ == "__main__":
    unittest.main()
