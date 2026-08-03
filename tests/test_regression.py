"""
End-to-end golden-file regression.

These pin the pipeline's output on fixed inputs so that a refactor which changes
results has to say so out loud. The synthetic case always runs; the real-DEM case
skips when input/ is absent, since it is gitignored.

Regenerate after an intended change with:

    UPDATE_GOLDEN=1 python -m unittest test_regression
"""

import json
import os
import shutil
import tempfile
import unittest

import _support
from _support import GOLDEN_DIR, run_pipeline, summarize
import synthetic

ORIGIN_LAT, ORIGIN_LON = -15.6, -72.3


def compare(case, got):
    """Compares against the stored golden file, or rewrites it when updating."""
    os.makedirs(GOLDEN_DIR, exist_ok=True)
    path = os.path.join(GOLDEN_DIR, f"{case}.json")
    if _support.updating_golden():
        with open(path, "w") as f:
            json.dump(got, f, indent=2, sort_keys=True)
        return None
    if not os.path.exists(path):
        raise AssertionError(f"missing golden file {path}; regenerate with UPDATE_GOLDEN=1")
    with open(path) as f:
        return json.load(f)


class GoldenCase(unittest.TestCase):
    def assert_matches_golden(self, case, got):
        expected = compare(case, got)
        if expected is None:
            self.skipTest(f"golden file for {case} regenerated")
        self.assertEqual(got["total_sites"], expected["total_sites"], "site count changed")
        self.assertEqual(got["total_capacity"], expected["total_capacity"], "total capacity changed")
        self.assertEqual(len(got["sites"]), len(expected["sites"]), "number of sites changed")
        for i, (g, e) in enumerate(zip(got["sites"], expected["sites"])):
            self.assertAlmostEqual(g["area_km2"], e["area_km2"], places=2, msg=f"site {i} area")
            self.assertEqual(g["capacity_exact"], e["capacity_exact"], f"site {i} capacity")
            self.assertEqual(g["facing_direction"], e["facing_direction"], f"site {i} facing")
        self.assertEqual(got["funnel"], expected["funnel"], "selection funnel changed")
        self.assertEqual(got["regions"], expected["regions"], "region accounting changed")


class TestSyntheticRegression(GoldenCase):
    """A deterministic ridge-and-slope map: candidates face west at the ridge."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="sitesearch_golden_")
        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(900, grid_x)
        cls.dem = synthetic.write_geotiff(os.path.join(cls.tmp, "ridge.tif"), z,
                                          ORIGIN_LAT, ORIGIN_LON)
        cls.results = run_pipeline(cls.dem, os.path.join(cls.tmp, "out"),
                                   ORIGIN_LAT, ORIGIN_LON)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_matches_golden(self):
        self.assert_matches_golden("synthetic_ridge", summarize(self.results))

    def test_finds_at_least_one_site(self):
        self.assertGreater(self.results["results"]["total_sites"], 0)

    def test_both_valley_flanks_are_found(self):
        """
        The scan finds targets on both sides of the valley: the slope east of it sees
        the ridge, and the ridge's own eastern flank sees the terrain rising beyond.

        A single ray cast along each pixel's aspect only ever caught the first of the
        two, which is the concrete form of the 2.2x under-acceptance measured in the
        review (roadmap 2.1).
        """
        facings = [s["facing_direction"] for s in self.results["results"]["sites"]]
        self.assertGreaterEqual(len(facings), 2)
        self.assertTrue(any(f in ("W", "SW", "NW") for f in facings),
                        f"expected a west-facing site, got {facings}")
        self.assertTrue(any(f in ("E", "SE", "NE") for f in facings),
                        f"expected an east-facing site, got {facings}")

    def test_the_larger_site_is_the_one_facing_the_ridge(self):
        sites = self.results["results"]["sites"]
        largest = max(sites, key=lambda s: s["area_km2"])
        self.assertIn(largest["facing_direction"], ("W", "SW", "NW"))

    def test_provenance_and_timings_are_recorded(self):
        self.assertIn("timings_sec", self.results)
        self.assertIn("morphology", self.results["timings_sec"])
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "out", "provenance.json")))

    def test_temporary_buffers_are_cleaned_up(self):
        out = os.path.join(self.tmp, "out")
        self.assertFalse(os.path.exists(os.path.join(out, "buffer_A.npy")))
        self.assertFalse(os.path.exists(os.path.join(out, "buffer_B.npy")))


class TestScanModePipeline(unittest.TestCase):
    """The arrival-scan physics mode, driven through the real entry point."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="sitesearch_scan_")
        grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
        z = synthetic.ridge_and_slope(700, grid_x)
        cls.dem = synthetic.write_geotiff(os.path.join(cls.tmp, "ridge.tif"), z,
                                          ORIGIN_LAT, ORIGIN_LON)
        cls.results = run_pipeline(
            cls.dem, os.path.join(cls.tmp, "out"), ORIGIN_LAT, ORIGIN_LON,
            n_azimuths=5, n_elev_bins=6,
            min_dist_km=2.0, max_dist_km=15.0, tile_size=350)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_scan_mode_finds_sites(self):
        self.assertGreater(self.results["results"]["total_sites"], 0)

    def test_funnel_reports_the_scan_stage(self):
        self.assertIn("directions accepted", self.results["funnel"])
        self.assertGreater(self.results["funnel"]["directions accepted"], 0)

    def test_sites_carry_their_scan_observables(self):
        """Distributions are stored per site so apertures can be folded in later."""
        site = self.results["results"]["sites"][0]
        self.assertIn("arrival_scan", site)
        obs = site["arrival_scan"]
        for key in ("solid_angle_sr_mean", "max_depth_gcm2_p50",
                    "mean_distance_m_mean", "horizon_deg_mean", "scanned_pixels"):
            self.assertIn(key, obs)
        self.assertGreater(obs["solid_angle_sr_mean"], 0.0)
        self.assertGreater(obs["max_depth_gcm2_p50"], 0.0)

    def test_recorded_distances_lie_inside_the_requested_window(self):
        obs = self.results["results"]["sites"][0]["arrival_scan"]
        self.assertGreaterEqual(obs["mean_distance_m_p50"], 2000.0)
        self.assertLessEqual(obs["mean_distance_m_p50"], 15000.0)

    def test_parameters_record_the_scan_configuration(self):
        params = self.results["parameters"]
        self.assertIsNotNone(params["scan"])
        self.assertEqual(params["scan"]["n_azimuths"], 5)


class TestEnergyDerivedWindow(unittest.TestCase):
    """An energy range replaces the hand-set distance window."""

    def test_energy_range_sets_the_distance_window(self):
        tmp = tempfile.mkdtemp(prefix="sitesearch_energy_")
        try:
            grid_x = synthetic.cell_sizes(ORIGIN_LAT)[1]
            z = synthetic.ridge_and_slope(500, grid_x)
            dem = synthetic.write_geotiff(os.path.join(tmp, "ridge.tif"), z,
                                          ORIGIN_LAT, ORIGIN_LON)
            results = run_pipeline(dem, os.path.join(tmp, "out"), ORIGIN_LAT, ORIGIN_LON,
                                   n_azimuths=3, n_elev_bins=4,
                                   energy_min_pev=1.0, energy_max_pev=100.0,
                                   tile_size=250)
            params = results["parameters"]
            # 1-100 PeV gives a 49 m - 4.9 km decay length, plus shower development
            self.assertLess(params["min_dist_km"], 0.1)
            self.assertGreater(params["max_dist_km"], 4.0)
            self.assertLess(params["max_dist_km"], 10.0)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@unittest.skipUnless(_support.have_real_dem(), "real DEM not present under input/dem/")
class TestRealDemRegression(GoldenCase):
    """A fixed crop of the Arequipa DEM, as a check against real terrain."""

    CROP = (2500, 4000, 900)   # row0, col0, size

    @classmethod
    def setUpClass(cls):
        try:
            import tifffile as tiff
        except ImportError as exc:                      # pragma: no cover
            raise unittest.SkipTest(f"tifffile unavailable: {exc}")
        cls.tmp = tempfile.mkdtemp(prefix="sitesearch_golden_real_")
        r0, c0, n = cls.CROP
        try:
            cls._load(tiff, r0, c0, n)
        except Exception as exc:
            # The shipped DEMs are LZW-compressed; without imagecodecs tifffile cannot
            # decode them. Skip rather than fail, since this case is environmental.
            shutil.rmtree(cls.tmp, ignore_errors=True)
            raise unittest.SkipTest(f"cannot read {_support.REAL_DEM}: {exc}")

    @classmethod
    def _load(cls, tiff, r0, c0, n):
        with tiff.TiffFile(_support.REAL_DEM) as tf:
            page = tf.pages[0]
            scale = page.tags["ModelPixelScaleTag"].value
            tie = page.tags["ModelTiepointTag"].value
            data = page.asarray()
        crop = data[r0:r0 + n, c0:c0 + n]
        cls.lat = tie[4] - r0 * scale[1]
        cls.lon = tie[3] + c0 * scale[0]
        cls.dem = synthetic.write_geotiff(os.path.join(cls.tmp, "crop.tif"), crop,
                                          cls.lat, cls.lon, cell_size_deg=scale[1])
        cls.results = run_pipeline(cls.dem, os.path.join(cls.tmp, "out"), cls.lat, cls.lon,
                                   min_dist_km=5.0, max_dist_km=25.0, min_sub_array_size=50)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_matches_golden(self):
        self.assert_matches_golden("arequipa_crop", summarize(self.results))


if __name__ == "__main__":
    unittest.main()
