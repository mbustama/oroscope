"""
Overlaying two runs: joint, union and co-location.

The arithmetic here produces headline numbers — "GRAND and TAMBO share 50 km²" — from
two rasters and a world file, and it is the kind of arithmetic that fails quietly. It
already did once: the loader took the wrong file and the report described a run that no
longer existed, with no error and a plausible number.

So the masks below are small enough to count by hand, and every expected value is
arithmetic rather than a previous run's output.
"""

import json
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
import tifffile as tiff

from _support import quiet  # noqa: F401  (also sets up sys.path)

import combine_experiments as ce
import site_searcher as ss


def write_mask(path, mask, cell_deg=1 / 3600, top_lat=-15.3, left_lon=-72.4):
    """A mask GeoTIFF and the world file beside it, as the searcher writes them."""
    tiff.imwrite(path, mask.astype(np.uint8))
    with open(os.path.splitext(path)[0] + ".tfw", "w") as f:
        f.write(f"{cell_deg:.10f}\n0.0\n0.0\n{-cell_deg:.10f}\n"
                f"{left_lon:.10f}\n{top_lat:.10f}\n")


class TestTheCurrentMaskIsTheOneLoaded(unittest.TestCase):
    """
    A directory re-run since the rename holds both ``oroscope_results_*.tif`` and a
    stale ``grand_search_results_*.tif``. The loader took the alphabetically first
    file, which is the legacy one, so it silently overlaid a superseded mask against a
    current one and reported an area for a run it was not describing.

    Measured when it was found: TAMBO at Colca came out as 44.5 km² from a stale
    48,663-pixel mask, against 83.6 km² for the run the report named. Nothing failed.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="oroscope_combine_")
        self.current = np.zeros((20, 20), dtype=np.uint8)
        self.current[:10, :] = 1                    # 200 pixels
        self.stale = np.zeros((20, 20), dtype=np.uint8)
        self.stale[:2, :] = 1                       # 40 pixels

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_the_current_prefix_wins_over_the_legacy_one(self):
        write_mask(os.path.join(self.tmp, ss.LEGACY_RESULTS_PREFIX + "colca.tif"),
                   self.stale)
        write_mask(os.path.join(self.tmp, ss.RESULTS_PREFIX + "colca.tif"),
                   self.current)
        run = ce.load_run(self.tmp)
        self.assertEqual(int(run["mask"].sum()), 200,
                         "loaded the stale legacy mask instead of the current one")

    def test_a_legacy_only_directory_still_loads(self):
        """The rename must not orphan output made before it."""
        write_mask(os.path.join(self.tmp, ss.LEGACY_RESULTS_PREFIX + "colca.tif"),
                   self.stale)
        self.assertEqual(int(ce.load_run(self.tmp)["mask"].sum()), 40)

    def test_an_unprefixed_mask_is_still_found(self):
        write_mask(os.path.join(self.tmp, "mask.tif"), self.current)
        self.assertEqual(int(ce.load_run(self.tmp)["mask"].sum()), 200)

    def test_a_directory_with_no_mask_is_refused(self):
        with self.assertRaises(SystemExit):
            ce.load_run(self.tmp)


class TestPixelArea(unittest.TestCase):
    """
    A pixel's ground area, from a world file in degrees. Everything downstream is this
    number times a pixel count, so an error here scales every reported area.
    """

    def test_a_pixel_shrinks_east_west_with_the_cosine_of_latitude(self):
        world = (1 / 3600, 0.0, 0.0, -1 / 3600, -72.4, -15.3)
        at_equator = ce.pixel_area_km2(world, 0.0)
        at_peru = ce.pixel_area_km2(world, -15.55)
        self.assertAlmostEqual(at_peru / at_equator, np.cos(np.radians(15.55)), places=6)

    def test_one_arcsecond_in_the_andes_is_about_900_square_metres(self):
        """~30.9 m north-south by ~29.8 m east-west at -15.5 degrees."""
        world = (1 / 3600, 0.0, 0.0, -1 / 3600, -72.4, -15.3)
        km2 = ce.pixel_area_km2(world, -15.5)
        self.assertAlmostEqual(km2 * 1e6, 920.0, delta=25.0)

    def test_area_scales_as_the_square_of_the_pixel_size(self):
        fine = ce.pixel_area_km2((1 / 3600, 0, 0, -1 / 3600, 0, 0), -15.5)
        coarse = ce.pixel_area_km2((4 / 3600, 0, 0, -4 / 3600, 0, 0), -15.5)
        self.assertAlmostEqual(coarse / fine, 16.0, places=6)


class TestCapacityOf(unittest.TestCase):
    """Reading a capacity out of a results file that may be absent or odd."""

    def test_it_reads_the_reported_capacity(self):
        self.assertEqual(ce.capacity_of({"results": {"total_capacity": 5317}}), 5317)

    def test_a_missing_results_file_is_not_an_error(self):
        self.assertIsNone(ce.capacity_of(None))

    def test_single_mode_reports_no_capacity_rather_than_crashing(self):
        """``search_mode: single`` writes the string 'N/A' there."""
        self.assertIsNone(ce.capacity_of({"results": {"total_capacity": "N/A"}}))


class TestAlignmentIsRefusedNotResampled(unittest.TestCase):
    """
    Two runs on differently-cropped DEMs would overlay cleanly and mean nothing.
    Comparing shapes alone is not enough, which is why the world file is compared too.
    """

    def make(self, shape=(10, 10), world=(1 / 3600, 0.0, 0.0, -1 / 3600, -72.4, -15.3)):
        return {"dir": "d", "mask": np.zeros(shape, bool), "world": world}

    def test_identical_runs_are_accepted(self):
        ce.check_alignment([self.make(), self.make()])       # must not raise

    def test_a_different_shape_is_refused(self):
        with self.assertRaises(SystemExit):
            ce.check_alignment([self.make(), self.make(shape=(10, 11))])

    def test_a_different_corner_is_refused(self):
        shifted = (1 / 3600, 0.0, 0.0, -1 / 3600, -72.1, -15.3)
        with self.assertRaises(SystemExit) as caught:
            ce.check_alignment([self.make(), self.make(world=shifted)])
        self.assertIn("same ground", str(caught.exception))

    def test_a_different_pixel_size_is_refused(self):
        coarser = (4 / 3600, 0.0, 0.0, -4 / 3600, -72.4, -15.3)
        with self.assertRaises(SystemExit):
            ce.check_alignment([self.make(), self.make(world=coarser)])

    def test_rounding_in_a_world_file_is_tolerated(self):
        """1e-9 degrees is well under a millimetre; refusing on that would be noise."""
        nudged = (1 / 3600, 0.0, 0.0, -1 / 3600, -72.4 + 1e-11, -15.3)
        ce.check_alignment([self.make(), self.make(world=nudged)])


class TestTheCombinedReport(unittest.TestCase):
    """
    End to end, on masks small enough to count.

    A: the left half of a 10x10 grid      -> 50 pixels
    B: the top half                       -> 50 pixels
    joint (both): the top-left quadrant   -> 25 pixels
    union (either)                        -> 75 pixels
    Jaccard = 25/75 = 1/3
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="oroscope_report_")
        a = np.zeros((10, 10), np.uint8)
        a[:, :5] = 1                         # left half
        b = np.zeros((10, 10), np.uint8)
        b[:5, :] = 1                         # top half

        for label, mask, capacity, sites in (("A", a, 1000, 3), ("B", b, 250, 7)):
            run_dir = os.path.join(cls.tmp, label)
            os.makedirs(run_dir)
            write_mask(os.path.join(run_dir, ss.RESULTS_PREFIX + "x.tif"), mask)
            with open(os.path.join(run_dir, ss.RESULTS_PREFIX + "x.json"), "w") as f:
                json.dump({"results": {"total_capacity": capacity,
                                       "total_sites": sites, "sites": []}}, f)

        cls.out = os.path.join(cls.tmp, "combined")
        argv = sys.argv
        sys.argv = ["combine_experiments.py",
                    os.path.join(cls.tmp, "A"), os.path.join(cls.tmp, "B"),
                    "--labels", "A", "B", "--out", cls.out, "--no_image"]
        try:
            with quiet():
                ce.main()
        finally:
            sys.argv = argv

        with open(os.path.join(cls.out, "combined_report.json")) as f:
            cls.report = json.load(f)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_the_joint_is_the_intersection(self):
        self.assertEqual(self.report["joint"]["pixels"], 25)

    def test_the_union_is_the_union(self):
        self.assertEqual(self.report["union"]["pixels"], 75)

    def test_each_run_reports_its_own_pixels(self):
        self.assertEqual([r["pixels"] for r in self.report["runs"]], [50, 50])

    def test_areas_are_pixels_times_the_pixel_area(self):
        px = self.report["pixel_area_km2"]
        self.assertAlmostEqual(self.report["joint"]["area_km2"], 25 * px, places=9)
        self.assertAlmostEqual(self.report["union"]["area_km2"], 75 * px, places=9)

    def test_the_fraction_of_each_in_the_joint(self):
        """Half of each mask is in the joint quadrant."""
        for run in self.report["runs"]:
            self.assertAlmostEqual(run["fraction_of_own_area_in_joint"], 0.5, places=9)

    def test_jaccard_is_intersection_over_union(self):
        stats = self.report["pairwise_overlap"]["A & B"]
        self.assertAlmostEqual(stats["jaccard"], 25 / 75, places=9)

    def test_the_reported_capacity_and_site_count_come_from_each_run(self):
        by_label = {r["label"]: r for r in self.report["runs"]}
        self.assertEqual(by_label["A"]["reported_capacity"], 1000)
        self.assertEqual(by_label["B"]["reported_sites"], 7)

    def test_the_masks_are_written_with_their_georeferencing(self):
        for name in ("combined_joint", "combined_union"):
            self.assertTrue(os.path.exists(os.path.join(self.out, name + ".tif")))
            self.assertTrue(os.path.exists(os.path.join(self.out, name + ".tfw")),
                            f"{name} needs a world file or it cannot be placed on a map")

    def test_the_membership_raster_encodes_one_bit_per_run(self):
        """Bit i set where run i's mask is: 1 = A only, 2 = B only, 3 = both."""
        membership = tiff.imread(os.path.join(self.out, "combined_membership.tif"))
        self.assertEqual(int((membership == 3).sum()), 25, "both")
        self.assertEqual(int((membership == 1).sum()), 25, "A alone")
        self.assertEqual(int((membership == 2).sum()), 25, "B alone")
        self.assertEqual(int((membership == 0).sum()), 25, "neither")


class TestDisjointAndIdenticalMasks(unittest.TestCase):
    """The two ends of the range, where the divisions are most likely to misbehave."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="oroscope_edge_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def report_for(self, a, b):
        for label, mask in (("A", a), ("B", b)):
            run_dir = os.path.join(self.tmp, label)
            os.makedirs(run_dir, exist_ok=True)
            write_mask(os.path.join(run_dir, ss.RESULTS_PREFIX + "x.tif"), mask)
        out = os.path.join(self.tmp, "out")
        argv = sys.argv
        sys.argv = ["combine_experiments.py",
                    os.path.join(self.tmp, "A"), os.path.join(self.tmp, "B"),
                    "--labels", "A", "B", "--out", out, "--no_image"]
        try:
            with quiet():
                ce.main()
        finally:
            sys.argv = argv
        with open(os.path.join(out, "combined_report.json")) as f:
            return json.load(f)

    def test_disjoint_masks_give_an_empty_joint_and_zero_jaccard(self):
        a = np.zeros((10, 10), np.uint8)
        a[:, :3] = 1
        b = np.zeros((10, 10), np.uint8)
        b[:, 7:] = 1
        report = self.report_for(a, b)
        self.assertEqual(report["joint"]["pixels"], 0)
        self.assertEqual(report["pairwise_overlap"]["A & B"]["jaccard"], 0.0)
        for run in report["runs"]:
            self.assertEqual(run["fraction_of_own_area_in_joint"], 0.0)

    def test_identical_masks_give_a_jaccard_of_one(self):
        a = np.zeros((10, 10), np.uint8)
        a[2:8, 2:8] = 1
        report = self.report_for(a, a.copy())
        self.assertEqual(report["pairwise_overlap"]["A & B"]["jaccard"], 1.0)
        self.assertEqual(report["joint"]["pixels"], report["union"]["pixels"])

    def test_two_empty_masks_do_not_divide_by_zero(self):
        empty = np.zeros((10, 10), np.uint8)
        report = self.report_for(empty, empty.copy())
        self.assertEqual(report["pairwise_overlap"]["A & B"]["jaccard"], 0.0)
        for run in report["runs"]:
            self.assertEqual(run["fraction_of_own_area_in_joint"], 0.0)


class TestCombineRefusesBadRequests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="oroscope_refuse_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_main(self, *args):
        argv = sys.argv
        sys.argv = ["combine_experiments.py", *args]
        try:
            with quiet():
                ce.main()
        finally:
            sys.argv = argv

    def test_one_run_is_not_a_combination(self):
        with self.assertRaises(SystemExit):
            self.run_main(self.tmp)

    def test_requiring_an_unknown_label_is_refused(self):
        for label in ("A", "B"):
            run_dir = os.path.join(self.tmp, label)
            os.makedirs(run_dir)
            write_mask(os.path.join(run_dir, ss.RESULTS_PREFIX + "x.tif"),
                       np.ones((4, 4), np.uint8))
        with self.assertRaises(SystemExit) as caught:
            self.run_main(os.path.join(self.tmp, "A"), os.path.join(self.tmp, "B"),
                          "--labels", "A", "B", "--require", "TAMBO",
                          "--out", os.path.join(self.tmp, "out"), "--no_image")
        self.assertIn("TAMBO", str(caught.exception))

    def test_a_mask_without_a_world_file_is_refused(self):
        """Without it, alignment cannot be confirmed and the overlay means nothing."""
        run_dir = os.path.join(self.tmp, "A")
        os.makedirs(run_dir)
        tiff.imwrite(os.path.join(run_dir, ss.RESULTS_PREFIX + "x.tif"),
                     np.ones((4, 4), np.uint8))
        with self.assertRaises(SystemExit) as caught:
            ce.load_run(run_dir)
        self.assertIn("world file", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
