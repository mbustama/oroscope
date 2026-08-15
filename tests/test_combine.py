"""
Overlaying two runs.

One test, for one bug worth never repeating: which mask a run directory hands over.
"""

import os
import shutil
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


if __name__ == "__main__":
    unittest.main()
