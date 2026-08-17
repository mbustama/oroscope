"""
The scripts under ``tools/``, which produce the committed results and figures.

They had no tests at all. They are not library code, but they decide two things the
library cannot check for itself: what the memory pre-flight is sized against, and what
lands in ``results/``. The first of those is the only thing standing between a run that
does not fit and the OOM killer, and it has been wrong twice.
"""

import doctest
import os
import sys
import unittest

import _support  # noqa: F401  (path setup)
from _support import ss

TOOLS_DIR = os.path.join(_support.REPO_ROOT, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

import run_full_dem  # noqa: E402


def load_tests(loader, tests, ignore):
    """
    unittest hook: runs this tool's own Examples.

    ``tests/test_doctests.py`` covers the library and deliberately excludes ``tools/``,
    whose scripts mostly take a DEM and write files. ``costliest_sampling`` is pure
    arithmetic on a pair of dictionaries, so its example is executable and -- per the
    project's rule that a doctest value is computed rather than predicted -- it is run.
    """
    tests.addTests(doctest.DocTestSuite(
        run_full_dem,
        optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS))
    return tests


class TestTheCostliestSamplingIsWhatIsPreflighted(unittest.TestCase):
    """
    One pre-flight covers several searches, so it must be sized against the dearest.

    Both knobs scale memory inversely, so the costliest configuration holds the
    *smallest* values. This was ``max`` on both, which sized the estimate against
    whichever search was cheaper and would have waved the other one through.
    """

    def test_the_smaller_of_a_mismatched_pair_wins(self):
        got = run_full_dem.costliest_sampling([
            {"downsample_factor": 4, "candidate_stride": 5},
            {"downsample_factor": 1, "candidate_stride": 1}])
        self.assertEqual(got, (1, 1))

    def test_each_knob_is_taken_independently(self):
        """A config need not be dearest on both axes to set the estimate on one."""
        got = run_full_dem.costliest_sampling([
            {"downsample_factor": 1, "candidate_stride": 15},
            {"downsample_factor": 4, "candidate_stride": 1}])
        self.assertEqual(got, (1, 1))

    def test_absent_and_null_values_default_to_one(self):
        self.assertEqual(run_full_dem.costliest_sampling([{}]), (1, 1))
        self.assertEqual(run_full_dem.costliest_sampling(
            [{"downsample_factor": None, "candidate_stride": None}]), (1, 1))
        self.assertEqual(run_full_dem.costliest_sampling([]), (1, 1))

    def test_the_choice_really_is_the_more_expensive_one(self):
        """
        The reason `min` is right, asserted against the estimator rather than assumed.

        If `estimate_peak_memory_gb` were ever to stop being monotonically decreasing in
        these two knobs, `costliest_sampling` would silently start choosing the cheaper
        configuration again.
        """
        pair = [{"downsample_factor": 4, "candidate_stride": 5},
                {"downsample_factor": 1, "candidate_stride": 1}]
        chosen = run_full_dem.costliest_sampling(pair)
        costs = {(int(c["downsample_factor"]), int(c["candidate_stride"])):
                 ss.estimate_peak_memory_gb(10204, 12603, int(c["downsample_factor"]),
                                            int(c["candidate_stride"]))
                 for c in pair}
        self.assertEqual(chosen, max(costs, key=costs.get))

    def test_every_shipped_config_pair_is_covered(self):
        """
        Whatever is in config/ today, the pre-flight must not under-read it.

        This is the property that matters, independent of whether any pair currently
        disagrees: no search may run at a sampling finer than the one estimated for.
        """
        import glob
        import json
        for region in ("arequipa", "ancash", "lima", "huaylas", "cajatambo"):
            paths = sorted(glob.glob(os.path.join(
                _support.REPO_ROOT, "config", f"*_{region}*.json")))
            configs = []
            for p in paths:
                with open(p) as f:
                    configs.append(json.load(f))
            configs = [c for c in configs if "downsample_factor" in c]
            if len(configs) < 2:
                continue
            ds, stride = run_full_dem.costliest_sampling(configs)
            for c in configs:
                self.assertLessEqual(ds, int(c.get("downsample_factor") or 1),
                                     msg=f"{region}: pre-flight downsample is coarser "
                                         f"than a run's own")
                self.assertLessEqual(stride, int(c.get("candidate_stride") or 1),
                                     msg=f"{region}: pre-flight stride is coarser "
                                         f"than a run's own")


class TestTheManifestDescribesTheStoreAndNotTheInvocation(unittest.TestCase):
    """
    A manifest that lists only what one run wrote is how stale data stays invisible.

    Both halves of this were live. ``--only tambo`` overwrote the manifest with three
    ``tambo_*`` names while nine files sat in ``results/arequipa_full/``, so the record
    said GRAND had never run. And ``results/huaylas_full/`` holds three committed
    ``*_control_*`` files that this tool never writes, which went stale at TAMBO's old
    100 m spacing with nothing naming them.
    """

    def _store(self, names):
        import tempfile
        d = tempfile.mkdtemp()
        for n in names:
            with open(os.path.join(d, n), "w") as f:
                f.write("x")
        return {"store": d, "dem": os.path.join(_support.REPO_ROOT, "dem.tif"),
                "configs": {"grand": "config/g.json", "tambo": "config/t.json"}}

    def _manifest(self, paths):
        import json
        with open(os.path.join(paths["store"], "manifest.json")) as f:
            return json.load(f)

    def test_files_is_what_the_store_holds_not_what_the_run_wrote(self):
        paths = self._store(["grand_results.json", "tambo_results.json"])
        run_full_dem.write_manifest(paths, "somewhere", ["tambo_results.json"])
        m = self._manifest(paths)
        self.assertIn("grand_results.json", m["files"],
                      "a partial run narrowed the manifest and hid the other half")
        self.assertEqual(m["written_by_this_run"], ["tambo_results.json"])

    def test_a_file_no_run_wrote_is_listed_and_dated(self):
        paths = self._store(["tambo_results.json", "tambo_control_ds4_stride5.json"])
        run_full_dem.write_manifest(paths, "somewhere", ["tambo_results.json"])
        m = self._manifest(paths)
        self.assertIn("tambo_control_ds4_stride5.json", m["also_present"])
        self.assertNotIn("tambo_results.json", m["also_present"])

    def test_a_fully_written_store_has_nothing_also_present(self):
        paths = self._store(["grand_results.json", "tambo_results.json"])
        run_full_dem.write_manifest(paths, "somewhere",
                                    ["grand_results.json", "tambo_results.json"])
        self.assertEqual(self._manifest(paths)["also_present"], {})

    def test_the_manifest_never_lists_itself(self):
        paths = self._store(["grand_results.json"])
        run_full_dem.write_manifest(paths, "somewhere", ["grand_results.json"])
        m = self._manifest(paths)
        self.assertNotIn("manifest.json", m["files"])
        self.assertNotIn("manifest.json", m["also_present"])


class TestTheStoreHoldsEveryArtefactItPromises(unittest.TestCase):
    """
    A store that is half-refreshed is worse than a stale one: it looks current.

    ``combine()`` copied ``combined_report.json`` back and left
    ``combination_explanation.txt`` behind, so a re-combine updated the numbers while
    the prose beside them kept the old ones. The region notebooks print that file
    verbatim, so notebooks 10 and 11 went on showing TAMBO's 100 m capacity after the
    150 m re-run -- a published number that no amount of re-executing would fix,
    because the input itself was stale.
    """

    # `colca` was added when the re-run gave it a store for the first time. It had been
    # missing here for the same reason it was missing from run_full_dem's region table
    # until 6.65: it was the crop most of the reasoning rested on and the one nothing
    # checked. A store left out of this tuple is not checked at all, silently.
    STORES = ("arequipa", "ancash", "lima", "huaylas", "cajatambo", "colca")

    def test_combine_copies_both_artefacts(self):
        """Asserted on the function, so the omission cannot come back quietly."""
        import inspect
        source = inspect.getsource(run_full_dem.combine)
        self.assertIn("combination_explanation.txt", source)
        self.assertIn("combined_report.json", source)

    def test_every_store_carries_its_combination_explanation(self):
        for region in self.STORES:
            store = os.path.join(_support.REPO_ROOT, "results", f"{region}_full")
            if not os.path.isdir(store):
                continue
            report = os.path.join(store, "combined_report.json")
            prose = os.path.join(store, "combination_explanation.txt")
            if not os.path.exists(report):
                continue
            self.assertTrue(os.path.exists(prose),
                            f"{region}: a combined report with no explanation beside it")

    def test_the_stored_prose_agrees_with_the_stored_numbers(self):
        """
        The check that would have caught it. Every capacity the explanation quotes must
        appear in the report it sits beside.
        """
        import json
        import re
        checked = []
        for region in self.STORES:
            store = os.path.join(_support.REPO_ROOT, "results", f"{region}_full")
            report = os.path.join(store, "combined_report.json")
            prose = os.path.join(store, "combination_explanation.txt")
            if not (os.path.exists(report) and os.path.exists(prose)):
                continue
            with open(report) as f:
                data = json.load(f)
            with open(prose) as f:
                text = f.read()
            capacities = {int(r["reported_capacity"]) for r in data.get("runs", [])
                          if r.get("reported_capacity")}
            self.assertTrue(capacities,
                            f"{region}: no capacities in the report, so this test would "
                            f"pass by checking nothing")
            checked.append(region)
            quoted = {int(m.replace(",", ""))
                      for m in re.findall(r"\b\d{1,3}(?:,\d{3})+\b", text)}
            missing = sorted(c for c in capacities if c not in quoted)
            self.assertEqual(
                missing, [],
                f"{region}: the stored explanation does not quote {missing}, so it was "
                f"written against a different run than the report beside it")
        self.assertTrue(checked, "no store was actually checked")


class TestTheEstimatorIsMonotonicInBothKnobs(unittest.TestCase):
    """
    ``costliest_sampling`` is only correct while these hold, so they are asserted.
    """

    def test_memory_falls_as_stride_rises(self):
        values = [ss.estimate_peak_memory_gb(10204, 12603, 4, s)
                  for s in (1, 2, 5, 10, 15)]
        self.assertEqual(values, sorted(values, reverse=True))

    def test_memory_falls_as_downsampling_rises(self):
        values = [ss.estimate_peak_memory_gb(10204, 12603, d, 5)
                  for d in (1, 2, 4, 8)]
        self.assertEqual(values, sorted(values, reverse=True))

    def test_the_map_estimate_falls_as_downsampling_rises(self):
        values = [ss.estimate_visualisation_memory_gb(10204, 12603, d)
                  for d in (1, 2, 4, 8)]
        self.assertEqual(values, sorted(values, reverse=True))


class TestEveryAnimationCanBuildItsFirstFrame(unittest.TestCase):
    """
    The frame nobody runs, on the notebook CI cannot execute.

    ``product_collapse`` opens on every weight at zero -- the empty product, 1
    everywhere, before any criterion is applied -- and 6.54 taught ``scoring.compose``
    to refuse an all-zero weighting. That is right for a search and wrong for this
    frame, so the animation broke. It stayed broken because notebook 07 needs an
    ffmpeg the runner does not have and is excluded from the notebook job, and the
    break is in frame 0 of a film nothing else builds.

    Rendering needs ffmpeg; *building a frame* does not. So this asks the builders for
    their first frame and nothing more, which is the part that regressed and the part
    CI can actually check.
    """

    def test_first_frame_builds(self):
        import warnings

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        import make_animations as ma

        for name, builder in sorted(ma.BUILDERS.items()):
            with self.subTest(animation=name):
                fig, anim = builder()
                try:
                    anim._func(0)
                finally:
                    # Not rendering is the point -- that is what needs the ffmpeg CI
                    # does not have -- so matplotlib's "deleted without rendering"
                    # warning is expected here rather than a symptom.
                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore", message=".*deleted without rendering.*")
                        del anim
                        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
