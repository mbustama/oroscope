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

    STORES = ("arequipa", "ancash", "lima", "huaylas", "cajatambo")

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


if __name__ == "__main__":
    unittest.main()
