"""
The command line's own logic: which value wins, and what reaches the pipeline.

``main()`` reconciles four sources — an explicitly typed option, a config file,
``config/fallbacks.json`` and the built-in default — and then translates config-file
spellings into pipeline parameters. Neither half was tested, and both have been wrong:
the command line used to lose to the config file silently, and one flag kept that old
behaviour after the rest was fixed.

The pipeline itself is not run here. It is replaced with a recorder, so these are fast
and test the resolution rather than the search.
"""

import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

from _support import quiet  # noqa: F401  (also sets up sys.path)
import synthetic

import site_searcher as ss


@contextlib.contextmanager
def cli(argv, cwd, fallbacks=None):
    """
    Runs ``main()`` with the given arguments and yields the pipeline's kwargs.

    ``main()`` resolves several paths relative to the working directory and replaces
    ``sys.stdout`` with a tee to its log file, restoring neither, so both are saved and
    put back here.
    """
    captured = {}

    def recorder(**kwargs):
        captured.update(kwargs)
        return {"results": {"total_sites": 0, "total_capacity": 0, "sites": []},
                "funnel": {}, "parameters": {}}

    config_dir = os.path.join(cwd, "config")
    os.makedirs(config_dir, exist_ok=True)
    if fallbacks is not None:
        with open(os.path.join(config_dir, "fallbacks.json"), "w") as f:
            json.dump(fallbacks, f)

    # main() reads ../config/fallbacks.json and writes ../output/, i.e. relative to a
    # working directory one level down. That is the `cd src` requirement.
    work = os.path.join(cwd, "src")
    os.makedirs(work, exist_ok=True)

    real_pipeline = ss.find_grand_regions_interactive
    real_argv, real_out, real_err, real_cwd = sys.argv, sys.stdout, sys.stderr, os.getcwd()
    ss.find_grand_regions_interactive = recorder
    sys.argv = ["site_searcher.py", *argv]
    os.chdir(work)
    try:
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            ss.main()
        yield captured
    finally:
        ss.find_grand_regions_interactive = real_pipeline
        sys.argv, sys.stdout, sys.stderr = real_argv, real_out, real_err
        os.chdir(real_cwd)


class CliCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="oroscope_cli_")
        # Real, but tiny: validate_parameters checks the DEM exists before anything
        # else, which is the fail-fast behaviour these tests also cover below.
        self.dem = synthetic.write_geotiff(
            os.path.join(self.tmp, "tiny.tif"),
            np.full((8, 8), 3000.0, dtype=np.float32), -15.6, -72.3)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write_config(self, **values):
        path = os.path.join(self.tmp, "config", "run.json")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        base = {"dem_path": self.dem, "origin_lat": -15.6, "origin_lon": -72.3}
        base.update(values)
        with open(path, "w") as f:
            json.dump(base, f)
        return path


class TestPrecedence(CliCase):
    """
    An explicitly typed option beats everything. This was the other way round, and
    silently: ``--generate_config`` writes every key, so a generated config made every
    flag on the command line a no-op with no warning.
    """

    def test_a_typed_option_beats_the_config_file(self):
        config = self.write_config(min_slope_deg=3.0)
        with cli(["--config_path", config, "--min_slope_deg", "9"], self.tmp) as kw:
            self.assertEqual(kw["min_slope_deg"], 9.0)

    def test_the_config_file_is_used_when_nothing_is_typed(self):
        config = self.write_config(min_slope_deg=7.0)
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertEqual(kw["min_slope_deg"], 7.0)

    def test_the_config_file_beats_the_fallbacks(self):
        config = self.write_config(max_slope_deg=44.0)
        with cli(["--config_path", config], self.tmp,
                 fallbacks={"max_slope_deg": 11.0}) as kw:
            self.assertEqual(kw["max_slope_deg"], 44.0)

    def test_the_fallbacks_are_used_when_the_config_is_silent(self):
        config = self.write_config()
        with cli(["--config_path", config], self.tmp,
                 fallbacks={"max_slope_deg": 11.0}) as kw:
            self.assertEqual(kw["max_slope_deg"], 11.0)

    def test_a_typed_option_beats_the_fallbacks_too(self):
        config = self.write_config()
        with cli(["--config_path", config, "--max_slope_deg", "33"], self.tmp,
                 fallbacks={"max_slope_deg": 11.0}) as kw:
            self.assertEqual(kw["max_slope_deg"], 33.0)

    def test_the_built_in_default_is_the_last_resort(self):
        config = self.write_config()
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertEqual(kw["min_slope_deg"], 3.0)

    def test_a_typed_option_that_equals_the_default_still_wins(self):
        """
        The reason ``explicitly_passed`` re-parses with SUPPRESS: argparse cannot
        otherwise tell ``--candidate_stride 5`` from the default of 5.
        """
        config = self.write_config(candidate_stride=9)
        with cli(["--config_path", config, "--candidate_stride", "5"], self.tmp) as kw:
            self.assertEqual(kw["candidate_stride"], 5)


class TestTheOutputDirectoryFollowsTheSameRule(CliCase):
    """
    It did not. This one resolved before the merge loop and so kept the old precedence:
    a config file beat an explicitly typed command line, for this flag only, while
    every other flag on the same line was honoured.
    """

    def test_a_typed_output_base_beats_the_config(self):
        base = os.path.join(self.tmp, "typed")
        config = self.write_config(output_directory_base_with_given_json="../output/")
        with cli(["--config_path", config,
                  "--output_directory_base_with_given_json", base], self.tmp) as kw:
            self.assertTrue(kw["run_output_dir"].startswith(base),
                            f"expected the typed base, got {kw['run_output_dir']}")

    def test_the_config_is_used_when_nothing_is_typed(self):
        base = os.path.join(self.tmp, "from_config")
        config = self.write_config(output_directory_base_with_given_json=base)
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertTrue(kw["run_output_dir"].startswith(base))

    def test_the_run_directory_is_named_after_the_config(self):
        config = self.write_config()
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertEqual(os.path.basename(kw["run_output_dir"]), "run")


class TestConfigSpellingsReachThePipeline(CliCase):
    """
    The translation layer: negative-form flags, presets and bands. Each of these is a
    place where a config key and a pipeline parameter have different names, which is
    exactly where a new parameter gets forgotten.
    """

    def test_require_sky_becomes_require_terrain(self):
        config = self.write_config()
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertTrue(kw["require_terrain"], "the default is to require rock")
        with cli(["--config_path", config, "--require_sky"], self.tmp) as kw:
            self.assertFalse(kw["require_terrain"])

    def test_nearest_sampling_becomes_bilinear_sampling(self):
        config = self.write_config()
        with cli(["--config_path", config, "--nearest_sampling"], self.tmp) as kw:
            self.assertFalse(kw["bilinear_sampling"])

    def test_no_geomagnetic_becomes_use_geomagnetic(self):
        config = self.write_config()
        with cli(["--config_path", config, "--no_geomagnetic"], self.tmp) as kw:
            self.assertFalse(kw["use_geomagnetic"])

    def test_no_explain_becomes_explain(self):
        config = self.write_config()
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertTrue(kw["explain"], "the summary is on by default")
        with cli(["--config_path", config, "--no_explain"], self.tmp) as kw:
            self.assertFalse(kw["explain"])

    def test_an_rfi_preset_name_becomes_a_zone_list(self):
        config = self.write_config(rfi_zones="arequipa")
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertEqual(kw["rfi_zones"], ss.AREQUIPA_RFI_ZONES)

    def test_rfi_none_becomes_no_zones(self):
        config = self.write_config(rfi_zones="none")
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertIsNone(kw["rfi_zones"])

    def test_a_custom_rfi_list_passes_through(self):
        zones = [["circle", -15.5, -72.3, 10.0, "test"]]
        config = self.write_config(rfi_zones=zones)
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertEqual(kw["rfi_zones"], zones)

    def test_bands_arrive_as_tuples(self):
        config = self.write_config(grammage_band_gcm2=[236.0, 1287.0],
                                   depth_band_gcm2=[1.0, 2.0],
                                   distance_band_m=[100.0, 200.0])
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertEqual(kw["grammage_band_gcm2"], (236.0, 1287.0))
            self.assertEqual(kw["depth_band_gcm2"], (1.0, 2.0))
            self.assertEqual(kw["distance_band_m"], (100.0, 200.0))

    def test_a_single_spectral_index_stays_a_number_and_a_pair_stays_a_pair(self):
        """One value pins the spectrum; two marginalise over the range."""
        with cli(["--config_path", self.write_config(decay_spectral_index=2.0)],
                 self.tmp) as kw:
            self.assertEqual(kw["decay_spectral_index"], 2.0)
        with cli(["--config_path", self.write_config(decay_spectral_index=[1.5, 2.7])],
                 self.tmp) as kw:
            self.assertEqual(tuple(kw["decay_spectral_index"]), (1.5, 2.7))

    def test_a_full_azimuth_sweep_is_requested_with_a_negative_half_width(self):
        config = self.write_config(azimuth_half_width_deg=-1)
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertIsNone(kw["azimuth_half_width_deg"],
                              "None means the scan sweeps all 360 degrees")

    def test_score_weights_are_parsed_from_their_string_form(self):
        config = self.write_config(score_weights="shower=2,depth=0.5")
        with cli(["--config_path", config], self.tmp) as kw:
            self.assertEqual(kw["score_weights"], {"shower": 2.0, "depth": 0.5})


class TestMainCleansUpAfterItself(CliCase):
    """
    ``main()`` tees stdout and stderr into the run's log. Both the swap and the open
    file used to outlive the call, so a process that ran it twice stacked a TeeLogger
    on the previous one and leaked a handle each time — which a test suite, a sweep, or
    anything driving the CLI in a loop does by definition.
    """

    def run_main(self, argv):
        def recorder(**kwargs):
            return {"results": {"total_sites": 0, "total_capacity": 0, "sites": []},
                    "funnel": {}, "parameters": {}}

        real_pipeline = ss.find_grand_regions_interactive
        real_argv, real_cwd = sys.argv, os.getcwd()
        work = os.path.join(self.tmp, "src")
        os.makedirs(work, exist_ok=True)
        ss.find_grand_regions_interactive = recorder
        sys.argv = ["site_searcher.py", *argv]
        os.chdir(work)
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                before = sys.stdout
                ss.main()
                after = sys.stdout
            return before, after
        finally:
            ss.find_grand_regions_interactive = real_pipeline
            sys.argv = real_argv
            os.chdir(real_cwd)

    def test_the_streams_are_restored(self):
        before, after = self.run_main(["--config_path", self.write_config()])
        self.assertIs(after, before, "main() left its TeeLogger installed on stdout")

    def test_the_log_file_is_closed(self):
        config = self.write_config()
        self.run_main(["--config_path", config])
        log = os.path.join(self.tmp, "output", "run", "log.txt")
        if not os.path.exists(log):                      # config-named run directory
            log = os.path.join(self.tmp, "src", "..", "output", "run", "log.txt")
        self.assertTrue(os.path.exists(log), "the run should have written a log")
        # A closed file is the only way this reads back cleanly on every platform.
        with open(log, encoding="utf-8") as f:
            self.assertIn("Execution started at", f.read())

    def test_running_twice_does_not_stack_interceptors(self):
        config = self.write_config()
        self.run_main(["--config_path", config])
        before, after = self.run_main(["--config_path", config])
        self.assertIs(after, before)


class TestGenerateConfig(CliCase):
    def test_it_writes_a_template_and_exits_without_searching(self):
        path = os.path.join(self.tmp, "new", "template.json")
        with self.assertRaises(SystemExit) as caught:
            with cli(["--generate_config", path, "--config_preset", "lima"], self.tmp):
                pass
        self.assertEqual(caught.exception.code, 0)
        self.assertTrue(os.path.exists(path))
        with open(path) as f:
            written = json.load(f)
        self.assertEqual(written["rfi_zones"], "lima")
        self.assertEqual(written, ss.default_config("lima"))


class TestFailFast(CliCase):
    """A missing essential should stop the run before anything expensive happens."""

    def test_no_dem_anywhere_is_refused(self):
        config = os.path.join(self.tmp, "config", "empty.json")
        os.makedirs(os.path.dirname(config), exist_ok=True)
        with open(config, "w") as f:
            json.dump({"origin_lat": -15.6, "origin_lon": -72.3}, f)
        with self.assertRaises(SystemExit):
            with cli(["--config_path", config], self.tmp):
                pass

    def test_an_impossible_slope_band_is_refused(self):
        config = self.write_config(min_slope_deg=40.0, max_slope_deg=10.0)
        with self.assertRaises(SystemExit):
            with cli(["--config_path", config], self.tmp):
                pass


if __name__ == "__main__":
    unittest.main()
