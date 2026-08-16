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

from oroscope import site_searcher as ss


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

    # main() reads ../config/fallbacks.json relative to a working directory one level
    # down. Config paths and the output base are resolved against the configuration
    # file now, but the fallbacks lookup is still cwd-relative, so these run from src/.
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


class TestTheThreeSourcesOfDefaultsAgree(unittest.TestCase):
    """
    A parameter has three places it can state a default: the pipeline's signature, the
    argparse parser, and ``default_config()``. They disagreed on **ten** of them, so
    omitting a parameter meant different things depending on which door you came in by
    — ``search_mode`` was ``single`` from Python and ``distributed`` from a shell, and
    ``min_dist_km`` was 30 km against 10.

    That is not the parity §6.24 closed, which was about capability. This is parity of
    *meaning*, and it is the one a user actually trips over.
    """

    # Placeholders, not defaults: a generated template is meant to be edited, and these
    # two say so by being obviously unreal. Everything else must agree.
    PLACEHOLDERS = {"dem_path", "region_name"}

    # Resolved by main() rather than passed through, or not a pipeline parameter.
    CLI_ONLY = {"config_path", "generate_config", "config_preset", "print_info",
                "output_directory_base_with_given_json", "require_sky", "resume_dir",
                "rfi_zones", "run_output_dir", "help"}

    @classmethod
    def setUpClass(cls):
        import argparse
        import inspect
        captured = {}
        real = argparse.ArgumentParser.parse_args

        def intercept(self, *args, **kwargs):
            captured["parser"] = self
            raise SystemExit(0)

        argparse.ArgumentParser.parse_args = intercept
        argv = sys.argv
        sys.argv = ["oroscope"]
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                ss.main()
        except SystemExit:
            pass
        finally:
            argparse.ArgumentParser.parse_args = real
            sys.argv = argv

        cls.signature = {k: v.default for k, v in inspect.signature(
            ss.find_grand_regions_interactive).parameters.items()
            if v.default is not inspect.Parameter.empty}
        cls.cli = {a.dest: a.default for a in captured["parser"]._actions}
        cls.template = ss.default_config()

    def disagreements(self, a, b, name_a, name_b):
        out = []
        for key in sorted(set(a) & set(b) - self.PLACEHOLDERS - self.CLI_ONLY):
            if a[key] != b[key]:
                out.append(f"{key}: {name_a}={a[key]!r} but {name_b}={b[key]!r}")
        return out

    def test_the_signature_agrees_with_the_command_line(self):
        bad = self.disagreements(self.signature, self.cli, "signature", "CLI")
        self.assertEqual(bad, [], "\n".join(bad))

    def test_the_signature_agrees_with_the_config_template(self):
        bad = self.disagreements(self.signature, self.template, "signature", "template")
        self.assertEqual(bad, [], "\n".join(bad))

    def test_the_command_line_agrees_with_the_config_template(self):
        bad = self.disagreements(self.cli, self.template, "CLI", "template")
        self.assertEqual(bad, [], "\n".join(bad))

    def test_the_template_does_not_offer_zero_as_an_origin(self):
        """
        Zero is a *valid* coordinate, in the Gulf of Guinea. A placeholder someone
        forgets to edit must not produce a run georeferenced to the wrong continent;
        null means "read the DEM's own tiepoint", which is the recommended use anyway.
        """
        self.assertIsNone(self.template["origin_lat"])
        self.assertIsNone(self.template["origin_lon"])


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


class TestConfigPathsAreRelativeToTheConfig(unittest.TestCase):
    """
    A configuration describes where its DEM sits relative to *itself*.

    Resolving against the working directory instead is what made the bundled configs
    run only from ``src/``. These pin the new rule and the fallback that keeps the old
    one working.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.cfg_dir = os.path.join(self.tmp, "config")
        self.dem_dir = os.path.join(self.tmp, "input", "dem")
        os.makedirs(self.cfg_dir)
        os.makedirs(self.dem_dir)
        self.dem = os.path.join(self.dem_dir, "x.tif")
        with open(self.dem, "w") as f:
            f.write("not really a tif")

    def write(self, **keys):
        path = os.path.join(self.cfg_dir, "run.json")
        with open(path, "w") as f:
            json.dump(keys, f)
        return path

    def test_a_relative_dem_resolves_against_the_config_from_any_cwd(self):
        path = self.write(dem_path=os.path.join("..", "input", "dem", "x.tif"))
        here = os.getcwd()
        try:
            for cwd in (self.tmp, self.cfg_dir, tempfile.gettempdir()):
                os.chdir(cwd)
                with self.subTest(cwd=cwd):
                    self.assertEqual(ss.load_config(path)["dem_path"], self.dem)
        finally:
            os.chdir(here)

    def test_an_absolute_path_is_left_alone(self):
        path = self.write(dem_path="/data/elsewhere.tif")
        self.assertEqual(ss.load_config(path)["dem_path"], "/data/elsewhere.tif")

    def test_a_working_directory_relative_path_still_works_with_a_warning(self):
        # The old behaviour. Breaking it silently to fix the wart would be a poor trade.
        path = self.write(dem_path=os.path.join("input", "dem", "x.tif"))
        here = os.getcwd()
        buf = io.StringIO()
        try:
            os.chdir(self.tmp)
            with contextlib.redirect_stdout(buf):
                loaded = ss.load_config(path)
        finally:
            os.chdir(here)
        self.assertEqual(loaded["dem_path"], os.path.join("input", "dem", "x.tif"))
        self.assertIn("working directory", buf.getvalue())

    def test_a_path_that_exists_nowhere_names_what_the_config_asked_for(self):
        path = self.write(dem_path=os.path.join("..", "input", "dem", "absent.tif"))
        resolved = ss.load_config(path)["dem_path"]
        self.assertTrue(resolved.endswith(os.path.join("input", "dem", "absent.tif")))
        self.assertTrue(os.path.isabs(resolved),
                        "an error should name the path the config meant, not a fragment")

    def test_the_repository_layout_reads_the_same_either_way(self):
        # config/ and src/ are both one level below the root, so no shipped config had
        # to change: ../input/... names the same file read from either.
        from_config = os.path.normpath(os.path.join("/repo/config", "../input/d.tif"))
        from_src = os.path.normpath(os.path.join("/repo/src", "../input/d.tif"))
        self.assertEqual(from_config, from_src)


class TestOneTranslationForEveryCaller(unittest.TestCase):
    """
    ``config_to_pipeline_kwargs`` is the single translation, and the reason it exists.

    The mapping was written out three times -- ``main()``, the child ``sensitivity``
    spawns, and ``tools/run_full_dem.py`` -- and the copies drifted. The sweep
    child splatted the config straight into the pipeline, so a preset name reached a
    function that iterates its argument, every character failed the ``item[0] ==
    'circle'`` test, and the point searched with no exclusion zones at all. No
    exception, no warning, a plausible wrong answer. These pin the translation itself,
    away from ``main()``, so every caller is covered by testing one function.
    """

    def test_an_rfi_preset_resolves_rather_than_being_iterated(self):
        kw = ss.config_to_pipeline_kwargs({"rfi_zones": "arequipa"}, quiet=True)
        self.assertEqual(len(kw["rfi_zones"]), 5)
        self.assertEqual(kw["rfi_zones"][0][0], "circle",
                         "a resolved zone list, not the letters of 'arequipa'")

    def test_the_silent_failure_this_prevents(self):
        # What the sweep child used to do: hand the raw string to a consumer that
        # iterates it. Every character is skipped, and the result is indistinguishable
        # from having asked for no zones.
        zones = [item for item in "arequipa" if item[0] == "circle"]
        self.assertEqual(zones, [], "the old path silently yielded nothing")
        self.assertEqual(len(ss.resolve_rfi_zones("arequipa")), 5,
                         "the translation yields the five real zones")

    def test_none_and_missing_both_mean_no_zones(self):
        for value in ("none", "NONE", None):
            self.assertIsNone(ss.resolve_rfi_zones(value))

    def test_require_sky_is_inverted_for_every_caller(self):
        self.assertTrue(ss.config_to_pipeline_kwargs({}, quiet=True)["require_terrain"])
        kw = ss.config_to_pipeline_kwargs({"require_sky": True}, quiet=True)
        self.assertFalse(kw["require_terrain"])
        self.assertNotIn("require_sky", kw, "the pipeline does not take that spelling")

    def test_command_line_only_keys_do_not_reach_the_pipeline(self):
        kw = ss.config_to_pipeline_kwargs(
            {"print_info": True, "output_directory_base_with_given_json": "../output/"},
            quiet=True)
        for key in ("print_info", "output_directory_base_with_given_json"):
            self.assertNotIn(key, kw)

    def test_bands_become_tuples_because_json_has_no_tuples(self):
        kw = ss.config_to_pipeline_kwargs(
            {"grammage_band_gcm2": [236.0, 1287.0], "depth_band_gcm2": [700.0, 2800.0]},
            quiet=True)
        self.assertEqual(kw["grammage_band_gcm2"], (236.0, 1287.0))
        self.assertEqual(kw["depth_band_gcm2"], (700.0, 2800.0))

    def test_a_misspelled_key_is_dropped_and_named(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            kw = ss.config_to_pipeline_kwargs({"min_slop_deg": 3.0})
        self.assertNotIn("min_slop_deg", kw, "an unknown key must not reach the pipeline")
        self.assertIn("min_slop_deg", buf.getvalue(),
                      "and it must be named, or a typo is silently ignored")

    def test_overrides_win_over_the_configuration(self):
        kw = ss.config_to_pipeline_kwargs({"dem_path": "from_config.tif"},
                                          quiet=True, dem_path="override.tif")
        self.assertEqual(kw["dem_path"], "override.tif")

    def test_every_real_config_translates(self):
        root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "config")
        expected_zones = {"grand_colca_config": 5, "tambo_colca_config": 0,
                          "grand_arequipa_full": 5, "tambo_arequipa_full": 0}
        for name, zones in expected_zones.items():
            path = os.path.join(root, f"{name}.json")
            if not os.path.exists(path):                 # pragma: no cover - layout
                continue
            with self.subTest(config=name):
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    kw = ss.config_to_pipeline_kwargs(ss.load_config(path))
                self.assertNotIn("Ignoring", buf.getvalue(),
                                 "a shipped config must translate without complaint")
                self.assertEqual(len(kw["rfi_zones"] or []), zones)

    def test_the_filter_is_bound_to_the_real_signature(self):
        # Read at call time, the filter followed whatever the pipeline name pointed at,
        # so a test double or a decorator presenting (*args, **kwargs) collapsed it and
        # every parameter was dropped.
        real = ss.find_grand_regions_interactive
        try:
            ss.find_grand_regions_interactive = lambda *a, **k: None
            kw = ss.config_to_pipeline_kwargs({"min_slope_deg": 11.0}, quiet=True)
        finally:
            ss.find_grand_regions_interactive = real
        self.assertEqual(kw["min_slope_deg"], 11.0)


if __name__ == "__main__":
    unittest.main()
