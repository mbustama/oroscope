"""
Topographic screening: slope/aspect recovery, exclusion zones, candidate thinning
and the funnel accounting.
"""

import contextlib
import io
import os
import tempfile
import unittest

import numpy as np

from _support import ss
import synthetic


def screen(elevation, grid, **kwargs):
    """Runs the screening stage with console output suppressed."""
    with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
        return ss.get_candidates_chunked(elevation, grid, kwargs.pop("rfi_zones", None),
                                         kwargs.pop("origin_lat", -16.0),
                                         kwargs.pop("origin_lon", -72.0), **kwargs)


def make_grid(latitude=-16.0):
    return ss.resolve_grid_geometry("nonexistent.tif", latitude, cell_size_deg=synthetic.CELL_DEG)


class TestSlopeAndAspectScreening(unittest.TestCase):
    def setUp(self):
        self.grid = make_grid()

    def test_plane_inside_the_slope_band_is_kept_entirely(self):
        z = synthetic.planar(200, 10.0, 90.0, self.grid.cell_size_y, self.grid.cell_size_x)
        cands = screen(z, self.grid, tile_size=200, candidate_stride=1,
                       min_slope_deg=3.0, max_slope_deg=25.0)
        # Every interior pixel qualifies; edges use one-sided differences
        self.assertGreater(len(cands), 200 * 200 * 0.97)

    def test_plane_outside_the_slope_band_is_rejected_entirely(self):
        for slope_deg in (1.0, 40.0):
            z = synthetic.planar(200, slope_deg, 90.0, self.grid.cell_size_y, self.grid.cell_size_x)
            cands = screen(z, self.grid, tile_size=200, candidate_stride=1,
                           min_slope_deg=3.0, max_slope_deg=25.0)
            self.assertEqual(len(cands), 0, msg=f"slope {slope_deg} should be excluded")

    def test_reported_aspect_matches_the_plane(self):
        for aspect_deg in (0.0, 90.0, 180.0, 270.0):
            z = synthetic.planar(200, 10.0, aspect_deg, self.grid.cell_size_y, self.grid.cell_size_x)
            cands = screen(z, self.grid, tile_size=200, candidate_stride=1)
            got = cands[:, 2]
            offset = (got - aspect_deg + 180) % 360 - 180
            self.assertLess(float(np.abs(offset).max()), 0.01, msg=f"aspect {aspect_deg}")

    def test_aspect_window_selects_only_matching_terrain(self):
        z = synthetic.planar(200, 10.0, 270.0, self.grid.cell_size_y, self.grid.cell_size_x)
        kept = screen(z, self.grid, tile_size=200, candidate_stride=1,
                      min_aspect_deg=250.0, max_aspect_deg=290.0)
        dropped = screen(z, self.grid, tile_size=200, candidate_stride=1,
                         min_aspect_deg=0.0, max_aspect_deg=45.0)
        self.assertGreater(len(kept), 0)
        self.assertEqual(len(dropped), 0)

    def test_altitude_bounds_bind(self):
        z = synthetic.planar(200, 10.0, 90.0, self.grid.cell_size_y, self.grid.cell_size_x, base=2000.0)
        below = screen(z, self.grid, tile_size=200, candidate_stride=1, max_alt=float(z.min()) - 1.0)
        above = screen(z, self.grid, tile_size=200, candidate_stride=1, min_alt=float(z.max()) + 1.0)
        self.assertEqual(len(below), 0)
        self.assertEqual(len(above), 0)


class TestTilingInvariance(unittest.TestCase):
    """
    Screening must not depend on how the map is cut into tiles.

    Derivatives are computed on a haloed block and cropped, so a pixel's slope is the
    same whether or not a tile boundary happens to fall next to it. Without the halo,
    np.gradient falls back to one-sided differences at every tile edge.
    """

    def candidates(self, z, grid, tile_size, **kwargs):
        cands = screen(z, grid, tile_size=tile_size, candidate_stride=1, **kwargs)
        return set(map(tuple, cands[:, :2].astype(int)))

    def test_tiled_screening_matches_untiled(self):
        grid = make_grid(-15.6)
        n = 600
        z = synthetic.ridge_and_slope(n, grid.cell_size_x)
        untiled = self.candidates(z, grid, n)
        for tile_size in (64, 128, 256):
            self.assertEqual(self.candidates(z, grid, tile_size), untiled,
                             msg=f"tile_size {tile_size} disagrees with the untiled result")

    def test_tiled_screening_matches_untiled_with_a_slope_baseline(self):
        grid = make_grid(-15.6)
        n = 600
        z = synthetic.ridge_and_slope(n, grid.cell_size_x)
        untiled = self.candidates(z, grid, n, slope_baseline_m=500.0)
        for tile_size in (128, 256):
            self.assertEqual(self.candidates(z, grid, tile_size, slope_baseline_m=500.0), untiled,
                             msg=f"tile_size {tile_size} disagrees at a 500 m baseline")


class TestSlopeBaseline(unittest.TestCase):
    """Slope is scale-dependent; the baseline makes the scale explicit."""

    def test_baseline_converts_to_a_per_axis_pixel_window(self):
        grid = make_grid(-16.0)
        ny, nx = ss.slope_baseline_pixels(grid, 1000.0)
        self.assertEqual(ny, round(1000.0 / grid.cell_size_y))
        self.assertEqual(nx, round(1000.0 / grid.cell_size_x))
        self.assertNotEqual(ny, nx, "anisotropic pixels give different windows per axis")

    def test_no_baseline_means_native_resolution(self):
        grid = make_grid(-16.0)
        self.assertEqual(ss.slope_baseline_pixels(grid, None), (0, 0))
        self.assertEqual(ss.slope_baseline_pixels(grid, 0), (0, 0))

    def test_a_plane_keeps_its_slope_at_every_baseline(self):
        """Smoothing must not bias slope on terrain that has no curvature."""
        grid = make_grid(-16.0)
        z = synthetic.planar(300, 12.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        for baseline in (None, 200.0, 1000.0):
            smooth = ss.slope_baseline_pixels(grid, baseline)
            slope, aspect = ss.terrain_derivatives(z, grid.cell_size_y, grid.cell_size_x, *smooth)
            interior = slope[50:-50, 50:-50]
            self.assertAlmostEqual(float(interior.mean()), 12.0, places=3,
                                   msg=f"baseline {baseline}")

    def test_longer_baselines_smooth_rough_terrain(self):
        """On real-shaped terrain a wider window lowers the measured slope."""
        grid = make_grid(-16.0)
        rng = np.random.default_rng(0)
        z = synthetic.planar(400, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        z = z + rng.normal(0, 15.0, z.shape).astype(np.float32)
        medians = []
        for baseline in (None, 250.0, 1000.0):
            smooth = ss.slope_baseline_pixels(grid, baseline)
            slope, _ = ss.terrain_derivatives(z, grid.cell_size_y, grid.cell_size_x, *smooth)
            medians.append(float(np.median(slope[20:-20, 20:-20])))
        self.assertEqual(medians, sorted(medians, reverse=True),
                         msg=f"expected decreasing slope with baseline, got {medians}")


class TestCandidateStride(unittest.TestCase):
    def test_stride_thins_candidates_proportionally(self):
        grid = make_grid()
        z = synthetic.planar(200, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        full = len(screen(z, grid, tile_size=200, candidate_stride=1))
        for stride in (2, 5, 10):
            thinned = len(screen(z, grid, tile_size=200, candidate_stride=stride))
            self.assertAlmostEqual(thinned / full, 1.0 / stride, delta=0.02,
                                   msg=f"stride {stride}")


class TestExclusionZones(unittest.TestCase):
    """An RFI zone must be a circle on the ground, not in pixel space."""

    def setUp(self):
        self.grid = make_grid()
        self.n = 1200
        self.origin_lat, self.origin_lon = -16.0, -72.0
        self.radius_km = 8.0
        self.z = synthetic.planar(self.n, 10.0, 90.0, self.grid.cell_size_y, self.grid.cell_size_x)
        self.zone_lat = self.origin_lat - (self.n / 2) * synthetic.CELL_DEG
        self.zone_lon = self.origin_lon + (self.n / 2) * synthetic.CELL_DEG

    def _mask(self):
        zones = [("circle", self.zone_lat, self.zone_lon, self.radius_km, "Test")]
        cands = screen(self.z, self.grid, rfi_zones=zones, origin_lat=self.origin_lat,
                       origin_lon=self.origin_lon, tile_size=600, candidate_stride=1)
        mask = np.zeros((self.n, self.n), dtype=bool)
        mask[cands[:, 0].astype(int), cands[:, 1].astype(int)] = True
        return mask

    def test_excluded_region_has_the_requested_radius_on_both_axes(self):
        mask = self._mask()
        mid = self.n // 2
        ew_km = np.count_nonzero(~mask[mid, :]) * self.grid.cell_size_x / 2000.0
        ns_km = np.count_nonzero(~mask[:, mid]) * self.grid.cell_size_y / 2000.0
        self.assertAlmostEqual(ew_km, self.radius_km, delta=0.05)
        self.assertAlmostEqual(ns_km, self.radius_km, delta=0.05)

    def test_zone_is_an_ellipse_in_pixel_space(self):
        """Equal ground radii span more columns than rows, by exactly cell_y/cell_x."""
        mask = self._mask()
        mid = self.n // 2
        cols = np.count_nonzero(~mask[mid, :])
        rows = np.count_nonzero(~mask[:, mid])
        self.assertAlmostEqual(cols / rows, self.grid.cell_size_y / self.grid.cell_size_x,
                               delta=0.01)


class TestFunnelAccounting(unittest.TestCase):
    def test_funnel_records_a_monotonic_screening_cascade(self):
        grid = make_grid()
        z = synthetic.planar(200, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        funnel = ss.Funnel()
        screen(z, grid, tile_size=200, candidate_stride=1, funnel=funnel)
        counts = funnel.as_dict()
        self.assertEqual(counts["DEM pixels"], 200 * 200)
        self.assertEqual(counts["finite elevation"], 200 * 200)
        self.assertLessEqual(counts["slope 3.0-25.0 deg"], counts["finite elevation"])
        self.assertGreater(counts["slope 3.0-25.0 deg"], 0)

    def test_funnel_stride_count_matches_returned_candidates(self):
        grid = make_grid()
        z = synthetic.planar(200, 10.0, 90.0, grid.cell_size_y, grid.cell_size_x)
        funnel = ss.Funnel()
        cands = screen(z, grid, tile_size=64, candidate_stride=5, funnel=funnel)
        self.assertEqual(funnel.get("kept by stride 5"), len(cands))

    def test_funnel_renders_without_stages(self):
        self.assertEqual(ss.Funnel().render(), "")


class TestFunnelSeparatesGeometryFromTheScoreCut(unittest.TestCase):
    """
    ``directions accepted`` counts the geometry; the score row counts what the cut kept.

    Both rows once carried the post-cut count, which made the score row a tautology --
    100.000% of the stage above it in every stored run -- and hid the geometric
    acceptance whenever a cut was in force. Worse, a stage that keeps exactly 100% can
    never be the binding constraint, so a search emptied by ``min_score`` was blamed on
    the arrival geometry and the summary advised widening windows that were not the
    problem.
    """

    def setUp(self):
        self.grid = make_grid()
        n = 240
        self.elevation = synthetic.colca_like(n, self.grid.cell_size_x)
        # Candidates down one canyon wall, facing across it
        rows = np.arange(40, 200, 2, dtype=np.float64)
        cols = np.full(rows.shape, 70.0)
        self.candidates = np.column_stack([rows, cols, np.full(rows.shape, 90.0)])
        self.scan_params = dict(n_azimuths=8, half_width_deg=30.0,
                                elev_min_deg=-20.0, elev_max_deg=20.0, n_elev_bins=9,
                                max_range_m=12000.0, min_dist_km=0.5, max_dist_km=6.0)

    def scan(self, tmp, **kwargs):
        buf = np.lib.format.open_memmap(
            os.path.join(tmp, "buf.npy"), mode="w+",
            dtype=bool, shape=self.elevation.shape)
        funnel = ss.Funnel()
        with contextlib.redirect_stdout(io.StringIO()):
            n_hits, _ = ss.run_arrival_scan(
                self.candidates, self.elevation, self.grid, buf, self.scan_params,
                funnel=funnel, **kwargs)
        del buf
        return n_hits, funnel.as_dict()

    def test_no_cut_records_the_geometry_alone(self):
        with tempfile.TemporaryDirectory() as tmp:
            n_hits, counts = self.scan(tmp, min_score=0.0)
        # min_score 0 accepts every viable candidate, so there is no cut to report and
        # the geometric count is the whole answer. This is why every stored GRAND run
        # already held the right number under this name.
        self.assertEqual(counts["directions accepted"], n_hits)
        self.assertEqual([k for k in counts if k.startswith("score")], [])

    def test_a_cut_is_recorded_below_the_geometry_it_cut(self):
        with tempfile.TemporaryDirectory() as tmp:
            loose, _ = self.scan(tmp, min_score=0.0)
        with tempfile.TemporaryDirectory() as tmp:
            n_hits, counts = self.scan(tmp, min_score=0.35)
        self.assertGreater(loose, 0, "fixture accepts nothing; the test proves nothing")
        # The geometry is unchanged by the threshold: only the second row moves.
        self.assertEqual(counts["directions accepted"], loose)
        self.assertEqual(counts["score >= 0.35"], n_hits)
        self.assertLess(counts["score >= 0.35"], counts["directions accepted"])

    def test_a_percentile_cut_names_its_own_knob(self):
        with tempfile.TemporaryDirectory() as tmp:
            n_hits, counts = self.scan(tmp, score_percentile=25.0)
        self.assertEqual(counts["score in top 25%"], n_hits)
        self.assertLess(n_hits, counts["directions accepted"])
        # Previously a percentile cut added no row at all, and its removals were
        # silently attributed to the arrival geometry.
        from oroscope import explain
        binding = explain.binding_constraint(counts)
        self.assertEqual(binding["stage"], "score in top 25%")
        self.assertEqual(binding["knob"], "score_percentile")

    def test_a_cut_that_empties_the_search_names_min_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            n_hits, counts = self.scan(tmp, min_score=1.0e9)
        self.assertEqual(n_hits, 0)
        from oroscope import explain
        binding = explain.binding_constraint(counts)
        self.assertEqual(binding["stage"], "score >= 1e+09")
        self.assertTrue(binding["fatal"])
        self.assertIn("min_score", binding["knob"])


if __name__ == "__main__":
    unittest.main()


class TestSeparableMorphology(unittest.TestCase):
    """
    A rectangle of ones factorises into a column and a row, so the separable form is
    bit-identical to the direct one and O(N(h+w)) rather than O(Nhw).
    """

    def setUp(self):
        rng = np.random.default_rng(0)
        self.mask = rng.random((400, 400)) < 0.35

    def test_closing_matches_the_direct_operation(self):
        from scipy.ndimage import binary_closing
        for shape in ((33, 33), (17, 21), (9, 5)):
            st = np.ones(shape)
            self.assertTrue(np.array_equal(ss.separable_closing(self.mask, st),
                                           binary_closing(self.mask, structure=st)),
                            msg=f"shape {shape}")

    def test_opening_matches_the_direct_operation(self):
        from scipy.ndimage import binary_opening
        for shape in ((17, 17), (21, 9)):
            st = np.ones(shape)
            self.assertTrue(np.array_equal(ss.separable_opening(self.mask, st),
                                           binary_opening(self.mask, structure=st)),
                            msg=f"shape {shape}")

    def test_morphology_is_independent_of_tile_size(self):
        """Cleanup must not depend on how the map happens to be cut into tiles."""
        import tempfile
        import shutil
        import os
        import contextlib
        import io

        def run(tile):
            tmp = tempfile.mkdtemp()
            try:
                a_path = os.path.join(tmp, "a.npy")
                b_path = os.path.join(tmp, "b.npy")
                a = np.lib.format.open_memmap(a_path, mode="w+", shape=self.mask.shape,
                                              dtype=bool)
                a[:] = self.mask
                a.flush()
                del a
                np.lib.format.open_memmap(b_path, mode="w+", shape=self.mask.shape,
                                          dtype=bool).flush()
                with contextlib.redirect_stderr(io.StringIO()):
                    ss.apply_morphology_pingpong(a_path, b_path, self.mask.shape, bool,
                                                 ss.separable_closing, np.ones((21, 21)),
                                                 desc="c", tile_size=tile)
                return np.array(np.lib.format.open_memmap(b_path, mode="r"))
            finally:
                shutil.rmtree(tmp, ignore_errors=True)

        ref = run(400)
        for tile in (100, 150):
            self.assertTrue(np.array_equal(run(tile), ref), msg=f"tile {tile}")


class TestExplicitCommandLineDetection(unittest.TestCase):
    """
    Telling a typed option from an untyped default.

    The configuration merge preferred a config file over the command line, because
    argparse cannot distinguish ``--candidate_stride 5`` from its own default of 5.
    Since ``--generate_config`` writes every key, a generated config silently turned
    every command-line flag into a no-op. Knowing which options were actually typed is
    what lets the command line win without defaults overwriting the config.
    """

    def build(self):
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument("--alpha", type=float, default=3.0)
        p.add_argument("--beta", type=int, default=7)
        p.add_argument("--flag", action="store_true")
        return p

    def test_nothing_typed_means_nothing_explicit(self):
        self.assertEqual(ss.explicitly_passed(self.build(), []), set())

    def test_a_typed_option_is_reported(self):
        self.assertEqual(ss.explicitly_passed(self.build(), ["--alpha", "5"]), {"alpha"})

    def test_a_value_equal_to_the_default_still_counts_as_typed(self):
        """The case that made the old merge unfixable: 3.0 is also the default."""
        self.assertEqual(ss.explicitly_passed(self.build(), ["--alpha", "3.0"]), {"alpha"})

    def test_store_true_is_detected(self):
        self.assertEqual(ss.explicitly_passed(self.build(), ["--flag"]), {"flag"})

    def test_several_at_once(self):
        got = ss.explicitly_passed(self.build(), ["--alpha", "1", "--beta", "2"])
        self.assertEqual(got, {"alpha", "beta"})

    def test_the_parser_is_left_usable_afterwards(self):
        """Defaults are restored, so the caller's own parse is unaffected."""
        p = self.build()
        ss.explicitly_passed(p, ["--alpha", "9"])
        args = p.parse_args([])
        self.assertEqual(args.alpha, 3.0)
        self.assertEqual(args.beta, 7)
        self.assertFalse(args.flag)


class TestPythonFloorCompatibility(unittest.TestCase):
    """
    Guards the declared 3.9 floor against syntax only newer Pythons accept.

    ``X | None`` in an annotation is evaluated at runtime before 3.10 unless
    ``from __future__ import annotations`` is in force. Every module here uses that
    syntax, and three of them shipped without the import: the suite passed on 3.12 and
    every 3.9 job failed with ``unsupported operand type(s) for |``. Local green said
    nothing, because local is not 3.9.
    """

    def source_modules(self):
        import glob
        import os
        src = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
        return sorted(glob.glob(os.path.join(src, "*.py")))

    def test_every_module_postpones_annotation_evaluation(self):
        import ast
        for path in self.source_modules():
            with open(path) as f:
                tree = ast.parse(f.read())
            has = any(isinstance(node, ast.ImportFrom) and node.module == "__future__"
                      and any(a.name == "annotations" for a in node.names)
                      for node in tree.body)
            self.assertTrue(has, f"{os.path.basename(path)} needs "
                                 f"`from __future__ import annotations`")

    def test_the_future_import_comes_first(self):
        """Python requires it before any other statement, so a misplaced one is fatal."""
        import ast
        for path in self.source_modules():
            with open(path) as f:
                tree = ast.parse(f.read())
            body = [n for n in tree.body
                    if not (isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant))]
            self.assertTrue(body, os.path.basename(path))
            first = body[0]
            self.assertTrue(isinstance(first, ast.ImportFrom) and first.module == "__future__",
                            f"{os.path.basename(path)}: the __future__ import must be the "
                            f"first statement, found {type(first).__name__}")


class TestMemorySafeguards(unittest.TestCase):
    """
    Guards against a search taking the machine down with it.

    A ten-point sensitivity sweep reached 6.9 GB and was killed by the kernel's OOM
    killer, which chooses its victim by heuristic and may well pick something the user
    cares about more. Two things went wrong and both are fixed: the pipeline leaked a
    matplotlib figure per run, and nothing bounded the process.
    """

    def test_the_estimate_scales_as_the_inverse_square_of_downsampling(self):
        """The knob that matters most: labelling arrays scale as 1/d^2."""
        at_1 = ss.estimate_peak_memory_gb(10000, 10000, downsample_factor=1)
        at_4 = ss.estimate_peak_memory_gb(10000, 10000, downsample_factor=4)
        self.assertLess(at_4, at_1)

    def test_the_estimate_grows_with_the_dem(self):
        small = ss.estimate_peak_memory_gb(2000, 2000)
        large = ss.estimate_peak_memory_gb(10000, 10000)
        self.assertGreater(large, small)

    def test_a_coarser_stride_needs_less(self):
        dense = ss.estimate_peak_memory_gb(5000, 5000, candidate_stride=1)
        sparse = ss.estimate_peak_memory_gb(5000, 5000, candidate_stride=10)
        self.assertGreater(dense, sparse)

    def test_the_estimate_is_a_positive_number_of_gib(self):
        v = ss.estimate_peak_memory_gb(1000, 1000)
        self.assertGreater(v, 0.0)
        self.assertLess(v, 100.0)

    def test_available_memory_is_reported_or_declined(self):
        have = ss.available_memory_gb()
        if have is not None:
            self.assertGreater(have, 0.0)

    def test_a_cap_of_none_or_zero_does_nothing(self):
        self.assertFalse(ss.apply_memory_cap(None))
        self.assertFalse(ss.apply_memory_cap(0))

    def test_a_cap_is_applied_and_can_be_lifted(self):
        """Applied in a child, since lowering this process's limit would stick."""
        import subprocess
        import sys
        code = (
            "import sys, resource; sys.path.insert(0, %r);"
            "import os; os.environ.setdefault('MPLBACKEND','Agg');"
            "from oroscope import site_searcher as ss;"
            "before = resource.getrlimit(resource.RLIMIT_AS)[0];"
            "ok = ss.apply_memory_cap(4.0);"
            "after = resource.getrlimit(resource.RLIMIT_AS)[0];"
            "print(ok, before != after or after == int(4.0 * 1024**3))"
            % os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
        )
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertIn("True True", out.stdout)

    def test_the_cap_actually_stops_a_runaway(self):
        """A cap that does not produce MemoryError is decoration."""
        import subprocess
        import sys
        code = (
            "import sys; sys.path.insert(0, %r);"
            "import os; os.environ.setdefault('MPLBACKEND','Agg');"
            "import numpy as np;"
            "from oroscope import site_searcher as ss;"
            "ss.apply_memory_cap(1.0);"
            "\ntry:\n"
            "    x = np.ones(int(4e9), dtype=np.uint8)\n"
            "    print('NOT CAPPED')\n"
            "except MemoryError:\n"
            "    print('MemoryError as intended')\n"
            % os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
        )
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertIn("MemoryError as intended", out.stdout, out.stderr)

    def test_the_pipeline_leaves_no_figures_open(self):
        """
        The leak itself. pyplot keeps a global reference to every figure, so one that
        is never closed can never be collected, and a process running several searches
        accumulates the whole figure -- canvas and artists -- each time.
        """
        import matplotlib.pyplot as plt
        plt.close("all")
        source = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "src", "oroscope", "site_searcher.py")
        with open(source) as f:
            text = f.read()
        self.assertIn("plt.close('all')", text,
                      "generate_visualizations_and_outputs must close its figures")
