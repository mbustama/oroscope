"""
Grid packing: does the antenna count match the geometry it claims to simulate?

Row and column strides are separate because an equal ground spacing is a different
number of pixels on each axis of a geographic grid.
"""

import math
import unittest

import numpy as np

from _support import ss
import synthetic

SQUARE, HEX = 0, 1


def full_mask(h, w):
    return np.ones((h, w), dtype=bool)


class TestSquareGrid(unittest.TestCase):
    def test_unit_spacing_fills_every_pixel(self):
        self.assertEqual(ss.count_grid_capacity(full_mask(20, 30), 1, 1, SQUARE), 600)

    def test_count_matches_the_number_of_grid_intersections(self):
        for h, w, sr, sc in [(100, 100, 10, 10), (100, 100, 7, 13), (55, 41, 8, 8)]:
            expected = math.ceil(h / sr) * math.ceil(w / sc)
            self.assertEqual(ss.count_grid_capacity(full_mask(h, w), sr, sc, SQUARE), expected,
                             msg=f"{h}x{w} stride {sr},{sc}")

    def test_separate_row_and_column_strides_are_honoured(self):
        mask = full_mask(100, 100)
        self.assertEqual(ss.count_grid_capacity(mask, 10, 20, SQUARE), 10 * 5)
        self.assertEqual(ss.count_grid_capacity(mask, 20, 10, SQUARE), 5 * 10)

    def test_empty_terrain_holds_nothing(self):
        self.assertEqual(ss.count_grid_capacity(np.zeros((50, 50), dtype=bool), 5, 5, SQUARE), 0)

    def test_only_valid_pixels_are_counted(self):
        mask = np.zeros((100, 100), dtype=bool)
        mask[:50, :] = True
        self.assertEqual(ss.count_grid_capacity(mask, 10, 10, SQUARE), 5 * 10)


class TestHexGrid(unittest.TestCase):
    def test_only_rows_on_the_sin60_step_are_sampled(self):
        """A single occupied row contributes only when it lands on int(0.866*spacing)."""
        spacing = 16
        v_step = int(spacing * 0.866025)          # 13
        self.assertEqual(v_step, 13)
        for row in range(0, 40):
            mask = np.zeros((200, 200), dtype=bool)
            mask[row, :] = True
            counted = ss.count_grid_capacity(mask, spacing, spacing, HEX) > 0
            self.assertEqual(counted, row % v_step == 0, msg=f"row {row}")

    def test_alternate_rows_are_staggered_by_half_the_column_spacing(self):
        """Row 0 starts at column 0; the next sampled row starts at spacing_c // 2."""
        spacing = 16
        v_step = int(spacing * 0.866025)
        unstaggered = np.zeros((200, 200), dtype=bool)
        unstaggered[0, :] = True
        staggered = np.zeros((200, 200), dtype=bool)
        staggered[v_step, :] = True
        # columns 0,16,...,192 -> 13 ;  columns 8,24,...,184 -> 12
        self.assertEqual(ss.count_grid_capacity(unstaggered, spacing, spacing, HEX), 13)
        self.assertEqual(ss.count_grid_capacity(staggered, spacing, spacing, HEX), 12)

    def test_hex_packs_more_densely_than_square(self):
        mask = full_mask(200, 200)
        square = ss.count_grid_capacity(mask, 10, 10, SQUARE)
        hexagonal = ss.count_grid_capacity(mask, 10, 10, HEX)
        self.assertGreater(hexagonal, square)
        # The gain is spacing/int(0.866*spacing), not the analytic 1/0.866 —
        # see TestHexQuantizationBias below.
        self.assertAlmostEqual(hexagonal / square, 10 / int(10 * 0.866025), delta=1e-6)

    def test_degenerate_spacing_does_not_divide_by_zero(self):
        self.assertGreater(ss.count_grid_capacity(full_mask(10, 10), 1, 1, HEX), 0)


class TestHexQuantizationBias(unittest.TestCase):
    """
    Characterization tests for a known defect: reported capacity exceeds what the
    requested ground spacing allows.

    Three independent truncations each shrink the grid — spacing_r, spacing_c and
    v_step are all produced with int() — so antennas end up closer together than
    asked. The error grows as the spacing approaches the pixel size:

        1000 m spacing (GRAND)  ->  +7%
         150 m spacing (TAMBO)  -> +42%
         100 m spacing (TAMBO)  -> +58%

    These tests pin the current behaviour so the fix (phase 1/2 layout models) shows
    up as a deliberate change rather than a silent one. When capacity is computed
    without integer stamping, they should be rewritten to assert the analytic density.
    """

    def setUp(self):
        self.grid = ss.resolve_grid_geometry("nonexistent.tif", -16.0,
                                             cell_size_deg=synthetic.CELL_DEG)
        self.h = self.w = 3000

    def measured_vs_analytic(self, spacing_m):
        spacing_r = max(1, int(spacing_m / self.grid.cell_size_y))
        spacing_c = max(1, int(spacing_m / self.grid.cell_size_x))
        count = ss.count_grid_capacity(full_mask(self.h, self.w), spacing_r, spacing_c, HEX)
        area_km2 = ((self.h * self.grid.cell_size_y / 1000.0)
                    * (self.w * self.grid.cell_size_x / 1000.0))
        analytic = area_km2 / (math.sqrt(3) / 2 * (spacing_m / 1000.0) ** 2)
        return count / analytic

    def test_overcounts_modestly_at_grand_spacing(self):
        ratio = self.measured_vs_analytic(1000.0)
        self.assertGreater(ratio, 1.0)
        self.assertAlmostEqual(ratio, 1.074, delta=0.01)

    def test_overcounts_severely_at_tambo_spacing(self):
        """At ~3 pixels per spacing the integer grid cannot represent the layout."""
        self.assertAlmostEqual(self.measured_vs_analytic(150.0), 1.423, delta=0.02)
        self.assertAlmostEqual(self.measured_vs_analytic(100.0), 1.581, delta=0.02)

    def test_bias_grows_as_spacing_approaches_the_pixel_size(self):
        ratios = [self.measured_vs_analytic(d) for d in (1000.0, 500.0, 200.0, 100.0)]
        self.assertEqual(ratios, sorted(ratios), "bias should increase at finer spacing")


if __name__ == "__main__":
    unittest.main()
