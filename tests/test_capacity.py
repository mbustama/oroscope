"""
Grid packing: does the detector count match the geometry it claims to simulate?

Positions are laid out in metres and then looked up in the pixel grid, so the count
follows the requested ground spacing at any resolution. The pixel sizes are passed
separately because they differ on each axis of a geographic grid.
"""

import math
import unittest

import numpy as np

from _support import ss
import synthetic

SQUARE, HEX = 0, 1

# A square metre-grid, so "spacing in metres" and "spacing in pixels" coincide and the
# geometry under test is easy to read off
UNIT = 1.0


def full_mask(h, w):
    return np.ones((h, w), dtype=bool)


def count(mask, spacing_m, code, cell_y=UNIT, cell_x=UNIT):
    return ss.count_grid_capacity(mask, cell_y, cell_x, spacing_m, code)


class TestSquareGrid(unittest.TestCase):
    def test_unit_spacing_fills_every_pixel(self):
        self.assertEqual(count(full_mask(20, 30), 1.0, SQUARE), 600)

    def test_count_matches_the_number_of_grid_intersections(self):
        for h, w, s in [(100, 100, 10.0), (100, 100, 7.0), (55, 41, 8.0)]:
            expected = math.ceil(h / s) * math.ceil(w / s)
            self.assertEqual(count(full_mask(h, w), s, SQUARE), expected,
                             msg=f"{h}x{w} spacing {s}")

    def test_pixel_sizes_are_honoured_separately(self):
        """100 rows of 100 m by 100 columns of 200 m is 10 km by 20 km: 10 rows of 20."""
        mask = full_mask(100, 100)
        self.assertEqual(count(mask, 1000.0, SQUARE, cell_y=100.0, cell_x=200.0), 10 * 20)

    def test_the_row_pitch_uses_the_row_pixel_size(self):
        """A count alone cannot catch a swapped axis, since it is symmetric; row
        placement can. At 100 m rows a 1 km lattice lands on every tenth pixel."""
        on_lattice = np.zeros((100, 1), dtype=bool)
        on_lattice[::10, 0] = True
        self.assertEqual(count(on_lattice, 1000.0, SQUARE, cell_y=100.0, cell_x=1000.0), 10)
        off_lattice = np.zeros((100, 1), dtype=bool)
        off_lattice[1::10, 0] = True
        self.assertEqual(count(off_lattice, 1000.0, SQUARE, cell_y=100.0, cell_x=1000.0), 0)

    def test_empty_terrain_holds_nothing(self):
        self.assertEqual(count(np.zeros((50, 50), dtype=bool), 5.0, SQUARE), 0)

    def test_only_valid_pixels_are_counted(self):
        mask = np.zeros((100, 100), dtype=bool)
        mask[:50, :] = True
        self.assertEqual(count(mask, 10.0, SQUARE), 5 * 10)


class TestHexGrid(unittest.TestCase):
    def test_rows_sit_on_the_sin60_pitch(self):
        """A single occupied row contributes only where a lattice row lands on it."""
        spacing = 16.0
        pitch = spacing * math.sin(math.radians(60.0))       # 13.856, not int() 13
        sampled = {int(k * pitch) for k in range(4)}
        for row in range(0, 40):
            mask = np.zeros((200, 200), dtype=bool)
            mask[row, :] = True
            self.assertEqual(count(mask, spacing, HEX) > 0, row in sampled,
                             msg=f"row {row}")

    def test_alternate_rows_are_staggered_by_half_the_spacing(self):
        """Row 0 starts at x=0; the next lattice row starts half a spacing along."""
        spacing = 16.0
        pitch = spacing * math.sin(math.radians(60.0))
        unstaggered = np.zeros((200, 200), dtype=bool)
        unstaggered[0, :] = True
        staggered = np.zeros((200, 200), dtype=bool)
        staggered[int(pitch), :] = True
        # x = 0,16,...,192 -> 13 ;  x = 8,24,...,184 -> 12
        self.assertEqual(count(unstaggered, spacing, HEX), 13)
        self.assertEqual(count(staggered, spacing, HEX), 12)

    def test_hex_packs_more_densely_than_square_by_the_analytic_ratio(self):
        # Large enough that the partial row at the far edge does not skew the ratio
        mask = full_mask(2000, 2000)
        square = count(mask, 10.0, SQUARE)
        hexagonal = count(mask, 10.0, HEX)
        self.assertGreater(hexagonal, square)
        # Now the true 1/sin(60) = 1.1547, not spacing/int(0.866*spacing)
        self.assertAlmostEqual(hexagonal / square, 1.0 / math.sin(math.radians(60.0)),
                               delta=0.01)

    def test_degenerate_spacing_returns_nothing_rather_than_dividing_by_zero(self):
        self.assertEqual(count(full_mask(10, 10), 0.0, HEX), 0)
        self.assertEqual(count(full_mask(10, 10), -5.0, HEX), 0)
        self.assertGreater(count(full_mask(10, 10), 1.0, HEX), 0)

    def test_spacing_finer_than_a_pixel_places_several_per_pixel(self):
        """The continuum limit: capacity is area over area-per-detector."""
        mask = full_mask(10, 10)
        self.assertGreater(count(mask, 0.5, SQUARE), mask.size)


class TestDensityMatchesTheRequestedSpacing(unittest.TestCase):
    """
    The count must follow the requested ground spacing, at any DEM resolution.

    This replaces a set of characterization tests that pinned a defect: the layout used
    to be stamped as an integer pixel stride, and the three separate int() truncations
    (spacing_r, spacing_c and the hex pitch) each pulled the detectors closer together
    than asked. Reported capacity came out at 1.074x the analytic density at GRAND's
    1 km, and 1.581x at TAMBO's 100 m, where one separation spans only about three
    pixels of a 30 m DEM.

    Placing positions in metres removes the strides, so the density is now correct on
    both scales and the old numbers are what a regression would look like.
    """

    def setUp(self):
        self.grid = ss.resolve_grid_geometry("nonexistent.tif", -16.0,
                                             cell_size_deg=synthetic.CELL_DEG)
        self.h = self.w = 3000

    def measured_vs_analytic(self, spacing_m, code=HEX):
        n = ss.count_grid_capacity(full_mask(self.h, self.w), self.grid.cell_size_y,
                                   self.grid.cell_size_x, spacing_m, code)
        area_km2 = ((self.h * self.grid.cell_size_y / 1000.0)
                    * (self.w * self.grid.cell_size_x / 1000.0))
        per_detector = (math.sqrt(3) / 2 if code == HEX else 1.0) * (spacing_m / 1000.0) ** 2
        return n / (area_km2 / per_detector)

    def test_matches_analytic_density_at_grand_spacing(self):
        self.assertAlmostEqual(self.measured_vs_analytic(1000.0), 1.0, delta=0.02)

    def test_matches_analytic_density_at_tambo_spacing(self):
        """The case the old integer stamping got wrong by 58%."""
        self.assertAlmostEqual(self.measured_vs_analytic(150.0), 1.0, delta=0.02)
        self.assertAlmostEqual(self.measured_vs_analytic(100.0), 1.0, delta=0.02)

    def test_square_grid_matches_its_own_analytic_density(self):
        for spacing in (1000.0, 150.0, 100.0):
            self.assertAlmostEqual(self.measured_vs_analytic(spacing, SQUARE), 1.0,
                                   delta=0.02, msg=f"{spacing} m")

    def test_accuracy_no_longer_degrades_as_spacing_approaches_the_pixel_size(self):
        for spacing in (1000.0, 500.0, 200.0, 100.0, 60.0):
            self.assertAlmostEqual(self.measured_vs_analytic(spacing), 1.0, delta=0.03,
                                   msg=f"{spacing} m")

    def test_capacity_scales_with_area_not_with_pixel_count(self):
        """Same ground area at half the DEM resolution must hold the same detectors."""
        fine = ss.count_grid_capacity(full_mask(1000, 1000), 30.0, 30.0, 300.0, HEX)
        coarse = ss.count_grid_capacity(full_mask(500, 500), 60.0, 60.0, 300.0, HEX)
        self.assertAlmostEqual(fine / coarse, 1.0, delta=0.02)


class TestManySites(unittest.TestCase):
    """
    A distributed search may select far more than 255 sites.

    The visualisation labelling used a uint8 array, so assigning the 256th site's
    colour raised OverflowError and took the whole run down after the physics had
    already been paid for. A layout of many small sub-arrays -- which is what phase 2
    needs for TAMBO's long strip -- reaches that ceiling routinely.
    """

    def setUp(self):
        import tempfile
        self.tmp = tempfile.mkdtemp(prefix="many_sites_")

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def blob_grid(self, n, per_side):
        """A grid of separated squares, each its own connected region."""
        mask = np.zeros((n, n), dtype=bool)
        pitch = n // per_side
        size = max(1, int(pitch * 0.5))
        for i in range(per_side):
            for j in range(per_side):
                mask[i * pitch:i * pitch + size, j * pitch:j * pitch + size] = True
        return mask

    def analyse(self, mask):
        import os
        path_A = os.path.join(self.tmp, "A.npy")
        m = np.lib.format.open_memmap(path_A, mode="w+", shape=mask.shape, dtype=bool)
        m[:] = mask
        m.flush()
        del m
        elevation = np.full(mask.shape, 4000.0, dtype=np.float32)
        rows, cols = mask.shape
        return ss.analyze_sites_and_capacity(
            path_A, elevation, rows, cols, 30.0, 30.0, 1, "distributed",
            target_antennas=1, min_sub_array_size=1, antenna_spacing_km=0.06,
            grid_type="hex")

    def test_selects_more_than_255_sites_without_overflowing(self):
        mask = self.blob_grid(600, 18)
        small_final, labeled_viz, site_details, capacity, count, _ = self.analyse(mask)
        self.assertGreater(count, 255, "test needs to cross the uint8 boundary to be meaningful")
        self.assertEqual(len(site_details), count)
        # Every selected site gets its own colour, and none of them wrapped to 0
        self.assertEqual(int(labeled_viz.max()), count)
        self.assertEqual(len(np.unique(labeled_viz)), count + 1)   # + background

    def test_the_selection_mask_agrees_with_the_labelling(self):
        mask = self.blob_grid(600, 18)
        small_final, labeled_viz, _, _, count, _ = self.analyse(mask)
        np.testing.assert_array_equal(small_final.astype(bool), labeled_viz > 0)


if __name__ == "__main__":
    unittest.main()
