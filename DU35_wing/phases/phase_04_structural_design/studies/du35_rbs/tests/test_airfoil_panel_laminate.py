"""Integration tests for the DU35 panel application case."""

import unittest
import sys
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(STUDY_ROOT / "code"))

from airfoil_panel_laminate import analyze_du35_skin_demo, analyze_du35_upper_panel


class AirfoilPanelLaminateTests(unittest.TestCase):
    def test_live_cad_snapshot_gives_expected_panel_geometry(self) -> None:
        result = analyze_du35_upper_panel()
        self.assertAlmostEqual(result.projected_width_m, 0.80853714, places=8)
        self.assertAlmostEqual(result.surface_width_m, 0.8108039732531884, places=10)
        self.assertGreater(result.surface_width_m, result.projected_width_m)
        self.assertAlmostEqual(result.skin_thickness_m, 0.003, places=10)
        self.assertAlmostEqual(result.cap_thickness_m, 0.065, places=10)
        self.assertAlmostEqual(result.total_thickness_m, 0.068, places=10)

    def test_full_panel_has_membrane_bending_coupling(self) -> None:
        result = analyze_du35_upper_panel()
        self.assertGreater(np.max(np.abs(result.stiffness.B)), 1.0)
        self.assertFalse(result.compatible_with_current_buckling_model)

    def test_educational_skin_case_remains_available(self) -> None:
        result = analyze_du35_skin_demo()
        np.testing.assert_allclose(result.stiffness.B, 0.0, atol=1e-9)
        self.assertGreater(abs(result.stiffness.D[0, 2]), 1.0)
        self.assertGreater(abs(result.stiffness.D[1, 2]), 1.0)


if __name__ == "__main__":
    unittest.main()
