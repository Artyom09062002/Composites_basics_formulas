"""Regression tests for the Day-2 RBS stability screening."""
import unittest
import sys
from pathlib import Path
import numpy as np

STUDY_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(STUDY_ROOT / "code"))

from airfoil_panel_laminate import analyze_du35_upper_panel
from composite_physics.buckling_rbs import critical_uniaxial_compression_rbs
from composite_physics.panel_loads import compressive_resultant_from_moment, compressive_resultant_from_stress
from composite_physics.reduced_stiffness import compute_reduced_D


class RbsTests(unittest.TestCase):
    def test_panel_rbs_is_symmetric_and_reduced_by_coupling(self):
        panel = analyze_du35_upper_panel().stiffness
        d_star = compute_reduced_D(panel.A, panel.B, panel.D)
        np.testing.assert_allclose(d_star, d_star.T, atol=1e-8)
        self.assertLess(d_star[0, 0], panel.D[0, 0])

    def test_stress_and_moment_load_conversions_have_expected_units(self):
        self.assertAlmostEqual(compressive_resultant_from_stress(100e6, 0.01).Nx_N_per_m, 1e6)
        by_moment = compressive_resultant_from_moment(10e6, 2.0, 0.5)
        self.assertAlmostEqual(by_moment.cap_force_N, 5e6)
        self.assertAlmostEqual(by_moment.Nx_N_per_m, 10e6)

    def test_buckling_margin_is_reported(self):
        panel = analyze_du35_upper_panel()
        d_star = compute_reduced_D(panel.stiffness.A, panel.stiffness.B, panel.stiffness.D)
        result = critical_uniaxial_compression_rbs(d_star, panel.surface_width_m, panel.surface_width_m, 1e6)
        self.assertGreater(result.Nx_cr_N_per_m, 0.0)
        self.assertIsNotNone(result.margin_of_safety)


if __name__ == "__main__":
    unittest.main()
