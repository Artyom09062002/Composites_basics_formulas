import unittest
import sys
from pathlib import Path
import numpy as np

STUDY_ROOT = Path(__file__).resolve().parents[1]
PHASE4_ROOT = Path(__file__).resolve().parents[3]
REPO_ROOT = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(PHASE4_ROOT / "studies" / "du35_rbs" / "code"))
sys.path.insert(0, str(STUDY_ROOT / "code"))

from du35_sandwich_assessment import MATERIAL, build_panels, load_cases
from du35_sandwich_geometry import LE_UPPER, TE_UPPER, normalized_upper_position
from composite_physics.sandwich_panel import SandwichMaterial, SandwichPanel, assess_sandwich_panel


class SandwichPanelTests(unittest.TestCase):
    def test_face_core_face_stack_is_symmetric(self):
        panel = build_panels()[LE_UPPER.name]
        np.testing.assert_allclose(panel.stiffness().B, 0.0, atol=1e-8)

    def test_thicker_core_increases_bending_stiffness(self):
        thin = SandwichPanel("thin", .8, .003, .010, MATERIAL)
        thick = SandwichPanel("thick", .8, .003, .025, MATERIAL)
        self.assertGreater(thick.stiffness().D[0, 0], thin.stiffness().D[0, 0])

    def test_shear_correction_is_below_clpt_result(self):
        panel = build_panels()[LE_UPPER.name]
        result = assess_sandwich_panel(panel, 1.0, 1e5)
        self.assertGreater(result.global_shear_corrected_Nx_N_per_m, 0)
        self.assertLess(result.global_shear_corrected_Nx_N_per_m, result.global_clpt_Nx_N_per_m)

    def test_hexcel_wrinkling_and_crimping_examples(self):
        # Values used in Hexcel's worked formula examples: Ef=70 GPa,
        # Ec=1 GPa, Gc=220 MPa, tc=25.4 mm, b=0.5 m.
        face = {"E1": 70e9, "E2": 10e9, "G12": 4e9, "v12": .3}
        core = {"E1": 1e9, "E2": 1e9, "G12": 220e6, "v12": .3}
        material = SandwichMaterial(face, core, 500e6, 1e9, 220e6)
        panel = SandwichPanel("published_example", .5, .0005, .0254, material)
        result = assess_sandwich_panel(panel, 1.0, 1e4)
        limits = {item.mode: item.critical_Nx_N_per_m for item in result.limits}
        wrinkle_stress = limits["face_wrinkling"] / (2 * panel.face_thickness_m)
        self.assertTrue(np.isclose(wrinkle_stress, 0.5 * (220e6 * 1e9 * 70e9) ** (1/3)))
        self.assertTrue(np.isclose(limits["core_shear_crimping"] * panel.width_m, 2.794e6))

    def test_two_load_paths_are_distinct_and_beam_is_lower(self):
        _, cases = load_cases()
        for panel_cases in cases.values():
            self.assertGreater(panel_cases["beam_mapped"], 0)
            self.assertLess(panel_cases["beam_mapped"], panel_cases["109MPa_reference"])

    def test_freecad_upper_panel_measurements_are_locked(self):
        self.assertTrue(np.isclose(LE_UPPER.midsurface_width_m, 0.8024399149355272))
        self.assertTrue(np.isclose(TE_UPPER.midsurface_width_m, 1.3879939402728796))
        self.assertTrue(0 < normalized_upper_position(LE_UPPER) < 1)
        self.assertTrue(0 < normalized_upper_position(TE_UPPER) < 1)


if __name__ == "__main__":
    unittest.main()
