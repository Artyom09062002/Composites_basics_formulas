"""Physics-based regression tests for laminate stiffness assembly."""

import unittest

import numpy as np

from clt import compute_ABD_matrix
from lamina_mechanics import compute_Q_matrix
from laminate import assemble_laminate_stiffness
from laminate_special_cases import is_specially_orthotropic
from reduced_stiffness import compute_reduced_D


MATERIAL = {"E1": 41.8e9, "E2": 14.0e9, "G12": 2.63e9, "v12": 0.28}
ANGLES = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]


def make_layup(thickness: float = 0.375e-3) -> list[dict]:
    return [{"theta": angle, "t": thickness, "mat": 0} for angle in ANGLES]


class LaminateAssemblyTests(unittest.TestCase):
    def test_symmetric_layup_has_zero_B_and_symmetric_A_D(self) -> None:
        result = assemble_laminate_stiffness(make_layup(), [MATERIAL])
        np.testing.assert_allclose(result.B, 0.0, atol=1e-9)
        np.testing.assert_allclose(result.A, result.A.T, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(result.D, result.D.T, rtol=0.0, atol=1e-12)

    def test_doubling_all_ply_thickness_scales_A_and_D(self) -> None:
        base = assemble_laminate_stiffness(make_layup(), [MATERIAL])
        doubled = assemble_laminate_stiffness(make_layup(0.750e-3), [MATERIAL])
        np.testing.assert_allclose(doubled.A, 2.0 * base.A, rtol=1e-12, atol=1e-8)
        np.testing.assert_allclose(doubled.D, 8.0 * base.D, rtol=1e-12, atol=1e-12)

    def test_reversing_laminate_keeps_A_D_and_flips_B(self) -> None:
        unsymmetric = [
            {"theta": 0.0, "t": 0.2e-3, "mat": 0},
            {"theta": 45.0, "t": 0.3e-3, "mat": 0},
            {"theta": 90.0, "t": 0.4e-3, "mat": 0},
        ]
        forward = assemble_laminate_stiffness(unsymmetric, [MATERIAL])
        reverse = assemble_laminate_stiffness(list(reversed(unsymmetric)), [MATERIAL])
        np.testing.assert_allclose(reverse.A, forward.A, rtol=1e-12, atol=1e-8)
        np.testing.assert_allclose(reverse.D, forward.D, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(reverse.B, -forward.B, rtol=1e-12, atol=1e-10)

    def test_all_zero_laminate_matches_closed_form(self) -> None:
        total_thickness = 2.0e-3
        layup = [{"theta": 0.0, "t": total_thickness / 4.0, "mat": 0}] * 4
        result = assemble_laminate_stiffness(layup, [MATERIAL])
        Q = compute_Q_matrix(**MATERIAL)
        np.testing.assert_allclose(result.A, Q * total_thickness, rtol=1e-12, atol=1e-8)
        np.testing.assert_allclose(result.B, 0.0, atol=1e-10)
        np.testing.assert_allclose(
            result.D, Q * total_thickness**3 / 12.0, rtol=1e-12, atol=1e-12
        )

    def test_clt_compatibility_api_uses_same_ABD(self) -> None:
        layup = make_layup()
        expected = assemble_laminate_stiffness(layup, [MATERIAL]).ABD
        actual = compute_ABD_matrix(layup, [MATERIAL])
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=0.0)

    def test_invalid_thickness_is_rejected(self) -> None:
        layup = [{"theta": 0.0, "t": 0.0, "mat": 0}]
        with self.assertRaisesRegex(ValueError, "thickness"):
            assemble_laminate_stiffness(layup, [MATERIAL])

    def test_plate_special_orthotropy_checks_D16_and_D26(self) -> None:
        result = assemble_laminate_stiffness(make_layup(), [MATERIAL])
        self.assertTrue(is_specially_orthotropic(result.A, result.B))
        self.assertFalse(is_specially_orthotropic(result.A, result.B, D=result.D))

    def test_reduced_stiffness_equals_D_for_symmetric_laminate(self) -> None:
        result = assemble_laminate_stiffness(make_layup(), [MATERIAL])
        np.testing.assert_allclose(
            compute_reduced_D(result.A, result.B, result.D), result.D, rtol=1e-12, atol=1e-12
        )


if __name__ == "__main__":
    unittest.main()
