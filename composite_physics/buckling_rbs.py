"""Simply-supported orthotropic plate buckling using reduced bending stiffness.

This is an RBS screening calculation, not a full ABD coupled-buckling solver.
It applies only when the resulting ``D*16`` and ``D*26`` are negligible.
"""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class BucklingResult:
    Nx_cr_N_per_m: float
    m: int
    n: int
    a_m: float
    b_m: float
    margin_of_safety: float | None


def critical_uniaxial_compression_rbs(
    D_star: np.ndarray, a_m: float, b_m: float, applied_Nx_N_per_m: float | None = None,
    max_mode: int = 12, coupling_tolerance: float = 1e-7,
) -> BucklingResult:
    """Minimise the Navier solution for simply-supported compression in x."""
    stiffness = np.asarray(D_star, dtype=float)
    if stiffness.shape != (3, 3) or not np.all(np.isfinite(stiffness)):
        raise ValueError("D_star must be a finite 3x3 matrix")
    if a_m <= 0 or b_m <= 0 or max_mode < 1:
        raise ValueError("a_m, b_m and max_mode must be positive")
    if abs(stiffness[0, 2]) > coupling_tolerance or abs(stiffness[1, 2]) > coupling_tolerance:
        raise ValueError("D*16/D*26 are non-zero: this orthotropic Navier formula is not applicable")

    best: tuple[float, int, int] | None = None
    for m in range(1, max_mode + 1):
        kx = m * np.pi / a_m
        for n in range(1, max_mode + 1):
            ky = n * np.pi / b_m
            numerator = (stiffness[0, 0] * kx**4
                         + 2.0 * (stiffness[0, 1] + 2.0 * stiffness[2, 2]) * kx**2 * ky**2
                         + stiffness[1, 1] * ky**4)
            nx_cr = numerator / kx**2
            if best is None or nx_cr < best[0]:
                best = (float(nx_cr), m, n)
    assert best is not None
    margin = None if applied_Nx_N_per_m is None else best[0] / applied_Nx_N_per_m - 1.0
    return BucklingResult(best[0], best[1], best[2], a_m, b_m, margin)
