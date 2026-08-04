"""Reduced bending stiffness (RBS) for extension--bending coupled laminates.

The Schur complement ``D* = D - B A^-1 B`` eliminates the in-plane strain
from the CLT constitutive relation.  It is an *equivalent, approximate* way to
use bending-only plate formulae for an unsymmetric laminate; it is not a full
coupled buckling eigenvalue solution.
"""

import numpy as np

__all__ = ["compute_reduced_D"]


def compute_reduced_D(A: np.ndarray, B: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Return the RBS matrix ``D*`` [N m].

    ``numpy.linalg.solve`` is used rather than forming ``A^-1`` explicitly.
    All inputs must be finite, symmetric 3x3 matrices; ``A`` must be
    nonsingular.
    """
    matrices = {"A": np.asarray(A, dtype=float), "B": np.asarray(B, dtype=float),
                "D": np.asarray(D, dtype=float)}
    for name, matrix in matrices.items():
        if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} must be a finite 3x3 matrix")
        if not np.allclose(matrix, matrix.T, rtol=1e-10, atol=1e-8):
            raise ValueError(f"{name} must be symmetric")

    reduced = matrices["D"] - matrices["B"] @ np.linalg.solve(matrices["A"], matrices["B"])
    return 0.5 * (reduced + reduced.T)
