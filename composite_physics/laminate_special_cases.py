"""Simple classification helpers for laminate plate formulae."""
import numpy as np


def is_specially_orthotropic(A: np.ndarray, B: np.ndarray, tol: float = 1e-8, *, D: np.ndarray | None = None) -> bool:
    """Return whether coupling terms excluded by the simple formula vanish."""
    if abs(A[0, 2]) > tol or abs(A[1, 2]) > tol or np.max(np.abs(B)) > tol:
        return False
    return D is None or (abs(D[0, 2]) <= tol and abs(D[1, 2]) <= tol)
