"""Laminate stiffness assembly for Classical Lamination Theory (CLT).

This module is intentionally limited to the step

    material -> Q -> ply angle -> Q-bar -> ply z-position -> A, B, D


References
----------
Todd Coburn, "Composite Strength", lecture L-08.
R. M. Jones, "Mechanics of Composite Materials", 2nd ed., Chapter 4.
"""

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from lamina_mechanics import compute_Q_matrix
from transformations import transform_Q

__all__ = [
    "LaminateStiffness",
    "assemble_laminate_stiffness",
    "compute_ABD_blocks",
    "compute_z_interfaces",
]


@dataclass(frozen=True)
class LaminateStiffness:
    """CLT stiffness result for one laminate.

    Attributes
    ----------
    A : ndarray, shape (3, 3)
        Extensional stiffness [N/m].
    B : ndarray, shape (3, 3)
        Extension-bending coupling stiffness [N].
    D : ndarray, shape (3, 3)
        Bending stiffness [N*m].
    ABD : ndarray, shape (6, 6)
        Combined laminate stiffness matrix.
    z : ndarray, shape (n_plies + 1,)
        Ply-interface coordinates [m], from ``-h/2`` to ``+h/2``.
    qbars : tuple of ndarray
        Transformed reduced stiffness matrix for each ply [Pa].
    """

    A: np.ndarray
    B: np.ndarray
    D: np.ndarray
    ABD: np.ndarray
    z: np.ndarray
    qbars: tuple[np.ndarray, ...]


def _validate_inputs(layup: list, materials: list) -> None:
    if not isinstance(layup, list) or not layup:
        raise ValueError("layup must be a non-empty list of ply dictionaries")
    if not isinstance(materials, list) or not materials:
        raise ValueError("materials must be a non-empty list of dictionaries")

    for index, material in enumerate(materials):
        missing = {"E1", "E2", "G12", "v12"} - material.keys()
        if missing:
            raise ValueError(f"material {index} is missing keys: {sorted(missing)}")
        for key in ("E1", "E2", "G12"):
            value = material[key]
            if not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"material {index} {key} must be finite and positive")
        v12 = material["v12"]
        if not isinstance(v12, Real) or not np.isfinite(v12):
            raise ValueError(f"material {index} v12 must be finite")
        v21 = v12 * material["E2"] / material["E1"]
        if 1.0 - v12 * v21 <= 0:
            raise ValueError(f"material {index} gives a non-positive plane-stress denominator")

    for index, ply in enumerate(layup):
        missing = {"theta", "t", "mat"} - ply.keys()
        if missing:
            raise ValueError(f"ply {index} is missing keys: {sorted(missing)}")
        theta = ply["theta"]
        thickness = ply["t"]
        material_index = ply["mat"]
        if not isinstance(theta, Real) or not np.isfinite(theta):
            raise ValueError(f"ply {index} theta must be finite")
        if not isinstance(thickness, Real) or not np.isfinite(thickness) or thickness <= 0:
            raise ValueError(f"ply {index} thickness must be finite and positive")
        if isinstance(material_index, bool) or not isinstance(material_index, Integral):
            raise ValueError(f"ply {index} mat must be an integer material index")
        if not 0 <= material_index < len(materials):
            raise ValueError(f"ply {index} mat index {material_index} is out of range")


def compute_z_interfaces(layup: list) -> np.ndarray:
    """Return ply-interface coordinates measured from the laminate mid-plane.

    Plies are traversed in the order supplied, from the ``-h/2`` face to the
    ``+h/2`` face.  Reversing the list therefore reverses the laminate through
    its thickness: A and D remain unchanged, while B changes sign.
    """

    if not isinstance(layup, list) or not layup:
        raise ValueError("layup must be a non-empty list of ply dictionaries")

    thicknesses = []
    for index, ply in enumerate(layup):
        if "t" not in ply:
            raise ValueError(f"ply {index} is missing key: 't'")
        thickness = ply["t"]
        if not isinstance(thickness, Real) or not np.isfinite(thickness) or thickness <= 0:
            raise ValueError(f"ply {index} thickness must be finite and positive")
        thicknesses.append(float(thickness))

    total_thickness = sum(thicknesses)
    z = np.empty(len(layup) + 1, dtype=float)
    z[0] = -total_thickness / 2.0
    z[1:] = z[0] + np.cumsum(thicknesses)
    z[-1] = total_thickness / 2.0  # eliminate accumulated round-off
    return z


def assemble_laminate_stiffness(layup: list, materials: list) -> LaminateStiffness:
    """Assemble A, B and D for a general multilayer laminate.

    Parameters
    ----------
    layup : list of dict
        One dictionary per ply with ``theta`` [degrees], ``t`` [m], and
        ``mat`` (index into ``materials``).  Order is from the ``-h/2`` face
        to the ``+h/2`` face.
    materials : list of dict
        Orthotropic plane-stress properties ``E1``, ``E2``, ``G12`` [Pa] and
        ``v12`` [-].

    Returns
    -------
    LaminateStiffness
        The A [N/m], B [N], D [N*m], combined ABD, interface coordinates and
        per-ply Q-bar matrices.

    Notes
    -----
    The through-thickness integrations are

    ``A_ij = sum(Qbar_ij * (z_k - z_{k-1}))``

    ``B_ij = 1/2 sum(Qbar_ij * (z_k^2 - z_{k-1}^2))``

    ``D_ij = 1/3 sum(Qbar_ij * (z_k^3 - z_{k-1}^3))``.
    """

    _validate_inputs(layup, materials)
    z = compute_z_interfaces(layup)

    A = np.zeros((3, 3), dtype=float)
    B = np.zeros((3, 3), dtype=float)
    D = np.zeros((3, 3), dtype=float)
    qbars = []

    for index, ply in enumerate(layup):
        material = materials[ply["mat"]]
        Q = compute_Q_matrix(
            material["E1"], material["E2"], material["G12"], material["v12"]
        )
        Qbar = transform_Q(Q, ply["theta"])
        z0, z1 = z[index], z[index + 1]
        A += Qbar * (z1 - z0)
        B += 0.5 * Qbar * (z1**2 - z0**2)
        D += (1.0 / 3.0) * Qbar * (z1**3 - z0**3)
        qbars.append(Qbar)

    # Floating-point trigonometry can leave tiny anti-symmetric terms.  The
    # constitutive matrices are symmetric by energy reciprocity, so remove
    # only that numerical noise explicitly.
    A = 0.5 * (A + A.T)
    B = 0.5 * (B + B.T)
    D = 0.5 * (D + D.T)

    ABD = np.block([[A, B], [B, D]])
    return LaminateStiffness(A=A, B=B, D=D, ABD=ABD, z=z, qbars=tuple(qbars))


def compute_ABD_blocks(layup: list, materials: list) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(A, B, D)`` as a compact interface for teaching examples."""

    stiffness = assemble_laminate_stiffness(layup, materials)
    return stiffness.A, stiffness.B, stiffness.D


if __name__ == "__main__":
    material = {"E1": 41.8e9, "E2": 14.0e9, "G12": 2.63e9, "v12": 0.28}
    ply_thickness = 0.375e-3
    angles = [0, 45, -45, 90, 90, -45, 45, 0]
    layup = [{"theta": angle, "t": ply_thickness, "mat": 0} for angle in angles]
    result = assemble_laminate_stiffness(layup, [material])

    print("=== laminate.py validation: [0/+45/-45/90]s ===")
    print("A [MN/m]:\n", result.A / 1e6)
    print("B [N]:\n", result.B)
    print("D [N*m]:\n", result.D)
    assert np.allclose(result.B, 0.0, atol=1e-9)
    assert np.allclose(result.A, result.A.T)
    assert np.allclose(result.D, result.D.T)
    print("Symmetry and B = 0 checks passed [OK]")