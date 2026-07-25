"""
lamina_mechanics.py - Stiffness and compliance matrices for composite plies.

Covers the full material symmetry hierarchy:
    Isotropic -> Transversely Isotropic -> Orthotropic -> Anisotropic (monoclinic)

Reference: Todd Coburn "Composite Strength" lectures L-02 - L-04
           Jones, R.M. "Mechanics of Composite Materials", 2nd ed., 2.2-2.4

Material symmetry classes (plane stress, 3 independent constants shown):
    Isotropic   : E, nu   -> Q11=Q22, Q66=(Q11-Q12)/2, Q16=Q26=0
    Orthotropic : E1,E2,G12,nu12 -> Q16=Q26=0 in material axes
    Anisotropic : 6 independent plane-stress constants -> Q16!=0, Q26!=0
"""

__all__ = ['compute_Q_isotropic', 'compute_S_isotropic', 'compute_Q_matrix', 'compute_S_matrix', 'engineering_constants_from_Q', 'compute_Q_anisotropic', 'compliance_from_anisotropic_Q', 'plane_stress_reduction']

import numpy as np


# -----------------------------------------------------------------------------
# ISOTROPIC MATERIAL (L-02, Jones 2.2)
# -----------------------------------------------------------------------------

def compute_Q_isotropic(E: float, v: float) -> np.ndarray:
    """
    Reduced stiffness matrix for an isotropic material (plane stress).

    Theory (L-02):
        Isotropic: E1=E2=E, G12=E/(2(1+nu)), nu12=nu.
        Q11 = Q22 = E / (1-nu^2)
        Q12 = nu*E / (1-nu^2)
        Q66 = G = E / (2(1+nu))   <- also equals (Q11-Q12)/2
        Q16 = Q26 = 0

    Parameters
    ----------
    E : float - Young's modulus [Pa]
    v : float - Poisson's ratio (dimensionless)

    Returns
    -------
    np.ndarray (3x3) - Reduced stiffness [Pa] in Voigt notation [sig1,sig2,tau12]
    """
    denom = 1.0 - v**2
    Q = np.zeros((3, 3))
    Q[0, 0] = Q[1, 1] = E / denom
    Q[0, 1] = Q[1, 0] = v * E / denom
    Q[2, 2] = E / (2.0 * (1.0 + v))
    return Q


def compute_S_isotropic(E: float, v: float) -> np.ndarray:
    """
    Compliance matrix for an isotropic material (plane stress).

    S11 = S22 = 1/E,  S12 = -nu/E,  S66 = 1/G = 2(1+nu)/E
    """
    S = np.zeros((3, 3))
    S[0, 0] = S[1, 1] = 1.0 / E
    S[0, 1] = S[1, 0] = -v / E
    S[2, 2] = 2.0 * (1.0 + v) / E
    return S


# -----------------------------------------------------------------------------
# ORTHOTROPIC MATERIAL - standard CLT case (L-03, L-04, Jones 2.4)
# -----------------------------------------------------------------------------

def compute_Q_matrix(E1: float, E2: float, G12: float, v12: float) -> np.ndarray:
    """
    Reduced stiffness matrix Q for an orthotropic ply (plane stress, material axes).

    Theory (L-04, Jones 2.4, Eq. 2.71-2.73):
        Reciprocal relation: nu21 = nu12 * E2/E1
        Denominator:         D = 1 - nu12*nu21
        Q11 = E1/D,  Q22 = E2/D,  Q12 = nu12*E2/D,  Q66 = G12
        Q16 = Q26 = 0  (no shear-extension coupling in material axes)

    Note on Voigt notation:
        Stress vector  : [sig1, sig2, tau12]
        Strain vector  : [eps1, eps2, gamma12]  (engineering shear strain)
        Relation       : {sig} = [Q]{eps}

    Parameters
    ----------
    E1   : float - Longitudinal (fibre-direction) modulus [Pa]
    E2   : float - Transverse modulus [Pa]
    G12  : float - In-plane shear modulus [Pa]
    v12  : float - Major Poisson's ratio (strain in 2 per unit strain in 1)

    Returns
    -------
    np.ndarray (3x3) - Reduced stiffness [Pa]

    Example (Jones Table 2.2, T300/5208):
        E1=181 GPa, E2=10.3 GPa, G12=7.17 GPa, nu12=0.28
        -> Q11~=181.6 GPa, Q22~=10.35 GPa, Q12~=2.897 GPa, Q66=7.17 GPa
        (Jones uses slightly different rounding; values match within <=0.3%)
    """
    v21 = v12 * E2 / E1
    denom = 1.0 - v12 * v21
    Q = np.zeros((3, 3))
    Q[0, 0] = E1 / denom
    Q[1, 1] = E2 / denom
    Q[0, 1] = Q[1, 0] = v12 * E2 / denom
    Q[2, 2] = G12
    return Q


def compute_S_matrix(E1: float, E2: float, G12: float, v12: float) -> np.ndarray:
    """
    Compliance matrix S for an orthotropic ply (plane stress, material axes).

    Theory (L-03, Jones 2.3):
        S11=1/E1, S22=1/E2, S12=-nu12/E1=-nu21/E2, S66=1/G12
        S*Q = I  (verified in validation below)

    Parameters - same as compute_Q_matrix.

    Returns
    -------
    np.ndarray (3x3) - Compliance [Pa^-^1]
    """
    v21 = v12 * E2 / E1
    S = np.zeros((3, 3))
    S[0, 0] = 1.0 / E1
    S[1, 1] = 1.0 / E2
    S[0, 1] = S[1, 0] = -v12 / E1   # = -v21/E2 by reciprocity
    S[2, 2] = 1.0 / G12
    return S


def engineering_constants_from_Q(Q: np.ndarray) -> dict:
    """
    Recover orthotropic engineering constants from a reduced stiffness matrix.

    Theory (Jones 2.4, inverse of compute_Q_matrix):
        S = Q^-^1
        E1  = 1/S11,  E2  = 1/S22
        G12 = 1/S66,  nu12 = -S12/S11

    Parameters
    ----------
    Q : np.ndarray (3x3) - Reduced stiffness matrix [Pa]

    Returns
    -------
    dict: E1, E2, G12, v12 [Pa and dimensionless]
    """
    S = np.linalg.inv(Q)
    return {
        "E1":  1.0 / S[0, 0],
        "E2":  1.0 / S[1, 1],
        "G12": 1.0 / S[2, 2],
        "v12": -S[0, 1] / S[0, 0],
    }


# -----------------------------------------------------------------------------
# ANISOTROPIC (MONOCLINIC) MATERIAL (L-04, Jones 2.5)
# -----------------------------------------------------------------------------

def compute_Q_anisotropic(
    Q11: float, Q12: float, Q22: float,
    Q16: float, Q26: float, Q66: float,
) -> np.ndarray:
    """
    Reduced stiffness for a generally anisotropic (monoclinic) ply,
    plane stress.

    Theory (L-04, Jones 2.5):
        For a material without any symmetry planes aligned with the reference
        axes (e.g., an orthotropic ply at an arbitrary angle already expressed
        in the rotated system), Q16 and Q26 are non-zero.  This is the most
        general plane-stress case - 6 independent constants.

        In CLT, Q (Q-bar) of any ply is always of this anisotropic form.
        This function lets you build that matrix directly from measured or
        pre-computed components.

    Parameters
    ----------
    Q11, Q12, Q22, Q16, Q26, Q66 : float - Stiffness components [Pa]
        Symmetry: Q21=Q12, Q61=Q16, Q62=Q26.

    Returns
    -------
    np.ndarray (3x3) - Full anisotropic reduced stiffness [Pa]
    """
    return np.array([
        [Q11, Q12, Q16],
        [Q12, Q22, Q26],
        [Q16, Q26, Q66],
    ], dtype=float)


def compliance_from_anisotropic_Q(Q: np.ndarray) -> dict:
    """
    Engineering constants for a generally anisotropic plane-stress material.

    Theory (Jones 2.5):
        S = Q^-^1
        Ex   = 1/S11,  Ey   = 1/S22,  Gxy  = 1/S66
        nuxy  = -S12/S11
        eta_xs = S16/S66  - shear-coupling coefficient (extension-shear coupling)
        eta_ys = S26/S66  - shear-coupling coefficient

    Returns
    -------
    dict: Ex, Ey, Gxy, v_xy, eta_xs, eta_ys
    """
    S = np.linalg.inv(Q)
    return {
        "Ex":    1.0 / S[0, 0],
        "Ey":    1.0 / S[1, 1],
        "Gxy":   1.0 / S[2, 2],
        "v_xy":  -S[0, 1] / S[0, 0],
        "eta_xs": S[0, 2] / S[2, 2],  # shear-extension coupling (x)
        "eta_ys": S[1, 2] / S[2, 2],  # shear-extension coupling (y)
    }


# -----------------------------------------------------------------------------
# 3D -> PLANE STRESS REDUCTION (L-03, Jones 2.3)
# -----------------------------------------------------------------------------

def plane_stress_reduction(C3d: np.ndarray) -> np.ndarray:
    """
    Reduce a 3D orthotropic stiffness matrix C (6x6, Voigt) to the plane-stress
    reduced stiffness Q (3x3).

    Theory (L-03, Jones 2.3):
        Plane-stress assumption: sig3=tau13=tau23=0.
        Solve for eps3 from row 3: eps3 = -(C13eps1 + C23eps2)/C33
        Substitute back -> reduced stiffness:
            Q11 = C11 - C13^2/C33
            Q22 = C22 - C23^2/C33
            Q12 = C12 - C13*C23/C33
            Q66 = C66  (shear unchanged)

    Parameters
    ----------
    C3d : np.ndarray (6x6) - 3D stiffness matrix [Pa] in Voigt order
                             [11,22,33,23,13,12]

    Returns
    -------
    np.ndarray (3x3) - Plane-stress Q [Pa]
    """
    C11, C12, C13 = C3d[0, 0], C3d[0, 1], C3d[0, 2]
    C22, C23      = C3d[1, 1], C3d[1, 2]
    C33, C66      = C3d[2, 2], C3d[5, 5]
    Q = np.zeros((3, 3))
    Q[0, 0] = C11 - C13**2 / C33
    Q[1, 1] = C22 - C23**2 / C33
    Q[0, 1] = Q[1, 0] = C12 - C13 * C23 / C33
    Q[2, 2] = C66
    return Q


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    print("=== lamina_mechanics.py validation ===\n")

    # 1) Orthotropic - T300/5208 (Jones Table 2.2)
    E1, E2, G12, v12 = 181e9, 10.3e9, 7.17e9, 0.28
    Q  = compute_Q_matrix(E1, E2, G12, v12)
    S  = compute_S_matrix(E1, E2, G12, v12)
    print("Orthotropic Q [GPa]:")
    print(np.round(Q / 1e9, 4))
    I  = S @ Q
    assert np.allclose(I, np.eye(3), atol=1e-10), "S*Q must equal I"
    print("S.Q = I  [OK]")

    # Recover constants from Q
    ec = engineering_constants_from_Q(Q)
    assert abs(ec["E1"] - E1) / E1 < 1e-10
    assert abs(ec["v12"] - v12) / v12 < 1e-10
    print("engineering_constants_from_Q [OK]")

    # 2) Isotropic - aluminium
    E_al, v_al = 70e9, 0.33
    Q_iso = compute_Q_isotropic(E_al, v_al)
    S_iso = compute_S_isotropic(E_al, v_al)
    print(f"\nIsotropic Q11={Q_iso[0,0]/1e9:.3f} GPa, Q66={Q_iso[2,2]/1e9:.3f} GPa")
    assert np.allclose(Q_iso[0, 0], Q_iso[1, 1]), "Isotropic: Q11=Q22"
    assert np.allclose(Q_iso[2, 2], (Q_iso[0,0]-Q_iso[0,1])/2, rtol=1e-10), "Q66=(Q11-Q12)/2"
    assert np.allclose(S_iso @ Q_iso, np.eye(3), atol=1e-10), "S*Q=I for isotropic"
    print("Isotropic symmetry checks [OK]")

    # 3) Anisotropic - build Q-bar of 30deg ply, treat as anisotropic
    import sys; sys.path.insert(0, ".")
    from transformations import transform_Q
    Qbar = transform_Q(Q, 30)
    ec2 = compliance_from_anisotropic_Q(Qbar)
    print(f"\nAnisotropic (Q-bar 30deg): Ex={ec2['Ex']/1e9:.2f} GPa, "
          f"eta_xs={ec2['eta_xs']:.4f}")
    assert abs(ec2["eta_xs"]) > 0, "30deg ply must have shear-extension coupling"
    print("Anisotropic coupling [OK]")

    # 4) Plane-stress reduction from 3D
    # Build 3D C matrix for orthotropic material (simplified)
    v13, v23 = 0.28, 0.40
    E3, G13, G23 = E2, G12, 3.4e9
    v31, v32 = v13 * E3 / E1, v23 * E3 / E2
    dD = 1 - v12*v12*E2/E1 - v23*v23*E3/E2 - v13*v13*E3/E1 - 2*v12*v23*v13*E3/E1
    C3d = np.zeros((6, 6))
    C3d[0,0]=(1-v23*v32)*E1/dD; C3d[1,1]=(1-v13*v31)*E2/dD; C3d[2,2]=(1-v12*v12*E2/E1)*E3/dD
    C3d[0,1]=C3d[1,0]=(v12+v13*v32)*E2/dD; C3d[0,2]=C3d[2,0]=(v13+v12*v23)*E3/dD
    v21_=v12*E2/E1; C3d[1,2]=C3d[2,1]=(v23+v21_*v13)*E3/dD
    C3d[3,3]=G23; C3d[4,4]=G13; C3d[5,5]=G12
    Q_ps = plane_stress_reduction(C3d)
    # Compare with direct Q: should match Q11,Q22,Q12,Q66 closely
    assert abs(Q_ps[0,0] - Q[0,0]) / Q[0,0] < 0.02, "Plane-stress Q11 mismatch"
    print("plane_stress_reduction [OK]")

    print("\nAll checks passed [OK]")
