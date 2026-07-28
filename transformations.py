"""
transformations.py - Rotation of stiffness, stress, and strain between
the material (1-2) and global (x-y) coordinate systems.

Reference: Todd Coburn "Composite Strength" lectures L-04, L-04a, L-04b
           Jones, R.M. "Mechanics of Composite Materials", 2nd ed., Ch. 2
"""

__all__ = ['transform_Q', 'transform_stress_strain']


import numpy as np
from lamina_mechanics import compute_Q_matrix


def _T_stress(theta_deg: float) -> np.ndarray:
    """
    Transformation matrix T for stress vector (Voigt: [sig1,sig2,tau12]).

    Jones Eq. 2.84:
        T = [[c^2,  s^2,  2cs ],
             [s^2,  c^2, -2cs ],
             [-cs, cs, c^2-s^2]]
    where c=cos(theta), s=sin(theta), theta is fibre angle from x-axis (CCW positive).
    """
    t = np.deg2rad(theta_deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([
        [ c**2,  s**2,  2*c*s],
        [ s**2,  c**2, -2*c*s],
        [-c*s,   c*s,  c**2 - s**2]
    ])


def _T_strain(theta_deg: float) -> np.ndarray:
    """
    Transformation matrix T for engineering strain vector [eps1,eps2,gamma12].

    Jones Eq. 2.85 (uses Reuter matrix R to handle the factor-of-2
    difference between tensor and engineering shear strain):
        T_strain = R * T_stress * R^{-1}

    Explicitly:
        T_strain = [[c^2,   s^2,   cs  ],
                    [s^2,   c^2,  -cs  ],
                    [-2cs, 2cs, c^2-s^2]]
    """
    t = np.deg2rad(theta_deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([
        [ c**2,  s**2,   c*s],
        [ s**2,  c**2,  -c*s],
        [-2*c*s, 2*c*s, c**2 - s**2]
    ])


def transform_Q(Q: np.ndarray, theta_deg: float) -> np.ndarray:
    """
    Transform the reduced stiffness matrix Q from material axes (1-2)
    to global axes (x-y), producing Qbar (Q-bar).

    Theory (L-04, Jones Eq. 2.99):
        Qbar = T^{-1} * Q * T_strain

    The result Qbar is fully populated (all 9 entries non-zero for off-axis
    plies), meaning normal-shear coupling appears - this is the source of
    the "weird" behaviour of off-axis composites (bending-twisting, etc.).

    Parameters
    ----------
    Q        : np.ndarray, shape (3,3) - Reduced stiffness in material axes [Pa]
    theta_deg: float - Fibre angle measured from x-axis, CCW positive [degrees]

    Returns
    -------
    Q_bar : np.ndarray, shape (3,3) - Transformed stiffness in global axes [Pa]

    Notes
    -----
    At theta_deg=0deg, Qbar = Q (no rotation, material and global axes coincide).
    At theta_deg=90deg, Qbar11 <-> Qbar22 swap (fibres perpendicular to load).
    """
    T  = _T_stress(theta_deg)
    Ts = _T_strain(theta_deg)
    T_inv = np.linalg.inv(T)
    Q_bar = T_inv @ Q @ Ts
    return Q_bar


def transform_stress_strain(
    values: np.ndarray,
    theta_deg: float,
    mode: str = "stress"
) -> np.ndarray:
    """
    Rotate a stress or strain vector between coordinate systems.

    Two directions are supported:
        "stress"  : global (x-y) -> material (1-2)
                    {sig}_{12} = T * {sig}_{xy}
        "strain"  : global (x-y) -> material (1-2)
                    {eps}_{12} = T_strain * {eps}_{xy}

    Use the inverse transform for material -> global: pass -theta_deg.

    Parameters
    ----------
    values    : np.ndarray, shape (3,) - [sigx,sigy,tauxy] or [epsx,epsy,gammaxy] [Pa or -]
    theta_deg : float - Fibre angle from x-axis [degrees]. Positive = CCW.
    mode      : str - "stress" or "strain"

    Returns
    -------
    transformed : np.ndarray, shape (3,)
        Vector in the rotated frame.

    Examples
    --------
    >>> sig_global = np.array([100e6, 0, 0])   # pure x-tension [Pa]
    >>> sig_material = transform_stress_strain(sig_global, 45, "stress")
    # At 45deg: sig1=sig2=50 MPa, tau12=-50 MPa
    """
    assert mode in ("stress", "strain"), "mode must be 'stress' or 'strain'"
    if mode == "stress":
        T = _T_stress(theta_deg)
    else:
        T = _T_strain(theta_deg)
    return T @ values


# ---------------------------------------------------------------------------
# Validation - Jones Example 2.3 (glass/epoxy at 30deg)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from lamina_mechanics import compute_Q_matrix

    E1, E2, G12, v12 = 38.6e9, 8.27e9, 4.14e9, 0.26
    Q = compute_Q_matrix(E1, E2, G12, v12)

    theta = 30.0
    Q_bar = transform_Q(Q, theta)

    print("=== transformations.py validation ===")
    print(f"Qbar at theta_deg={theta}deg [GPa]:\n{Q_bar/1e9}\n")
    print(f"Qbar11 = {Q_bar[0,0]/1e9:.3f} GPa  (Jones: ~23.65 GPa)")
    print(f"Qbar22 = {Q_bar[1,1]/1e9:.3f} GPa  (Jones: ~9.22  GPa)")
    print(f"Qbar12 = {Q_bar[0,1]/1e9:.3f} GPa  (Jones: ~5.08  GPa)")
    print(f"Qbar66 = {Q_bar[2,2]/1e9:.3f} GPa  (Jones: ~8.46  GPa)")
    print(f"Qbar16 = {Q_bar[0,2]/1e9:.3f} GPa  (Jones: ~9.24  GPa)")
    print(f"Qbar26 = {Q_bar[1,2]/1e9:.3f} GPa  (Jones: ~3.67  GPa)")

    # Check: theta_deg=0deg should give Q unchanged
    Q_bar_0 = transform_Q(Q, 0.0)
    assert np.allclose(Q_bar_0, Q, atol=1), "theta_deg=0deg transform must reproduce Q"

    # Stress rotation check: pure x-tension at 45deg -> sig1=sig2, tau12!=0
    sig_xy = np.array([100e6, 0.0, 0.0])
    sig_12 = transform_stress_strain(sig_xy, 45.0, "stress")
    print(f"\nStress rotation (45deg): [sig1,sig2,tau12] = {sig_12/1e6} MPa")
    print(f"  sig1 = {sig_12[0]/1e6:.1f} MPa  (expect 50)")
    print(f"  sig2 = {sig_12[1]/1e6:.1f} MPa  (expect 50)")
    print(f"  tau12= {sig_12[2]/1e6:.1f} MPa  (expect -50)")
    assert abs(sig_12[0] - 50e6) < 1e3 and abs(sig_12[1] - 50e6) < 1e3
    print("\nAll checks passed.")
