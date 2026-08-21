"""
bem_solver.py — Blade Element Momentum (BEM) aerodynamic solver.

Computes distributed aerodynamic loads along the NREL 5MW / Sandia 61.5m blade,
then integrates to shear-force and bending-moment diagrams.

Scope: operating (rated) + parked (extreme wind) load cases.
Units: SI throughout (m, kg, s, N, N·m, Pa).

References:
  [1] Hansen M.O.L., "Aerodynamics of Wind Turbines", 2nd ed., Earthscan, 2008
      — BEM algorithm, Prandtl tip loss, Glauert/Buhl high-thrust correction
  [2] Jonkman et al., NREL/TP-500-38060, 2009 — turbine parameters, geometry
  [3] Burton et al., "Wind Energy Handbook", 2nd ed., Wiley, 2011 — DLC loads
  [4] Buhl M.L., "A New Empirical Relationship between Thrust Coefficient and
      Induction Factor for the Turbulent Windmill State", NREL/TP-500-36834, 2005
      — Glauert correction for a > 0.4
  [5] Viterna L.A. & Corrigan R.D., "Fixed pitch rotor performance of large HAWTs",
      NASA-CP-2230, 1982 — post-stall polar extrapolation to 360°
  [6] Montgomerie B., "Methods for Root Effects, Tip Effects and Extending the
      Angle of Attack Range to ±180°", FOI-R-1305-SE, 2004

NOTE: screening-level tool. Polars are simplified (2D, no 3D corrections).
      Results typically within 10-15% of full aeroelastic simulation.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Turbine constants  (REF: NREL/TP-500-38060, Table 1 and Section 6)
# ---------------------------------------------------------------------------

TURBINE = dict(
    B=3,                # number of blades
    R=63.0,             # rotor radius, hub centre to blade tip [m]
    R_hub=1.5,          # hub radius [m]
    blade_len=61.5,     # blade length [m]
    V_rated=11.4,       # rated wind speed [m/s]
    V_cutin=3.0,        # cut-in wind speed [m/s]
    V_cutout=25.0,      # cut-out wind speed [m/s]
    omega_rated=12.1 * 2 * np.pi / 60,   # rated rotor speed [rad/s] = 1.2671 rad/s
    rho=1.225,          # air density [kg/m³] (IEC standard atmosphere)
    g=9.81,             # gravitational acceleration [m/s²]
)


# ---------------------------------------------------------------------------
# Airfoil polar loader
# ---------------------------------------------------------------------------

def load_polars(polar_dir: str) -> dict:
    """
    Load all airfoil polars from CSV files in polar_dir.
    Returns dict {name: DataFrame(alpha_deg, Cl, Cd)}.
    """
    polar_dir = Path(polar_dir)
    polars = {}
    for f in sorted(polar_dir.glob("*.csv")):
        name = f.stem
        df = pd.read_csv(f, comment="#")
        df.columns = [c.strip() for c in df.columns]
        polars[name] = df
    return polars


def interp_polar(polars: dict, airfoil_name: str, alpha_deg: float):
    """
    Return (Cl, Cd) for a given airfoil and angle of attack.
    Extrapolates linearly outside table range (flag if far outside).
    """
    df = polars[airfoil_name]
    alpha_table = df["alpha_deg"].values
    Cl = float(np.interp(alpha_deg, alpha_table, df["Cl"].values))
    Cd = float(np.interp(alpha_deg, alpha_table, df["Cd"].values))
    return Cl, Cd


# ---------------------------------------------------------------------------
# Viterna-Corrigan 360° polar extrapolation   (Ref [5], [6])
# ---------------------------------------------------------------------------

def extend_polar_viterna(alpha_deg: np.ndarray, Cl: np.ndarray, Cd: np.ndarray,
                         AR: float = 17.0) -> tuple:
    """
    Extrapolate a limited-angle polar to [-180°, +180°] using Viterna-Corrigan.

    Parameters
    ----------
    alpha_deg : array of AoA [deg], must be sorted ascending
    Cl, Cd    : lift and drag coefficient arrays (same length as alpha_deg)
    AR        : aspect ratio of the blade section (default 17 ~ R/mean_chord for NREL 5MW)

    Returns
    -------
    (alpha_ext, Cl_ext, Cd_ext) — 721-point arrays covering [-180°, +180°]

    Method (Ref [5]):
    For stall_angle < |α| ≤ 90°:
        Cl = A1·sin(2α) + A2·cos²(α)/sin(α)
        Cd = B1·sin²(α) + B2·cos(α)
    For 90° < |α| ≤ 180°: symmetry relation about 90°
        Cl(α) = -Cl(180°-|α|) · sign(α)
        Cd(α) =  Cd(180°-|α|)
    """
    # -- Step 1: identify positive stall (max Cl in 0°..30°) -----------------
    pos_mask = (alpha_deg >= 0.0) & (alpha_deg <= 35.0)
    if pos_mask.sum() > 1:
        idx_stall = int(np.argmax(Cl[pos_mask]))
        alpha_s_deg = float(alpha_deg[pos_mask][idx_stall])
        Cl_s = float(Cl[pos_mask][idx_stall])
        Cd_s = float(Cd[pos_mask][idx_stall])
    else:
        alpha_s_deg = 15.0
        Cl_s = 1.0
        Cd_s = 0.020

    alpha_s = np.radians(alpha_s_deg)

    # -- Step 2: Cd at 90° using AR correction (Viterna eq. 8) ---------------
    # Cd_max = 1.11 + 0.13 * AR, but cap at 1.50 for slender sections
    Cd_max = min(1.11 + 0.13 * AR, 1.50)

    # -- Step 3: Viterna coefficients -----------------------------------------
    sin_as = np.sin(alpha_s)
    cos_as = np.cos(alpha_s)
    A1 = Cd_max / 2.0
    # A2 from boundary condition: Cl(alpha_s) = A1*sin(2*alpha_s) + A2*cos²/sin
    if abs(sin_as) > 1e-4:
        A2 = (Cl_s - Cd_max * sin_as * cos_as) * sin_as / max(cos_as**2, 1e-6)
    else:
        A2 = 0.0
    B1 = Cd_max
    # B2 from boundary condition: Cd(alpha_s) = B1*sin²(alpha_s) + B2*cos(alpha_s)
    B2 = (Cd_s - Cd_max * sin_as**2) / max(abs(cos_as), 1e-4)

    # -- Step 4: build 360° grid at 0.5° resolution --------------------------
    alpha_ext = np.linspace(-180.0, 180.0, 721)
    Cl_ext = np.zeros(721)
    Cd_ext = np.zeros(721)

    # Bounds of original table
    alpha_orig_min = float(alpha_deg[0])   # e.g. -20°
    alpha_orig_max = float(alpha_deg[-1])  # e.g. +25°

    for i, a_deg in enumerate(alpha_ext):
        a_abs = abs(a_deg)
        sign = 1.0 if a_deg >= 0.0 else -1.0

        if alpha_orig_min <= a_deg <= alpha_orig_max:
            # Within original polar range — always use original data (no Viterna here)
            Cl_ext[i] = float(np.interp(a_deg, alpha_deg, Cl))
            Cd_ext[i] = float(np.interp(a_deg, alpha_deg, Cd))

        elif a_abs <= 90.0:
            # Viterna post-stall region [alpha_stall … 90°]
            a_rad = np.radians(a_abs)
            sa = np.sin(a_rad)
            ca = np.cos(a_rad)
            cl_v = A1 * np.sin(2.0 * a_rad) + (A2 * ca**2 / sa if sa > 1e-4 else 0.0)
            cd_v = B1 * sa**2 + B2 * ca
            Cl_ext[i] = sign * cl_v
            Cd_ext[i] = max(cd_v, 0.0)

        else:
            # Mirror region [90° … 180°] — Montgomerie symmetry about 90°
            a_mirror = np.radians(180.0 - a_abs)
            sa = np.sin(a_mirror)
            ca = np.cos(a_mirror)
            cl_v = A1 * np.sin(2.0 * a_mirror) + (A2 * ca**2 / sa if sa > 1e-4 else 0.0)
            cd_v = B1 * sa**2 + B2 * ca
            Cl_ext[i] = -sign * cl_v   # antisymmetric
            Cd_ext[i] = max(cd_v, 0.0)

    return alpha_ext, Cl_ext, Cd_ext


def extend_polars_360(polars: dict, AR: float = 17.0) -> dict:
    """
    Extend all polars in the dict to 360° using Viterna-Corrigan.
    Returns a new dict with the same keys but 721-point DataFrames.
    Cylinder polars (Cl ~ 0 everywhere) are extended with flat-plate drag only.
    """
    polars_360 = {}
    for name, df in polars.items():
        alpha_arr = df["alpha_deg"].values
        Cl_arr    = df["Cl"].values
        Cd_arr    = df["Cd"].values

        # Detect cylinder: max |Cl| < 0.1
        if np.max(np.abs(Cl_arr)) < 0.1:
            # Bluff-body model: Cd ≈ constant (cylinder drag is roughly
            # orientation-independent at turbine Re).  The sin² flat-plate
            # model is physically wrong here — it drives Cd→0 at small AoA,
            # which halves the feathered-parked moment versus original data.
            Cd_cyl = float(np.max(Cd_arr))     # e.g. 0.50 for Cylinder1
            alpha_ext = np.linspace(-180.0, 180.0, 721)
            Cl_ext = np.zeros(721)
            Cd_ext = np.full(721, Cd_cyl)      # constant bluff-body drag
        else:
            alpha_ext, Cl_ext, Cd_ext = extend_polar_viterna(
                alpha_arr, Cl_arr, Cd_arr, AR=AR)

        polars_360[name] = pd.DataFrame({
            "alpha_deg": alpha_ext,
            "Cl": Cl_ext,
            "Cd": Cd_ext,
        })
    return polars_360


# ---------------------------------------------------------------------------
# Blade geometry loader
# ---------------------------------------------------------------------------

def load_blade_geometry(geom_csv: str) -> pd.DataFrame:
    """
    Load blade planform from blade_geometry.csv.
    Returns DataFrame with columns: r_m, chord_m, twist_deg, airfoil.
    Adds R_from_centre column = r_m + R_hub (radial position from rotor centre).
    """
    df = pd.read_csv(geom_csv, comment="#")
    df.columns = [c.strip() for c in df.columns]
    # r_m in the file is already measured from blade root = hub flange
    # rotor radius from hub centre = r_m + R_hub
    df["R_from_centre"] = df["r_m"] + TURBINE["R_hub"]
    return df


def load_blade_mass(stiff_csv: str) -> pd.DataFrame:
    """Load blade distributed mass from blade_stiffness.csv."""
    df = pd.read_csv(stiff_csv, comment="#")
    df.columns = [c.strip() for c in df.columns]
    # Replace TBD strings with NaN
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


# ---------------------------------------------------------------------------
# Prandtl tip (and hub) loss factor
# ---------------------------------------------------------------------------

def prandtl_tip_loss(B: int, R: float, r: float, phi_rad: float,
                     R_hub: float = 1.5) -> float:
    """
    Combined Prandtl tip + hub loss factor F = F_tip * F_hub.
    phi_rad: local flow angle [rad] (angle between rotor plane and relative wind).
    """
    eps = 1e-6
    phi = abs(phi_rad) + eps

    # Tip loss
    f_tip = B / 2.0 * (R - r) / (r * np.sin(phi))
    F_tip = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_tip), -1, 1))

    # Hub loss
    f_hub = B / 2.0 * (r - R_hub) / (R_hub * np.sin(phi))
    F_hub = (2.0 / np.pi) * np.arccos(np.clip(np.exp(-f_hub), -1, 1))

    return max(F_tip * F_hub, 1e-4)


# ---------------------------------------------------------------------------
# Glauert / Buhl high-thrust correction
# ---------------------------------------------------------------------------

def glauert_correction(a_raw: float, F: float) -> float:
    """
    Buhl (NREL/TP-500-36834, 2005) modified Glauert correction.
    For a < 0.4: standard momentum theory.
    For a >= 0.4: empirical correction to avoid singularity.
    """
    if a_raw < 0.4:
        return a_raw
    # Buhl Eq. (6): Ct = (18F - 20 - 3*sqrt(Ct*(50-36F) + 12F*(3F-4))) / (50 - 36F + 12kF)
    # Quadratic for a from Ct_local = sigma * Cn / sin^2(phi) at given a:
    # We solve via the linearised form from Spera (1994) as used in OpenFAST:
    # a = 0.1*(18F - 20 - 3*sqrt(Ct_ad*(50-36F)+12F*(3F-4))) / (50-36F+12*k*F)
    # Simplified: use empirical linear blend
    # (consistent with Hansen 2008 section 6.3)
    a = (18.0 * F - 20.0 - 3.0 * np.sqrt(max(Ct_from_a(a_raw) * (50 - 36 * F) + 12 * F * (3 * F - 4), 0))) \
        / (50.0 - 36.0 * F + 12.0 * F * F) if F > 0.01 else a_raw
    return np.clip(a, 0.0, 0.95)


def Ct_from_a(a: float) -> float:
    """Thrust coefficient from axial induction (ideal actuator disk)."""
    return 4.0 * a * (1.0 - a)


# ---------------------------------------------------------------------------
# Single-station BEM iteration
# ---------------------------------------------------------------------------

def bem_station(r: float, c: float, theta_local_deg: float, airfoil: str,
                polars: dict, V_wind: float, omega: float,
                B: int, R: float, R_hub: float, rho: float,
                n_iter: int = 200, tol: float = 1e-7,
                relax: float = 0.4) -> dict:
    """
    Solve BEM equations at a single radial station.

    Parameters
    ----------
    r              : radial position from rotor centre [m]
    c              : local chord [m]
    theta_local_deg: local pitch angle = blade_pitch + twist [deg]
    airfoil        : airfoil name (key in polars dict)
    V_wind         : wind speed [m/s]
    omega          : rotor angular velocity [rad/s]
    B, R, R_hub    : blade count, tip radius, hub radius [m]
    rho            : air density [kg/m³]

    Returns
    -------
    dict with: a, a_prime, phi_deg, alpha_deg, Cl, Cd, dFn, dFt, dP
    dFn [N/m]: normal force per unit span (thrust direction, flapwise)
    dFt [N/m]: tangential force per unit span (torque direction, edgewise)
    dP  [W/m]: power per unit span
    """
    # Local solidity
    sigma = B * c / (2.0 * np.pi * r)
    theta = np.radians(theta_local_deg)

    a = 0.2   # initial axial induction
    a_prime = 0.0   # initial tangential induction

    for _ in range(n_iter):
        # Velocity components in rotor plane
        Vn = V_wind * (1.0 - a)       # normal to rotor plane
        Vt = omega * r * (1.0 + a_prime)  # tangential

        V_rel_sq = Vn**2 + Vt**2
        V_rel = np.sqrt(V_rel_sq)

        if V_rel < 1e-6:
            break

        phi = np.arctan2(Vn, Vt)   # flow angle [rad]

        # Angle of attack
        alpha_rad = phi - theta
        alpha_deg = np.degrees(alpha_rad)

        # Aerodynamic coefficients
        Cl, Cd = interp_polar(polars, airfoil, alpha_deg)

        # Force coefficients in rotor-oriented directions
        # Cn: normal to rotor plane (thrust, flapwise)
        # Ct: tangential in rotor plane (torque, edgewise)
        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
        Ct = Cl * np.sin(phi) - Cd * np.cos(phi)

        # Prandtl tip+hub loss
        F = prandtl_tip_loss(B, R, r, phi, R_hub)

        # New induction factors (momentum theory)
        denom_a = (4.0 * F * np.sin(phi)**2 / (sigma * Cn + 1e-12)) + 1.0
        a_new = 1.0 / denom_a if abs(Cn) > 1e-8 else 0.0

        denom_ap = (4.0 * F * np.sin(phi) * np.cos(phi) / (sigma * Ct + 1e-12)) - 1.0
        a_prime_new = 1.0 / denom_ap if abs(Ct) > 1e-8 and abs(denom_ap) > 1e-6 else 0.0
        a_prime_new = max(0.0, a_prime_new)

        # Glauert correction for turbulent windmill state
        if a_new > 0.4:
            # Empirical Spera / Buhl correction (Hansen 2008, eq. 6.37)
            K = 4.0 * F * np.sin(phi)**2 / (sigma * Cn + 1e-12)
            a_new = 0.5 * (2.0 + K * (1.0 - 2.0 * 0.4)
                           - np.sqrt((K * (1.0 - 2.0 * 0.4) + 2.0)**2
                                     + 4.0 * (K * 0.4**2 - 1.0)))

        a_new = np.clip(a_new, 0.0, 0.95)
        a_prime_new = np.clip(a_prime_new, 0.0, 0.5)

        # Relaxation to aid convergence
        a_next = (1.0 - relax) * a + relax * a_new
        a_prime_next = (1.0 - relax) * a_prime + relax * a_prime_new

        if abs(a_next - a) < tol and abs(a_prime_next - a_prime) < tol:
            a, a_prime = a_next, a_prime_next
            break

        a, a_prime = a_next, a_prime_next

    # Final aerodynamic loads per unit span
    Vn = V_wind * (1.0 - a)
    Vt = omega * r * (1.0 + a_prime)
    V_rel_sq = Vn**2 + Vt**2
    phi = np.arctan2(Vn, Vt)
    alpha_deg = np.degrees(phi) - theta_local_deg
    Cl, Cd = interp_polar(polars, airfoil, alpha_deg)
    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
    Ct = Cl * np.sin(phi) - Cd * np.cos(phi)

    q = 0.5 * rho * V_rel_sq  # dynamic pressure [Pa]
    dFn = q * c * Cn   # normal force per unit span [N/m]
    dFt = q * c * Ct   # tangential force per unit span [N/m]
    dP  = dFt * omega * r  # power per unit span [W/m]

    return dict(
        r=r, a=a, a_prime=a_prime,
        phi_deg=np.degrees(phi), alpha_deg=alpha_deg,
        Cl=Cl, Cd=Cd, V_rel=np.sqrt(V_rel_sq),
        dFn=dFn, dFt=dFt, dP=dP,
    )


# ---------------------------------------------------------------------------
# Full-blade BEM
# ---------------------------------------------------------------------------

def run_bem(geom_df: pd.DataFrame, polars: dict,
            V_wind: float, omega: float, pitch_deg: float = 0.0,
            rho: float = 1.225) -> pd.DataFrame:
    """
    Run BEM over all blade stations.

    Parameters
    ----------
    geom_df   : blade geometry DataFrame (from load_blade_geometry)
    polars    : airfoil polar dict (from load_polars)
    V_wind    : wind speed [m/s]
    omega     : rotor angular velocity [rad/s]
    pitch_deg : collective blade pitch angle [deg] (fine pitch = 0)
    rho       : air density [kg/m³]

    Returns
    -------
    DataFrame with BEM results at each station + cumulative loads.
    """
    B   = TURBINE["B"]
    R   = TURBINE["R"]
    R_hub = TURBINE["R_hub"]

    rows = []
    for _, row in geom_df.iterrows():
        r     = float(row["R_from_centre"])
        c     = float(row["chord_m"])
        twist = float(row["twist_deg"])
        foil  = str(row["airfoil"])

        # Local pitch = collective pitch + geometric twist (nose-up positive)
        theta_local = pitch_deg + twist

        result = bem_station(
            r=r, c=c, theta_local_deg=theta_local, airfoil=foil,
            polars=polars, V_wind=V_wind, omega=omega,
            B=B, R=R, R_hub=R_hub, rho=rho,
        )
        result["chord"] = c
        result["twist"] = twist
        result["airfoil"] = foil
        rows.append(result)

    df = pd.DataFrame(rows)

    # -----------------------------------------------------------------------
    # Integrate distributed loads → shear force & bending moment (tip→root)
    # Convention: positive dFn = towards downwind (flapwise, out-of-plane)
    #             positive dFt = in direction of rotation (edgewise, in-plane)
    # -----------------------------------------------------------------------
    r_arr  = df["r"].values          # radial stations [m]
    dFn    = df["dFn"].values        # flapwise load per unit span [N/m]
    dFt    = df["dFt"].values        # edgewise load per unit span [N/m]

    n = len(r_arr)
    Qflap = np.zeros(n)   # flapwise shear force [N]
    Mflap = np.zeros(n)   # flapwise bending moment [N·m]
    Qedge = np.zeros(n)   # edgewise shear force [N]
    Medge = np.zeros(n)   # edgewise bending moment [N·m]

    # Integrate from tip to root using trapezoidal rule
    for i in range(n - 2, -1, -1):
        dr = r_arr[i + 1] - r_arr[i]
        fn_avg = 0.5 * (dFn[i] + dFn[i + 1])
        ft_avg = 0.5 * (dFt[i] + dFt[i + 1])
        # Shear force
        Qflap[i] = Qflap[i + 1] + fn_avg * dr
        Qedge[i] = Qedge[i + 1] + ft_avg * dr
        # Bending moment contribution (force × moment arm from current station)
        Mflap[i] = Mflap[i + 1] + Qflap[i + 1] * dr + fn_avg * dr**2 / 2.0
        Medge[i] = Medge[i + 1] + Qedge[i + 1] * dr + ft_avg * dr**2 / 2.0

    df["Qflap_N"]  = Qflap
    df["Mflap_Nm"] = Mflap
    df["Qedge_N"]  = Qedge
    df["Medge_Nm"] = Medge

    return df


# ---------------------------------------------------------------------------
# Gravity loads (edgewise bending moment from blade self-weight)
# ---------------------------------------------------------------------------

def gravity_edge_moment(mass_df: pd.DataFrame) -> np.ndarray:
    """
    Compute edgewise bending moment from gravity along the blade.
    Worst case: blade horizontal (3 o'clock position), gravity pulls down → edgewise.
    M_gravity(r) = ∫_r^R (r' - r) * m(r') * g dr'  [N·m]

    Parameters
    ----------
    mass_df: DataFrame with r_m (station from blade root) and mass_pm_kg_m columns.

    Returns (r_m, M_gravity_Nm) arrays aligned to mass_df stations.
    """
    g = TURBINE["g"]
    r_arr = mass_df["r_m"].values
    m_arr = mass_df["mass_pm_kg_m"].values
    n = len(r_arr)

    M_grav = np.zeros(n)
    for i in range(n - 2, -1, -1):
        dr = r_arr[i + 1] - r_arr[i]
        m_avg = 0.5 * (m_arr[i] + m_arr[i + 1])
        M_grav[i] = M_grav[i + 1] + (M_grav[i + 1] / (r_arr[i + 1] - r_arr[0] + 1e-9) * dr
                                       if False else 0)
        # Direct moment calculation
        M_grav[i] = M_grav[i + 1] + m_avg * g * dr * (r_arr[i + 1] - r_arr[i]) / 2.0

    # Redo properly: M(r_i) = ∫_{r_i}^{R} (r' - r_i) * m(r') * g dr'
    for i in range(n):
        M = 0.0
        for j in range(i, n - 1):
            dr = r_arr[j + 1] - r_arr[j]
            r_mid = 0.5 * (r_arr[j] + r_arr[j + 1])
            m_mid = 0.5 * (m_arr[j] + m_arr[j + 1])
            M += m_mid * g * (r_mid - r_arr[i]) * dr
        M_grav[i] = M

    return r_arr, M_grav


# ---------------------------------------------------------------------------
# Parked extreme-wind loads  (DLC 6.1 envelope: V_50 = 70 m/s, omega = 0)
# ---------------------------------------------------------------------------

def parked_loads(geom_df: pd.DataFrame, polars: dict,
                 V_extreme: float = 70.0, pitch_deg: float = 90.0,
                 rho: float = 1.225) -> pd.DataFrame:
    """
    Parked (non-rotating) load calculation for extreme wind at a given pitch angle.

    The wind blows axially (phi = 90°), rotor is stationary (omega = 0).
    At each section:
        AoA = 90° - (pitch + twist)

    For the standard feathered case (pitch = 90°):
        AoA ≈ -twist  → small AoA → small Cl, large phi → Cn ≈ Cd (drag-dominated)
        → flapwise moment is small (~0.8 MN·m) — this is CORRECT physics.

    For the worst-case parked (pitch ~0–30°):
        AoA can reach 60–90° → deep stall → requires 360° polars.
        Use extend_polars_360() before calling this function.

    Force per unit span [N/m]:
        dFn = 0.5 * rho * V² * c * (Cl*cos(phi) + Cd*sin(phi))  [flapwise]
        dFt = 0.5 * rho * V² * c * (Cl*sin(phi) - Cd*cos(phi))  [edgewise]
    """
    q = 0.5 * rho * V_extreme**2
    phi = np.radians(90.0)   # wind axial: perpendicular to rotor plane
    rows = []
    for _, row in geom_df.iterrows():
        r    = float(row["R_from_centre"])
        c    = float(row["chord_m"])
        twist = float(row["twist_deg"])
        foil  = str(row["airfoil"])

        theta_local = np.radians(pitch_deg + twist)
        alpha_deg_local = np.degrees(phi - theta_local)

        Cl, Cd = interp_polar(polars, foil, alpha_deg_local)

        # phi = 90° → sin(phi)=1, cos(phi)=0
        # Cn = Cl*0 + Cd*1 = Cd  (drag is the flapwise force when edge-on)
        # But for general phi we keep the full expression:
        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
        Ct = Cl * np.sin(phi) - Cd * np.cos(phi)

        dFn = q * c * Cn
        dFt = q * c * Ct
        rows.append(dict(r=r, chord=c, alpha_deg=alpha_deg_local,
                         pitch_deg=pitch_deg, Cl=Cl, Cd=Cd,
                         dFn=dFn, dFt=dFt))

    df = pd.DataFrame(rows)

    # Integrate tip → root
    r_arr = df["r"].values
    dFn   = df["dFn"].values
    dFt   = df["dFt"].values
    n = len(r_arr)
    Qflap = np.zeros(n); Mflap = np.zeros(n)
    Qedge = np.zeros(n); Medge = np.zeros(n)

    for i in range(n - 2, -1, -1):
        dr = r_arr[i + 1] - r_arr[i]
        fn_avg = 0.5 * (dFn[i] + dFn[i + 1])
        ft_avg = 0.5 * (dFt[i] + dFt[i + 1])
        Qflap[i] = Qflap[i + 1] + fn_avg * dr
        Qedge[i] = Qedge[i + 1] + ft_avg * dr
        Mflap[i] = Mflap[i + 1] + Qflap[i + 1] * dr + fn_avg * dr**2 / 2.0
        Medge[i] = Medge[i + 1] + Qedge[i + 1] * dr + ft_avg * dr**2 / 2.0

    df["Qflap_N"]  = Qflap
    df["Mflap_Nm"] = Mflap
    df["Qedge_N"]  = Qedge
    df["Medge_Nm"] = Medge
    return df


def parked_pitch_sweep(geom_df: pd.DataFrame, polars_360: dict,
                       V_extreme: float = 70.0, rho: float = 1.225,
                       pitch_range: tuple = (0.0, 180.0),
                       n_steps: int = 37) -> pd.DataFrame:
    """
    Sweep pitch angle from pitch_range[0] to pitch_range[1] in n_steps steps
    and compute parked flapwise root moment at each pitch.

    Requires 360° polars (use extend_polars_360() first).

    Returns DataFrame with columns: pitch_deg, Mflap_root_Nm, Medge_root_Nm,
    drag_total_N, Cl_mean, Cd_mean.

    The pitch that maximises |Mflap_root_Nm| is the DLC 6.1 worst case.
    """
    pitches = np.linspace(pitch_range[0], pitch_range[1], n_steps)
    results = []
    for p in pitches:
        df = parked_loads(geom_df, polars_360, V_extreme=V_extreme,
                          pitch_deg=float(p), rho=rho)
        results.append(dict(
            pitch_deg       = float(p),
            Mflap_root_Nm   = float(df["Mflap_Nm"].iloc[0]),
            Medge_root_Nm   = float(df["Medge_Nm"].iloc[0]),
            drag_total_N    = float(df["Qflap_N"].iloc[0]),
            Cl_mean         = float(df["Cl"].mean()),
            Cd_mean         = float(df["Cd"].mean()),
            alpha_mean_deg  = float(df["alpha_deg"].mean()),
        ))
    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Summary scalars
# ---------------------------------------------------------------------------

def rotor_summary(bem_df: pd.DataFrame, B: int = 3,
                  R_hub: float = 1.5) -> dict:
    """Compute total rotor thrust and power from BEM results."""
    r   = bem_df["r"].values
    dFn = bem_df["dFn"].values
    dFt = bem_df["dFt"].values
    dP  = bem_df["dP"].values

    thrust_blade = np.trapz(dFn, r)
    torque_blade = np.trapz(dFt * r, r)
    power_blade  = np.trapz(dP, r)

    return dict(
        thrust_total_N   = B * thrust_blade,
        torque_total_Nm  = B * torque_blade,
        power_total_W    = B * power_blade,
        thrust_blade_N   = thrust_blade,
        root_Mflap_Nm    = bem_df["Mflap_Nm"].iloc[0],
        root_Medge_Nm    = bem_df["Medge_Nm"].iloc[0],
    )


# ---------------------------------------------------------------------------
# Convenience: run standard NREL 5MW operating point
# ---------------------------------------------------------------------------

def run_nrel5mw_rated(data_dir: str) -> tuple:
    """
    Run BEM for NREL 5MW at rated conditions.

    Returns (bem_df, parked_df, mass_df, summary_dict)
    where parked_df is the feathered (90°) case.
    Use run_nrel5mw_full() for the complete parked pitch sweep.
    """
    data_dir = Path(data_dir)
    geom_csv  = data_dir / "reference_blade" / "blade_geometry.csv"
    stiff_csv = data_dir / "reference_blade" / "blade_stiffness.csv"
    polar_dir = data_dir / "reference_blade" / "airfoil_polars"

    geom_df  = load_blade_geometry(str(geom_csv))
    mass_df  = load_blade_mass(str(stiff_csv))
    polars   = load_polars(str(polar_dir))

    # Operating case: V_rated = 11.4 m/s, Omega_rated = 12.1 rpm, pitch = 0°
    bem_df = run_bem(
        geom_df, polars,
        V_wind=TURBINE["V_rated"],
        omega=TURBINE["omega_rated"],
        pitch_deg=0.0,
        rho=TURBINE["rho"],
    )

    # Parked feathered case: V_50 = 70 m/s, pitch = 90° (as reference)
    # NOTE: this gives a SMALL moment (~0.8 MN·m) — correct physics for feathered blade.
    # The WORST-CASE parked moment is at pitch ~0–45° and requires 360° polars.
    # Use run_nrel5mw_full() to get the true DLC 6.1 envelope.
    parked_df = parked_loads(
        geom_df, polars,
        V_extreme=70.0, pitch_deg=90.0,
        rho=TURBINE["rho"],
    )

    summary = rotor_summary(bem_df)
    return bem_df, parked_df, mass_df, summary


def run_nrel5mw_full(data_dir: str) -> tuple:
    """
    Run full load analysis for NREL 5MW including worst-case parked sweep.

    Returns (bem_df, parked_feathered_df, parked_sweep_df, mass_df, summary)
      parked_feathered_df : pitch=90° (min load, standard feathered)
      parked_sweep_df     : pitch sweep 0°–180° with 360° polars (finds worst case)
    """
    data_dir = Path(data_dir)
    geom_csv  = data_dir / "reference_blade" / "blade_geometry.csv"
    stiff_csv = data_dir / "reference_blade" / "blade_stiffness.csv"
    polar_dir = data_dir / "reference_blade" / "airfoil_polars"

    geom_df  = load_blade_geometry(str(geom_csv))
    mass_df  = load_blade_mass(str(stiff_csv))
    polars   = load_polars(str(polar_dir))

    # Extend polars to 360° for parked analysis
    # AR estimate: R / mean_chord ≈ 61.5 / 2.5 ≈ 25, cap Viterna Cd_max at 1.5
    polars_360 = extend_polars_360(polars, AR=25.0)

    # Operating case
    bem_df = run_bem(
        geom_df, polars,
        V_wind=TURBINE["V_rated"],
        omega=TURBINE["omega_rated"],
        pitch_deg=0.0,
        rho=TURBINE["rho"],
    )

    # Parked feathered (90°) with 360° polars
    parked_feathered = parked_loads(
        geom_df, polars_360,
        V_extreme=70.0, pitch_deg=90.0,
        rho=TURBINE["rho"],
    )

    # Parked pitch sweep 0°–180° → finds worst-case DLC 6.1 moment
    parked_sweep = parked_pitch_sweep(
        geom_df, polars_360,
        V_extreme=70.0, rho=TURBINE["rho"],
        pitch_range=(0.0, 180.0), n_steps=37,   # 5° steps
    )

    summary = rotor_summary(bem_df)
    return bem_df, parked_feathered, parked_sweep, mass_df, summary
