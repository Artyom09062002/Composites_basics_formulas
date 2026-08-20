"""Phase 4 analytical structural model for the NREL/Sandia 61.5 m blade.

The model is deliberately screening-level.  Published distributed stiffness
and mass are the calibration backbone; spar-cap and web layups are explicit so
that glass and hybrid alternatives can be compared with CLT failure checks.
All units are SI.
"""

from dataclasses import dataclass, replace
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PHASE_ROOT = Path(__file__).resolve().parents[1]
PHASES_ROOT = PHASE_ROOT.parent
for code_dir in (
    PHASE_ROOT / "code",
    PHASES_ROOT / "phase_03_materials" / "code",
    PHASES_ROOT / "phase_02_aerodynamics" / "code",
):
    sys.path.insert(0, str(code_dir))

AERO_DATA = PHASES_ROOT / "phase_02_aerodynamics" / "data"

from bem_solver import (extend_polars_360, gravity_edge_moment,
                        load_blade_geometry, load_blade_mass, load_polars,
                        parked_loads, run_nrel5mw_full)
from clt_engine import Laminate, Material, Ply, tsai_wu_R
from materials_db import get_ply, load_materials


GAMMA_F = 1.35
GAMMA_M_NORMAL = 2.0
GAMMA_M_SHEAR = 2.5
TIP_LIMIT_M = 5.5
REFERENCE_MASS_KG = 17_740.0
REFERENCE_CAP_MASS_KG = 7_000.0
CAP_EI_FRACTION = 0.25
HYBRID_CARBON_EI_FRACTION = 0.70
WEB_SHEAR_FRACTION = 0.50
SPAR_START_M = 6.2
SPAR_END_M = 52.3
WEB_START_M = 3.0
WEB_END_M = 55.0


@dataclass
class DesignResult:
    name: str
    stations: pd.DataFrame
    mass_kg: float
    tip_uop_m: float
    tip_parked_m: float
    min_cap_rf: float
    min_web_rf: float
    min_shell_rf: float
    worst_pitch_deg: float


def _material(name: str, gamma_m: float) -> Material:
    p = get_ply(name, load_materials())
    return Material(name=name, E1=p["E1"], E2=p["E2"], G12=p["G12"],
                    v12=p["nu12"], rho=p["rho"], Xt=p["Xt"] / gamma_m,
                    Xc=p["Xc"] / gamma_m, Yt=p["Yt"] / gamma_m,
                    Yc=p["Yc"] / gamma_m, S=p["S"] / gamma_m,
                    source="data/materials/materials_db.csv")


def _round_even(values: np.ndarray) -> np.ndarray:
    return (2.0 * np.ceil(np.maximum(values, 0.0) / 2.0)).astype(int)


def _smooth_outboard(counts: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Prevent ply-count increases after the inner spar build-up."""
    out = counts.copy()
    ids = np.where(mask)[0]
    if len(ids):
        out[ids] = np.maximum.accumulate(out[ids][::-1])[::-1]
    return out


def integrate_beam(x: np.ndarray, moment: np.ndarray,
                   bending_stiffness: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Integrate curvature M/EI twice for a clamped Euler-Bernoulli beam."""
    curvature = np.asarray(moment) / np.asarray(bending_stiffness)
    dx = np.diff(x)
    slope = np.r_[0.0, np.cumsum(0.5 * (curvature[1:] + curvature[:-1]) * dx)]
    deflection = np.r_[0.0, np.cumsum(0.5 * (slope[1:] + slope[:-1]) * dx)]
    return slope, deflection


def _geometry_at(x: np.ndarray, geom: pd.DataFrame) -> tuple[np.ndarray, ...]:
    gx = geom["r_m"].to_numpy() - 1.5
    chord_src = geom["chord_m"].to_numpy()
    tc_map = {"Cylinder1": 1.00, "Cylinder2": 0.55, "DU40_A17": 0.40,
              "DU35_A17": 0.35, "DU30_A17": 0.30, "DU25_A17": 0.25,
              "DU21_A17": 0.21, "NACA64_A17": 0.18}
    tc_src = np.array([tc_map[a] for a in geom["airfoil"]])
    chord = np.interp(x, gx, chord_src, left=chord_src[0], right=chord_src[-1])
    tc = np.interp(x, gx, tc_src, left=tc_src[0], right=tc_src[-1])
    cap_separation = 0.75 * tc * chord
    cap_width = np.clip(0.18 * chord, 0.25, 0.85)
    web_height = 0.65 * tc * chord
    return chord, tc, cap_separation, cap_width, web_height


def _scaled_distribution(df: pd.DataFrame, column: str, x: np.ndarray,
                         target_root: float | None = None) -> np.ndarray:
    local_x = df["r"].to_numpy() - float(df["r"].iloc[0])
    values = np.abs(df[column].to_numpy())
    out = np.interp(x, local_x, values, left=values[0], right=0.0)
    if target_root is not None:
        out *= target_root / out[0]
    return out


def load_envelopes(x: np.ndarray) -> tuple[dict, float]:
    op, feathered, sweep, _, _ = run_nrel5mw_full(str(AERO_DATA))
    worst_row = sweep.loc[sweep["Mflap_root_Nm"].abs().idxmax()]
    worst_pitch = float(worst_row["pitch_deg"])
    geom = load_blade_geometry(str(AERO_DATA / "reference_blade" / "blade_geometry.csv"))
    polars = extend_polars_360(load_polars(str(AERO_DATA / "reference_blade" / "airfoil_polars")), AR=25.0)
    parked_worst = parked_loads(geom, polars, V_extreme=70.0, pitch_deg=worst_pitch)
    cases = {
        "U_op": {"M": _scaled_distribution(op, "Mflap_Nm", x, 7.530700295e6),
                 "Q": _scaled_distribution(op, "Qflap_N", x)},
        "U_park_feather": {"M": _scaled_distribution(feathered, "Mflap_Nm", x),
                           "Q": _scaled_distribution(feathered, "Qflap_N", x)},
        "U_park_fault": {"M": _scaled_distribution(parked_worst, "Mflap_Nm", x,
                                                     abs(float(worst_row["Mflap_root_Nm"]))),
                         "Q": _scaled_distribution(parked_worst, "Qflap_N", x)},
    }
    for case in cases.values():
        case["M_design"] = GAMMA_F * case["M"]
        case["Q_design"] = GAMMA_F * case["Q"]
    return cases, worst_pitch


def _cap_laminate(n_glass: int, n_carbon: int) -> Laminate:
    glass = _material("ELT5500_UD", GAMMA_M_NORMAL)
    carbon = _material("Newport307_CarbonUD", GAMMA_M_NORMAL)
    plies = ([Ply(glass, 0.0, 0.00091)] * (n_glass // 2) +
             [Ply(carbon, 0.0, 0.00036)] * n_carbon +
             [Ply(glass, 0.0, 0.00091)] * (n_glass // 2))
    return Laminate(plies)


def _laminate_rf(laminate: Laminate, axial_strain: float) -> float:
    A, _, _ = laminate.abd()
    eps = np.array([axial_strain, -A[1, 0] / A[1, 1] * axial_strain, 0.0])
    rfs = []
    for sign in (1.0, -1.0):
        for _, _, stress in laminate.ply_stresses(sign * eps, np.zeros(3)):
            material = laminate.plies[len(rfs) % len(laminate.plies)].material
            rfs.append(tsai_wu_R(stress, material))
    return float(np.min(rfs))


def _zone(x: float) -> str:
    if x < 3.0: return "root"
    if x < 20.0: return "inner"
    if x < 40.0: return "mid"
    if x < 55.0: return "outer"
    return "tip"


def build_design(kind: str = "glass") -> DesignResult:
    if kind not in {"glass", "hybrid"}:
        raise ValueError("kind must be 'glass' or 'hybrid'")
    stiff = pd.read_csv(AERO_DATA / "reference_blade" / "blade_stiffness.csv", comment="#")
    geom = load_blade_geometry(str(AERO_DATA / "reference_blade" / "blade_geometry.csv"))
    x = stiff["r_m"].to_numpy()
    ei_ref = stiff["EI_flap_Nm2"].to_numpy()
    chord, tc, dcap, bcap, hweb = _geometry_at(x, geom)
    spar = (x >= SPAR_START_M) & (x <= SPAR_END_M)
    web = (x >= WEB_START_M) & (x <= WEB_END_M)
    cases, worst_pitch = load_envelopes(x)

    glass = get_ply("ELT5500_UD", load_materials())
    carbon = get_ply("Newport307_CarbonUD", load_materials())
    db = get_ply("Saertex_DB", load_materials())
    t_ref = np.where(spar, 2 * CAP_EI_FRACTION * ei_ref /
                     (glass["E1"] * bcap * dcap**2), 0.0)
    if kind == "glass":
        n_glass = _round_even(t_ref / glass["t_ply"])
        n_carbon = np.zeros_like(n_glass)
    else:
        n_glass = _round_even((1 - HYBRID_CARBON_EI_FRACTION) * t_ref / glass["t_ply"])
        n_carbon = _round_even(HYBRID_CARBON_EI_FRACTION * t_ref * glass["E1"] /
                               (carbon["E1"] * carbon["t_ply"]))
    n_glass = _smooth_outboard(n_glass, spar)
    n_carbon = _smooth_outboard(n_carbon, spar)

    q_max = np.max(np.vstack([c["Q_design"] for c in cases.values()]), axis=0)
    tau_allow = db["S"] / GAMMA_M_SHEAR
    t_web_req = WEB_SHEAR_FRACTION * q_max / (4.0 * hweb * tau_allow)
    n_web = np.where(web, np.maximum(4, _round_even(t_web_req / db["t_ply"])), 0)
    n_web = _smooth_outboard(n_web, (x >= SPAR_START_M) & web)

    ei_noncap = ei_ref * (1.0 - np.where(spar, CAP_EI_FRACTION, 0.0))
    ea_cap = (glass["E1"] * bcap * n_glass * glass["t_ply"] +
              carbon["E1"] * bcap * n_carbon * carbon["t_ply"])
    ei_design = ei_noncap + 0.5 * ea_cap * dcap**2

    cap_mass = 2.0 * np.trapz(bcap * (n_glass * glass["t_ply"] * glass["rho"] +
                                     n_carbon * carbon["t_ply"] * carbon["rho"]), x)
    nominal_web_mass = np.trapz(4 * hweb * np.where(web, 4 * db["t_ply"], 0.0) * db["rho"], x)
    actual_web_mass = np.trapz(4 * hweb * n_web * db["t_ply"] * db["rho"], x)
    mass = REFERENCE_MASS_KG - REFERENCE_CAP_MASS_KG + cap_mass + actual_web_mass - nominal_web_mass

    cap_rf = np.full(len(x), np.inf)
    shell_rf = np.full(len(x), np.inf)
    triax = _material("SNL_Triax", GAMMA_M_NORMAL)
    for i in range(len(x)):
        max_strain = max(c["M_design"][i] / ei_design[i] * dcap[i] / 2 for c in cases.values())
        if spar[i]:
            cap_rf[i] = _laminate_rf(_cap_laminate(int(n_glass[i]), int(n_carbon[i])), max_strain)
        else:
            shell_rf[i] = min(tsai_wu_R(np.array([triax.E1 * max_strain, 0.0, 0.0]), triax),
                              tsai_wu_R(np.array([-triax.E1 * max_strain, 0.0, 0.0]), triax))

    db_mat = _material("Saertex_DB", GAMMA_M_SHEAR)
    tau = np.divide(WEB_SHEAR_FRACTION * q_max, 4 * hweb * n_web * db["t_ply"],
                    out=np.zeros_like(q_max), where=n_web > 0)
    web_rf = np.array([tsai_wu_R(np.array([0.0, 0.0, v]), db_mat) if web[i] else np.inf
                       for i, v in enumerate(tau)])

    _, uop_defl = integrate_beam(x, cases["U_op"]["M_design"], ei_design)
    _, park_defl = integrate_beam(x, cases["U_park_fault"]["M_design"], ei_design)
    _, gravity_m = gravity_edge_moment(load_blade_mass(str(AERO_DATA / "reference_blade" / "blade_stiffness.csv")))
    edge_strain = GAMMA_F * gravity_m / stiff["EI_edge_Nm2"].to_numpy() * 0.45 * chord
    edge_rf = np.array([min(tsai_wu_R(np.array([triax.E1 * e, 0.0, 0.0]), triax),
                                tsai_wu_R(np.array([-triax.E1 * e, 0.0, 0.0]), triax)) for e in edge_strain])
    shell_rf = np.minimum(shell_rf, edge_rf)

    table = pd.DataFrame({"r_m": x, "zone": [_zone(v) for v in x], "chord_m": chord,
                          "tc": tc, "cap_width_m": bcap, "cap_separation_m": dcap,
                          "n_cap_glass": n_glass, "n_cap_carbon": n_carbon,
                          "n_web_db_per_skin": n_web, "EI_flap_Nm2": ei_design,
                          "M_uop_design_Nm": cases["U_op"]["M_design"],
                          "M_park_fault_design_Nm": cases["U_park_fault"]["M_design"],
                          "deflection_uop_m": uop_defl,
                          "deflection_park_fault_m": park_defl,
                          "cap_RF": cap_rf, "web_RF": web_rf, "shell_RF": shell_rf})
    return DesignResult(kind, table, float(mass), float(uop_defl[-1]), float(park_defl[-1]),
                        float(np.min(cap_rf)), float(np.min(web_rf)),
                        float(np.min(shell_rf)), worst_pitch)


def write_phase4_outputs() -> tuple[DesignResult, DesignResult]:
    out = PHASE_ROOT / "results"
    out.mkdir(parents=True, exist_ok=True)
    glass = build_design("glass")
    hybrid = build_design("hybrid")
    glass.stations.to_csv(out / "glass_station_results.csv", index=False)
    hybrid.stations.to_csv(out / "hybrid_station_results.csv", index=False)
    schedule = glass.stations[["r_m", "zone", "n_cap_glass", "n_web_db_per_skin"]].copy()
    schedule["hybrid_n_cap_glass"] = hybrid.stations["n_cap_glass"]
    schedule["hybrid_n_cap_carbon"] = hybrid.stations["n_cap_carbon"]
    schedule["shell_schedule"] = schedule["zone"].map({
        "root": "SNL Triax effective layers to 80-100 mm; no foam",
        "inner": "8 DB + 16 Triax effective layers; PVC H100 50 mm",
        "mid": "4 DB + 8 Triax effective layers; PVC H100 40 mm",
        "outer": "4 DB + 4 Triax effective layers; PVC H100 30 mm",
        "tip": "6 SNL Triax effective layers; no spar cap"})
    schedule.to_csv(PHASE_ROOT / "data" / "layup_schedule.csv", index=False)
    summary = pd.DataFrame([{
        "design": d.name, "mass_kg": d.mass_kg,
        "mass_delta_vs_reference_pct": 100 * (d.mass_kg / REFERENCE_MASS_KG - 1),
        "tip_uop_m": d.tip_uop_m, "tip_parked_fault_m": d.tip_parked_m,
        "min_cap_RF": d.min_cap_rf, "min_web_RF": d.min_web_rf,
        "min_shell_RF": d.min_shell_rf, "worst_pitch_deg": d.worst_pitch_deg}
        for d in (glass, hybrid)])
    summary.to_csv(out / "phase4_summary.csv", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(glass.stations.r_m, glass.stations.M_uop_design_Nm / 1e6,
                 label="U_op design", lw=2.2)
    axes[0].plot(glass.stations.r_m, glass.stations.M_park_fault_design_Nm / 1e6,
                 label="U_park fault design", lw=2.2)
    axes[0].set(xlabel="Span r [m]", ylabel="Flapwise moment [MN·m]",
                title="Phase 4 design load envelopes")
    axes[0].grid(alpha=.3); axes[0].legend()
    axes[1].plot(glass.stations.r_m, glass.stations.deflection_uop_m,
                 label="Glass, U_op", lw=2.2)
    axes[1].plot(glass.stations.r_m, glass.stations.deflection_park_fault_m,
                 label="Glass, parked fault", lw=2.2)
    axes[1].plot(hybrid.stations.r_m, hybrid.stations.deflection_park_fault_m,
                 "--", label="Hybrid, parked fault", lw=2.2)
    axes[1].axhline(TIP_LIMIT_M, color="red", ls=":", label="5.5 m limit")
    axes[1].set(xlabel="Span r [m]", ylabel="Deflection [m]",
                title="Analytical beam deflection")
    axes[1].grid(alpha=.3); axes[1].legend()
    fig.tight_layout(); fig.savefig(out / "structural_performance.png", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].step(glass.stations.r_m, glass.stations.n_cap_glass, where="mid",
                 label="Glass baseline", lw=2.2)
    axes[0].step(hybrid.stations.r_m, hybrid.stations.n_cap_glass, where="mid",
                 label="Hybrid: glass outer", lw=2.2)
    axes[0].step(hybrid.stations.r_m, hybrid.stations.n_cap_carbon, where="mid",
                 label="Hybrid: carbon core", lw=2.2)
    axes[0].set(xlabel="Span r [m]", ylabel="Ply count per cap",
                title="Symmetric spar-cap schedule")
    axes[0].grid(alpha=.3); axes[0].legend()
    axes[1].step(glass.stations.r_m, glass.stations.n_web_db_per_skin,
                 where="mid", color="tab:purple", lw=2.2)
    axes[1].set(xlabel="Span r [m]", ylabel="DB effective plies per web skin",
                title="Two-web shear schedule")
    axes[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(out / "layup_schedule.png", dpi=170)
    plt.close(fig)
    return glass, hybrid


if __name__ == "__main__":
    g, h = write_phase4_outputs()
    for d in (g, h):
        print(f"{d.name}: mass={d.mass_kg:.0f} kg, U_op tip={d.tip_uop_m:.2f} m, "
              f"parked tip={d.tip_parked_m:.2f} m, RF cap/web/shell="
              f"{d.min_cap_rf:.2f}/{d.min_web_rf:.2f}/{d.min_shell_rf:.2f}")
