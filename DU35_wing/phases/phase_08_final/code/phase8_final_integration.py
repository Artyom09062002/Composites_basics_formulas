"""Phase 8 final whole-blade integration for the NREL 5 MW screening blade.

The authoritative structural definition is the glass-cap baseline with the
FreeCAD v5 robust sandwich webs: 5 DB layers/side forward, 6 DB layers/side
aft, and 60 mm H100 core.  This module reconciles mass, reruns beam static and
modal response, checks static strength and sandwich buckling, and carries the
closed 20-year fatigue cases into one auditable result.

The model is a reproducible beam/shell-equivalent engineering screen.  It is
not an IEC 61400-5 certification analysis or a detailed shell/solid joint FE
model.  All units are SI.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PHASE_ROOT = Path(__file__).resolve().parents[1]
PHASES_ROOT = PHASE_ROOT.parent
for code_dir in (
    PHASE_ROOT / "code",
    PHASES_ROOT / "phase_07_fatigue" / "code",
    PHASES_ROOT / "phase_05_fea" / "code",
    PHASES_ROOT / "phase_04_structural_design" / "code",
    PHASES_ROOT / "phase_03_materials" / "code",
    PHASES_ROOT / "phase_02_aerodynamics" / "code",
):
    sys.path.insert(0, str(code_dir))

from clt_engine import Laminate, Material, Ply, tsai_wu_R
from fatigue_analysis import ROBUST_DAMAGE_TARGET, assess_phase7, robust_case_summary
from fea_model import BladeBeamFE, distributed_load_from_shear, panel_buckling
from materials_db import get_ply, load_materials
from structural_model import (GAMMA_F, GAMMA_M_NORMAL, REFERENCE_CAP_MASS_KG,
                              REFERENCE_MASS_KG, SPAR_END_M, SPAR_START_M,
                              _material, build_design, load_envelopes)


STATIC_RF_TARGET = 1.0
TIP_DEFLECTION_LIMIT_M = 5.5
THREE_P_HZ = 0.606
FREQUENCY_TARGET_HZ = 1.2 * THREE_P_HZ
FACE_LAYER_T_M = 0.001
FACE_DENSITY_KG_M3 = 1830.0
CORE_DENSITY_KG_M3 = 100.0
FACE_G_PA = 11.8e9
CORE_G_PA = 22e6
FACE_SHEAR_STRENGTH_PA = 62e6
FACE_COMPRESSION_STRENGTH_PA = 213e6
# NASA/TM-2012-217694 Table 3 gives Divinycell H100 shear yield = 1.13 MPa.
CORE_SHEAR_STRENGTH_PA = 1.13e6
FINAL_WEB_LAYERS = {"forward": 5, "aft": 6}
FINAL_CORE_T_M = 0.060
REFERENCE_WEB_LAYERS = 2
REFERENCE_CORE_T_M = 0.050
WEB_SPLIT_ERRORS_PCT = (-10.0, 0.0, 10.0)
LONG_PANEL_ASPECT = 10.0

STIFFNESS_CSV = PHASES_ROOT / "phase_02_aerodynamics" / "data" / "reference_blade" / "blade_stiffness.csv"
WEB_GEOMETRY_CSV = (PHASES_ROOT / "phase_04_structural_design" / "studies"
                    / "full_blade_webs" / "data" / "cad_web_geometry.csv")
OUT_DIR = PHASE_ROOT / "results"
OUT_JSON = OUT_DIR / "phase8_results.json"


@dataclass(frozen=True)
class AcceptanceCheck:
    category: str
    check_id: str
    value: float
    criterion: str
    limit: float | None
    margin: float | None
    status: str
    governing: str
    basis: str


def _finite_min(values: pd.Series | np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.min(array[np.isfinite(array)]))


def _integral(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.trapz(np.asarray(y, float), np.asarray(x, float)))


def final_mass_distribution(stations: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    """Reconcile the source mass shape and replace caps/webs with v6 values."""
    stiffness = pd.read_csv(STIFFNESS_CSV, comment="#")
    web_geometry = pd.read_csv(WEB_GEOMETRY_CSV)
    materials = load_materials()
    glass = get_ply("ELT5500_UD", materials)
    x = stations.r_m.to_numpy(float)

    raw_reference_pm = np.interp(x, stiffness.r_m, stiffness.mass_pm_kg_m)
    raw_reference_mass = _integral(raw_reference_pm, x)
    reference_scale = REFERENCE_MASS_KG / raw_reference_mass
    reference_pm = raw_reference_pm * reference_scale

    cap_pm = (2.0 * stations.cap_width_m.to_numpy(float)
              * stations.n_cap_glass.to_numpy(float)
              * glass["t_ply"] * glass["rho"])
    cap_mass = _integral(cap_pm, x)
    cap_delta_mass = cap_mass - REFERENCE_CAP_MASS_KG
    cap_delta_pm = cap_pm / cap_mass * cap_delta_mass

    web_delta_pm = np.zeros_like(x)
    reference_web_pm = np.zeros_like(x)
    final_web_pm = np.zeros_like(x)
    for web in ("forward", "aft"):
        rows = web_geometry[web_geometry.web == web].sort_values("blade_span_m")
        span = rows.blade_span_m.to_numpy(float)
        source_height = rows.web_height_m.to_numpy(float)
        height = np.interp(x, span, source_height)
        active = (x >= span.min() - 1e-3) & (x <= span.max() + 1e-3)
        height = np.where(active, height, 0.0)
        reference_areal = (2.0 * REFERENCE_WEB_LAYERS * FACE_LAYER_T_M
                           * FACE_DENSITY_KG_M3 + REFERENCE_CORE_T_M
                           * CORE_DENSITY_KG_M3)
        final_areal = (2.0 * FINAL_WEB_LAYERS[web] * FACE_LAYER_T_M
                       * FACE_DENSITY_KG_M3 + FINAL_CORE_T_M
                       * CORE_DENSITY_KG_M3)
        # Preserve the exact CAD-station area integral; the beam stations are
        # coarser and otherwise overestimate the two web areas by about 3%.
        coarse_area = _integral(height, x)
        exact_area = _integral(source_height, span)
        mapped_height = height * exact_area / coarse_area
        reference_web_pm += mapped_height * reference_areal
        final_web_pm += mapped_height * final_areal
        web_delta_pm += mapped_height * (final_areal - reference_areal)

    final_pm = reference_pm + cap_delta_pm + web_delta_pm
    table = pd.DataFrame({
        "r_m": x,
        "source_mass_pm_raw_kg_m": raw_reference_pm,
        "source_mass_pm_normalized_kg_m": reference_pm,
        "cap_mass_delta_pm_kg_m": cap_delta_pm,
        "reference_web_mass_pm_kg_m": reference_web_pm,
        "final_v6_web_mass_pm_kg_m": final_web_pm,
        "web_mass_delta_pm_kg_m": web_delta_pm,
        "final_mass_pm_kg_m": final_pm,
    })
    breakdown = {
        "reference_reported_mass_kg": REFERENCE_MASS_KG,
        "raw_mass_table_integral_kg": raw_reference_mass,
        "mass_shape_normalization_factor": reference_scale,
        "reference_cap_budget_kg": REFERENCE_CAP_MASS_KG,
        "final_glass_cap_mass_kg": cap_mass,
        "cap_mass_delta_kg": cap_delta_mass,
        "reference_2DB_50mm_web_mass_kg": _integral(reference_web_pm, x),
        "final_v6_web_mass_kg": _integral(final_web_pm, x),
        "web_mass_delta_kg": _integral(web_delta_pm, x),
        "final_v6_blade_mass_kg": _integral(final_pm, x),
    }
    breakdown["mass_delta_vs_reference_pct"] = 100.0 * (
        breakdown["final_v6_blade_mass_kg"] / REFERENCE_MASS_KG - 1.0
    )
    return table, breakdown


def _line_load_moment(x: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Internal bending-moment magnitude for a cantilever line load."""
    return np.array([
        _integral(q[i:] * (x[i:] - x[i]), x[i:]) if i < len(x) - 1 else 0.0
        for i in range(len(x))
    ])


def _shell_strength_with_final_gravity(
    stations: pd.DataFrame, final_mass_pm: np.ndarray, cases: dict,
) -> tuple[float, float]:
    """Re-evaluate shell/edge Tsai-Wu RF with the reconciled v6 mass."""
    x = stations.r_m.to_numpy(float)
    triax = _material("SNL_Triax", GAMMA_M_NORMAL)
    dcap = stations.cap_separation_m.to_numpy(float)
    ei_flap = stations.EI_flap_Nm2.to_numpy(float)
    chord = stations.chord_m.to_numpy(float)
    stiffness = pd.read_csv(STIFFNESS_CSV, comment="#")
    ei_edge = np.interp(x, stiffness.r_m, stiffness.EI_edge_Nm2)

    shell_flap_rf = np.full(len(x), np.inf)
    for i, span in enumerate(x):
        if not (SPAR_START_M <= span <= SPAR_END_M):
            strain = max(c["M_design"][i] / ei_flap[i] * dcap[i] / 2.0
                         for c in cases.values())
            shell_flap_rf[i] = min(
                tsai_wu_R(np.array([triax.E1 * strain, 0.0, 0.0]), triax),
                tsai_wu_R(np.array([-triax.E1 * strain, 0.0, 0.0]), triax),
            )

    gravity_q = GAMMA_F * 9.80665 * final_mass_pm
    gravity_moment = _line_load_moment(x, gravity_q)
    edge_strain = np.divide(
        gravity_moment * 0.45 * chord, ei_edge,
        out=np.zeros_like(gravity_moment), where=ei_edge > 0,
    )
    shell_edge_rf = np.array([
        min(tsai_wu_R(np.array([triax.E1 * e, 0.0, 0.0]), triax),
            tsai_wu_R(np.array([-triax.E1 * e, 0.0, 0.0]), triax))
        for e in edge_strain
    ])
    return float(np.min(np.minimum(shell_flap_rf, shell_edge_rf))), float(gravity_moment[0])


def _sandwich_D(face_layers: int, core_t_m: float) -> np.ndarray:
    face = Material("Saertex DB face", 13.6e9, 13.3e9, FACE_G_PA, 0.49)
    core = Material("H100 equivalent core", 256e6, 256e6, CORE_G_PA, 0.30)
    face_t = face_layers * FACE_LAYER_T_M
    return Laminate([
        Ply(face, 0.0, face_t),
        Ply(core, 0.0, core_t_m),
        Ply(face, 0.0, face_t),
    ]).abd()[2]


def final_web_static_screen(stations: pd.DataFrame) -> pd.DataFrame:
    """Check v6 web strength/buckling across the accepted split uncertainty."""
    geometry = pd.read_csv(WEB_GEOMETRY_CSV)
    load_r = stations.r_m.to_numpy(float)
    moment = stations.M_park_fault_design_Nm.to_numpy(float)
    design_shear = np.abs(np.gradient(moment, load_r, edge_order=2))
    rows: list[dict[str, object]] = []
    wrinkle = 0.5 * (CORE_G_PA * 256e6 * 13.6e9) ** (1.0 / 3.0)

    for span, pair_rows in geometry.groupby("blade_span_m", sort=True):
        pair = pair_rows.set_index("web")
        total_v = float(np.interp(float(span), load_r, design_shear))
        heights = {web: float(pair.loc[web, "web_height_m"])
                   for web in ("forward", "aft")}
        a66 = {web: (2.0 * FACE_G_PA * FINAL_WEB_LAYERS[web] * FACE_LAYER_T_M
                     + CORE_G_PA * FINAL_CORE_T_M)
               for web in ("forward", "aft")}
        kf = a66["forward"] * heights["forward"]
        ka = a66["aft"] * heights["aft"]
        nominal_forward = kf / (kf + ka)

        for split_error in WEB_SPLIT_ERRORS_PCT:
            forward_fraction = float(np.clip(
                nominal_forward * (1.0 + split_error / 100.0), 1e-9, 1.0 - 1e-9
            ))
            fractions = {"forward": forward_fraction, "aft": 1.0 - forward_fraction}
            for web in ("forward", "aft"):
                height = heights[web]
                nxy = total_v * fractions[web] / height
                gamma = nxy / a66[web]
                face_tau = FACE_G_PA * gamma
                core_tau = CORE_G_PA * gamma
                d = _sandwich_D(FINAL_WEB_LAYERS[web], FINAL_CORE_T_M)
                ks = 5.34 + 4.0 / LONG_PANEL_ASPECT**2
                nxy_cr = ks * np.pi**2 * np.sqrt(d[0, 0] * d[1, 1]) / height**2
                modes = {
                    "global_shear_buckling": nxy_cr / nxy,
                    "face_shear": FACE_SHEAR_STRENGTH_PA / face_tau,
                    "core_shear_NASA_yield": CORE_SHEAR_STRENGTH_PA / core_tau,
                    "face_wrinkling": wrinkle / face_tau,
                    "face_compression": FACE_COMPRESSION_STRENGTH_PA / face_tau,
                }
                governing_mode = min(modes, key=modes.get)
                rows.append({
                    "r_m": float(span), "web": web,
                    "split_error_percent": split_error,
                    "nominal_forward_fraction": nominal_forward,
                    "load_fraction": fractions[web],
                    "layers_each_side": FINAL_WEB_LAYERS[web],
                    "core_thickness_m": FINAL_CORE_T_M,
                    "design_shear_total_N": total_v,
                    "applied_Nxy_N_per_m": nxy,
                    "global_buckling_RF": modes["global_shear_buckling"],
                    "face_shear_RF": modes["face_shear"],
                    "core_shear_RF": modes["core_shear_NASA_yield"],
                    "face_wrinkling_RF": modes["face_wrinkling"],
                    "face_compression_RF": modes["face_compression"],
                    "governing_mode": governing_mode,
                    "governing_RF": modes[governing_mode],
                    "status": "PASS" if modes[governing_mode] >= STATIC_RF_TARGET else "FAIL",
                })
    return pd.DataFrame(rows)


def run_phase8() -> dict[str, object]:
    design = build_design("glass")
    stations = design.stations.copy()
    x = stations.r_m.to_numpy(float)
    cases, _ = load_envelopes(x)
    mass_table, mass = final_mass_distribution(stations)

    stiffness = pd.read_csv(STIFFNESS_CSV, comment="#")
    ei_edge = np.interp(x, stiffness.r_m, stiffness.EI_edge_Nm2)
    model = BladeBeamFE(x, stations.EI_flap_Nm2.to_numpy(float), ei_edge,
                        mass_table.final_mass_pm_kg_m.to_numpy(float))
    static_results = []
    for case_name in ("U_op", "U_park_feather", "U_park_fault"):
        q = distributed_load_from_shear(
            x, cases[case_name]["Q_design"], cases[case_name]["M_design"][0]
        )
        static_results.append(model.static("flap", q, case_name))
    gravity_q = GAMMA_F * 9.80665 * mass_table.final_mass_pm_kg_m.to_numpy(float)
    static_results.append(model.static("edge", gravity_q, "U_gravity_v6_mass"))
    flap_modal = model.modal("flap")
    edge_modal = model.modal("edge")

    shell_rf, gravity_root_moment = _shell_strength_with_final_gravity(
        stations, mass_table.final_mass_pm_kg_m.to_numpy(float), cases
    )
    cap_rf = _finite_min(stations.cap_RF)
    cap_buckling = panel_buckling(stations, "phase8_v6")
    cap_buckling_rf = _finite_min(cap_buckling.cap_buckling_RF)
    web_static = final_web_static_screen(stations)
    web_strength_rf = float(web_static[["face_shear_RF", "core_shear_RF",
                                        "face_wrinkling_RF", "face_compression_RF"]].min().min())
    web_buckling_rf = float(web_static.global_buckling_RF.min())

    components, webs, _ = assess_phase7()
    fatigue_cases = robust_case_summary(webs)
    shell_d20 = float(components[components.component == "shell_and_TE_SNL_Triax"].damage_20y.max())
    cap_d20 = float(components[components.component == "spar_cap_glass_UD"].damage_20y.max())
    web_d20 = float(fatigue_cases.max_face_D20.max())
    core_d20 = float(webs.attrs["robust_detail"].core_damage_20y.max())

    statics = pd.DataFrame([{
        "case": result.case, "direction": result.direction,
        "root_shear_N": result.root_shear_N,
        "root_moment_Nm": result.root_moment_Nm,
        "tip_displacement_m": result.tip_displacement_m,
    } for result in static_results])
    modal = pd.DataFrame([
        {"direction": direction, "mode": i + 1, "frequency_hz": float(frequency)}
        for direction, result in (("flap", flap_modal), ("edge", edge_modal))
        for i, frequency in enumerate(result.frequencies_hz)
    ])
    op = statics[statics.case == "U_op"].iloc[0]
    parked = statics[statics.case == "U_park_fault"].iloc[0]
    gravity = statics[statics.case == "U_gravity_v6_mass"].iloc[0]
    first_flap = float(modal[(modal.direction == "flap") & (modal["mode"] == 1)].iloc[0].frequency_hz)

    checks = [
        AcceptanceCheck("loads", "operating_root_moment", abs(float(op.root_moment_Nm)),
                        "included factored load case", None, None, "INFO", "root",
                        "NREL BEM operating envelope with gamma_F=1.35"),
        AcceptanceCheck("loads", "parked_fault_root_moment", abs(float(parked.root_moment_Nm)),
                        "included factored load case", None, None, "INFO", "root",
                        "70 m/s Viterna parked-fault envelope with gamma_F=1.35"),
        AcceptanceCheck("loads", "gravity_root_moment", abs(float(gravity.root_moment_Nm)),
                        "included factored load case", None, None, "INFO", "root",
                        "v6 reconciled mass distribution with gamma_F=1.35"),
        AcceptanceCheck("mass", "final_blade_mass", mass["final_v6_blade_mass_kg"],
                        "reference comparison only", REFERENCE_MASS_KG,
                        REFERENCE_MASS_KG - mass["final_v6_blade_mass_kg"], "INFO", "whole blade",
                        "Sandia reported mass plus explicit cap/web replacement deltas"),
        AcceptanceCheck("static", "cap_strength_RF", cap_rf, "RF >= 1.0",
                        STATIC_RF_TARGET, cap_rf - STATIC_RF_TARGET,
                        "PASS" if cap_rf >= STATIC_RF_TARGET else "FAIL", "spanwise minimum",
                        "Phase 4 CLT under factored flap loads"),
        AcceptanceCheck("static", "shell_strength_RF", shell_rf, "RF >= 1.0",
                        STATIC_RF_TARGET, shell_rf - STATIC_RF_TARGET,
                        "PASS" if shell_rf >= STATIC_RF_TARGET else "FAIL", "spanwise minimum",
                        "CLT with v6 gravity mass and factored flap loads"),
        AcceptanceCheck("static", "web_strength_RF", web_strength_rf, "RF >= 1.0",
                        STATIC_RF_TARGET, web_strength_rf - STATIC_RF_TARGET,
                        "PASS" if web_strength_rf >= STATIC_RF_TARGET else "FAIL", "all split cases",
                        "v6 sandwich face/core strength screen"),
        AcceptanceCheck("buckling", "cap_buckling_RF", cap_buckling_rf, "RF >= 1.0",
                        STATIC_RF_TARGET, cap_buckling_rf - STATIC_RF_TARGET,
                        "PASS" if cap_buckling_rf >= STATIC_RF_TARGET else "FAIL", "spanwise minimum",
                        "Phase 5 classical cap panel screen"),
        AcceptanceCheck("buckling", "web_buckling_RF", web_buckling_rf, "RF >= 1.0",
                        STATIC_RF_TARGET, web_buckling_rf - STATIC_RF_TARGET,
                        "PASS" if web_buckling_rf >= STATIC_RF_TARGET else "FAIL", "a/h=10 bound",
                        "CAD-height v6 sandwich D matrix; split error -10/0/+10%"),
        AcceptanceCheck("deflection", "operating_tip_deflection", abs(float(op.tip_displacement_m)),
                        "tip <= 5.5 m", TIP_DEFLECTION_LIMIT_M,
                        TIP_DEFLECTION_LIMIT_M - abs(float(op.tip_displacement_m)),
                        "PASS" if abs(float(op.tip_displacement_m)) <= TIP_DEFLECTION_LIMIT_M else "FAIL",
                        "blade tip", "beam FE operating envelope"),
        AcceptanceCheck("deflection", "parked_fault_tip_deflection", abs(float(parked.tip_displacement_m)),
                        "tip <= 5.5 m", TIP_DEFLECTION_LIMIT_M,
                        TIP_DEFLECTION_LIMIT_M - abs(float(parked.tip_displacement_m)),
                        "PASS" if abs(float(parked.tip_displacement_m)) <= TIP_DEFLECTION_LIMIT_M else "FAIL",
                        "blade tip", "beam FE parked-fault envelope"),
        AcceptanceCheck("modal", "first_flap_frequency", first_flap,
                        "f1 >= 1.2 x 3P", FREQUENCY_TARGET_HZ,
                        first_flap - FREQUENCY_TARGET_HZ,
                        "PASS" if first_flap >= FREQUENCY_TARGET_HZ else "FAIL", "mode 1 flap",
                        "beam FE with reconciled v6 mass"),
        AcceptanceCheck("fatigue", "shell_TE_D20", shell_d20, "D20 <= 0.70",
                        ROBUST_DAMAGE_TARGET, ROBUST_DAMAGE_TARGET - shell_d20,
                        "PASS" if shell_d20 <= ROBUST_DAMAGE_TARGET else "FAIL", "r=10.25 m",
                        "SAND reference Miner damage transfer"),
        AcceptanceCheck("fatigue", "spar_cap_D20", cap_d20, "D20 <= 0.70",
                        ROBUST_DAMAGE_TARGET, ROBUST_DAMAGE_TARGET - cap_d20,
                        "PASS" if cap_d20 <= ROBUST_DAMAGE_TARGET else "FAIL", "r=10.25 m",
                        "glass UD equivalent-strain transfer"),
        AcceptanceCheck("fatigue", "web_face_worst_D20", web_d20, "D20 <= 0.70",
                        ROBUST_DAMAGE_TARGET, ROBUST_DAMAGE_TARGET - web_d20,
                        "PASS" if web_d20 <= ROBUST_DAMAGE_TARGET else "FAIL", "r=26.65 m",
                        "combined b=9 and A66h +10% robust case"),
        AcceptanceCheck("fatigue", "H100_core_D20", core_d20, "D20 <= 0.70",
                        ROBUST_DAMAGE_TARGET, ROBUST_DAMAGE_TARGET - core_d20,
                        "PASS" if core_d20 <= ROBUST_DAMAGE_TARGET else "FAIL", "spanwise maximum",
                        "NASA/TM-2012-217694 H100 shear-fatigue curve"),
    ]
    check_table = pd.DataFrame([asdict(check) for check in checks])

    return {
        "metadata": {
            "phase": 8,
            "scope": "full-scale NREL 5 MW beam/shell-equivalent research screening",
            "design": "FreeCAD v6; glass caps; forward/aft 5/6 DB per side; 60 mm H100",
            "acceptance": "RF>=1.0; D20<=0.70; tip<=5.5 m; f1>=1.2x3P; no FAIL",
        },
        "final_design": {
            "blade_length_m": 61.5,
            "spar_cap_material": "E-LT-5500 UD glass",
            "forward_web_DB_layers_each_side": 5,
            "aft_web_DB_layers_each_side": 6,
            "web_core": "Divinycell H100 equivalent",
            "web_core_thickness_m": FINAL_CORE_T_M,
            "web_zones_per_web": 1,
            "external_aero_geometry_changed": False,
            "web_positions_changed": False,
        },
        "mass_breakdown": mass,
        "checks": check_table.to_dict("records"),
        "static_cases": statics.to_dict("records"),
        "modal_results": modal.to_dict("records"),
        "mass_distribution": mass_table.to_dict("records"),
        "web_static_detail": web_static.to_dict("records"),
        "fatigue_cases": fatigue_cases.to_dict("records"),
        "station_response": pd.DataFrame({
            "r_m": x,
            "final_mass_pm_kg_m": mass_table.final_mass_pm_kg_m,
            "operating_moment_design_Nm": cases["U_op"]["M_design"],
            "parked_fault_moment_design_Nm": cases["U_park_fault"]["M_design"],
            "operating_deflection_m": static_results[0].displacement_m,
            "parked_fault_deflection_m": static_results[2].displacement_m,
        }).to_dict("records"),
        "verification_summary": {
            "open_fail_count": int((check_table.status == "FAIL").sum()),
            "engineering_pass_count": int((check_table.status == "PASS").sum()),
            "information_count": int((check_table.status == "INFO").sum()),
        },
    }


def _json_safe(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def write_phase8_outputs() -> dict[str, object]:
    result = run_phase8()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(_json_safe(result), indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    result = write_phase8_outputs()
    checks = pd.DataFrame(result["checks"])
    print(checks[["category", "check_id", "value", "criterion", "status"]].to_string(index=False))
    print(f"\nFinal v6 mass: {result['mass_breakdown']['final_v6_blade_mass_kg']:.1f} kg")
    print(f"Open FAIL: {result['verification_summary']['open_fail_count']}")
