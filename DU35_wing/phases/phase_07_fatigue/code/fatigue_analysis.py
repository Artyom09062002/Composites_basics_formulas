"""Phase 7 fatigue screening for the 61.5 m glass blade baseline.

The calculation deliberately does not invent a time history.  Published
20-year Miner damage from SAND2013-2569 Table 26 is transferred to the current
glass/CAD baseline.  The published NREL 60 m fatigue-target distribution is
used only to extend the reference stations outboard and to derive an
equivalent flapwise shear range for the two CAD webs.

This is a screening model, not a certification fatigue substantiation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


PHASE_ROOT = Path(__file__).resolve().parents[1]
PHASES_ROOT = PHASE_ROOT.parent
STATIONS = PHASES_ROOT / "phase_04_structural_design" / "results" / "glass_station_results.csv"
STIFFNESS = PHASES_ROOT / "phase_02_aerodynamics" / "data" / "reference_blade" / "blade_stiffness.csv"
REFERENCE_DAMAGE = PHASE_ROOT / "data" / "sandia_reference_damage.csv"
TARGETS = PHASE_ROOT / "data" / "nrel_60m_fatigue_targets.csv"
WEB_GEOMETRY = (PHASES_ROOT / "phase_04_structural_design" / "studies"
                / "full_blade_webs" / "data" / "cad_web_geometry.csv")

DESIGN_LIFE_YEARS = 20.0
FATIGUE_FACTOR = 1.38
FLAP_EQUIVALENT_CYCLES = 1_000_000.0
EDGE_EQUIVALENT_CYCLES = 2_000_000.0
ROOT_FLAP_RANGE_NM = 8_741_000.0
ROOT_EDGE_RANGE_NM = 10_500_000.0


@dataclass(frozen=True)
class SNCurve:
    name: str
    b: float
    C_MPa: float
    E_GPa: float
    source: str


GLASS_UD = SNCurve("E-LT-5500 UD", 10.0, 1000.0, 41.8, "SAND2013-2569 Table 24")
CARBON_UD = SNCurve("Newport 307 UD", 14.0, 1546.0, 114.5, "SAND2013-2569 Tables 5/24")
TRIAX = SNCurve("SNL Triax", 10.0, 700.0, 27.7, "SAND2013-2569 Table 24")
# No published Saertex DB shear S-N curve is available in the project inputs.
# Static shear strength plus the GL glass slope is a deliberately conservative
# replaceable screening proxy; it is never labelled as validated material data.
DB_SHEAR_PROXY = SNCurve("Saertex DB shear proxy", 10.0, 62.0, 13.1,
                         "Conservative design envelope: C=62 MPa from MAT-002 "
                         "static shear strength and b=10")
# Mandell, Miller & Samborsky (AIAA 2013-1630), Figure 7, Fabric L:
# The Saertex VU 90079 R=-1 trend was digitized from Figure 7 as
# applied axial S=156*N^-0.103 MPa; it is not a stated numeric source fit.
# Their Table 2 / ASTM D3518 relation gives ply shear = applied stress / 2,
# hence the digitized reference shear trend is S6=78*N^-0.103 MPa.
SAERTEX_VU90079_MEASURED_RM1 = SNCurve(
    "Saertex VU 90079 digitized R=-1 shear reference",
    1.0 / 0.103,
    78.0,
    13.1,
    "Digitized from Mandell, Miller & Samborsky 2013, AIAA-2013-1630, "
    "Fig. 7 Fabric L; Table 2 gives ASTM D3518 shear = axial/2",
)
# The robust design envelope deliberately remains more severe than that
# digitized curve: C=62 MPa and b=9. It is not fitted to make the design pass.
DB_ROBUST_B9 = SNCurve(
    "Saertex DB conservative robust envelope",
    9.0,
    62.0,
    13.1,
    "Estimated conservative envelope: digitized Saertex VU 90079 R=-1 "
    "C_shear=78 MPa, b=9.709 rounded down to C=62 MPa, b=9",
)
# Bednarcyk et al., NASA/TM-2012-217694, Sec. 4 / Fig. 6 report the
# Divinycell H100 shear-fatigue fit Smax=2.34*N^(-1/12.08) MPa at R=0.1.
# This replaces the earlier estimated bound.  The project fatigue factor is
# still applied, so the implementation remains conservative relative to fit.
H100_CORE_CURVE = SNCurve(
    "Divinycell H100 measured shear-fatigue curve",
    12.08,
    2.34,
    0.040,
    "NASA/TM-2012-217694 Section 4, Table 3 and Figure 6; R=0.1",
)

ROBUST_DAMAGE_TARGET = 0.70
CORE_THICKNESS_CANDIDATES_M = (0.050, 0.060, 0.080, 0.100)
ROBUST_CASES = (
    ("nominal_b10_split0", DB_SHEAR_PROXY, 0.0),
    ("b9_split0", DB_ROBUST_B9, 0.0),
    ("b10_split_minus10", DB_SHEAR_PROXY, -10.0),
    ("b10_split_plus10", DB_SHEAR_PROXY, 10.0),
    ("b9_split_minus10", DB_ROBUST_B9, -10.0),
    ("b9_split_plus10", DB_ROBUST_B9, 10.0),
)


def miner_damage(stress_ranges_MPa: np.ndarray, cycle_counts: np.ndarray,
                 curve: SNCurve, factor: float = FATIGUE_FACTOR) -> float:
    """Palmgren-Miner sum for N=(C/(factor*S))**b."""
    stress = np.asarray(stress_ranges_MPa, dtype=float)
    cycles = np.asarray(cycle_counts, dtype=float)
    if np.any(stress < 0) or np.any(cycles < 0):
        raise ValueError("stress ranges and cycle counts must be non-negative")
    return float(np.sum(cycles * np.power(factor * stress / curve.C_MPa, curve.b)))


def equivalent_stress_MPa(damage: float, cycles: float, curve: SNCurve,
                          factor: float = FATIGUE_FACTOR) -> float:
    if damage < 0 or cycles <= 0:
        raise ValueError("damage must be non-negative and cycles positive")
    if damage == 0:
        return 0.0
    return float(curve.C_MPa / factor * (damage / cycles) ** (1.0 / curve.b))


def transfer_damage_by_strain(damage: float, cycles: float,
                              source_curve: SNCurve, target_curve: SNCurve) -> float:
    """Transfer a published equivalent damage through the same strain range."""
    source_stress = equivalent_stress_MPa(damage, cycles, source_curve)
    strain = source_stress / (source_curve.E_GPa * 1000.0)
    target_stress = strain * target_curve.E_GPa * 1000.0
    return miner_damage(np.array([target_stress]), np.array([cycles]), target_curve)


def stiffness_based_web_split(a66_forward: float, h_forward: float,
                              a66_aft: float, h_aft: float) -> tuple[float, float]:
    """Compatible-shear split, V_i proportional to A66_i*h_i."""
    kf = a66_forward * h_forward
    ka = a66_aft * h_aft
    if min(kf, ka) <= 0:
        raise ValueError("web shear stiffness terms must be positive")
    total = kf + ka
    return kf / total, ka / total


def _log_interp(x: np.ndarray, xp: np.ndarray, damage: np.ndarray) -> np.ndarray:
    return np.power(10.0, np.interp(x, xp, np.log10(damage)))


def _target_range(r_m: np.ndarray, column: str, root_range_Nm: float) -> np.ndarray:
    targets = pd.read_csv(TARGETS)
    fraction = np.clip(np.asarray(r_m) / 60.0, 0.0, 1.0)
    normalized = np.interp(fraction, targets.span_fraction, targets[column])
    return root_range_Nm * normalized


def _evaluate_web_pair(
    span: float,
    total_v: float,
    pair: pd.DataFrame,
    forward_layers: int,
    aft_layers: int,
    core_thickness_m: float,
    design_variant: str,
    curve: SNCurve = DB_SHEAR_PROXY,
    split_error_percent: float = 0.0,
    split_method: str = "stiffness_based",
    case_id: str = "nominal_b10_split0",
) -> list[dict[str, object]]:
    """Evaluate one CAD web pair with independent face layups.

    The split perturbation scales the nominal forward A66h fraction and gives
    the remainder to the aft web, so total section shear is conserved.
    """
    face_G = 11.8e9
    core_G = 22e6  # Sandia project material, retained from the accepted CAD basis.
    layers = {"forward": int(forward_layers), "aft": int(aft_layers)}
    heights = {web: float(pair.loc[web, "web_height_m"])
               for web in ("forward", "aft")}
    a66 = {
        web: 2.0 * face_G * layers[web] * 0.001 + core_G * core_thickness_m
        for web in ("forward", "aft")
    }
    nominal_f, nominal_a = stiffness_based_web_split(
        a66["forward"], heights["forward"], a66["aft"], heights["aft"]
    )
    if split_method == "equal_50_50":
        nominal_f, nominal_a = 0.5, 0.5
    fraction_f = float(np.clip(
        nominal_f * (1.0 + split_error_percent / 100.0), 1e-9, 1.0 - 1e-9
    ))
    fractions = {"forward": fraction_f, "aft": 1.0 - fraction_f}
    nominal_fractions = {"forward": nominal_f, "aft": nominal_a}

    rows: list[dict[str, object]] = []
    for web in ("forward", "aft"):
        v_web = total_v * fractions[web]
        nxy = v_web / heights[web]
        gamma = nxy / a66[web]
        face_tau_MPa = face_G * gamma / 1e6
        core_tau_MPa = core_G * gamma / 1e6
        face_damage = miner_damage(
            np.array([face_tau_MPa]), np.array([FLAP_EQUIVALENT_CYCLES]), curve
        )
        core_damage = miner_damage(
            np.array([core_tau_MPa]), np.array([FLAP_EQUIVALENT_CYCLES]), H100_CORE_CURVE
        )
        rows.append({
            "case_id": case_id,
            "r_m": span,
            "web": web,
            "design_variant": design_variant,
            "split_method": split_method,
            "split_error_percent": split_error_percent,
            "nominal_load_fraction": nominal_fractions[web],
            "load_fraction": fractions[web],
            "face_layers_each_side": layers[web],
            "core_thickness_m": core_thickness_m,
            "A66_N_per_m": a66[web],
            "A66h_N": a66[web] * heights[web],
            "equivalent_total_shear_range_N": total_v,
            "face_shear_range_MPa": face_tau_MPa,
            "curve_C_MPa": curve.C_MPa,
            "curve_b": curve.b,
            "curve_source": curve.source,
            "face_damage_20y_proxy": face_damage,
            "face_status": "PASS" if face_damage <= ROBUST_DAMAGE_TARGET else "FAIL",
            "core_shear_range_MPa": core_tau_MPa,
            "core_damage_20y": core_damage,
            "core_fatigue_status": (
                "SOURCED_CURVE_PASS" if core_damage <= ROBUST_DAMAGE_TARGET
                else "SOURCED_CURVE_FAIL"
            ),
        })
    return rows


def _candidate_mass_kg(web_geom: pd.DataFrame, forward_layers: int,
                       aft_layers: int, core_thickness_m: float) -> float:
    """Approximate two-web face+core material mass for candidate ranking."""
    face_density = 1830.0
    core_density = 100.0
    mass = 0.0
    for web, count in (("forward", forward_layers), ("aft", aft_layers)):
        part = web_geom[web_geom.web == web].sort_values("blade_span_m")
        area = float(np.trapz(part.web_height_m, part.blade_span_m))
        areal_mass = 2.0 * count * 0.001 * face_density + core_thickness_m * core_density
        mass += area * areal_mass
    return mass


def optimize_robust_web_design(
    web_geom: pd.DataFrame,
    unique_r: np.ndarray,
    total_shear_range: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Select the lightest one-zone-per-web design that closes every robust case."""
    candidates: list[dict[str, object]] = []
    detail_by_candidate: dict[tuple[int, int, float], pd.DataFrame] = {}
    for forward_layers in range(2, 13):
        for aft_layers in range(2, 13):
            for core_t in CORE_THICKNESS_CANDIDATES_M:
                rows: list[dict[str, object]] = []
                for span, total_v in zip(unique_r, total_shear_range):
                    pair = web_geom[np.isclose(web_geom.blade_span_m, span)].set_index("web")
                    for case_id, curve, split_error in ROBUST_CASES:
                        rows.extend(_evaluate_web_pair(
                            float(span), float(total_v), pair,
                            forward_layers, aft_layers, core_t,
                            "fatigue_robust_final", curve, split_error,
                            "stiffness_based", case_id,
                        ))
                detail = pd.DataFrame(rows)
                max_face = float(detail.face_damage_20y_proxy.max())
                max_core = float(detail.core_damage_20y.max())
                max_split_imbalance = float(
                    (detail.nominal_load_fraction - 0.5).abs().max()
                )
                key = (forward_layers, aft_layers, core_t)
                detail_by_candidate[key] = detail
                candidates.append({
                    "forward_layers_each_side": forward_layers,
                    "aft_layers_each_side": aft_layers,
                    "core_thickness_m": core_t,
                    "spanwise_zones_per_web": 1,
                    "max_adjacent_layer_step": 0,
                    "estimated_web_material_mass_kg": _candidate_mass_kg(
                        web_geom, forward_layers, aft_layers, core_t
                    ),
                    "max_nominal_split_imbalance": max_split_imbalance,
                    "max_face_D20_all_cases": max_face,
                    "max_core_D20_all_cases": max_core,
                    "status": "PASS" if max(max_face, max_core) <= ROBUST_DAMAGE_TARGET else "FAIL",
                })
    candidate_table = pd.DataFrame(candidates)
    feasible = candidate_table[candidate_table.status == "PASS"].sort_values(
        ["spanwise_zones_per_web", "estimated_web_material_mass_kg",
         "max_nominal_split_imbalance", "max_face_D20_all_cases"]
    )
    if feasible.empty:
        raise RuntimeError("No one-zone DB/core scheme closes the robust D20 target")
    selected = feasible.iloc[0]
    key = (int(selected.forward_layers_each_side),
           int(selected.aft_layers_each_side), float(selected.core_thickness_m))
    detail = detail_by_candidate[key].copy()
    schedule = pd.DataFrame([
        {"web": "forward", "start_m": float(unique_r.min()), "end_m": float(unique_r.max()),
         "face_layers_each_side": key[0], "core_thickness_m": key[2]},
        {"web": "aft", "start_m": float(unique_r.min()), "end_m": float(unique_r.max()),
         "face_layers_each_side": key[1], "core_thickness_m": key[2]},
    ])
    candidate_table["selected"] = False
    selected_mask = (
        (candidate_table.forward_layers_each_side == key[0]) &
        (candidate_table.aft_layers_each_side == key[1]) &
        np.isclose(candidate_table.core_thickness_m, key[2])
    )
    candidate_table.loc[selected_mask, "selected"] = True
    return schedule, detail, candidate_table


def _extend_damage(r: np.ndarray, reference_r: np.ndarray, reference_damage: np.ndarray,
                   stress_shape: np.ndarray, b: float) -> np.ndarray:
    """Interpolate published points; extend outboard by stress-range ratio^b."""
    out = _log_interp(np.minimum(r, reference_r[-1]), reference_r, reference_damage)
    anchor_shape = float(np.interp(reference_r[-1], r, stress_shape))
    mask = r > reference_r[-1]
    if anchor_shape > 0:
        out[mask] = reference_damage[-1] * np.power(
            np.maximum(stress_shape[mask], 0.0) / anchor_shape, b
        )
    return out


def assess_phase7() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    stations = pd.read_csv(STATIONS)
    stiffness = pd.read_csv(STIFFNESS, comment="#")
    reference = pd.read_csv(REFERENCE_DAMAGE)
    r = stations.r_m.to_numpy(float)
    chord = stations.chord_m.to_numpy(float)
    ei_flap = stations.EI_flap_Nm2.to_numpy(float)
    ei_edge = stiffness.EI_edge_Nm2.to_numpy(float)
    flap_range = _target_range(r, "dynamic_flap", ROOT_FLAP_RANGE_NM)
    edge_range = _target_range(r, "dynamic_lead_lag", ROOT_EDGE_RANGE_NM)
    flap_strain_shape = np.divide(flap_range * stations.cap_separation_m.to_numpy(float) / 2.0,
                                  ei_flap, out=np.zeros_like(r), where=ei_flap > 0)
    edge_strain_shape = np.divide(edge_range * 0.45 * chord, ei_edge,
                                  out=np.zeros_like(r), where=ei_edge > 0)

    rr = reference.r_m.to_numpy(float)
    cap_ref_glass = np.array([
        transfer_damage_by_strain(d, FLAP_EQUIVALENT_CYCLES, CARBON_UD, GLASS_UD)
        for d in reference.primary_ud_flap_damage.to_numpy(float)
    ])
    cap_damage = _extend_damage(r, rr, cap_ref_glass, flap_strain_shape, GLASS_UD.b)
    shell_flap = _extend_damage(r, rr, reference.snl_triax_flap_damage.to_numpy(float),
                                flap_strain_shape, TRIAX.b)
    shell_edge = _extend_damage(r, rr, reference.snl_triax_edge_damage.to_numpy(float),
                                edge_strain_shape, TRIAX.b)

    component_rows: list[dict[str, object]] = []
    for i, span in enumerate(r):
        if 6.2 <= span <= 52.3:
            component_rows.append({
                "r_m": span, "component": "spar_cap_glass_UD",
                "damage_20y": cap_damage[i], "life_years_linear": DESIGN_LIFE_YEARS / cap_damage[i],
                "status": "PASS" if cap_damage[i] <= 1.0 else "FAIL",
                "basis": "SAND carbon-cap damage transferred by equivalent strain to E-LT-5500",
            })
        shell_total = shell_flap[i] + shell_edge[i]
        component_rows.append({
            "r_m": span, "component": "shell_and_TE_SNL_Triax",
            "damage_20y": shell_total,
            "life_years_linear": DESIGN_LIFE_YEARS / shell_total if shell_total > 0 else float("inf"),
            "status": "PASS" if shell_total <= 1.0 else "FAIL",
            "basis": "SAND Table 26 flap+edge Miner damage; outboard target-shape extension",
        })
    components = pd.DataFrame(component_rows)

    web_geom = pd.read_csv(WEB_GEOMETRY)
    web_rows: list[dict[str, object]] = []
    unique_r = np.sort(web_geom.blade_span_m.unique())
    web_moment = _target_range(unique_r, "dynamic_flap", ROOT_FLAP_RANGE_NM)
    total_shear_range = np.abs(np.gradient(web_moment, unique_r, edge_order=2))
    for span, total_v in zip(unique_r, total_shear_range):
        pair = web_geom[np.isclose(web_geom.blade_span_m, span)].set_index("web")
        baseline_layers = 3 if 10.25 <= span <= 38.95 else 2
        web_rows.extend(_evaluate_web_pair(
            float(span), float(total_v), pair, baseline_layers, baseline_layers,
            0.050, "day4_baseline", DB_SHEAR_PROXY, 0.0,
            "stiffness_based", "nominal_b10_split0",
        ))
        web_rows.extend(_evaluate_web_pair(
            float(span), float(total_v), pair, baseline_layers, baseline_layers,
            0.050, "day4_baseline", DB_SHEAR_PROXY, 0.0,
            "equal_50_50", "nominal_b10_equal50",
        ))
        recommended_layers = baseline_layers
        while recommended_layers < 12:
            candidate = _evaluate_web_pair(
                float(span), float(total_v), pair,
                recommended_layers, recommended_layers, 0.050,
                "fatigue_reinforced_legacy", DB_SHEAR_PROXY, 0.0,
                "stiffness_based", "nominal_b10_split0",
            )
            if max(float(row["face_damage_20y_proxy"]) for row in candidate) <= 1.0:
                web_rows.extend(candidate)
                break
            recommended_layers += 1
        else:
            raise RuntimeError(f"No web-face fatigue solution found at r={span:g} m")
    schedule, robust_detail, candidate_table = optimize_robust_web_design(
        web_geom, unique_r, total_shear_range
    )
    robust_nominal = robust_detail[robust_detail.case_id == "nominal_b10_split0"]
    web_rows.extend(robust_nominal.to_dict("records"))
    webs = pd.DataFrame(web_rows)

    all_quantified = pd.concat([
        components[["r_m", "component", "damage_20y", "status"]],
        webs[(webs.split_method == "stiffness_based") &
             (webs.design_variant == "fatigue_robust_final")].rename(
            columns={"face_damage_20y_proxy": "damage_20y", "face_status": "status"}
        ).assign(component=lambda d: "web_face_" + d.web)[["r_m", "component", "damage_20y", "status"]]
    ], ignore_index=True)
    webs.attrs["robust_schedule"] = schedule
    webs.attrs["robust_detail"] = robust_detail
    webs.attrs["robust_candidates"] = candidate_table
    return components, webs, all_quantified


def db_proxy_sensitivity(webs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate the fixed reinforced design across C/b proxy variations."""
    design = webs[(webs.design_variant == "fatigue_robust_final") &
                  (webs.split_method == "stiffness_based")].copy()
    rows: list[dict[str, object]] = []
    for c_change in (-20.0, -10.0, 0.0, 10.0, 20.0):
        c_value = DB_SHEAR_PROXY.C_MPa * (1.0 + c_change / 100.0)
        for b_value in (8.0, 9.0, 10.0, 11.0, 12.0):
            curve = SNCurve("DB sensitivity", b_value, c_value,
                            DB_SHEAR_PROXY.E_GPa, "parametric sensitivity")
            damages = design.face_shear_range_MPa.map(
                lambda stress: miner_damage(np.array([stress]),
                                             np.array([FLAP_EQUIVALENT_CYCLES]), curve)
            )
            idx = damages.idxmax()
            rows.append({
                "C_change_percent": c_change,
                "C_MPa": c_value,
                "b": b_value,
                "max_D20": float(damages.loc[idx]),
                "status": "PASS" if float(damages.loc[idx]) <= 1.0 else "FAIL",
                "governing_r_m": float(design.loc[idx, "r_m"]),
                "governing_web": str(design.loc[idx, "web"]),
            })
    sweep = pd.DataFrame(rows)
    max_stress = float(design.face_shear_range_MPa.max())
    thresholds: list[dict[str, object]] = []
    for b_value in (8.0, 9.0, 10.0, 11.0, 12.0):
        critical_c = FATIGUE_FACTOR * max_stress * FLAP_EQUIVALENT_CYCLES ** (1.0 / b_value)
        degradation = 100.0 * (1.0 - critical_c / DB_SHEAR_PROXY.C_MPa)
        thresholds.append({
            "b": b_value,
            "critical_C_MPa_at_D20_1": critical_c,
            "allowable_C_degradation_percent": degradation,
            "interpretation": ("current C already fails" if degradation < 0
                               else "maximum C decrease before failure"),
        })
    return sweep, pd.DataFrame(thresholds)


def load_split_sensitivity(webs: pd.DataFrame) -> pd.DataFrame:
    """Perturb the forward A66h fraction while conserving total section shear."""
    design = webs[(webs.design_variant == "fatigue_robust_final") &
                  (webs.split_method == "stiffness_based")].copy()
    rows: list[dict[str, object]] = []
    for error in (-10.0, 0.0, 10.0):
        candidate_rows: list[dict[str, object]] = []
        for span, pair in design.groupby("r_m"):
            by_web = pair.set_index("web")
            f0 = float(by_web.loc["forward", "load_fraction"])
            f = float(np.clip(f0 * (1.0 + error / 100.0), 1e-9, 1.0 - 1e-9))
            fractions = {"forward": f, "aft": 1.0 - f}
            for web in ("forward", "aft"):
                nominal_fraction = float(by_web.loc[web, "load_fraction"])
                nominal_stress = float(by_web.loc[web, "face_shear_range_MPa"])
                stress = nominal_stress * fractions[web] / nominal_fraction
                damage = miner_damage(np.array([stress]),
                                      np.array([FLAP_EQUIVALENT_CYCLES]), DB_SHEAR_PROXY)
                candidate_rows.append({
                    "forward_split_error_percent": error,
                    "r_m": float(span), "web": web,
                    "load_fraction": fractions[web],
                    "face_shear_range_MPa": stress, "D20": damage,
                })
        candidate = pd.DataFrame(candidate_rows)
        governing = candidate.loc[candidate.D20.idxmax()]
        rows.append({
            "forward_split_error_percent": error,
            "governing_r_m": float(governing.r_m),
            "governing_web": str(governing.web),
            "governing_load_fraction": float(governing.load_fraction),
            "max_D20": float(governing.D20),
            "status": "PASS" if float(governing.D20) <= 1.0 else "FAIL",
        })
    return pd.DataFrame(rows)


def robust_case_summary(webs: pd.DataFrame) -> pd.DataFrame:
    detail = webs.attrs["robust_detail"]
    rows: list[dict[str, object]] = []
    for case_id, group in detail.groupby("case_id", sort=False):
        governing = group.loc[group.face_damage_20y_proxy.idxmax()]
        rows.append({
            "case_id": case_id,
            "curve_C_MPa": float(governing.curve_C_MPa),
            "curve_b": float(governing.curve_b),
            "split_error_percent": float(governing.split_error_percent),
            "governing_r_m": float(governing.r_m),
            "governing_web": str(governing.web),
            "max_face_D20": float(governing.face_damage_20y_proxy),
            "target_D20": ROBUST_DAMAGE_TARGET,
            "margin_to_target": ROBUST_DAMAGE_TARGET - float(governing.face_damage_20y_proxy),
            "status": "PASS" if float(governing.face_damage_20y_proxy) <= ROBUST_DAMAGE_TARGET else "FAIL",
        })
    return pd.DataFrame(rows)


def write_phase7_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = PHASE_ROOT / "results"
    out.mkdir(parents=True, exist_ok=True)
    components, webs, all_quantified = assess_phase7()
    components.to_csv(out / "component_damage.csv", index=False)
    webs.to_csv(out / "web_fatigue_comparison.csv", index=False)
    all_quantified.to_csv(out / "phase7_summary.csv", index=False)
    db_sweep, db_thresholds = db_proxy_sensitivity(webs)
    db_sweep.to_csv(out / "db_proxy_sensitivity.csv", index=False)
    db_thresholds.to_csv(out / "db_proxy_thresholds.csv", index=False)
    load_split_sensitivity(webs).to_csv(out / "a66h_split_sensitivity.csv", index=False)
    robust_case_summary(webs).to_csv(out / "robust_design_cases.csv", index=False)
    webs.attrs["robust_schedule"].to_csv(out / "robust_web_schedule.csv", index=False)
    webs.attrs["robust_detail"].to_csv(out / "robust_design_station_results.csv", index=False)
    webs.attrs["robust_candidates"].to_csv(out / "robust_design_candidates.csv", index=False)
    webs[["r_m", "web", "design_variant", "split_method", "face_layers_each_side",
          "core_thickness_m", "core_shear_range_MPa", "core_damage_20y",
          "core_fatigue_status"]].to_csv(out / "core_fatigue_curve_results.csv", index=False)
    pd.DataFrame([
        {"curve": DB_SHEAR_PROXY.name, "C_MPa": DB_SHEAR_PROXY.C_MPa,
         "b": DB_SHEAR_PROXY.b, "status": "estimated conservative design envelope",
         "source": DB_SHEAR_PROXY.source},
        {"curve": DB_ROBUST_B9.name, "C_MPa": DB_ROBUST_B9.C_MPa,
         "b": DB_ROBUST_B9.b, "status": "estimated conservative robust envelope",
         "source": DB_ROBUST_B9.source},
        {"curve": SAERTEX_VU90079_MEASURED_RM1.name,
         "C_MPa": SAERTEX_VU90079_MEASURED_RM1.C_MPa,
         "b": SAERTEX_VU90079_MEASURED_RM1.b, "status": "digitized reference",
         "source": SAERTEX_VU90079_MEASURED_RM1.source},
        {"curve": H100_CORE_CURVE.name, "C_MPa": H100_CORE_CURVE.C_MPa,
         "b": H100_CORE_CURVE.b, "status": "measured curve",
         "source": H100_CORE_CURVE.source},
    ]).to_csv(out / "fatigue_material_curves.csv", index=False)
    return components, webs, all_quantified


if __name__ == "__main__":
    components, webs, summary = write_phase7_outputs()
    governing = summary.loc[summary.damage_20y.idxmax()]
    print(f"Governing quantified result: {governing.component} at r={governing.r_m:.3f} m, "
          f"D20={governing.damage_20y:.6g}, {governing.status}")
    baseline = webs[webs.design_variant == "day4_baseline"]
    legacy = webs[webs.design_variant == "fatigue_reinforced_legacy"]
    reinforced = webs[webs.design_variant == "fatigue_robust_final"]
    stiff = baseline[baseline.split_method == "stiffness_based"]
    equal = baseline[baseline.split_method == "equal_50_50"]
    print(f"Web max D20 proxy: baseline stiffness={stiff.face_damage_20y_proxy.max():.6g}, "
          f"baseline 50/50={equal.face_damage_20y_proxy.max():.6g}, "
          f"legacy nominal reinforcement={legacy.face_damage_20y_proxy.max():.6g}, "
          f"robust nominal={reinforced.face_damage_20y_proxy.max():.6g}; "
          f"core sourced={reinforced.core_damage_20y.max():.6g}")
    robust = robust_case_summary(webs)
    print(f"Robust all-case maximum D20={robust.max_face_D20.max():.6g} "
          f"against target {ROBUST_DAMAGE_TARGET:.2f}")
    print(webs.attrs["robust_schedule"].to_string(index=False))
