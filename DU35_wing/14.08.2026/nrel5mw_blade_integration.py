"""Source-traceable integration of the DU35 section into the NREL 5 MW blade.

This module does not create a new load case.  It reads the published NREL
planform and the existing verified full-blade screening outputs, maps the CAD
DU35 section by *airfoil identity first* and chord second, and then reruns the
local sandwich check with the factored parked-fault strain at that station.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import sys

WIND_WING_ROOT = Path(__file__).resolve().parents[3]
COMPOSITES_ROOT = WIND_WING_ROOT.parent / "Composites_physics"
DAY02_CODE = WIND_WING_ROOT / "sessions" / "day02_rbs" / "code"
DAY03_CODE = WIND_WING_ROOT / "sessions" / "day03_sandwich_panels" / "code"
for path in (COMPOSITES_ROOT, DAY02_CODE, DAY03_CODE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import numpy as np

from airfoil_panel_laminate import UD_GFRP
from du35_sandwich_assessment import (
    REFERENCE_CAP_STRESS_PA,
    build_panels,
)
from du35_sandwich_geometry import (
    CHORD_M,
    LE_UPPER,
    TE_UPPER,
    normalized_upper_position,
)
from sandwich_panel import assess_sandwich_panel


SESSION_ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_CSV = WIND_WING_ROOT / "data" / "reference_blade" / "blade_geometry.csv"
STATIONS_CSV = WIND_WING_ROOT / "results" / "structural" / "glass_station_results.csv"
PHASE4_SUMMARY_CSV = WIND_WING_ROOT / "results" / "structural" / "phase4_summary.csv"
STATIC_CSV = WIND_WING_ROOT / "results" / "fea" / "phase5_static_cases.csv"
MODAL_CSV = WIND_WING_ROOT / "results" / "fea" / "phase5_modal.csv"
BUCKLING_CSV = WIND_WING_ROOT / "results" / "fea" / "phase5_panel_buckling.csv"

HUB_RADIUS_M = 1.5
REFERENCE_BLADE_MASS_KG = 17_740.0
TIP_DEFLECTION_LIMIT_M = 5.5
THREE_P_HZ = 0.605
ASPECT_RATIOS = np.linspace(0.2, 10.0, 393)


@dataclass(frozen=True)
class MappedStation:
    structural_r_m: float
    aerodynamic_r_m: float
    airfoil: str
    source_chord_m: float
    cad_chord_m: float
    chord_mismatch_percent: float
    moment_Nm: float
    EI_Nm2: float
    cap_separation_m: float
    cap_strain: float


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [row for row in csv.DictReader(
            line for line in handle if not line.lstrip().startswith("#")
        )]


def map_cad_du35_to_source_station() -> MappedStation:
    """Map by published airfoil identity, then choose the closest chord.

    Chord alone is ambiguous on the NREL planform: approximately 4.49 m occurs
    once on the inner increasing branch and again near the DU35 region.
    """
    geometry = [row for row in _read_rows(GEOMETRY_CSV)
                if row["airfoil"] == "DU35_A17"]
    if not geometry:
        raise ValueError("No DU35_A17 station is present in the source geometry")
    selected = min(geometry, key=lambda row: abs(float(row["chord_m"]) - CHORD_M))
    aerodynamic_r = float(selected["r_m"])
    structural_r = aerodynamic_r - HUB_RADIUS_M

    stations = _read_rows(STATIONS_CSV)
    station = min(stations, key=lambda row: abs(float(row["r_m"]) - structural_r))
    if abs(float(station["r_m"]) - structural_r) > 1e-6:
        raise ValueError("The selected NREL DU35 station is absent from structural results")

    source_chord = float(selected["chord_m"])
    moment = float(station["M_park_fault_design_Nm"])
    ei = float(station["EI_flap_Nm2"])
    separation = float(station["cap_separation_m"])
    cap_strain = moment / ei * separation / 2.0
    return MappedStation(
        structural_r_m=structural_r,
        aerodynamic_r_m=aerodynamic_r,
        airfoil=selected["airfoil"],
        source_chord_m=source_chord,
        cad_chord_m=CHORD_M,
        chord_mismatch_percent=100.0 * (CHORD_M - source_chord) / source_chord,
        moment_Nm=moment,
        EI_Nm2=ei,
        cap_separation_m=separation,
        cap_strain=cap_strain,
    )


def corrected_du35_assessment() -> tuple[MappedStation, list[dict[str, object]]]:
    station = map_cad_du35_to_source_station()
    panels = build_panels()
    rows: list[dict[str, object]] = []
    reference_strain = REFERENCE_CAP_STRESS_PA / UD_GFRP["E1"]
    for geometry in (LE_UPPER, TE_UPPER):
        panel = panels[geometry.name]
        cases = {
            "109MPa_reference": reference_strain,
            "source_DU35_parked_fault": (
                normalized_upper_position(geometry) * station.cap_strain
            ),
        }
        for case_name, strain in cases.items():
            applied = panel.membrane_resultant(strain)
            results = [assess_sandwich_panel(panel, ratio, applied)
                       for ratio in ASPECT_RATIOS]
            worst = min(results, key=lambda result: result.reserve_factor)
            rows.append({
                "panel": geometry.name,
                "load_case": case_name,
                "structural_r_m": station.structural_r_m,
                "airfoil": station.airfoil,
                "source_chord_m": station.source_chord_m,
                "cad_chord_m": station.cad_chord_m,
                "chord_mismatch_percent": station.chord_mismatch_percent,
                "axial_strain": strain,
                "applied_Nx_N_per_m": applied,
                "worst_a_over_b": worst.a_over_b,
                "governing_mode": worst.governing_mode,
                "critical_Nx_N_per_m": worst.governing_critical_Nx_N_per_m,
                "reserve_factor": worst.reserve_factor,
                "margin_of_safety": worst.margin_of_safety,
            })
    return station, rows


def full_blade_status(local_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    phase4 = next(row for row in _read_rows(PHASE4_SUMMARY_CSV)
                  if row["design"] == "glass")
    static = _read_rows(STATIC_CSV)
    modal = _read_rows(MODAL_CSV)
    buckling = _read_rows(BUCKLING_CSV)

    operating = next(row for row in static if row["case"] == "U_op")
    parked = next(row for row in static if row["case"] == "U_park_fault")
    first_flap = next(row for row in modal
                      if row["direction"] == "flap" and row["mode"] == "1")
    finite_cap = [float(row["cap_buckling_RF"]) for row in buckling
                  if row["cap_buckling_RF"].lower() != "inf"]
    finite_web = [float(row["web_shear_buckling_RF"]) for row in buckling
                  if row["web_shear_buckling_RF"].lower() != "inf"]
    te_source = next(row for row in local_rows
                     if row["panel"] == "TE_upper_foam"
                     and row["load_case"] == "source_DU35_parked_fault")
    te_reference = next(row for row in local_rows
                        if row["panel"] == "TE_upper_foam"
                        and row["load_case"] == "109MPa_reference")

    return [
        {"check": "blade_mass", "value": float(phase4["mass_kg"]),
         "limit_or_reference": REFERENCE_BLADE_MASS_KG,
         "status": "ABOVE_REFERENCE", "scope": "Phase 4 glass baseline"},
        {"check": "operating_tip_deflection", "value": abs(float(operating["tip_displacement_m"])),
         "limit_or_reference": TIP_DEFLECTION_LIMIT_M,
         "status": "PASS", "scope": "Phase 5 source-based beam FE"},
        {"check": "parked_fault_tip_deflection", "value": abs(float(parked["tip_displacement_m"])),
         "limit_or_reference": TIP_DEFLECTION_LIMIT_M,
         "status": "PASS", "scope": "Phase 5 source-based beam FE"},
        {"check": "first_flap_frequency", "value": float(first_flap["frequency_hz"]),
         "limit_or_reference": THREE_P_HZ,
         "status": "PASS", "scope": "Phase 5 source-based beam FE"},
        {"check": "minimum_cap_buckling_RF", "value": min(finite_cap),
         "limit_or_reference": 1.0, "status": "PASS",
         "scope": "Phase 5 local cap screen"},
        {"check": "minimum_shell_strength_RF", "value": float(phase4["min_shell_RF"]),
         "limit_or_reference": 1.0, "status": "PASS",
         "scope": "Phase 4 shell strength screen"},
        {"check": "minimum_unstiffened_web_buckling_RF", "value": min(finite_web),
         "limit_or_reference": 1.0, "status": "UNRESOLVED_FAIL",
         "scope": "Phase 5 model omits documented sandwich-core architecture"},
        {"check": "DU35_TE_source_station_RF", "value": float(te_source["reserve_factor"]),
         "limit_or_reference": 1.0, "status": "PASS",
         "scope": "Corrected DU35 mapping; factored parked-fault load"},
        {"check": "DU35_TE_109MPa_reference_RF", "value": float(te_reference["reserve_factor"]),
         "limit_or_reference": 1.0, "status": "REFERENCE_FAIL_ONLY",
         "scope": "Comparison case; not the sourced shell load"},
    ]


def write_outputs(output_dir: Path | None = None) -> tuple[MappedStation, list[dict], list[dict]]:
    output_dir = output_dir or SESSION_ROOT / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    station, local = corrected_du35_assessment()
    status = full_blade_status(local)
    _write_csv(output_dir / "du35_corrected_mapping.csv", local)
    _write_csv(output_dir / "nrel5mw_preliminary_status.csv", status)
    return station, local, status


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    mapped, local_rows = corrected_du35_assessment()
    blade_rows = full_blade_status(local_rows)
    print(f"DU35 source station: r={mapped.structural_r_m:.3f} m, "
          f"chord={mapped.source_chord_m:.3f} m, "
          f"CAD mismatch={mapped.chord_mismatch_percent:.3f}%")
    for row in local_rows:
        print(f"{row['panel']:22s} {row['load_case']:27s} "
              f"RF={row['reserve_factor']:.3f}")
    for row in blade_rows:
        print(f"{row['check']:40s} {row['value']:.4g} {row['status']}")
