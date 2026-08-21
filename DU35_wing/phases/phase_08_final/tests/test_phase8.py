"""Objective checks for Phase 8 final whole-blade integration."""

from pathlib import Path
import json
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

from phase8_final_integration import (FREQUENCY_TARGET_HZ, FINAL_CORE_T_M,
                                      FINAL_WEB_LAYERS, REFERENCE_MASS_KG,
                                      ROBUST_DAMAGE_TARGET, STATIC_RF_TARGET,
                                      TIP_DEFLECTION_LIMIT_M, OUT_JSON,
                                      run_phase8, write_phase8_outputs)


passed = 0


def check(name: str, condition: bool) -> None:
    global passed
    if not condition:
        raise AssertionError(name)
    passed += 1
    print(f"PASS {passed:02d}: {name}")


result = run_phase8()
checks = pd.DataFrame(result["checks"])
mass = result["mass_breakdown"]
statics = pd.DataFrame(result["static_cases"])
modal = pd.DataFrame(result["modal_results"])
web = pd.DataFrame(result["web_static_detail"])
fatigue = pd.DataFrame(result["fatigue_cases"])

check("authoritative v6 web layups are forward 5 and aft 6 DB per side",
      result["final_design"]["forward_web_DB_layers_each_side"] == FINAL_WEB_LAYERS["forward"]
      and result["final_design"]["aft_web_DB_layers_each_side"] == FINAL_WEB_LAYERS["aft"])
check("authoritative v6 web core is 60 mm H100",
      np.isclose(result["final_design"]["web_core_thickness_m"], FINAL_CORE_T_M))
check("aerodynamic geometry and web positions remain unchanged",
      not result["final_design"]["external_aero_geometry_changed"]
      and not result["final_design"]["web_positions_changed"])
check("raw mass table discrepancy is explicitly retained",
      mass["raw_mass_table_integral_kg"] > 1.05 * REFERENCE_MASS_KG)
check("normalized source mass equals the published 17740 kg basis",
      np.isclose(mass["reference_reported_mass_kg"], REFERENCE_MASS_KG))
check("final mass reconciles reference plus cap and web deltas",
      np.isclose(mass["final_v6_blade_mass_kg"],
                 REFERENCE_MASS_KG + mass["cap_mass_delta_kg"] + mass["web_mass_delta_kg"]))
check("v6 web mass matches the exact CAD-station Phase 7 candidate integral",
      np.isclose(mass["final_v6_web_mass_kg"], 3721.421707486302))
check("final mass is finite and within 25 percent of reference",
      np.isfinite(mass["final_v6_blade_mass_kg"])
      and mass["final_v6_blade_mass_kg"] < 1.25 * REFERENCE_MASS_KG)
check("operating parked fault and gravity cases are integrated",
      {"U_op", "U_park_feather", "U_park_fault", "U_gravity_v6_mass"}
      == set(statics.case))
check("parked-fault factored root moment remains 31.2 MN m within 1 percent",
      np.isclose(abs(float(statics[statics.case == "U_park_fault"].iloc[0].root_moment_Nm)),
                 31.204e6, rtol=0.01))
check("operating tip deflection passes 5.5 m limit",
      abs(float(statics[statics.case == "U_op"].iloc[0].tip_displacement_m))
      <= TIP_DEFLECTION_LIMIT_M)
check("parked-fault tip deflection passes 5.5 m limit",
      abs(float(statics[statics.case == "U_park_fault"].iloc[0].tip_displacement_m))
      <= TIP_DEFLECTION_LIMIT_M)
first_flap = float(modal[(modal.direction == "flap") & (modal["mode"] == 1)].iloc[0].frequency_hz)
check("first flap frequency clears 1.2 times 3P", first_flap >= FREQUENCY_TARGET_HZ)
check("all modal frequencies are positive and ordered",
      all((g.frequency_hz.to_numpy() > 0).all()
          and (np.diff(g.frequency_hz.to_numpy()) > 0).all()
          for _, g in modal.groupby("direction")))
check("cap and shell static strength RF pass",
      checks[checks.check_id.isin(["cap_strength_RF", "shell_strength_RF"])].value.min()
      >= STATIC_RF_TARGET)
check("v6 web strength passes for every station and split case",
      web.status.eq("PASS").all()
      and web[["face_shear_RF", "core_shear_RF", "face_wrinkling_RF",
               "face_compression_RF"]].min().min() >= STATIC_RF_TARGET)
check("v6 sandwich web buckling passes conservative a over h bound",
      web.global_buckling_RF.min() >= STATIC_RF_TARGET)
check("web static screen covers minus nominal and plus ten percent split",
      set(web.split_error_percent) == {-10.0, 0.0, 10.0})
check("each web pair conserves total shear fraction",
      all(np.isclose(group.load_fraction.sum(), 1.0)
          for _, group in web.groupby(["r_m", "split_error_percent"])))
check("all six fatigue sensitivity combinations pass D20 target",
      len(fatigue) == 6 and fatigue.status.eq("PASS").all()
      and fatigue.max_face_D20.max() <= ROBUST_DAMAGE_TARGET)
check("combined b9 plus ten percent remains governing at 0.692872",
      np.isclose(float(fatigue.loc[fatigue.max_face_D20.idxmax(), "max_face_D20"]),
                 0.6928719334586165))
check("every engineering acceptance check passes with no open FAIL",
      not checks.status.eq("FAIL").any()
      and result["verification_summary"]["open_fail_count"] == 0)
check("mass remains an information metric rather than a hidden pass criterion",
      checks[checks.check_id == "final_blade_mass"].iloc[0].status == "INFO")

written = write_phase8_outputs()
check("single compact Phase 8 JSON result is written", OUT_JSON.exists())
reloaded = json.loads(OUT_JSON.read_text(encoding="utf-8"))
check("written result preserves zero open failures",
      reloaded["verification_summary"]["open_fail_count"] == 0
      and written["verification_summary"]["open_fail_count"] == 0)

print(f"\nPhase 8 verification: {passed}/{passed} PASS")
