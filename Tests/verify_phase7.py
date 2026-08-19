"""Objective checks for Phase 7 fatigue screening."""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fatigue_analysis import (CARBON_UD, DB_ROBUST_B9, DB_SHEAR_PROXY,
                              GLASS_UD, H100_CORE_CURVE, ROBUST_DAMAGE_TARGET,
                              SAERTEX_VU90079_MEASURED_RM1, SNCurve, assess_phase7,
                              db_proxy_sensitivity, load_split_sensitivity,
                              equivalent_stress_MPa, miner_damage,
                              robust_case_summary, stiffness_based_web_split,
                              transfer_damage_by_strain)


passed = 0


def check(name: str, condition: bool) -> None:
    global passed
    if not condition:
        raise AssertionError(name)
    passed += 1
    print(f"PASS {passed:02d}: {name}")


curve = SNCurve("benchmark", 10.0, 1000.0, 40.0, "analytical")
check("one-cycle intercept gives unit Miner damage",
      np.isclose(miner_damage(np.array([1000 / 1.38]), np.array([1.0]), curve), 1.0))
check("ten identical half-range cycles follow analytical power law",
      np.isclose(miner_damage(np.array([500 / 1.38]), np.array([10.0]), curve), 10 / 2**10))
stress = equivalent_stress_MPa(0.25, 2e6, curve)
check("equivalent stress exactly reconstructs damage",
      np.isclose(miner_damage(np.array([stress]), np.array([2e6]), curve), 0.25))
converted = transfer_damage_by_strain(3.31e-4, 1e6, CARBON_UD, GLASS_UD)
check("carbon-to-glass transfer is positive and finite", np.isfinite(converted) and converted > 0)
f, a = stiffness_based_web_split(10.0, 2.0, 10.0, 1.0)
check("stiffness split conserves total load", np.isclose(f + a, 1.0))
check("taller equal-layup web takes two thirds", np.isclose(f, 2 / 3))

components, webs, summary = assess_phase7()
check("whole-blade component table is populated", len(components) >= 30)
check("both load-sharing methods are present", set(webs.split_method) == {"stiffness_based", "equal_50_50"})
check("each web pair conserves load under stiffness split",
      all(np.isclose(g.load_fraction.sum(), 1.0)
          for _, g in webs[(webs.split_method == "stiffness_based") &
                           (webs.design_variant == "day4_baseline")].groupby("r_m")))
check("Day 4 reinforcement zone has three DB layers",
      set(webs[(webs.design_variant == "day4_baseline") &
               (webs.r_m >= 10.25) & (webs.r_m <= 38.95)].face_layers_each_side) == {3})
check("outside reinforcement zone has two DB layers",
      set(webs[(webs.design_variant == "day4_baseline") &
               ((webs.r_m < 10.25) | (webs.r_m > 38.95))].face_layers_each_side) == {2})
check("all quantified damage values are finite and non-negative",
      np.isfinite(summary.damage_20y).all() and (summary.damage_20y >= 0).all())
check("published governing shell region remains near 10.25 m",
      abs(float(components.loc[components.damage_20y.idxmax(), "r_m"]) - 10.25) < 0.01)
check("core fatigue uses an explicitly sourced curve status",
      set(webs.core_fatigue_status) <= {"SOURCED_CURVE_PASS", "SOURCED_CURVE_FAIL"})
check("nominal quantified Phase 7 case passes D<=0.70",
      float(summary.damage_20y.max()) <= ROBUST_DAMAGE_TARGET)
check("Day 4 web baseline fails the conservative shear-fatigue proxy",
      float(webs[(webs.design_variant == "day4_baseline") &
                 (webs.split_method == "stiffness_based")].face_damage_20y_proxy.max()) > 1.0)
robust_web = webs[webs.design_variant == "fatigue_robust_final"]
check("robust web passes the nominal proxy target",
      float(robust_web.face_damage_20y_proxy.max()) <= ROBUST_DAMAGE_TARGET)

db_sweep, db_thresholds = db_proxy_sensitivity(webs)
nominal_threshold = db_thresholds[np.isclose(db_thresholds.b, 10.0)].iloc[0]
check("robust design has more than 20 percent C margin to D20=1 at b=10",
      float(nominal_threshold.allowable_C_degradation_percent) > 20.0)
check("b=9 passes the D20=0.70 target at the nominal C value",
      float(db_sweep[(db_sweep.C_change_percent == 0.0) &
                     (db_sweep.b == 9.0)].iloc[0].max_D20) <= ROBUST_DAMAGE_TARGET)
split = load_split_sensitivity(webs)
check("nominal split sensitivity reproduces governing D20",
      np.isclose(float(split[split.forward_split_error_percent == 0.0].iloc[0].max_D20),
                 float(robust_web.face_damage_20y_proxy.max())))
check("plus ten percent forward split error passes D20 target",
      float(split[split.forward_split_error_percent == 10.0].iloc[0].max_D20)
      <= ROBUST_DAMAGE_TARGET)
check("minus ten percent forward split error passes D20 target",
      float(split[split.forward_split_error_percent == -10.0].iloc[0].max_D20)
      <= ROBUST_DAMAGE_TARGET)
check("H100 sourced curve result is finite and below D20 target",
      np.isfinite(webs.core_damage_20y).all() and
      float(webs.core_damage_20y.max()) < ROBUST_DAMAGE_TARGET)

robust_cases = robust_case_summary(webs)
check("all six individual and combined robust cases pass",
      len(robust_cases) == 6 and set(robust_cases.status) == {"PASS"})
check("combined b=9 plus ten percent split is below D20 target",
      float(robust_cases[robust_cases.case_id == "b9_split_plus10"].iloc[0].max_face_D20)
      <= ROBUST_DAMAGE_TARGET)
check("combined b=9 minus ten percent split is below D20 target",
      float(robust_cases[robust_cases.case_id == "b9_split_minus10"].iloc[0].max_face_D20)
      <= ROBUST_DAMAGE_TARGET)
schedule = webs.attrs["robust_schedule"]
check("final scheme has one spanwise zone per web", len(schedule) == 2)
check("selected constant layups are forward 5 DB and aft 6 DB per side",
      int(schedule[schedule.web == "forward"].iloc[0].face_layers_each_side) == 5 and
      int(schedule[schedule.web == "aft"].iloc[0].face_layers_each_side) == 6)
check("selected H100 core thickness is 60 mm",
      np.allclose(schedule.core_thickness_m, 0.060))
candidates = webs.attrs["robust_candidates"]
selected = candidates[candidates.selected].iloc[0]
feasible = candidates[candidates.status == "PASS"]
check("final scheme respects the at-most-one-layer transition rule",
      int(selected.max_adjacent_layer_step) <= 1)
check("selected one-zone scheme has minimum estimated mass among feasible candidates",
      np.isclose(float(selected.estimated_web_material_mass_kg),
                 float(feasible.estimated_web_material_mass_kg.min())))
governing_pair = webs.attrs["robust_detail"]
governing_pair = governing_pair[(governing_pair.case_id == "nominal_b10_split0") &
                                np.isclose(governing_pair.r_m, 26.65)]
check("DB asymmetry balances governing nominal split within two percent",
      float((governing_pair.nominal_load_fraction - 0.5).abs().max()) < 0.02)
check("measured Saertex reference curve is stronger than robust envelope",
      SAERTEX_VU90079_MEASURED_RM1.C_MPa > DB_ROBUST_B9.C_MPa and
      SAERTEX_VU90079_MEASURED_RM1.b > DB_ROBUST_B9.b)
check("H100 curve uses NASA published B and gamma",
      np.isclose(H100_CORE_CURVE.C_MPa, 2.34) and np.isclose(H100_CORE_CURVE.b, 12.08))

print(f"\nPhase 7 verification: {passed}/{passed} PASS")
