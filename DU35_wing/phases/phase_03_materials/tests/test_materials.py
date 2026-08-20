"""
verify/verify_materials.py — Phase 3 gate verification.

Checks:
  1. All 5 materials load without error.
  2. Physical range checks for every property.
  3. Cross-source divergence < 15% for key properties.
  4. CLT Q-matrix is positive-definite (sanity: E1*E2*(1-nu12*nu21) > 0).
  5. Specific stiffness ranking: Carbon > Glass (expected physics).

Run from repo root:
    python verify/verify_materials.py
"""

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "code")))

import numpy as np
from materials_db import load_materials, get_ply, specific_props

# ── Cross-source reference values (secondary sources) ──────────────────────
# Format: {material_name: {property: (secondary_value, tolerance_frac, "source")}}
# All properties in the same units as materials_db.csv (GPa or MPa).
#
# Secondary sources:
#  ELT5500_UD   → MAT-003 (DOE/MSU Mandell & Samborsky 2010)
#  Saertex_DB   → MAT-001 (DOE/MSU fatigue database effective properties)
#  SNL_Triax    → MAT-001 (cross-check from biax+UD mix-law)
#  Carbon UD    → MAT-005 (Hexcel IM7/8552 + Soutis 2005)
#  PVC Foam H100→ MAT-006 (DIAB Divinycell H100 datasheet)
#
SECONDARY = {
    "ELT5500_UD": {
        "E1_GPa" : (41.4, 0.15, "MAT-003 / Mandell 2010"),
        "Xt_MPa" : (980.0, 0.15, "MAT-003 / Mandell 2010"),
        "Xc_MPa" : (554.0, 0.15, "MAT-003 / Mandell 2010"),
        "rho_kgm3": (1920.0, 0.05, "MAT-003"),
    },
    "Saertex_DB": {
        "E1_GPa" : (9.65, 0.15, "MAT-001 biax effective"),
        "G12_GPa": (12.5, 0.15, "MAT-001 biax shear"),
    },
    "SNL_Triax": {
        "E1_GPa" : (26.0, 0.15, "MAT-001 mix-law estimate"),
    },
    "Newport307_CarbonUD": {
        "E1_GPa" : (110.0, 0.15, "MAT-005 / Soutis 2005 T300-class"),
        "Xt_MPa" : (1400.0, 0.15, "MAT-005 / Hexcel IM7 ref"),
        "rho_kgm3": (1560.0, 0.10, "MAT-005 / Hexcel datasheet"),
    },
    "PVC_Foam_H100": {
        "E1_GPa" : (0.125, 0.15, "MAT-006 / DIAB H100 datasheet"),
        "Xt_MPa" : (3.10,  0.15, "MAT-006 / DIAB H100 datasheet"),
        "rho_kgm3": (100.0, 0.05, "MAT-006 / DIAB H100 datasheet"),
    },
}

# ── Physical range bounds ────────────────────────────────────────────────────
BOUNDS = {
    # (property, min, max, unit)
    "ELT5500_UD"      : [("E1_GPa", 30, 55, "GPa"),  ("Xt_MPa", 700, 1300, "MPa"),
                          ("rho_kgm3", 1800, 2100, "kg/m³"), ("nu12", 0.15, 0.45, "-")],
    "Saertex_DB"      : [("E1_GPa", 5, 15, "GPa"), ("G12_GPa", 8, 20, "GPa"),
                          ("rho_kgm3", 1600, 2000, "kg/m³")],
    "SNL_Triax"       : [("E1_GPa", 18, 40, "GPa"), ("rho_kgm3", 1700, 2100, "kg/m³")],
    "Newport307_CarbonUD": [("E1_GPa", 90, 160, "GPa"), ("Xt_MPa", 1000, 2000, "MPa"),
                             ("rho_kgm3", 1400, 1700, "kg/m³")],
    "PVC_Foam_H100"   : [("E1_GPa", 0.05, 0.30, "GPa"), ("rho_kgm3", 80, 130, "kg/m³")],
}

# ── Helpers ─────────────────────────────────────────────────────────────────
_PASS = 0
_FAIL = 0

def check(label, cond, note=""):
    global _PASS, _FAIL
    status = "PASS" if cond else "FAIL"
    if cond:
        _PASS += 1
        print(f"  PASS  {label}")
    else:
        _FAIL += 1
        print(f"  FAIL  {label}  ({note})")


def main():
    global _PASS, _FAIL
    print("=" * 65)
    print("  MATERIALS DB VERIFICATION — Phase 3")
    print("=" * 65)

    # ── 1. Load database ────────────────────────────────────────────
    try:
        db = load_materials()
        check("Database loads without error", True)
    except Exception as e:
        check("Database loads without error", False, str(e))
        print("Cannot continue — aborting."); sys.exit(1)

    expected_mats = ["ELT5500_UD", "Saertex_DB", "SNL_Triax",
                     "Newport307_CarbonUD", "PVC_Foam_H100"]
    for m in expected_mats:
        check(f"Material '{m}' present", m in db.index)

    # ── 2. Physical range checks ─────────────────────────────────────
    print("\n  — Physical range checks —")
    for mat, bounds in BOUNDS.items():
        if mat not in db.index:
            continue
        for (prop, lo, hi, unit) in bounds:
            val = db.loc[mat, prop]
            check(f"  {mat} {prop} in [{lo}, {hi}] {unit}  (got {val:.4g})",
                  lo <= val <= hi)

    # ── 3. Cross-source divergence < 15% ─────────────────────────────
    print("\n  — Cross-source divergence checks (tol ≤ 15%) —")
    for mat, props in SECONDARY.items():
        if mat not in db.index:
            continue
        for prop, (ref_val, tol, src) in props.items():
            primary = db.loc[mat, prop]
            div = abs(primary - ref_val) / ref_val
            check(f"  {mat} {prop}: primary={primary:.4g}, "
                  f"secondary={ref_val:.4g} (err={100*div:.1f}%, tol={100*tol:.0f}%)",
                  div <= tol,
                  f"secondary source: {src}")

    # ── 4. Q-matrix positive definite ───────────────────────────────
    print("\n  — Q-matrix definiteness (per ply) —")
    for mat in expected_mats:
        if mat not in db.index:
            continue
        ply = get_ply(mat, db)
        E1, E2, nu12 = ply["E1"], ply["E2"], ply["nu12"]
        nu21 = nu12 * E2 / E1
        denom = 1.0 - nu12 * nu21
        check(f"  {mat} Q11 denominator > 0  (1-nu12*nu21={denom:.4f})", denom > 0)
        Q11 = E1 / denom
        Q22 = E2 / denom
        Q12 = nu12 * E2 / denom
        check(f"  {mat} Q11*Q22 > Q12^2  (positive definite)",
              Q11 * Q22 > Q12 ** 2)

    # ── 5. Specific stiffness ranking ───────────────────────────────
    print("\n  — Specific stiffness ranking —")
    sp = specific_props(db)
    e_carbon = sp.loc["Newport307_CarbonUD", "E1_specific_MNm_kg"]
    e_glass  = sp.loc["ELT5500_UD", "E1_specific_MNm_kg"]
    check(f"Carbon E1/rho ({e_carbon:.1f}) > Glass E1/rho ({e_glass:.1f}) [MN·m/kg]",
          e_carbon > e_glass)
    check(f"Carbon E1/rho ratio ≥ 2.5× glass  (got {e_carbon/e_glass:.2f}×)",
          e_carbon / e_glass >= 2.5)

    xt_carbon = sp.loc["Newport307_CarbonUD", "Xt_specific_kNm_kg"]
    xt_glass  = sp.loc["ELT5500_UD", "Xt_specific_kNm_kg"]
    check(f"Carbon Xt/rho ({xt_carbon:.1f}) > Glass Xt/rho ({xt_glass:.1f}) [kN·m/kg]",
          xt_carbon > xt_glass)

    # ── Summary ──────────────────────────────────────────────────────
    total = _PASS + _FAIL
    print("\n" + "=" * 65)
    print(f"  Result: {_PASS}/{total} PASS")
    print("=" * 65)
    if _FAIL > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
