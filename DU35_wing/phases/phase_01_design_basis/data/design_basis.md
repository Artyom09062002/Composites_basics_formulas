# Design Basis — Composite Wind Turbine Blade Study
**Version:** 1.0 (Phase 1 gate)
**Date:** 2026-07-20
**Status:** Research/screening level — not for certification

---

## 1. Turbine and Blade

| Parameter | Value | Source |
|-----------|-------|--------|
| Reference turbine | NREL 5MW Baseline | REF-001 |
| Rated power | 5.0 MW | REF-001 |
| Rotor diameter | 126 m | REF-001 |
| Blade length | 61.5 m | REF-001, REF-002 |
| Hub radius | 1.5 m | REF-001 |
| IEC wind class | **Class IB** (high wind, high turbulence) | REF-001, STD-001 |
| Hub height | 90 m | REF-001 |
| Rated wind speed V_rated | 11.4 m/s | REF-001 |
| Cut-in / cut-out | 3 m/s / 25 m/s | REF-001 |
| Rated rotor speed | 12.1 rpm | REF-001 |
| 1P frequency | 0.202 Hz | REF-001 |
| 3P frequency | 0.605 Hz | REF-001 |
| Reference blade model | Sandia SNL 61.5m (SAND2013-2569) | REF-002 |

---

## 2. Wind Class Parameters (IEC Class IB)

| Parameter | Symbol | Value | Source |
|-----------|--------|-------|--------|
| Reference wind speed (50-yr) | V_ref | 50 m/s | STD-001 Table 1 |
| Annual average wind speed | V_ave = 0.2 V_ref | 10 m/s | STD-001 Table 1 |
| 50-yr extreme wind (hub ht) | V_50 = 1.4 V_ref | 70 m/s | STD-001 §6.3.2 |
| 1-yr extreme wind | V_1 = 0.8 V_50 | 56 m/s | STD-001 §6.3.2 |
| Turbulence intensity (Class B) | I_ref | 0.14 | STD-001 Table 1 |
| Turbulence characteristic | Category B | — | STD-001 |

---

## 3. Design Objectives (this study)

| Objective | Target | Priority | Notes |
|-----------|--------|----------|-------|
| Minimum mass | ≤ Sandia reference (17,740 kg) | Primary | REF-002; hybrid carbon/glass expected −10–20% |
| Tip deflection | Δ_tip ≤ 5.5 m at U_op | Constraint | ≤ 8.9% of blade length; STD-002 §7.6 |
| Strength reserve factor | R ≥ 1.0 (all DLC, Tsai-Wu/Hashin with γ_d) | Constraint | STD-002 §7.5, STD-003 §7.4 |
| Buckling reserve | RF_buck ≥ 1.0 (with γ_f=1.35) | Constraint | STD-003 §7.5 |
| Fatigue life | Miner's sum D ≤ 1.0 over 20 years | Constraint | STD-001 §7.4; S-N from MAT-003 |
| 1P/3P margin | |f_blade − f_1P/3P| ≥ 10% | Constraint | STD-001 §9; avoids resonance |
| Manufacturability | Infusion (VARIM) compatible | Preference | Company intent |

---

## 4. Structural Architecture

Per Sandia SAND2013-2569 reference model (REF-002):

| Zone | Structural element | Material set |
|------|-------------------|--------------|
| Root (r/R = 0–0.05) | Thick glass laminate, cylindrical | SNL Triax + UD glass |
| Transition (0.05–0.2) | Tapering panels + shear webs | SNL Triax + DB glass |
| Spar cap (0.1–0.85) | Main load-bearing UD strips | E-LT-5500 UD glass (baseline); UD carbon (hybrid option) |
| Shear webs (2×) | Sandwich core + glass skins | Saertex DB + foam core |
| Shell skins | Biax glass + foam sandwich | Saertex DB + foam |
| Tip (0.85–1.0) | Light glass panels | SNL Triax |

---

## 5. Materials (Baseline — Sandia set)

Full properties are registered in `DU35_wing/SOURCES.md` (MAT-001 to MAT-004)
and stored in `phase_03_materials/data/`.

| Material | Role | E1 [GPa] | ρ [kg/m³] | Source |
|----------|------|----------|-----------|--------|
| E-LT-5500 UD glass (0°) | Spar cap | ~41 | ~2000 | MAT-002 |
| Saertex ±45° biax (DB) | Shear web / shell | ~13 | ~1830 | MAT-002 |
| SNL Triax (0/±45) | Root + transitions | ~27 | ~1850 | MAT-002 |
| Newport 307 UD carbon | Hybrid spar cap | ~121 | ~1540 | MAT-004 |
| PVC/PET foam (core) | Sandwich panels | ~0.06–0.26 | ~60–200 | MAT-002 |

*Detailed elastic + strength constants and S-N curves: Phase 3.*

---

## 6. Design Load Cases (Screening subset)

The full list is in `DU35_wing/SOURCES.md`. For Phases 2–4 the following
screening envelopes are used:

| Case | Label | Description | V [m/s] | γ_f | Type |
|------|-------|-------------|---------|-----|------|
| DLC 1.3 env | U_op | Extreme turbulence at rated; flap-wise bending dominant | 11.4 | 1.35 | Ultimate |
| DLC 6.1 env | U_park_feather | Parked, 50-yr extreme wind; safely feathered at 90° | 70 | 1.35 | Ultimate |
| Screening fault envelope | U_park_fault | Parked extreme wind with unavailable/failed safe pitch; 360° pitch sweep | 70 | 1.35* | Ultimate screening |
| DLC 1.1 | Fatigue | NTM, full life spectrum | — | 1.0 | Fatigue |

Verified Phase 2 root targets before material factor: U_op = 7.53 MN·m,
U_park_feather = 0.80 MN·m, U_park_fault = 23.11 MN·m, and gravity edgewise
= 3.07 MN·m. With the provisional load factor, the Phase 4 design-load targets
are 10.2, 1.1, 31.2, and 4.1 MN·m respectively. *The U_park_fault DLC mapping
and partial factor require confirmation in the later IEC/OpenFAST validation;
it is retained now as a conservative screening envelope.*

---

## 7. Safety Factors (Screening values)

| Factor | Symbol | Value | Source |
|--------|--------|-------|--------|
| Load partial factor (ultimate) | γ_f | 1.35 | STD-001 Table 3 |
| Material factor — UD glass tension/compression | γ_m | 2.0 | STD-002 Table 3 (screening) |
| Material factor — shear | γ_m | 2.5 | STD-002 Table 3 (screening) |
| **Combined design factor (strength)** | **γ_d = γ_f × γ_m** | **2.70 / 3.375** | Derived |
| Stiffness factor (deflection) | γ_f,stiff | 1.1 | STD-003 §7.4 |

Required reserve factor: **R = σ_strength / (γ_d × σ_applied) ≥ 1.0**
Equivalently (CLT output): **Tsai-Wu R ≥ γ_f × γ_m = 2.70** at applied design load without partial factors.

---

## 8. Validation Reference

All simulation outputs in Phases 2–5 will be compared to the Sandia SAND2013-2569 reference blade:

| Quantity | Reference value | Source |
|----------|----------------|--------|
| Blade mass | 17,740 kg | REF-002 §2.2 |
| Flap-wise stiffness at root | ≈ 1.5×10¹⁰ N·m² | REF-002 §3 |
| 1st flapwise eigenfrequency | ~0.88 Hz | REF-002 §4 |
| 1st edgewise eigenfrequency | ~1.07 Hz | REF-002 §4 |
| Max tip deflection (rated, steady) | ~5.5–6.0 m | REF-001/REF-002 derived |

*Tolerance: ±10% on distributed quantities (stiffness, loads); ±5% on global metrics (mass, eigenfrequencies).*

---

## 9. Scope and Limitations

This study is **research/screening level**:
- Materials: literature values from open databases (MAT-001–004). Not company-specific coupon data.
- Geometry: Sandia 61.5m reference model; manufacturing geometry not included.
- Analysis: CLT (Phases 0–4), parametric shell FEA (Phase 5), surrogate (Phase 6).
- Standards: IEC 61400-1/5 and DNV-ST-0376 for factor selection; full code compliance not assessed.
- Certification requires: validated commercial FEA (Abaqus/ANSYS), physical testing, licensed engineer sign-off.

**All numbers in this document require Artem's review before being reported to company.**
