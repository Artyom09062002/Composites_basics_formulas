# DU35 / NREL 5 MW Blade — Source Registry

All material properties, loads, geometry, and standard requirements used in this project must be logged here.
Format: `[ID] Author/Org, Title, Year, URL/DOI/Report-number, relevant section/table.`

Unsourced numbers are marked **TBD** in code and must not silently drive design decisions (rule 4, CLAUDE.md).

---

## Standards

| ID | Reference | Relevance |
|----|-----------|-----------|
| STD-001 | IEC 61400-1 Ed.4, "Wind energy generation systems – Part 1: Design requirements", IEC, 2019 | Wind classes, DLC list, load partial factors |
| STD-002 | IEC 61400-5, "Wind energy generation systems – Part 5: Wind turbine rotor blades", IEC, 2020 | Blade-specific design requirements, safety factors, test requirements |
| STD-003 | DNV-ST-0376, "Rotor Blades for Wind Turbines", DNV, 2022 (rev. 2023) | Structural safety factors, laminate design, fatigue |
| STD-004 | GL Guidelines for the Certification of Wind Turbines, Ed. 2010 (Germanischer Lloyd) | Legacy safety factors; superseded by IEC 61400-1 Ed.4 / DNV-ST-0376 but still referenced |

---

## Turbine / Blade Reference

| ID | Reference | Relevance |
|----|-----------|-----------|
| REF-001 | Jonkman J., Butterfield S., Musial W., Scott G., "Definition of a 5-MW Reference Wind Turbine for Offshore System Development", NREL/TP-500-38060, NREL, 2009. https://www.nrel.gov/docs/fy09osti/38060.pdf | Rotor diameter 126 m, blade length 61.5 m, rated power, IEC Class IB, hub height, rated RPM, cut-in/out speeds |
| REF-002 | Resor B.R., "Definition of a 5MW/61.5m Wind Turbine Blade Reference Model", SAND2013-2569, Sandia National Laboratories, 2013. https://energy.sandia.gov/wp-content/gallery/uploads/SAND2013-2569.pdf | Blade geometry (chord, twist, airfoils by span), structural layup zones, material set, created explicitly for optimization and materials studies |
| REF-003 | Jonkman J.M., Buhl M.L. Jr., "FAST User's Guide", NREL/EL-500-38230, NREL, 2005. | Aerodynamic and structural simulation methodology; AeroDyn blade data tables |

---

## Materials (SNL/MSU/DOE composite database)

| ID | Reference | Material | Relevance |
|----|-----------|----------|-----------|
| MAT-001 | Sandia/MSU/DOE Composite Material Fatigue Database, Montana State University, continuously updated. https://www.montana.edu/composites/ | General | Master reference for glass and carbon composite fatigue data |
| MAT-002 | SAND2013-2569 (=REF-002), Table 4 and Appendix A | E-LT-5500 UD glass, Saertex biax/DB, SNL Triax, core foam | Elastic constants, strength, density as used in Sandia 5MW blade model |
| MAT-003 | Mandell J.F., Samborsky D.D., "DOE/MSU Composite Material Fatigue Database", latest release | E-LT-5500, other glass UD | S-N curves, static properties, multi-source validation |
| MAT-004 | SAND2013-2569, Section 3.3 | Newport 307-based UD carbon (spar cap hybrid option) | Elastic + strength properties for carbon UD |

---

## Aerodynamic Profiles

| ID | Reference | Relevance |
|----|-----------|-----------|
| AERO-001 | Jonkman et al. NREL/TP-500-38060 (=REF-001), Table 6-1 | Airfoil distribution along span: Cylinder1, Cylinder2, DU40, DU35, DU30, DU25, DU21, NACA64-618 |
| AERO-002 | Lindenburg C., "Aeroelastic Modelling of the LMH64-5 Blade", ECN-C-02-016, 2002 | DU-series airfoil Cl/Cd polars at Re relevant for 5MW blade |
| AERO-003 | Rfoil / XFOIL computed polars from NREL 5MW AeroDyn input deck (public repository: https://github.com/NREL/openfast) | Cl, Cd, Cm tables at 0°–360° for all airfoils |

---

## Load Cases (DLC) — Phase 1 definition

Sources: STD-001 (IEC 61400-1 Ed.4, Table 1), STD-002 (IEC 61400-5, Section 7), STD-003 (DNV-ST-0376, Section 5).

Selected subset for screening analysis (Phases 2–4):

| DLC | Description | Wind model | Analysis type | Ref standard |
|-----|-------------|------------|---------------|--------------|
| DLC 1.1 | Normal power production, statistical loads | NTM (Normal Turbulence Model) | Fatigue (F) | IEC 61400-1 §7.4 |
| DLC 1.3 | Normal power production, extreme turbulence | ETM (Extreme Turbulence Model) | Ultimate (U) | IEC 61400-1 §7.4 |
| DLC 1.4 | Normal power production, extreme coherent gust with direction change | ECD | Ultimate (U) | IEC 61400-1 §7.4 |
| DLC 2.1 | Fault during power production | NTM | Ultimate (U) | IEC 61400-1 §7.5 |
| DLC 6.1 | Parked (idling), extreme wind | EWM 50-year | Ultimate (U) | IEC 61400-1 §7.9 |
| DLC 6.4 | Parked, normal turbulence | NTM | Fatigue (F) | IEC 61400-1 §7.9 |

**For analytical screening (Phases 2–4) two envelope cases are used:**
- **U_op**: rated operation at V_rated = 11.4 m/s, extreme turbulence (≈DLC 1.3 envelope)
- **U_park**: parked extreme wind 50-year, V_50 = 70 m/s for IEC Class IB (≈DLC 6.1)

---

## Safety Factors

Sources: STD-002 (IEC 61400-5, §7.5, Table 3), STD-003 (DNV-ST-0376, §7.4).

### Partial load factors γ_f (IEC 61400-1 Ed.4, Table 3)
| Load case type | γ_f (unfavourable) |
|---------------|-------------------|
| Ultimate (normal/abnormal) | 1.35 |
| Fatigue | 1.0 (damage sum uses γ_f=1) |

### Material / consequence factor γ_m (IEC 61400-5 / DNV-ST-0376)
| Property | γ_m (IEC 61400-5 Table 3) | γ_m (DNV-ST-0376) | Note |
|----------|--------------------------|-------------------|------|
| Tensile / compressive strength (UD glass) | 1.65–2.0 | 1.7–2.1 | Depends on laminate quality level; screening uses 2.0 |
| Shear strength | 2.0–2.5 | 2.0 | Screening uses 2.5 |
| Stiffness (tip deflection) | 1.1 | 1.1 | |

**Combined design factor for strength:** γ_d = γ_f × γ_m ≥ 2.73 (=1.35 × 2.0) for ultimate static.

### Deflection limit
Tip clearance to tower: Δ_tip ≤ L_blade × 0.10 = 6.15 m (screening rule; actual tower-clearance geometry sets tighter limit — typically ~5% for onshore upwind, but NREL 5MW is a reference machine with relaxed limits). Conservative screening: **Δ_tip ≤ 5.5 m** (≈8.9% of 61.5 m).

---

## Phase 3 — Additional Material Sources

| ID | Reference | Material | Relevance |
|----|-----------|----------|-----------|
| MAT-005 | Hexcel HexTow IM7/8552 datasheet, Hexcel Corp., 2016; and Soutis C., "Fibre Reinforced Composites in Aircraft Construction", Progress in Aerospace Sciences 41 (2005) 143–151 | UD carbon generic | Secondary check for Newport 307 carbon UD properties (E1, Xt, rho) |
| MAT-006 | DIAB Group, "Divinycell H Technical Data", 2020. https://www.diabgroup.com/products/h-grade | PVC foam H100 | Manufacturer datasheet: E, shear modulus, tensile/compressive strength, density |
| MAT-007 | Mandell J.F., Samborsky D.D., Cairns D.S., "Fatigue of Composite Materials and Substructures for Wind Turbine Blades", SAND2002-0771, Sandia NL, 2002. https://www.osti.gov/biblio/801442 | E-glass UD, biax | Extensive fatigue S-N data; static properties confirm E-LT-5500 values |
| MAT-008 | Berggreen C., Branner K., Jensen J.F., Schultz J.P., "Application and Analysis of Sandwich Elements in the Primary Structure of Large Wind Turbine Blades", J. Sandwich Structures & Materials 9 (2007) 525–552. DOI:10.1177/1099636207065529 | Sandwich core | Secondary properties for PVC foam cores in wind blade context |

---

## Revision log

| Date | Change |
|------|--------|
| 2026-07-20 | Initial population — Phase 1 (design basis) |
| 2026-07-21 | Phase 3: Added MAT-005 through MAT-008 for materials database cross-check |
| 2026-08-17 | Phase 7: Added published fatigue damage and normalized moment-range spectrum sources |
| 2026-08-19 | Phase 7 closure: Added H100 shear-fatigue evidence and DB-slope sensitivity basis |
| 2026-08-19 | Phase 7 robust closure: Added FAT-006 through FAT-008 and replaced the H100 estimate with the NASA measured curve |

## Phase 7 — Fatigue

| ID | Reference | Relevance |
|----|-----------|-----------|
| FAT-001 | Resor B.R., *Definition of a 5MW/61.5m Wind Turbine Blade Reference Model*, SAND2013-2569, 2013. https://digital.library.unt.edu/ark:/67531/metadc836285/m2/1/high_res_d/1095962.pdf | Tables 18–19: 20-year fatigue setup and stations; Tables 24–26: two-parameter S–N constants and published Miner damage. |
| FAT-002 | NREL/TP-5000-65227, *Wind Turbine Blade Fatigue Test Development*, 2016. https://docs.nrel.gov/docs/fy16osti/65227.pdf | Published normalized flapwise and lead-lag fatigue moment-range distributions for a representative 60 m blade; used as the spanwise reference spectrum shape. |
| FAT-003 | Sandia/MSU/DOE Composite Material Fatigue Database, Montana State University. https://fdn.montana.edu/composites/index.html | Public fatigue-data repository and future source for replacement of the DB-face screening proxy. |
| FAT-004 | Burman M., Zenkert D., *Fatigue of foam core sandwich beams—1: undamaged specimens*, International Journal of Fatigue 19(7), 1997. https://doi.org/10.1016/S0142-1123(97)00069-8 | H100 PVC foam tested in core-shear-dominated sandwich fatigue. Used only to define an explicitly estimated lower-bound screen: Basquin exponent about 12 and approximately 10^4 cycles at the core shear-yield level. |
| FAT-005 | Burman M., *Fatigue Crack Initiation and Propagation in Sandwich Structures*, supporting full text. https://www.diva-portal.org/smash/get/diva2%3A472688/FULLTEXT01.pdf | Accessible supporting evidence for the H100 high-cycle slope and approximate cycle count at shear yield. Its DBLT-850 slope near 7.8 supports including b=8 in sensitivity, but is not substituted as a Saertex DB shear curve. |
| FAT-006 | Mandell J.F., Miller D.A., Samborsky D.D., *Creep/Fatigue Behavior of Resin Infused Biaxial Glass Fabric Laminates*, AIAA 2013-1630. https://doi.org/10.2514/6.2013-1630 | Table 1 identifies Saertex VU 90079; the Figure 7 R=-1 Fabric L trend was digitized as `S=156 N^-0.103 MPa`; Table 2 supplies the ASTM D3518 axial-to-ply-shear relation. Digitized transformed reference: `C_shear=78 MPa`, `b=9.709` (not a stated numeric fit in the source). |
| FAT-007 | Bednarcyk B.A., Yarrington P.W., Arnold S.M., *Multiscale Fatigue Life Prediction for Composite Panels*, NASA/TM-2012-217694, 2012. https://ntrs.nasa.gov/citations/20120015395 | Section 4, Table 3 and Figure 6: measured H100 shear fatigue at R=0.1, C=2.34 MPa, b=12.08; G=40 MPa; shear yield=1.13 MPa. |
| FAT-008 | Project robust-envelope derivation documented in `phases/phase_07_fatigue/data/db_h100_sources.md` | Traceable estimate: FAT-006 reference reduced from C=78 MPa, b=9.709 to C=62 MPa, b=9; crossed with A66h split errors -10/0/+10% and D20 target 0.70. |
