# Structural Layup Zones — Sandia 61.5m Blade (SAND2013-2569)

Source: Resor, SAND2013-2569, Sandia NL, 2013, Sections 3–4.

## Zone map (r/R)

```
Root       |  Transition  |  Spar region (aerodynamic blade)  |  Tip
0          0.02          0.10                                  0.90   1.0
|--Cyl----|-transition--|---------aero sections---------------|-tip--|
```

## Layup schedule by region

### Root (r/R = 0–0.02, ~0–1.2 m)
- Shell: [0/±45/0]_s SNL Triax, t ≈ 80–100 mm (thick plug)
- Purpose: load introduction to hub bolts
- Ply angles: [0°, +45°, -45°, 0°] symmetric

### Root transition (r/R = 0.02–0.10, 1.2–6.2 m)
- Shell outer: SNL Triax [0/±45] dropping plies
- Shell inner: Saertex DB (±45°) fabric
- Core: 50mm foam between inner/outer skins
- Spar cap starts: E-LT-5500 UD [0°] building up from ~4 plies

### Spar cap (r/R = 0.10–0.85, 6.2–52.3 m)
- Two caps (suction + pressure side), symmetric about chord plane
- Material: E-LT-5500 UD glass [0°] (baseline); carbon UD (hybrid option)
- Layup: [0°]_n symmetric, n varies with span (max at ~r/R=0.15–0.20)
- Max cap width: ~0.35–0.40 m; thickness peaks ~60–70mm
- **Key design variable in optimization**

### Shear webs (2×, at ~15% and ~50% chord)
- Web skins: Saertex DB (±45°), 2–4 plies each side
- Web core: PVC foam 50–100 mm
- Runs r/R ≈ 0.05–0.90

### Shell skins (suction + pressure, outside spar cap)
- Outer laminate: 2 × Saertex DB (±45°) + SNL Triax (leading/trailing edge reinforcement)
- Core: PVC/PET foam 20–50 mm
- Leading edge: additional UD band for buckling resistance

### Tip (r/R = 0.85–1.0)
- Thin shell only: SNL Triax [0/±45], t ≈ 4–6 mm
- No spar cap; tapers to ~0 at tip

## Ply properties (baseline)

| Material | t_ply [mm] | E1 [GPa] | ρ [kg/m³] |
|----------|-----------|----------|-----------|
| E-LT-5500 UD (0°) | 0.91 | 41.8 | 1920 |
| Saertex DB (±45°) | 0.47 | 13.6 (Ex) | 1830 |
| SNL Triax (0/±45) | 0.94 | 27.7 (Ex) | 1850 |
| PVC foam H100 | — | 0.13 (E3) | 100 |

Source: SAND2013-2569, Table 4 and Appendix A (= MAT-002 in sources.md).

## Reference blade mass budget

| Region | Mass [kg] | % total |
|--------|-----------|---------|
| Root (0–6m) | ~2,500 | ~14% |
| Spar caps | ~7,000 | ~40% |
| Shear webs | ~1,800 | ~10% |
| Shell skins | ~5,000 | ~28% |
| Tip + misc | ~1,440 | ~8% |
| **Total** | **~17,740** | 100% |

Source: SAND2013-2569, Section 2.2 (REF-002). Individual zone breakdown is approximate.

## Phase 3/4 action
Load ply properties from the Phase 3 material database.
Implement the layup schedule in `phase_04_structural_design/code/`.
Run CLT at each station with DLC loads → check R_TW ≥ 2.70.
