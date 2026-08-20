# Composite Mechanics and NREL 5 MW Blade Screening

This repository contains two connected parts:

1. `composite_physics/` — reusable classical lamination, reduced-stiffness,
   plate-buckling, load-mapping and sandwich-panel formulas;
2. `DU35_wing/` — a phase-organised, full-scale NREL 5 MW blade example based
   on public NREL/Sandia data and CAD-measured structural geometry.

The project is an engineering screening and educational model. It is not an
IEC certification package or a substitute for a validated shell/solid FE
model and physical blade tests.

## Final screening result

The Phase 8 configuration uses a 61.5 m blade, glass spar caps, two sandwich
shear webs, 5/6 DB layers per side on the forward/aft webs and a 60 mm H100
core. The integrated screening checks report:

- minimum static reserve factor: `RF = 2.055`;
- parked-fault tip deflection: `4.816 m` against a `5.5 m` limit;
- first flap frequency: `1.422 Hz` against a `0.727 Hz` target;
- worst combined fatigue sensitivity: `D20 = 0.693` against `0.70`;
- estimated blade mass: `19,614 kg`, 10.6% above the 17,740 kg reference;
- no open numerical FAIL in the published final acceptance table.

The result is reproducible at beam/shell-equivalent research-screening level.

## Repository layout

```text
composite_physics/       reusable composite mechanics
tests/core/              reusable-physics tests
DU35_wing/phases/        blade code, data, tests and results by phase
scripts/run_all_tests.py one-command public verification
```

Each blade phase owns its calculation code, source data, tests and selected
results. Earlier DU35 RBS, sandwich-panel and full-blade-web studies are kept
under `phase_04_structural_design/studies/` because they support the accepted
structural architecture rather than representing separate formal phases.

## Run the public verification

Requirements: Python 3.12 or newer.

```bash
python -m pip install -r requirements.txt
python scripts/run_all_tests.py
```

The tests demonstrate mathematical and computational reproducibility. They do
not constitute experimental validation.

## Final CAD model

The final FreeCAD v6 file is stored in the project Google Drive. Its public
manifest and SHA-256 are in
`DU35_wing/phases/phase_08_final/cad/cad_manifest.json`; a compact rotation GIF
is tracked for repository preview.

[Open the project Google Drive](https://drive.google.com/drive/folders/1568kY9HVPt-y1rDeaz4n72Q8LZ079cSK)

## Primary references

- Jones, *Mechanics of Composite Materials*, 2nd ed.
- Jonkman et al., NREL/TP-500-38060.
- Griffith and Ashwill, SAND2011-3779.
- Sandia National Laboratories, SAND2013-2569.
- Bednarcyk et al., NASA/TM-2012-217694.

Detailed traceability is provided in `DU35_wing/SOURCES.md` and the phase data
folders.
