# Composite Mechanics — Core Formulas and Wind-Blade Example

A compact Python implementation of classical lamination theory and composite-panel screening methods.

The repository is divided into:

- reusable composite-mechanics modules in the repository root;
- `DU35_wing/` — a worked wind-turbine blade example based on measured CAD geometry and source-based NREL 5 MW data.

The calculations are intended for engineering screening and education. They are not certification analyses and have not been validated experimentally.

---

## Core calculation chain

```text
lamina_mechanics.py     engineering constants       → Q
transformations.py      Q and ply angle             → Qbar
laminate.py             through-thickness integration → A, B, D
reduced_stiffness.py    A, B, D                     → D*
buckling_rbs.py         D*, panel dimensions        → critical plate load
panel_loads.py          stress or bending moment    → applied Nx
sandwich_panel.py       faces/core geometry         → sandwich limit states
```

The reduced bending stiffness used for extension–bending coupled laminates is:

```text
D* = D - B A^-1 B
```

Each stage is kept separate so that assumptions, material properties and calculation methods can be checked or replaced independently.

---

## Repository structure

```text
repository/
├── Tests/                         core physics tests
├── buckling_rbs.py
├── lamina_mechanics.py
├── laminate.py
├── panel_loads.py
├── reduced_stiffness.py
├── sandwich_panel.py
├── transformations.py
│
└── DU35_wing/
    ├── Day 2: spar-cap RBS screening
    ├── Day 3: DU35 sandwich-panel screening
    └── day04_full_blade/
        ├── README.md
        ├── code/
        ├── data/
        ├── results/
        └── SOURCES.md
```

Files specific to the blade geometry, NREL/Sandia inputs and individual work sessions belong inside `DU35_wing`. Only generally reusable mechanics remain in the repository root.

---

## Worked example

### Day 2 — DU35 upper spar cap

A measured DU35 upper spar-cap panel was evaluated using reduced bending stiffness because its equivalent stack has non-zero extension–bending coupling.

Within the support-spacing sweep:

- applied compression resultant: `7.085 MN/m`;
- minimum critical load: `27.49 MN/m`;
- reserve factor: `RF = 3.880`.

This is a simply supported equivalent-panel screening result, not a manufacturing laminate definition.

### Day 3 — DU35 upper sandwich panels

The leading-edge and trailing-edge upper panels were assessed for:

1. global buckling with core-shear correction;
2. face wrinkling;
3. core shear crimping;
4. face compression.

The trailing-edge foam panel governed:

- severe `109 MPa` comparison case: `RF = 0.848`;
- beam-mapped case: `RF = 2.998`.

The `109 MPa` case is retained for comparison but is not the source-based shell load at the current blade station.

### Day 4 — full-scale NREL 5 MW blade

The DU35 section was connected to a full-scale 61.5 m blade screening model with two source-based sandwich shear webs.

The earlier web result `RF = 0.0058` came from a superseded monolithic-wall model without the documented 50 mm core and separated DB faces. It is therefore not used as a design-decision basis.

For the Sandia reference web construction:

```text
2 DB layers / 50 mm foam / 2 DB layers
```

the checks gave:

- global shear buckling: `RF = 1.285` — PASS;
- governing face shear: `RF = 0.784` — FAIL at the aft web, `r = 26.65 m`.

The minimum tested reinforcement was then applied:

```text
3 DB layers per side on both webs from 10.25 to 38.95 m
```

The reinforced design gave:

- minimum reserve factor: `RF = 1.078` — PASS;
- new governing point: aft web at approximately `r = 43.05 m`, outside the reinforced zone;
- corrected source-based DU35 trailing-edge result: `RF = 2.227` — PASS.

“Minimum tested reinforcement” means the smallest change evaluated in this study. It is not a formal mass-optimized solution.

---

## Day 4 files

Only files required to understand or reproduce the numerical assessment are published:

```text
DU35_wing/day04_full_blade/
├── README.md
├── code/
│   ├── nrel5mw_blade_integration.py
│   └── nrel5mw_sandwich_web_assessment.py
├── data/
│   └── cad_web_geometry.csv
├── results/
│   └── Day4_results.xlsx
└── SOURCES.md
```

The final FreeCAD model is stored separately in the project Google Drive:

```text
NREL5MW_61p5m_CADBasedBlade_v3_ReinforcedWebs.FCStd
```

[Open the project Google Drive](https://drive.google.com/drive/folders/1568kY9HVPt-y1rDeaz4n72Q8LZ079cSK)

Generated screenshots, FreeCAD backup files and internal CAD-building utilities are not tracked because they can be recreated from the final model.

---

## Assumptions and limitations

The Day 4 shear-web calculation uses:

- factored parked-fault shear derived as `V = |dM/dr|`;
- equal load sharing between the two webs;
- a support-spacing sweep of `a/h = 0.5–10`;
- `RF ≥ 1.000` as the screening PASS threshold;
- a conservative H100 core-strength proxy where the Sandia source does not report foam strength.

The model does not yet include:

- measured load sharing between the two webs;
- confirmed transverse-support spacing;
- validated shell FEA;
- adhesive and debonding failure;
- fatigue and manufacturing tolerances;
- certification safety factors.

The next stage is a mass-and-stiffness update based on the reinforced CAD model, followed by shell FEA using the real support arrangement.

---

## Running the core tests

Requirements:

- Python 3.12 or newer;
- NumPy.

```bash
python -m unittest discover -s Tests -v
```

The tests check the implemented mathematics and stored reference behaviour. Passing tests demonstrate reproducibility of the code, not experimental validation of the blade.

---

## References

- R. M. Jones, *Mechanics of Composite Materials*, 2nd ed.
- Hexcel, *Honeycomb Sandwich Design Technology*.
- J. Jonkman et al., *Definition of a 5-MW Reference Wind Turbine for Offshore System Development*, NREL/TP-500-38060.
- D. T. Griffith and T. D. Ashwill, *The Sandia 100-Meter All-Glass Baseline Wind Turbine Blade*, SAND2011-3779.
- Sandia National Laboratories, *SAND2013-2569*.
