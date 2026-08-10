# Composites basics formulas

Classical lamination theory built up from first principles, followed by two
worked screening applications for the upper structure of a DU35 wind-turbine
blade section:

1. reduced-bending-stiffness buckling of the upper spar-cap panel;
2. compression screening of the upper leading-edge and trailing-edge sandwich
   panels.

The library is deliberately written as a visible calculation chain: material
properties to `Q`, rotation to `Qbar`, through-thickness integration to `A`,
`B`, `D`, reduction to `D*`, and finally plate or sandwich-panel buckling. Each
stage can therefore be checked or replaced independently.

The applications use geometry measured from a FreeCAD model. They are
engineering screening calculations, not certified analyses, and they have not
been validated experimentally.

---

## Quick start

Requirements: Python 3.12+ and NumPy.

```bash
uv sync

# Day 2: upper spar-cap reduced-bending-stiffness screening
python du35_rbs_assessment.py

# Day 3: upper sandwich-panel screening and CSV export
python du35_sandwich_assessment.py

# Complete validation suite
python -m unittest discover -s tests -v
```

The current validation suite contains **20 checks**. Six of them are specific
to the Day-3 sandwich extension; the other fourteen cover the laminate and
spar-cap calculation chain.

The three base modules can also run independently and print reference cases:

```bash
python lamina_mechanics.py
python transformations.py
python laminate.py
```

---

## Calculation chain

```text
lamina_mechanics.py     E1, E2, G12, v12     -> Q
transformations.py      Q, ply angle         -> Qbar
laminate.py             Qbar, ply positions  -> A, B, D
reduced_stiffness.py    A, B, D              -> D*
buckling_rbs.py         D*, a, b             -> critical plate load
panel_loads.py          stress or moment      -> applied Nx
```

The two DU35 applications reuse that library through separate, traceable
paths:

```text
Day 2
airfoil_panel_laminate.py
        -> du35_rbs_assessment.py
        -> upper spar-cap critical load and margin across a/b

Day 3
du35_sandwich_geometry.py     measured panel widths + beam-load mapping
sandwich_panel.py             four sandwich compression limit states
        -> du35_sandwich_assessment.py
        -> load cases + a/b sweep + panel comparison CSV files
```

---

## Current engineering results

### Day 2 — upper spar cap

The upper spar-cap panel is represented by a 3 mm triaxial skin bonded to a
65 mm unidirectional cap. Its measured midsurface width between the two web
centrelines is 0.810804 m.

The stack is not mirrored about its mid-plane, so the standard plate equation
written directly in `D` is not admissible. The calculation uses reduced bending
stiffness instead:

```text
D* = D - B A^-1 B
```

Within the swept range `a/b = 0.5...4.0`, the minimum critical load is
**27.490 MN/m** against **7.085 MN/m** applied. The corresponding reserve factor
is **3.880** and the margin of safety is **+2.880**.

This is a screening result under simply supported edge conditions and an
assumed equivalent-material stack. It is not a manufacturing laminate or a
certification result.

### Day 3 — upper sandwich panels

Two additional thin upper-shell regions are evaluated as symmetric
face/core/face sandwich panels:

| Panel | Chordwise zone | Midsurface width | Stack |
|---|---:|---:|---:|
| LE upper sandwich | `0.01c-0.15c` | 0.8024 m | 3/25/3 mm |
| TE upper foam | `0.50c-0.80c` | 1.3880 m | 3/25/3 mm |

The width is the developed arc along the foam midsurface, not its horizontal
projection. This is the length of material represented by the equivalent flat
plate. The actual spanwise support spacing is not available from the
cross-section, so `a/b` is swept from 0.5 to 4.0 in 141 steps rather than set to
one invented value.

Four compression limit states are evaluated independently:

1. global orthotropic buckling with a core-shear correction;
2. face wrinkling;
3. core shear crimping;
4. face material compression.

The smallest critical load governs. For both sandwich panels, the governing
mode is global buckling with the core-shear correction.

#### Two load mappings

The earlier 109 MPa result belongs to the UD spar cap and cannot be assigned
directly to the triaxial shell. Two load paths are therefore retained:

- `109MPa_reference`: the cap stress is converted to cap strain, and the same
  strain is imposed on both sandwich faces. This gives 0.433 MN/m and is a
  deliberately severe comparison case.
- `beam_mapped`: the shell strain is mapped from the existing beam-model
  curvature and each panel's measured vertical location. This gives
  0.118 MN/m for LE and 0.123 MN/m for TE.

Neither value is a verified shell load. Their difference represents the
unresolved load sharing between the shell and the spar cap.

#### Governing comparison

The minimum reserve factor within the swept support-spacing range is:

| Item | Load case | Applied Nx | Critical Nx | Reserve factor | Governing mode |
|---|---|---:|---:|---:|---|
| LE upper sandwich | 109 MPa reference | 0.433 MN/m | 0.717 MN/m | **1.654** | global buckling |
| LE upper sandwich | beam mapped | 0.118 MN/m | 0.717 MN/m | **6.088** | global buckling |
| TE upper foam | 109 MPa reference | 0.433 MN/m | 0.368 MN/m | **0.848** | global buckling |
| TE upper foam | beam mapped | 0.123 MN/m | 0.368 MN/m | **2.998** | global buckling |
| Upper spar cap | 109 MPa reference | 7.085 MN/m | 27.490 MN/m | **3.880** | RBS buckling |
| Upper spar cap | beam mapped | 2.875 MN/m | 27.490 MN/m | **9.561** | RBS buckling |

**The TE upper foam panel governs the current screening.** It does not pass the
severe reference case, but it passes the beam-mapped case. The final design
verdict therefore requires the actual transverse-support spacing and validated
load sharing between the shell and the spar cap.

---

## Day-3 result files

The Day-3 assessment writes three tracked, English-language result tables used
by the report:

| File | Contents |
|---|---|
| `results/sandwich/du35_load_cases.csv` | beam station, moment, `EI`, curvature, cap strain, and reference stress |
| `results/sandwich/du35_sandwich_sweep.csv` | all 141 `a/b` points for both panels and both load cases, including all four critical loads and RF |
| `results/sandwich/du35_panel_comparison.csv` | compact governing comparison of LE, TE, and the upper spar cap |

These files are generated outputs. They should not be edited manually. Change
the input or calculation code and run `du35_sandwich_assessment.py` again.

The additional file `results/sandwich/du35_sandwich_geometry.csv` may be
generated locally as a geometry audit, but it is not required for the reported
result because the measured values are stored and tested in
`du35_sandwich_geometry.py`.

---

## Files

### Core library

| Module | Responsibility |
|---|---|
| `lamina_mechanics.py` | Builds reduced ply stiffness `Q` and compliance `S` from engineering constants and recovers constants from stiffness. |
| `transformations.py` | Rotates stiffness between material and global axes and produces `Qbar`. |
| `laminate.py` | Integrates `Qbar` through the thickness to assemble `A`, `B`, `D`, and `ABD`. |
| `clt.py` | Compatibility wrapper for the earlier `compute_ABD_matrix()` entry point. |
| `laminate_special_cases.py` | Checks whether the coupling terms ignored by simplified plate formulae vanish. |
| `reduced_stiffness.py` | Computes `D* = D - B A^-1 B` without explicitly forming `A^-1`. |
| `buckling_rbs.py` | Simply supported orthotropic plate buckling under uniaxial compression, minimized over integer modes. |
| `panel_loads.py` | Keeps the stress-based and beam-moment-based conversions to `Nx` separate. |

### DU35 applications

| Module | Responsibility |
|---|---|
| `airfoil_panel_laminate.py` | Stores the measured spar-cap geometry and equivalent triax/UD material model. |
| `du35_rbs_assessment.py` | Runs the Day-2 spar-cap RBS applicability and buckling sweep. |
| `du35_sandwich_geometry.py` | Stores measured LE/TE sandwich geometry and maps the existing beam result to each panel location. |
| `sandwich_panel.py` | Implements the symmetric sandwich stack, applied membrane resultant, four compression limit states, core-shear correction, and governing RF. |
| `du35_sandwich_assessment.py` | Runs both Day-3 load paths, the 141-point spacing sweep, the spar-cap comparison, and CSV export. |

### Tests

| File | Scope |
|---|---|
| `tests/test_laminate.py` | laminate assembly invariants, thickness scaling, stack reversal, closed-form limits, compatibility wrapper, and input rejection |
| `tests/test_rbs_buckling.py` | reduced stiffness, both load conversions, buckling output, and margin reporting |
| `tests/test_airfoil_panel_laminate.py` | stored spar-cap CAD geometry, coupling, and the symmetric teaching case |
| `tests/test_sandwich_panel.py` | symmetric sandwich stack, measured LE/TE widths, published wrinkling and crimping expressions, shear correction, core-thickness trend, and distinct load paths |

The tests check the implemented mathematics, stored geometry, and selected
reference behaviours. They do not demonstrate that the simplified model
reproduces a real blade.

---

## Required external input

The beam-mapped Day-3 path reads:

```text
../Wind_wing/results/structural/glass_station_results.csv
```

The file supplies the existing blade-model station results from which moment,
`EI`, cap separation, curvature, and cap strain are interpolated. The
`Composites_physics` and `Wind_wing` projects must therefore preserve this
relative folder arrangement, or an explicit CSV path must be passed to
`load_inner_beam_snapshot()`.

The annotated CAD model used to measure the section is stored separately in
the [project Google Drive](https://drive.google.com/drive/folders/1568kY9HVPt-y1rDeaz4n72Q8LZ079cSK):

```text
DU35_BoxSpar_2Web_v2_Sandwich_Annotated.FCStd
```

---

## Input status and limits

| Input or assumption | Current value | Evidence state |
|---|---|---|
| DU35 chord | 4.491873 m | measured |
| Spar-cap panel width | 0.810804 m midsurface arc | measured |
| LE sandwich width | 0.802440 m foam midsurface arc | measured |
| TE sandwich width | 1.387994 m foam midsurface arc | measured |
| Spar-cap skin/cap thickness | 3/65 mm | measured geometry |
| Sandwich face/core/face thickness | 3/25/3 mm | measured geometry used by the screening model |
| Triax, UD, and H100 properties | values stored in the application modules | assumed equivalent engineering properties |
| 109 MPa cap stress | carried over from earlier beam work | reference case, not a verified shell load |
| Beam-mapped shell load | calculated from the existing beam model | calculated screening load, not locally validated |
| Support spacing `a` | unknown; `a/b = 0.5...4.0` swept | assumed range |
| Panel edges | simply supported | assumed |

Not evaluated here:

- measured shell-to-cap load sharing;
- actual rib or transverse-support spacing;
- shell curvature and coupled shell/web/cap behaviour;
- adhesive failure, debonding, ply-level failure, and fatigue;
- manufacturing tolerances and certification load factors;
- experimental validation.

If the TE panel remains close to `RF = 1` after the real spacing and load split
are supplied, the next step is a coupled shell finite-element model.

---

## References

- R. M. Jones, *Mechanics of Composite Materials*, 2nd ed., Chapters 2 and 4.
- T. Coburn, *Composite Strength*, lectures L-04 and L-08.
- J. Schilling and C. Mittelstedt, "Studies on the validity of the reduced
  bending stiffness method for eccentrically stacked laminates", *PAMM* (2021),
  DOI [10.1002/pamm.202000199](https://doi.org/10.1002/pamm.202000199).
- Hexcel, *Honeycomb Sandwich Design Technology*.
- D. T. Griffith and T. D. Ashwill, *The Sandia 100-Meter All-Glass Baseline
  Wind Turbine Blade: SNL100-00*, SAND2011-3779.
- J. Jonkman et al., *Definition of a 5-MW Reference Wind Turbine for Offshore
  System Development*, NREL/TP-500-38060.
