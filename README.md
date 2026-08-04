# Composites basics formulas

![Python](https://img.shields.io/badge/python-3.12%2B-blue) ![numpy](https://img.shields.io/badge/numpy-2.5%2B-blue) ![tests](https://img.shields.io/badge/tests-14-lightgrey)

Classical lamination theory built up from first principles, plus one worked
application: a local buckling screening of the upper spar-cap panel of a DU35
wind-turbine blade section.

The library half is deliberately written the long way round — material
properties to `Q`, rotation to `Q̄`, through-thickness integration to `A`, `B`,
`D` — so each step can be checked on its own against a textbook. The
application half runs it on geometry measured from a FreeCAD model.

Everything here is a screening calculation. It is not certified analysis, and
there is no experimental validation.

<p align="center">
  <img src="docs/panel-section.png" width="720" alt="DU35 section with the screened panel between the two web centrelines">
</p>

---

## Quick start

Python 3.12+ and numpy.

```bash
uv sync                          # install
python du35_rbs_assessment.py    # reproduce the application case
pytest                           # 14 checks
```

The three library modules also run standalone and print a validation case
against a published reference:

```bash
python lamina_mechanics.py       # Q and S for a single ply
python transformations.py        # Q̄ at 30°, against Jones Example 2.3
python laminate.py               # A, B, D for a symmetric [0/±45/90]s layup
```

Figures under `docs/` are regenerated from live results, not stored numbers:

```bash
python figures/make_figures.py   # needs matplotlib
```

---

## The calculation chain

```
lamina_mechanics.py     E1, E2, G12, v12   ->  Q          ply in its own axes
transformations.py      Q, theta           ->  Q̄          rotated to global axes
laminate.py             Q̄, ply positions   ->  A, B, D    integrated through thickness
reduced_stiffness.py    A, B, D            ->  D*         coupling folded into bending
buckling_rbs.py         D*, a, b           ->  Nx,cr      critical load and mode
panel_loads.py          sigma, t   or   M  ->  Nx         applied resultant
                                              ------
du35_rbs_assessment.py                     ->  margin over the span sweep
```

Each stage is a separate module so it can be tested and replaced on its own.

---

## Expected output

`du35_rbs_assessment.py` prints the following. If these numbers move, either an
input changed or something broke.

```text
DU35 upper cap — RBS buckling screening
b = 0.810804 m; h_cap = 65.0 mm
Applied Nx = 7.085 MN/m (screening conversion Nx = representative cap stress × loaded thickness)
D* [kN m]:
[[1082.99315067  111.03516314    0.        ]
 [ 111.03516314  377.91941593    0.        ]
 [   0.            0.           82.37338969]]

--- method applicability ---
max |B|        = 1,265,824.2 N
B11            = 1,265,824.2 N   (mid-plane offset B11/A11 = 0.4395 mm)
D16,  D26      = 0, 0 N m
D*16, D*26     = 0, 0 N m
stack mirrored about mid-plane : False
reduced panel orthotropic      : True
-> REJECT D-only: the stack is extension-bending coupled, so a
   formula written in D alone does not apply to this panel.
-> ACCEPT RBS: D* carries the coupling and leaves no bend-twist
   terms, so the orthotropic Navier solution applies to D*.

a/b     a [m]    Nx,cr [MN/m]   m,n   MS
0.5    0.4054        74.735  (1,1)  +9.548
1.0    0.8108        30.213  (1,1)  +3.264
2.0    1.6216        30.213  (2,1)  +3.264
4.0    3.2432        27.513  (3,1)  +2.883
```

`pytest` should report 14 passed.

---

## The DU35 case

The upper spar cap is a 3 mm triaxial skin bonded to a 65 mm unidirectional
cap, spanning 0.8108 m between the two shear-web centrelines. When the blade
bends, that panel goes into compression and can buckle long before the material
itself is near its limit.

### Why the obvious formula does not apply

The standard route is the simply-supported plate formula written in the bending
stiffness `D` alone. It assumes the stack is mirrored about its mid-plane. This
one is not — the skin sits on one face only, which puts the stiffness centre
0.44 mm off the mid-plane and produces `B11 = 1.27 MN`.

<p align="center">
  <img src="docs/stack-eccentricity.png" width="720" alt="Mid-plane and stiffness centre 0.44 mm apart, with a magnified inset">
</p>

That is not a small inaccuracy to be absorbed — it is a violated assumption, so
`D` alone is inadmissible here.

### What replaces it

Reduced bending stiffness:

$$D^* = D - B\,A^{-1}B$$

`D*` is the bending stiffness a mirrored panel would need in order to behave
like this one, so the coupling is folded in rather than discarded. The
resulting `D*16` and `D*26` come out at exactly zero, which is what makes the
orthotropic Navier solution legitimate on `D*`; the critical load is then
minimised over buckling modes in the usual way.

Keeping `B` rather than dropping it changed `D11` by 0.060% and `D66` by 1.24%.
The coupling is real but weak in this particular stack — which could only be
established by computing it, not before.

### The unknown span

The spanwise support spacing `a` is not known: the CAD model is a 100 mm slice
and does not contain it. So `a` is swept from `a/b = 0.5` to `4` rather than
assumed.

<p align="center">
  <img src="docs/span-sweep.png" width="720" alt="Critical load against assumed span ratio, worst case 27.49 MN/m against 7.085 MN/m applied">
</p>

The curve is festooned because the panel buckles into a whole number of
half-waves, and switches to one more as it lengthens. The worst case anywhere
in the swept range is **27.49 MN/m against 7.085 MN/m applied — a ratio of 3.88
and a margin of safety of +2.88.**

---

## Files

### Library

| Module | What it does |
|---|---|
| `lamina_mechanics.py` | A single ply in its own axes. Builds the reduced stiffness `Q` and compliance `S` from `E1`, `E2`, `G12`, `v12`; isotropic and fully anisotropic variants; recovery of engineering constants from `Q`; plane-stress reduction of a 3D stiffness. Validated against Jones. |
| `transformations.py` | Rotation between material (1-2) and global (x-y) axes. Stress and strain transformation matrices, and `Q̄ = T⁻¹ Q T_strain`. The off-diagonal `Q̄16`, `Q̄26` terms that appear for off-axis plies are where normal–shear coupling comes from. Validated against Jones Example 2.3. |
| `laminate.py` | Assembly of `A`, `B`, `D` and the 6×6 `ABD` by integrating `Q̄` through the thickness. Ply interfaces run from the `−h/2` face to the `+h/2` face, so reversing the list reverses the laminate. Input validation is deliberately strict. |
| `clt.py` | Thin compatibility wrapper exposing the older `compute_ABD_matrix(layup, materials)` entry point. |
| `laminate_special_cases.py` | `is_specially_orthotropic()` — whether the coupling terms a simplified plate formula ignores actually vanish. |
| `reduced_stiffness.py` | `D* = D − B A⁻¹ B`, computed with `linalg.solve` rather than by forming `A⁻¹`. Symmetry of the inputs is checked. |
| `buckling_rbs.py` | Simply-supported orthotropic plate buckling under uniaxial compression, minimised over modes `m, n = 1…12`. Refuses to run when `D*16` or `D*26` are not negligible, so it cannot be applied outside its own assumptions by accident. |
| `panel_loads.py` | Stress or beam moment to panel resultant `Nx`. Two independent paths — `σ·t`, and `M/z_caps` then `/b` — kept separate so the load can be traced or replaced without touching the buckling solver. |

### Application

| Module | What it does |
|---|---|
| `airfoil_panel_laminate.py` | The DU35 case as data: measured widths, thicknesses, chord, and the two-zone equivalent material model. Every geometric number in the assessment originates here. Also keeps a symmetric quasi-isotropic skin case for teaching. |
| `du35_rbs_assessment.py` | Runs the whole thing: geometry → `ABD` → `D*` → applicability verdict → critical load and margin across the span sweep. |
| `figures/make_figures.py` | Regenerates `docs/stack-eccentricity.png` and `docs/span-sweep.png` from live results. |

### Tests

| File | What it covers |
|---|---|
| `tests/test_laminate.py` | Physical invariants of the assembly: a symmetric layup gives `B = 0` and `D* = D`; doubling every ply thickness scales `A` by 2 and `D` by 8; reversing the stack leaves `A` and `D` alone and flips the sign of `B`; a pure 0° laminate matches the closed forms `A = Qh` and `D = Qh³/12`; the `clt.py` wrapper returns the same `ABD`; special-orthotropy classification; rejection of invalid input. |
| `tests/test_rbs_buckling.py` | `D*` stays symmetric and is reduced by the coupling; both load conversions give the expected units; the margin is reported. |
| `tests/test_airfoil_panel_laminate.py` | The stored CAD geometry is what it should be, the panel really is coupled, and the teaching case still has `B = 0` with non-zero `D16`, `D26`. |

These test the mathematics of the assembly against closed-form limits. They say
nothing about whether the model resembles a real blade.

---

## Inputs, and how far they can be trusted

| Input | Value | Status |
|---|---|---|
| Panel width `b` | 0.810804 m | Measured — midsurface arc between the two web centrelines. The straight projection, 0.808537 m, is also stored; using it gives 27.67 instead of 27.49 MN/m. |
| Section chord | 4.491873 m | Measured. |
| Skin / cap thickness | 3 mm / 65 mm | Measured geometry. |
| Triax and UD properties | see `airfoil_panel_laminate.py` | Assumed — an equivalent two-zone engineering model, not a manufacturing ply book. |
| Cap stress | 109 MPa | A screening value carried over from earlier beam work. Not an IEC load case. |
| Spanwise spacing `a` | — | Unknown. Swept over `a/b = 0.5…4`; nothing outside that range has been calculated. |
| Edge conditions | simply supported | Assumed. Real restraint is stiffer and would raise the critical load. |

RBS is an approximation. It carries the effect of `B` into a bending-only
formula; it does not replace a fully coupled eigenvalue solution, and none has
been run.

Only the stiffest panel in the section has been screened. The thin leading- and
trailing-edge sandwich panels — the likelier buckling drivers — have not been
evaluated.

---

## References

- R. M. Jones, *Mechanics of Composite Materials*, 2nd ed., Ch. 2 and 4.
- T. Coburn, *Composite Strength*, lectures L-04 and L-08.
- J. Schilling and C. Mittelstedt, "Studies on the validity of the reduced
  bending stiffness method for eccentrically stacked laminates", *PAMM* (2021),
  DOI [10.1002/pamm.202000199](https://doi.org/10.1002/pamm.202000199).
- D. T. Griffith and T. D. Ashwill, *The Sandia 100-Meter All-Glass Baseline
  Wind Turbine Blade: SNL100-00*, SAND2011-3779.
- J. Jonkman et al., *Definition of a 5-MW Reference Wind Turbine for Offshore
  System Development*, NREL/TP-500-38060.
