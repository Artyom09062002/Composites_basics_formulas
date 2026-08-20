# Airfoil List — NREL 5MW Blade

Source: REF-001 (NREL/TP-500-38060), Table 6-1; AERO-001, AERO-003.

## Airfoils used

| Name | r/R range | t/c [%] | Notes |
|------|-----------|---------|-------|
| Cylinder1 | 0–0.044 | 100 | Root cylinder, structural only, Cd=0.50 |
| Cylinder2 | 0.044–0.093 | 100 | Root transition cylinder, Cd=0.35 |
| DU40_A17 | 0.093–0.130 | 40 | Delft University thick root section |
| DU35_A17 | 0.130–0.220 | 35 | Delft University |
| DU30_A17 | 0.220–0.255 | 30 | Delft University |
| DU25_A17 | 0.255–0.360 | 25 | Delft University |
| DU21_A17 | 0.360–0.440 | 21 | Delft University |
| NACA64_A17 | 0.440–1.000 | 18 | NACA 64-618; standard thin outboard section |

## Where to get polar data

- **OpenFAST repository:** https://github.com/NREL/openfast — path: `reg_tests/r-test/glue-codes/fast/5MW_Baseline/AeroData/`
  Files: `Cylinder1.dat`, `DU40_A17.dat`, `DU35_A17.dat`, `DU30_A17.dat`, `DU25_A17.dat`, `DU21_A17.dat`, `NACA64_A17.dat`
  Format: AeroDyn v15, columns = [Alpha, Cl, Cd, Cm] at multiple Re.
- **ECN report (AERO-002):** DU-series polars at Re = 3–6×10⁶ (representative for 5MW).
- Note: polars include 360° extrapolation (Viterna method) needed for BEM parked DLC.

## Phase 2 action
Download or copy polar files into `data/reference_blade/airfoil_polars/` directory.
Parse into a dict {airfoil_name: DataFrame(alpha, Cl, Cd)} for BEM code.
