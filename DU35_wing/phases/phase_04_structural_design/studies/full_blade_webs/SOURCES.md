# Full-Blade Web Study Sources

- `SAND2013-2569`, Table 4 — blade length, stations, chord, twist and offsets.
- `SAND2013-2569`, Table 5 — Saertex DB and foam thickness and elastic data.
- `SAND2013-2569`, Table 11 — 600 mm spar-cap region width.
- `SAND2013-2569`, Table 13 — two DB layers per side and 50 mm foam core.
- NREL archived FAST `CertTest/5MW_Baseline/Airfoils` — official coordinates
  for the eight reference airfoils.
- Phase 4 factored parked-fault bending-moment envelope — applied web shear
  derived as `V = |dM/dr|`.

The study uses the primary Sandia elastic data: Saertex DB `t=1 mm`,
`Ex=13.6 GPa`, `Ey=13.3 GPa`, `Gxy=11.8 GPa`; foam `Ex=Ey=256 MPa`,
`Gxy=22 MPa`, `rho=200 kg/m3`. The `1.8 MPa` foam shear strength remains
explicitly identified as a conservative H100 proxy because Table 5 does not
state foam strength.
