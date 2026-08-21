# Phase 7 DB and H100 fatigue source basis

## FAT-006 — Saertex VU 90079 digitized DB reference

Mandell, J. F., Miller, D. A., and Samborsky, D. D., *Creep/Fatigue
Behavior of Resin Infused Biaxial Glass Fabric Laminates*, AIAA 2013-1630,
2013.

Source: https://www.montana.edu/composites/documents/SDM2013%20CreepFatigue%20Behavior%20of%20Resin%20Infused%20Biaxial%20Glass%20Fabric%20Laminates.pdf

- Table 1 identifies Fabric L as Saertex VU 90079, 96.8% ±45-degree fibre.
- The R=-1 Fabric L trend was digitized from Figure 7 as
  `S_axial = 156 N^-0.103 MPa`. This is a digitized trend-line fit, not a
  numeric fit stated in the paper text or a source table.
- Table 2 states the ASTM D3518 ply shear relation `tau_6 = S_axial/2`.
- Therefore the transformed digitized Figure 7 trend-line reference is
  `tau_6 = 78 N^-0.103 MPa`, equivalent to `C=78 MPa`, `b=9.709`
  in the project form `N=(C/tau)^b`. Both constants are classified as
  **digitized**, not as a stated numeric source fit.

The final design does not use the digitized curve to weaken the earlier
screen. It uses an explicitly estimated conservative envelope, `C=62 MPa`,
`b=9`. Relative to FAT-006, C is reduced by 20.5% and b is rounded downward.

## FAT-007 — Divinycell H100 measured shear-fatigue curve

Bednarcyk, B. A., Yarrington, P. W., and Arnold, S. M., *Multiscale Fatigue
Life Prediction for Composite Panels*, NASA/TM-2012-217694, 2012.

Source: https://ntrs.nasa.gov/citations/20120015395

- Section 4 and Figure 6 give the H100 shear-fatigue Basquin fit at R=0.1:
  `tau_max = 2.34 N^(-1/12.08) MPa`.
- Project form: `C=2.34 MPa`, `b=12.08`.
- Table 3 reports `G=40 MPa` and shear yield stress `1.13 MPa` for the tested
  H100. The accepted Sandia CAD model retains its project stiffness value;
  the measured NASA curve is used for fatigue damage.
- The project fatigue factor 1.38 is still applied, making the implementation
  more conservative than the plotted best-fit curve.

## FAT-008 — robust design envelope and load-split matrix

This is a traceable project extrapolation, not a measured material curve.

- DB robust envelope: FAT-006 transformed shear curve, then reduced from
  `C=78 MPa`, `b=9.709` to `C=62 MPa`, `b=9`.
- Load split: `V_i = V_total (A66_i h_i) / sum(A66 h)`, perturbed by ±10% on
  the forward fraction while conserving total shear.
- Required matrix: `b in {10,9}` crossed with forward split error
  `{-10%,0,+10%}` — six cases.
- Acceptance: maximum 20-year Miner damage `D20 <= 0.70` for every case.

The selected one-zone scheme is forward 5 DB layers/side, aft 6 DB
layers/side and 60 mm H100 core. The additional aft layer is structural: it
brings the governing nominal forward split to 0.487, close to equal sharing.
