# Phase 2 Reference Data — NREL 5 MW / Sandia 61.5 m Blade

**Primary sources:**
- REF-001: Jonkman et al., NREL/TP-500-38060, 2009 (turbine definition, geometry Table 6-1)
- REF-002: Resor, SAND2013-2569, Sandia, 2013 (structural layup, material zones, stiffness distribution)

## Files

| File | Content | Source |
|------|---------|--------|
| blade_geometry.csv | Chord, twist, airfoil name vs span station | REF-001 Table 6-1 |
| blade_stiffness.csv | Distributed EI_flap, EI_edge, GJ, EA, mass per unit span | REF-002 |
| blade_airfoils.md | Airfoil list and where to get polar data | REF-001, AERO-001 |
| blade_layup_zones.md | Structural zone definitions and layup schedule (Sandia model) | REF-002 |

## Coordinate system
- r = radial station from blade root [m]; r=0 at hub flange, r=61.5m at tip
- Chord c [m], twist θ [deg] (positive nose-up / towards feather)
- Sections follow IEC 61400-5 convention: flap-wise = out-of-rotor-plane, edge-wise = in-plane
