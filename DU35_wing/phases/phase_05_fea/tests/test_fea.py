"""Gate 5 verification for the shell-equivalent finite-element model."""

import sys
from pathlib import Path

import numpy as np

PHASE_ROOT = Path(__file__).resolve().parents[1]
PHASES_ROOT = PHASE_ROOT.parent
for code_dir in (
    PHASE_ROOT / "code",
    PHASES_ROOT / "phase_04_structural_design" / "code",
    PHASES_ROOT / "phase_03_materials" / "code",
    PHASES_ROOT / "phase_02_aerodynamics" / "code",
):
    sys.path.insert(0, str(code_dir))

from fea_model import BladeBeamFE, build_phase5_model, write_phase5_outputs
from structural_model import build_design

passed = 0
failed = 0


def check(label, condition, note=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"PASS  {label}")
    else:
        failed += 1
        print(f"FAIL  {label}: {note}")


print("=" * 70)
print("PHASE 5 SHELL-EQUIVALENT FEA VERIFICATION")
print("=" * 70)

# Independent finite-element benchmark: clamped beam under uniform line load.
L, q, EI, mass = 20.0, 1500.0, 8.0e8, 120.0
x = np.linspace(0.0, L, 81)
benchmark = BladeBeamFE(x, np.full_like(x, EI), np.full_like(x, EI),
                        np.full_like(x, mass))
static = benchmark.static("flap", np.full_like(x, q), "uniform")
exact_tip = q * L**4 / (8.0 * EI)
exact_root = q * L**2 / 2.0
check("FE cantilever tip displacement matches qL^4/(8EI) within 0.01%",
      abs(static.tip_displacement_m / exact_tip - 1.0) < 1e-4,
      f"FE={static.tip_displacement_m:.8f}, exact={exact_tip:.8f}")
check("FE cantilever root moment matches qL^2/2 within 0.01%",
      abs(abs(static.root_moment_Nm) / exact_root - 1.0) < 1e-4,
      f"FE={static.root_moment_Nm:.2f}, exact={exact_root:.2f}")

result = write_phase5_outputs()
model, stations, cases = build_phase5_model("glass")
phase4 = build_design("glass")
static_by_case = {r.case: r for r in result["static"]}

for case, analytic_tip in (("U_op", phase4.tip_uop_m),
                           ("U_park_fault", phase4.tip_parked_m)):
    fe = static_by_case[case]
    target_moment = cases[case]["M_design"][0]
    check(f"{case}: FE root moment matches Phase 4 load envelope within 1%",
          abs(abs(fe.root_moment_Nm) / target_moment - 1.0) < 0.01,
          f"FE={fe.root_moment_Nm:.2f}, target={target_moment:.2f}")
    check(f"{case}: FE tip deflection matches Phase 4 integration within 5%",
          abs(abs(fe.tip_displacement_m) / analytic_tip - 1.0) < 0.05,
          f"FE={fe.tip_displacement_m:.4f}, analytical={analytic_tip:.4f}")

flap = result["flap_modal"].frequencies_hz
edge = result["edge_modal"].frequencies_hz
check("Flap modal frequencies are positive and ascending", np.all(flap > 0.0) and np.all(np.diff(flap) > 0.0))
check("Edge modal frequencies are positive and ascending", np.all(edge > 0.0) and np.all(np.diff(edge) > 0.0))
check("First flap mode clears 3P frequency (0.606 Hz) by at least 20%",
      flap[0] > 1.2 * 0.606, f"f1={flap[0]:.3f} Hz")

buckling = result["buckling"]
cap_rf = buckling["cap_buckling_RF"].replace(np.inf, np.nan).min()
web_rf = buckling["web_shear_buckling_RF"].replace(np.inf, np.nan).min()
check("Spar-cap compression buckling screen has RF >= 1", cap_rf >= 1.0,
      f"RF={cap_rf:.3f}")
check("Unstiffened-web screen detects the expected local buckling vulnerability",
      web_rf < 1.0, f"RF={web_rf:.4f}")

print("-" * 70)
print(f"Result: {passed}/{passed + failed} PASS")
if failed:
    raise SystemExit(1)
