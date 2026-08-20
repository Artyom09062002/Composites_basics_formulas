"""Gate 4 verification for the analytical blade structural model."""

import sys
from pathlib import Path
import numpy as np

PHASE_ROOT = Path(__file__).resolve().parents[1]
PHASES_ROOT = PHASE_ROOT.parent
for code_dir in (
    PHASE_ROOT / "code",
    PHASES_ROOT / "phase_03_materials" / "code",
    PHASES_ROOT / "phase_02_aerodynamics" / "code",
):
    sys.path.insert(0, str(code_dir))

from structural_model import (REFERENCE_MASS_KG, TIP_LIMIT_M, build_design,
                              integrate_beam, write_phase4_outputs)

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
print("PHASE 4 ANALYTICAL STRUCTURAL VERIFICATION")
print("=" * 70)

# Independent Euler-Bernoulli benchmark: cantilever under uniform load.
L, q, EI = 20.0, 1500.0, 8.0e8
x = np.linspace(0, L, 1001)
moment = 0.5 * q * (L - x) ** 2
_, displacement = integrate_beam(x, moment, np.full_like(x, EI))
exact = q * L**4 / (8 * EI)
check("Beam integration matches qL^4/(8EI) within 0.01%",
      abs(displacement[-1] / exact - 1) < 1e-4,
      f"numeric={displacement[-1]:.6f}, exact={exact:.6f}")

glass, hybrid = write_phase4_outputs()
for design in (glass, hybrid):
    check(f"{design.name}: cap CLT reserve factor >= 1", design.min_cap_rf >= 1,
          f"RF={design.min_cap_rf:.3f}")
    check(f"{design.name}: web CLT reserve factor >= 1", design.min_web_rf >= 1,
          f"RF={design.min_web_rf:.3f}")
    check(f"{design.name}: shell/edge reserve factor >= 1", design.min_shell_rf >= 1,
          f"RF={design.min_shell_rf:.3f}")
    check(f"{design.name}: U_op tip deflection <= {TIP_LIMIT_M:.1f} m",
          design.tip_uop_m <= TIP_LIMIT_M, f"tip={design.tip_uop_m:.3f} m")
    check(f"{design.name}: parked-fault tip deflection <= {TIP_LIMIT_M:.1f} m",
          design.tip_parked_m <= TIP_LIMIT_M, f"tip={design.tip_parked_m:.3f} m")
    counts = design.stations.loc[design.stations.n_cap_glass > 0, "n_cap_glass"].to_numpy()
    check(f"{design.name}: cap ply counts are even", np.all(counts % 2 == 0))

check("Glass baseline mass is within +10% of Sandia reference",
      glass.mass_kg <= 1.10 * REFERENCE_MASS_KG,
      f"mass={glass.mass_kg:.0f} kg")
check("Hybrid alternative is at least 10% lighter than glass baseline",
      hybrid.mass_kg <= 0.90 * glass.mass_kg,
      f"glass={glass.mass_kg:.0f}, hybrid={hybrid.mass_kg:.0f} kg")
check("Worst parked pitch reproduced at 175 deg", abs(glass.worst_pitch_deg - 175) < 0.1)

print("-" * 70)
print(f"Result: {passed}/{passed + failed} PASS")
if failed:
    raise SystemExit(1)
