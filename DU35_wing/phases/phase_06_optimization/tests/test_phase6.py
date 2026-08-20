"""Gate 6 verification for the web-support DoE and screening surrogate."""

import sys
from pathlib import Path

PHASE_ROOT = Path(__file__).resolve().parents[1]
PHASES_ROOT = PHASE_ROOT.parent
for code_dir in (
    PHASE_ROOT / "code",
    PHASES_ROOT / "phase_05_fea" / "code",
    PHASES_ROOT / "phase_04_structural_design" / "code",
    PHASES_ROOT / "phase_03_materials" / "code",
    PHASES_ROOT / "phase_02_aerodynamics" / "code",
):
    sys.path.insert(0, str(code_dir))

from phase6_optimization import (BUCKLING_RF_TARGET, PLY_MULTIPLIERS,
                                  SUPPORT_SPACINGS_M, run_web_support_doe,
                                  validate_surrogate_and_select)

passed = failed = 0


def check(label, condition, note=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"PASS  {label}")
    else:
        failed += 1
        print(f"FAIL  {label}: {note}")


print("=" * 70)
print("PHASE 6 DOE AND SURROGATE VERIFICATION")
print("=" * 70)
doe = run_web_support_doe()
validation, selected = validate_surrogate_and_select(doe)

check("DoE contains every planned web/support combination",
      len(doe) == len(PLY_MULTIPLIERS) * len(SUPPORT_SPACINGS_M),
      f"rows={len(doe)}")
check("DoE includes feasible buckling candidates", doe.feasible_screening.any())
check("All selected candidates pass direct FE buckling verification",
      selected.direct_FE_pass.all(), selected.to_string(index=False))
check("Selected direct FE candidates meet RF target",
      (selected.direct_FE_web_buckling_RF >= BUCKLING_RF_TARGET).all())
check("Surrogate mean held-out buckling error is <= 10%",
      validation.relative_error_pct.abs().mean() <= 10.0,
      f"mean={validation.relative_error_pct.abs().mean():.2f}%")
check("Surrogate maximum held-out buckling error is <= 15%",
      validation.relative_error_pct.abs().max() <= 15.0,
      f"max={validation.relative_error_pct.abs().max():.2f}%")
check("Selected candidate table is finite and non-empty",
      not selected.empty and selected.select_dtypes(include="number").notna().all().all())

print("-" * 70)
print(f"Result: {passed}/{passed + failed} PASS")
if failed:
    raise SystemExit(1)
