"""Run every public verification suite in dependency order."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
PHASES = ROOT / "DU35_wing" / "phases"

CODE_DIRS = [
    ROOT,
    PHASES / "phase_02_aerodynamics" / "code",
    PHASES / "phase_03_materials" / "code",
    PHASES / "phase_04_structural_design" / "code",
    PHASES / "phase_04_structural_design" / "studies" / "du35_rbs" / "code",
    PHASES / "phase_04_structural_design" / "studies" / "du35_sandwich" / "code",
    PHASES / "phase_05_fea" / "code",
    PHASES / "phase_06_optimization" / "code",
    PHASES / "phase_07_fatigue" / "code",
    PHASES / "phase_08_final" / "code",
]

TESTS = [
    ROOT / "tests" / "core" / "test_laminate.py",
    PHASES / "phase_04_structural_design" / "studies" / "du35_rbs" / "tests" / "test_airfoil_panel_laminate.py",
    PHASES / "phase_04_structural_design" / "studies" / "du35_rbs" / "tests" / "test_rbs_buckling.py",
    PHASES / "phase_04_structural_design" / "studies" / "du35_sandwich" / "tests" / "test_sandwich_panel.py",
    PHASES / "phase_02_aerodynamics" / "tests" / "test_bem.py",
    PHASES / "phase_03_materials" / "tests" / "test_materials.py",
    PHASES / "phase_04_structural_design" / "tests" / "test_clt.py",
    PHASES / "phase_04_structural_design" / "tests" / "test_structural.py",
    PHASES / "phase_05_fea" / "tests" / "test_fea.py",
    PHASES / "phase_06_optimization" / "tests" / "test_phase6.py",
    PHASES / "phase_07_fatigue" / "tests" / "test_phase7.py",
    PHASES / "phase_08_final" / "tests" / "test_phase8.py",
]


def main() -> int:
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(
        [str(path) for path in CODE_DIRS] + ([existing] if existing else [])
    )
    env["PYTHONUTF8"] = "1"

    for index, test in enumerate(TESTS, start=1):
        relative = test.relative_to(ROOT)
        print(f"\n=== [{index}/{len(TESTS)}] {relative} ===", flush=True)
        result = subprocess.run(
            [sys.executable, "-X", "utf8", str(test)],
            cwd=ROOT,
            env=env,
            check=False,
        )
        if result.returncode:
            print(f"FAILED: {relative}", file=sys.stderr)
            return result.returncode

    print(f"\nAll {len(TESTS)} public verification suites passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
