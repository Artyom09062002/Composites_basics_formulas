"""
verify_clt.py — Phase 0 gate: CLT engine vs. textbook benchmark and
analytical invariants. All checks must PASS before the engine is used
in the blade project.

Benchmark 1: reduced stiffnesses of T300/5208 graphite/epoxy
  (Jones, "Mechanics of Composite Materials", 2nd ed.:
   E1=181 GPa, E2=10.3 GPa, G12=7.17 GPa, v12=0.28
   -> Q11=181.8, Q22=10.35, Q12=2.897, Q66=7.17 GPa)

Benchmarks 2-6: exact analytical properties any correct CLT must satisfy.
"""

import sys
from pathlib import Path

PHASE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PHASE_ROOT / "code"))
import numpy as np
from clt_engine import (Material, Laminate, Ply, make_laminate,
                        tsai_wu_R, tsai_wu_index)

GPa = 1e9
MPa = 1e6
results = []

def check(name, ok, detail=""):
    results.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


# --- Benchmark material: Jones T300/5208 --------------------------------
t300 = Material("T300/5208", E1=181*GPa, E2=10.3*GPa, G12=7.17*GPa, v12=0.28,
                Xt=1500*MPa, Xc=1500*MPa, Yt=40*MPa, Yc=246*MPa, S=68*MPa,
                source="Jones 2nd ed., typical graphite/epoxy")

# 1. Q-matrix vs textbook values
Q = t300.Q / GPa
ref = {"Q11": 181.8, "Q22": 10.35, "Q12": 2.897, "Q66": 7.17}
ok = (abs(Q[0,0]-ref["Q11"])/ref["Q11"] < 5e-3 and
      abs(Q[1,1]-ref["Q22"])/ref["Q22"] < 5e-3 and
      abs(Q[0,1]-ref["Q12"])/ref["Q12"] < 5e-3 and
      abs(Q[2,2]-ref["Q66"])/ref["Q66"] < 5e-3)
check("Q matrix matches Jones T300/5208 benchmark", ok,
      f"Q11={Q[0,0]:.1f}, Q22={Q[1,1]:.2f}, Q12={Q[0,1]:.3f}, Q66={Q[2,2]:.2f} GPa")

# 2. Single 0-deg ply: effective constants must equal lamina constants
single = Laminate([Ply(t300, 0.0, 1e-3)])
eff = single.effective_constants()
ok = (abs(eff["Ex"]-t300.E1)/t300.E1 < 1e-9 and
      abs(eff["Ey"]-t300.E2)/t300.E2 < 1e-9 and
      abs(eff["Gxy"]-t300.G12)/t300.G12 < 1e-9 and
      abs(eff["vxy"]-t300.v12) < 1e-9)
check("Single 0-deg ply: Ex=E1, Ey=E2, Gxy=G12, vxy=v12", ok)

# 3. Symmetric laminate: B = 0
lam_sym = make_laminate([0, 45, -45, 90], t300, 0.125e-3, symmetric=True)
A, B, D = lam_sym.abd()
ok = np.max(np.abs(B)) < 1e-6 * np.max(np.abs(A)) * lam_sym.h
check("Symmetric laminate: B = 0", ok, f"|B|max/|A|max*h = {np.max(np.abs(B))/(np.max(np.abs(A))*lam_sym.h):.1e}")

# 4. Balanced laminate: A16 = A26 = 0
lam_bal = make_laminate([45, -45, 45, -45], t300, 0.125e-3, symmetric=False)
A, _, _ = lam_bal.abd()
ok = abs(A[0,2]) < 1e-9*abs(A[0,0]) and abs(A[1,2]) < 1e-9*abs(A[0,0])
check("Balanced +-45 laminate: A16 = A26 = 0", ok)

# 5. Quasi-isotropic [0/45/-45/90]s: in-plane isotropy of A
lam_qi = make_laminate([0, 45, -45, 90], t300, 0.125e-3, symmetric=True)
A, _, _ = lam_qi.abd()
iso1 = abs(A[0,0]-A[1,1])/A[0,0]
iso2 = abs(A[2,2]-(A[0,0]-A[0,1])/2)/A[2,2]
ok = iso1 < 1e-9 and iso2 < 1e-9
check("Quasi-isotropic: A11=A22 and A66=(A11-A12)/2", ok,
      f"dev1={iso1:.1e}, dev2={iso2:.1e}")

# 6. Load-solve roundtrip: apply Nx, recompute Nx from strains
lam = lam_qi
Nx = 1e5  # N/m
eps0, kappa = lam.solve(N=(Nx, 0, 0))
A, B, D = lam.abd()
N_back = A @ eps0 + B @ kappa
ok = abs(N_back[0]-Nx)/Nx < 1e-9 and abs(N_back[1]) < 1e-6*Nx
check("ABD solve roundtrip: N -> (eps,kappa) -> N", ok)

# 7. Tsai-Wu closed-form limits: pure s1=Xt and pure s2=-Yc give R=1
r1 = tsai_wu_R(np.array([t300.Xt, 0, 0]), t300)
r2 = tsai_wu_R(np.array([0, -t300.Yc, 0]), t300)
i1 = tsai_wu_index(np.array([t300.Xt, 0, 0]), t300)
ok = abs(r1-1) < 1e-9 and abs(r2-1) < 1e-9 and abs(i1-1) < 1e-9
check("Tsai-Wu: R=1 exactly at uniaxial strength limits", ok,
      f"R(s1=Xt)={r1:.6f}, R(s2=-Yc)={r2:.6f}")

# --- Summary ------------------------------------------------------------
n_pass = sum(1 for _, ok, _ in results if ok)
print(f"\n{n_pass}/{len(results)} checks passed.")
sys.exit(0 if n_pass == len(results) else 1)
