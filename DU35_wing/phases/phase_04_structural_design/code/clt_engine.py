"""
clt_engine.py — Classical Laminate Theory engine (Phase 0 of blade project).

Units: SI (Pa, m, N/m for line loads, N for moments per unit width).
Conventions follow Jones, "Mechanics of Composite Materials", 2nd ed.

Verified against:
  - Jones T300/5208 reduced stiffness benchmark (see verify_clt.py)
  - Analytical invariants: B=0 (symmetric), A16=A26=0 (balanced),
    in-plane isotropy of quasi-isotropic laminates
  - Tsai-Wu closed-form limit cases

NOTE: screening-level tool. Not a substitute for validated FEA or
certified software for final design.
"""

from dataclasses import dataclass, field
import numpy as np


# ----------------------------------------------------------------------
# Material
# ----------------------------------------------------------------------

@dataclass
class Material:
    """Unidirectional lamina, plane stress."""
    name: str
    E1: float          # longitudinal modulus [Pa]
    E2: float          # transverse modulus [Pa]
    G12: float         # in-plane shear modulus [Pa]
    v12: float         # major Poisson ratio [-]
    rho: float = 0.0   # density [kg/m^3]
    # Strengths (positive magnitudes) [Pa]
    Xt: float = 0.0    # longitudinal tension
    Xc: float = 0.0    # longitudinal compression
    Yt: float = 0.0    # transverse tension
    Yc: float = 0.0    # transverse compression
    S: float = 0.0     # in-plane shear
    source: str = ""   # provenance of the numbers

    @property
    def v21(self) -> float:
        return self.v12 * self.E2 / self.E1

    @property
    def Q(self) -> np.ndarray:
        """Reduced stiffness matrix in material axes [Pa]."""
        d = 1.0 - self.v12 * self.v21
        return np.array([
            [self.E1 / d,             self.v12 * self.E2 / d, 0.0],
            [self.v12 * self.E2 / d,  self.E2 / d,            0.0],
            [0.0,                     0.0,                    self.G12],
        ])


# ----------------------------------------------------------------------
# Transformations
# ----------------------------------------------------------------------

def qbar(Q: np.ndarray, theta_deg: float) -> np.ndarray:
    """Transformed reduced stiffness Q̄ for a ply at angle theta (deg)."""
    th = np.radians(theta_deg)
    c, s = np.cos(th), np.sin(th)
    Q11, Q12, Q22, Q66 = Q[0, 0], Q[0, 1], Q[1, 1], Q[2, 2]
    return np.array([
        [Q11*c**4 + 2*(Q12 + 2*Q66)*s**2*c**2 + Q22*s**4,
         (Q11 + Q22 - 4*Q66)*s**2*c**2 + Q12*(s**4 + c**4),
         (Q11 - Q12 - 2*Q66)*s*c**3 + (Q12 - Q22 + 2*Q66)*s**3*c],
        [(Q11 + Q22 - 4*Q66)*s**2*c**2 + Q12*(s**4 + c**4),
         Q11*s**4 + 2*(Q12 + 2*Q66)*s**2*c**2 + Q22*c**4,
         (Q11 - Q12 - 2*Q66)*s**3*c + (Q12 - Q22 + 2*Q66)*s*c**3],
        [(Q11 - Q12 - 2*Q66)*s*c**3 + (Q12 - Q22 + 2*Q66)*s**3*c,
         (Q11 - Q12 - 2*Q66)*s**3*c + (Q12 - Q22 + 2*Q66)*s*c**3,
         (Q11 + Q22 - 2*Q12 - 2*Q66)*s**2*c**2 + Q66*(s**4 + c**4)],
    ])


def stress_to_material(sig_xy: np.ndarray, theta_deg: float) -> np.ndarray:
    """Rotate stress vector [sx, sy, txy] from laminate to material axes."""
    th = np.radians(theta_deg)
    c, s = np.cos(th), np.sin(th)
    T = np.array([
        [c*c, s*s,  2*s*c],
        [s*s, c*c, -2*s*c],
        [-s*c, s*c, c*c - s*s],
    ])
    return T @ sig_xy


# ----------------------------------------------------------------------
# Laminate
# ----------------------------------------------------------------------

@dataclass
class Ply:
    material: Material
    theta: float   # deg
    t: float       # thickness [m]


@dataclass
class Laminate:
    plies: list = field(default_factory=list)  # list[Ply], bottom -> top

    @property
    def h(self) -> float:
        return sum(p.t for p in self.plies)

    def z_interfaces(self) -> np.ndarray:
        z = np.zeros(len(self.plies) + 1)
        z[0] = -self.h / 2.0
        for k, p in enumerate(self.plies):
            z[k + 1] = z[k] + p.t
        return z

    def abd(self):
        """Return A [N/m], B [N], D [N*m] matrices."""
        A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))
        z = self.z_interfaces()
        for k, p in enumerate(self.plies):
            Qb = qbar(p.material.Q, p.theta)
            A += Qb * (z[k+1] - z[k])
            B += Qb * (z[k+1]**2 - z[k]**2) / 2.0
            D += Qb * (z[k+1]**3 - z[k]**3) / 3.0
        return A, B, D

    def abd_full(self) -> np.ndarray:
        A, B, D = self.abd()
        return np.block([[A, B], [B, D]])

    def solve(self, N=(0, 0, 0), M=(0, 0, 0)):
        """Midplane strains eps0 and curvatures kappa for loads N [N/m], M [N]."""
        rhs = np.concatenate([np.asarray(N, float), np.asarray(M, float)])
        sol = np.linalg.solve(self.abd_full(), rhs)
        return sol[:3], sol[3:]

    def effective_constants(self):
        """In-plane engineering constants of a symmetric laminate."""
        A, B, _ = self.abd()
        if np.max(np.abs(B)) > 1e-3 * np.max(np.abs(A)) * self.h:
            raise ValueError("Effective constants defined for symmetric laminates (B ~ 0).")
        a = np.linalg.inv(A)
        h = self.h
        return {
            "Ex":  1.0 / (h * a[0, 0]),
            "Ey":  1.0 / (h * a[1, 1]),
            "Gxy": 1.0 / (h * a[2, 2]),
            "vxy": -a[0, 1] / a[0, 0],
        }

    def ply_stresses(self, eps0, kappa, where="mid"):
        """
        Per-ply stresses in MATERIAL axes [Pa].
        where: 'mid', 'bot', 'top' evaluation point within each ply.
        Returns list of (ply_index, theta, sigma_123).
        """
        z = self.z_interfaces()
        out = []
        for k, p in enumerate(self.plies):
            zk = {"bot": z[k], "top": z[k+1], "mid": 0.5*(z[k]+z[k+1])}[where]
            eps_xy = np.asarray(eps0) + zk * np.asarray(kappa)
            sig_xy = qbar(p.material.Q, p.theta) @ eps_xy
            out.append((k, p.theta, stress_to_material(sig_xy, p.theta)))
        return out


# ----------------------------------------------------------------------
# Failure criteria (plane stress). Return reserve factor R:
#   R > 1 -> load can be multiplied by R before failure; R < 1 -> failed.
# ----------------------------------------------------------------------

def tsai_wu_R(sig, m: Material) -> float:
    """Tsai-Wu strength ratio for stress state sig=[s1, s2, t12]."""
    s1, s2, t12 = sig
    F1  = 1.0/m.Xt - 1.0/m.Xc
    F2  = 1.0/m.Yt - 1.0/m.Yc
    F11 = 1.0/(m.Xt*m.Xc)
    F22 = 1.0/(m.Yt*m.Yc)
    F66 = 1.0/m.S**2
    F12 = -0.5*np.sqrt(F11*F22)
    a = (F11*s1*s1 + F22*s2*s2 + F66*t12*t12 + 2.0*F12*s1*s2)
    b = F1*s1 + F2*s2
    if a < 1e-30:
        return np.inf if b <= 0 else 1.0/b
    return (-b + np.sqrt(b*b + 4.0*a)) / (2.0*a)


def tsai_wu_index(sig, m: Material) -> float:
    """Tsai-Wu failure index (=1 at failure) for the given stress state."""
    s1, s2, t12 = sig
    F1  = 1.0/m.Xt - 1.0/m.Xc
    F2  = 1.0/m.Yt - 1.0/m.Yc
    F11 = 1.0/(m.Xt*m.Xc)
    F22 = 1.0/(m.Yt*m.Yc)
    F66 = 1.0/m.S**2
    F12 = -0.5*np.sqrt(F11*F22)
    return (F1*s1 + F2*s2 + F11*s1*s1 + F22*s2*s2
            + F66*t12*t12 + 2.0*F12*s1*s2)


def hashin_R(sig, m: Material) -> dict:
    """
    Hashin 2D reserve factors per mode. Transverse shear strength
    approximated as Yc/2 (common screening assumption).
    """
    s1, s2, t12 = sig
    St = m.Yc / 2.0
    out = {}
    # Fiber modes
    if s1 >= 0:
        f = (s1/m.Xt)**2 + (t12/m.S)**2
        out["fiber_tension"] = 1.0/np.sqrt(f) if f > 0 else np.inf
    else:
        out["fiber_compression"] = m.Xc/abs(s1) if s1 != 0 else np.inf
    # Matrix modes
    if s2 >= 0:
        f = (s2/m.Yt)**2 + (t12/m.S)**2
        out["matrix_tension"] = 1.0/np.sqrt(f) if f > 0 else np.inf
    else:
        # quadratic in R: a R^2 + b R = 1
        a = (s2/(2.0*St))**2 + (t12/m.S)**2
        b = ((m.Yc/(2.0*St))**2 - 1.0) * (s2/m.Yc)
        if a < 1e-30:
            out["matrix_compression"] = np.inf if b <= 0 else 1.0/b
        else:
            out["matrix_compression"] = (-b + np.sqrt(b*b + 4.0*a)) / (2.0*a)
    return out


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def make_laminate(stack, material: Material, t_ply: float, symmetric=True) -> Laminate:
    """
    stack: list of angles for the half-laminate (if symmetric=True) or
           the full laminate (symmetric=False), bottom -> mid.
    """
    angles = list(stack) + (list(reversed(stack)) if symmetric else [])
    return Laminate([Ply(material, a, t_ply) for a in angles])
