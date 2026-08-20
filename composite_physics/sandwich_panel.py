"""Screening mechanics for a symmetric faces/core/faces sandwich panel.

The module keeps four end-compression limit states separate:

* global orthotropic plate buckling with a core-shear correction;
* face wrinkling;
* core shear crimping;
* face material compression.

The closed-form wrinkling and crimping expressions follow Hexcel's
``Honeycomb Sandwich Design Technology`` guide.  They are order-of-magnitude
screening equations, not substitutes for validated shell/solid FEA or tests.
"""

from dataclasses import dataclass
import numpy as np

from .buckling_rbs import critical_uniaxial_compression_rbs
from .laminate import LaminateStiffness, assemble_laminate_stiffness

__all__ = [
    "SandwichMaterial",
    "SandwichPanel",
    "SandwichLimitState",
    "SandwichAssessment",
    "assess_sandwich_panel",
]


@dataclass(frozen=True)
class SandwichMaterial:
    face: dict
    core: dict
    face_compression_Pa: float
    core_compression_modulus_Pa: float
    core_shear_modulus_Pa: float


@dataclass(frozen=True)
class SandwichPanel:
    name: str
    width_m: float
    face_thickness_m: float
    core_thickness_m: float
    material: SandwichMaterial

    @property
    def face_centroid_separation_m(self) -> float:
        return self.core_thickness_m + self.face_thickness_m

    @property
    def total_thickness_m(self) -> float:
        return 2.0 * self.face_thickness_m + self.core_thickness_m

    def stiffness(self) -> LaminateStiffness:
        layup = [
            {"theta": 0.0, "t": self.face_thickness_m, "mat": 0},
            {"theta": 0.0, "t": self.core_thickness_m, "mat": 1},
            {"theta": 0.0, "t": self.face_thickness_m, "mat": 0},
        ]
        return assemble_laminate_stiffness(layup, [self.material.face, self.material.core])

    def membrane_resultant(self, axial_strain: float) -> float:
        """Return applied compression resultant assuming both faces share strain."""
        if axial_strain < 0:
            raise ValueError("axial_strain is a positive compression magnitude")
        return 2.0 * self.material.face["E1"] * self.face_thickness_m * axial_strain

    def face_stress(self, axial_strain: float) -> float:
        return self.material.face["E1"] * axial_strain


@dataclass(frozen=True)
class SandwichLimitState:
    mode: str
    critical_Nx_N_per_m: float


@dataclass(frozen=True)
class SandwichAssessment:
    panel: str
    a_over_b: float
    a_m: float
    applied_Nx_N_per_m: float
    global_clpt_Nx_N_per_m: float
    global_shear_corrected_Nx_N_per_m: float
    global_mode_m: int
    global_mode_n: int
    limits: tuple[SandwichLimitState, ...]
    governing_mode: str
    governing_critical_Nx_N_per_m: float
    reserve_factor: float
    margin_of_safety: float


def _face_wrinkling_stress(panel: SandwichPanel) -> float:
    mat = panel.material
    return 0.5 * (
        mat.core_shear_modulus_Pa
        * mat.core_compression_modulus_Pa
        * mat.face["E1"]
    ) ** (1.0 / 3.0)


def assess_sandwich_panel(
    panel: SandwichPanel,
    a_over_b: float,
    applied_Nx_N_per_m: float,
    max_mode: int = 12,
) -> SandwichAssessment:
    """Evaluate all four screening limit states for one assumed span ratio."""
    if a_over_b <= 0 or applied_Nx_N_per_m <= 0:
        raise ValueError("a_over_b and applied_Nx_N_per_m must be positive")

    stiffness = panel.stiffness()
    if not np.allclose(stiffness.B, 0.0, atol=1e-8):
        raise ValueError("the compact sandwich model requires a symmetric stack (B=0)")
    a_m = a_over_b * panel.width_m
    clpt = critical_uniaxial_compression_rbs(
        stiffness.D, a_m=a_m, b_m=panel.width_m, max_mode=max_mode
    )

    # Hexcel end-load correction written per unit width:
    # Ncr = N_E / (1 + N_E/S), S = Gc*h, h = face-centre separation.
    shear_rigidity = (
        panel.material.core_shear_modulus_Pa * panel.face_centroid_separation_m
    )
    global_corrected = clpt.Nx_cr_N_per_m / (
        1.0 + clpt.Nx_cr_N_per_m / shear_rigidity
    )

    wrinkling = 2.0 * panel.face_thickness_m * _face_wrinkling_stress(panel)
    crimping = panel.core_thickness_m * panel.material.core_shear_modulus_Pa
    face_compression = (
        2.0 * panel.face_thickness_m * panel.material.face_compression_Pa
    )
    limits = (
        SandwichLimitState("global_buckling", float(global_corrected)),
        SandwichLimitState("face_wrinkling", float(wrinkling)),
        SandwichLimitState("core_shear_crimping", float(crimping)),
        SandwichLimitState("face_compression", float(face_compression)),
    )
    governing = min(limits, key=lambda item: item.critical_Nx_N_per_m)
    rf = governing.critical_Nx_N_per_m / applied_Nx_N_per_m
    return SandwichAssessment(
        panel=panel.name,
        a_over_b=a_over_b,
        a_m=a_m,
        applied_Nx_N_per_m=applied_Nx_N_per_m,
        global_clpt_Nx_N_per_m=clpt.Nx_cr_N_per_m,
        global_shear_corrected_Nx_N_per_m=global_corrected,
        global_mode_m=clpt.m,
        global_mode_n=clpt.n,
        limits=limits,
        governing_mode=governing.mode,
        governing_critical_Nx_N_per_m=governing.critical_Nx_N_per_m,
        reserve_factor=rf,
        margin_of_safety=rf - 1.0,
    )
