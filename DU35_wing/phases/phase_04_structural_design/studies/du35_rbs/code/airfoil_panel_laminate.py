"""Measured DU35 upper spar-cap panel reduced to a two-zone CLT laminate.

The CAD provides geometry, not a manufacturing ply book.  Consequently this
module is deliberately an equivalent-material model: 3 mm triaxial skin plus
65 mm longitudinal UD cap.  It must not be read as a certified layup.
"""

from dataclasses import dataclass
import numpy as np

from composite_physics.laminate import LaminateStiffness, assemble_laminate_stiffness

CHORD_M = 4.491873
PANEL_WIDTH_M = 0.8108039732531884  # FreeCAD midsurface arc, web centreline to centreline
PROJECTED_PANEL_WIDTH_M = 0.80853714
SKIN_THICKNESS_M = 0.003
CAP_THICKNESS_M = 0.065

TRIAX_GFRP = {"E1": 27.7e9, "E2": 13.65e9, "G12": 7.20e9, "v12": 0.39}
UD_GFRP = {"E1": 41.8e9, "E2": 14.0e9, "G12": 2.63e9, "v12": 0.28}


@dataclass(frozen=True)
class PanelLaminateResult:
    stiffness: LaminateStiffness
    chord_m: float
    projected_width_m: float
    surface_width_m: float
    skin_thickness_m: float
    cap_thickness_m: float
    total_thickness_m: float
    compatible_with_current_buckling_model: bool


def analyze_du35_upper_panel() -> PanelLaminateResult:
    """Return the actual CAD-geometry / equivalent-material panel snapshot."""
    stack = [
        {"theta": 0.0, "t": SKIN_THICKNESS_M, "mat": 0},
        {"theta": 0.0, "t": CAP_THICKNESS_M, "mat": 1},
    ]
    stiffness = assemble_laminate_stiffness(stack, [TRIAX_GFRP, UD_GFRP])
    return PanelLaminateResult(
        stiffness=stiffness, chord_m=CHORD_M, projected_width_m=PROJECTED_PANEL_WIDTH_M,
        surface_width_m=PANEL_WIDTH_M, skin_thickness_m=SKIN_THICKNESS_M,
        cap_thickness_m=CAP_THICKNESS_M, total_thickness_m=SKIN_THICKNESS_M + CAP_THICKNESS_M,
        compatible_with_current_buckling_model=bool(np.max(np.abs(stiffness.B)) < 1.0),
    )


def analyze_du35_skin_demo() -> PanelLaminateResult:
    """Keep a symmetric quasi-isotropic skin example for teaching CLT."""
    stack = [{"theta": theta, "t": 0.375e-3, "mat": 0}
             for theta in (0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0)]
    stiffness = assemble_laminate_stiffness(stack, [UD_GFRP])
    return PanelLaminateResult(stiffness, CHORD_M, PROJECTED_PANEL_WIDTH_M, PANEL_WIDTH_M,
                               0.003, 0.0, 0.003, True)
