"""FreeCAD-measured geometry and traceable load mapping for upper DU35 panels."""

from dataclasses import dataclass
from pathlib import Path
import csv
import numpy as np

CHORD_M = 4.491873
FACE_THICKNESS_M = 0.003
CORE_THICKNESS_M = 0.025


@dataclass(frozen=True)
class UpperPanelGeometry:
    name: str
    x_over_c_start: float
    x_over_c_end: float
    midsurface_width_m: float
    centroid_y_m: float


LE_UPPER = UpperPanelGeometry("LE_upper_sandwich", 0.01, 0.15,
                              0.8024399149355272, 0.4027816968651594)
TE_UPPER = UpperPanelGeometry("TE_upper_foam", 0.50, 0.80,
                              1.3879939402728796, 0.4222627226531165)

# B-Rep centres of the top and bottom spar-cap solids.  Their midpoint is used
# only to map the 2D CAD location into the normalized through-depth coordinate
# of the existing beam model.
CAD_TOP_CAP_Y_M = 0.6366488741693934
CAD_BOTTOM_CAP_Y_M = -0.778913209758633
CAD_REFERENCE_AXIS_Y_M = 0.5 * (CAD_TOP_CAP_Y_M + CAD_BOTTOM_CAP_Y_M)


@dataclass(frozen=True)
class BeamLoadSnapshot:
    r_m: float
    chord_m: float
    moment_Nm: float
    EI_Nm2: float
    cap_separation_m: float
    curvature_per_m: float
    cap_strain: float


def load_inner_beam_snapshot(
    csv_path: Path | None = None,
    target_chord_m: float = CHORD_M,
) -> BeamLoadSnapshot:
    """Interpolate the increasing-chord inner branch of the existing beam model."""
    if csv_path is None:
        phase4_root = Path(__file__).resolve().parents[3]
        csv_path = phase4_root / "results" / "glass_station_results.csv"
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    inner = [row for row in rows if float(row["r_m"]) <= 14.35]
    inner.sort(key=lambda row: float(row["chord_m"]))
    for lower, upper in zip(inner, inner[1:]):
        c0, c1 = float(lower["chord_m"]), float(upper["chord_m"])
        if c0 <= target_chord_m <= c1:
            f = (target_chord_m - c0) / (c1 - c0)
            interp = lambda key: float(lower[key]) + f * (float(upper[key]) - float(lower[key]))
            moment = interp("M_park_fault_design_Nm")
            ei = interp("EI_flap_Nm2")
            separation = interp("cap_separation_m")
            curvature = moment / ei
            return BeamLoadSnapshot(
                r_m=interp("r_m"), chord_m=target_chord_m, moment_Nm=moment,
                EI_Nm2=ei, cap_separation_m=separation, curvature_per_m=curvature,
                cap_strain=curvature * separation / 2.0,
            )
    raise ValueError("target chord is not bracketed on the inner blade branch")


def normalized_upper_position(panel: UpperPanelGeometry) -> float:
    """Map a CAD panel centroid to 0 at the cap midpoint and 1 at the top cap."""
    return ((panel.centroid_y_m - CAD_REFERENCE_AXIS_Y_M)
            / (CAD_TOP_CAP_Y_M - CAD_REFERENCE_AXIS_Y_M))


def beam_panel_strain(panel: UpperPanelGeometry, beam: BeamLoadSnapshot) -> float:
    return normalized_upper_position(panel) * beam.cap_strain
