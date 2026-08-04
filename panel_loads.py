"""Traceable conversions from spar-cap stress or beam moment to panel ``N_x``."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PanelLoad:
    """Applied axial compression resultant for plate buckling.

    ``Nx_N_per_m`` is positive in magnitude; the buckling module treats it as
    compression by convention.
    """

    Nx_N_per_m: float
    cap_force_N: float
    basis: str


def compressive_resultant_from_stress(stress_Pa: float, loaded_thickness_m: float) -> PanelLoad:
    """Convert a representative cap stress to a plate resultant ``N_x=sigma*t``.

    This assumes stress is uniform through the stated loaded thickness.  It is
    appropriate only as a screening conversion, not as a recovered ply stress.
    """
    if stress_Pa <= 0 or loaded_thickness_m <= 0:
        raise ValueError("stress and loaded thickness must be positive")
    nx = stress_Pa * loaded_thickness_m
    return PanelLoad(nx, nx, "screening conversion Nx = representative cap stress × loaded thickness")


def compressive_resultant_from_moment(
    moment_Nm: float, cap_separation_m: float, panel_width_m: float
) -> PanelLoad:
    """Convert beam moment to one-cap force then to panel resultant.

    Uses ``F_cap=M/z_caps`` and ``Nx=F_cap/b``.  It requires the *actual*
    vertical separation of the tension and compression cap resultants.
    """
    if moment_Nm <= 0 or cap_separation_m <= 0 or panel_width_m <= 0:
        raise ValueError("moment, cap separation and panel width must be positive")
    cap_force = moment_Nm / cap_separation_m
    return PanelLoad(cap_force / panel_width_m, cap_force,
                     "beam equilibrium: F_cap = M/z_caps; Nx = F_cap/b_panel")
