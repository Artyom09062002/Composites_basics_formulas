"""Phase 5 shell-equivalent finite-element model for the NREL 5 MW blade.

The available public inputs contain distributed beam properties rather than
section-coordinate and laminate-shell meshes.  This module therefore uses a
three-dimensional beam finite-element idealisation: independent flapwise and
edgewise Euler--Bernoulli bending fields share the Phase 4 spanwise mass and
stiffness.  Local web and spar-cap panel buckling is evaluated from the same
parametric layup schedule.  It is a screening FEA model, not a replacement for
a production shell/solid model with adhesive joints, offsets and nonlinear
contact.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.linalg import eigh

PHASE_ROOT = Path(__file__).resolve().parents[1]
PHASES_ROOT = PHASE_ROOT.parent
for code_dir in (
    PHASE_ROOT / "code",
    PHASES_ROOT / "phase_04_structural_design" / "code",
    PHASES_ROOT / "phase_03_materials" / "code",
    PHASES_ROOT / "phase_02_aerodynamics" / "code",
):
    sys.path.insert(0, str(code_dir))

AERO_DATA = PHASES_ROOT / "phase_02_aerodynamics" / "data"

from bem_solver import gravity_edge_moment, load_blade_mass
from materials_db import get_ply, load_materials
from structural_model import (CAP_EI_FRACTION, GAMMA_F, _geometry_at,
                              build_design, load_envelopes)


@dataclass
class StaticResult:
    case: str
    direction: str
    displacement_m: np.ndarray
    rotation_rad: np.ndarray
    root_shear_N: float
    root_moment_Nm: float
    tip_displacement_m: float


@dataclass
class ModalResult:
    direction: str
    frequencies_hz: np.ndarray
    mode_shapes: np.ndarray


def beam_element_stiffness(ei: float, length: float) -> np.ndarray:
    """Cubic Euler--Bernoulli bending element, DOFs [w1, theta1, w2, theta2]."""
    c = ei / length**3
    return c * np.array([
        [12.0, 6.0 * length, -12.0, 6.0 * length],
        [6.0 * length, 4.0 * length**2, -6.0 * length, 2.0 * length**2],
        [-12.0, -6.0 * length, 12.0, -6.0 * length],
        [6.0 * length, 2.0 * length**2, -6.0 * length, 4.0 * length**2],
    ])


def beam_element_mass(mass_per_length: float, length: float) -> np.ndarray:
    """Consistent Euler--Bernoulli mass matrix for a uniform element."""
    c = mass_per_length * length / 420.0
    return c * np.array([
        [156.0, 22.0 * length, 54.0, -13.0 * length],
        [22.0 * length, 4.0 * length**2, 13.0 * length, -3.0 * length**2],
        [54.0, 13.0 * length, 156.0, -22.0 * length],
        [-13.0 * length, -3.0 * length**2, -22.0 * length, 4.0 * length**2],
    ])


def distributed_load_from_shear(x: np.ndarray, shear: np.ndarray,
                                target_root_moment: float) -> np.ndarray:
    """Convert a BEM shear envelope to a non-negative FE line load.

    The BEM root/hub region has discrete interpolation artefacts.  The shape
    comes from -dV/dx; one scalar correction preserves the factored root
    moment used by the Phase 4 analytical model.
    """
    raw = np.maximum(-np.gradient(np.abs(shear), x), 0.0)
    raw_moment = np.trapz(raw * x, x)
    if raw_moment <= 0.0:
        raise ValueError("Cannot construct a distributed load from the shear envelope.")
    return raw * (target_root_moment / raw_moment)


class BladeBeamFE:
    """Clamped, variable-property blade finite-element model."""

    def __init__(self, x: np.ndarray, ei_flap: np.ndarray, ei_edge: np.ndarray,
                 mass_per_length: np.ndarray):
        self.x = np.asarray(x, dtype=float)
        self.ei_flap = np.asarray(ei_flap, dtype=float)
        self.ei_edge = np.asarray(ei_edge, dtype=float)
        self.mass_per_length = np.asarray(mass_per_length, dtype=float)
        if not (np.all(np.diff(self.x) > 0) and len(self.x) >= 3):
            raise ValueError("Stations must be strictly increasing and contain at least 3 nodes.")

    @property
    def n_dof(self) -> int:
        return 2 * len(self.x)

    def assemble(self, direction: str) -> tuple[np.ndarray, np.ndarray]:
        if direction not in {"flap", "edge"}:
            raise ValueError("direction must be 'flap' or 'edge'")
        ei = self.ei_flap if direction == "flap" else self.ei_edge
        k = np.zeros((self.n_dof, self.n_dof))
        m = np.zeros_like(k)
        for i, length in enumerate(np.diff(self.x)):
            ids = np.array([2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3])
            k[np.ix_(ids, ids)] += beam_element_stiffness(0.5 * (ei[i] + ei[i + 1]), length)
            m[np.ix_(ids, ids)] += beam_element_mass(
                0.5 * (self.mass_per_length[i] + self.mass_per_length[i + 1]), length)
        return k, m

    def static(self, direction: str, line_load: np.ndarray, case: str) -> StaticResult:
        line_load = np.asarray(line_load, dtype=float)
        if line_load.shape != self.x.shape:
            raise ValueError("line_load must be defined at every FE station")
        k, _ = self.assemble(direction)
        f = np.zeros(self.n_dof)
        for i, length in enumerate(np.diff(self.x)):
            q = 0.5 * (line_load[i] + line_load[i + 1])
            ids = np.array([2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3])
            f[ids] += q * np.array([length / 2.0, length**2 / 12.0,
                                    length / 2.0, -length**2 / 12.0])
        free = np.arange(2, self.n_dof)
        u = np.zeros(self.n_dof)
        u[free] = np.linalg.solve(k[np.ix_(free, free)], f[free])
        reactions = k @ u - f
        return StaticResult(case, direction, u[0::2], u[1::2],
                            float(reactions[0]), float(reactions[1]),
                            float(u[-2]))

    def modal(self, direction: str, n_modes: int = 4) -> ModalResult:
        k, m = self.assemble(direction)
        free = np.arange(2, self.n_dof)
        values, vectors = eigh(k[np.ix_(free, free)], m[np.ix_(free, free)],
                               subset_by_index=[0, n_modes - 1])
        frequencies = np.sqrt(np.maximum(values, 0.0)) / (2.0 * np.pi)
        modes = np.zeros((len(self.x), n_modes))
        modes[1:, :] = vectors[0::2, :]
        for i in range(n_modes):
            scale = np.max(np.abs(modes[:, i]))
            if scale > 0:
                modes[:, i] /= scale
        return ModalResult(direction, frequencies, modes)


def panel_buckling(stations: pd.DataFrame, design: str,
                   web_ply_multiplier: float = 1.0,
                   web_support_spacing_m: float | None = None) -> pd.DataFrame:
    """Classical local panel checks driven by Phase 4 cap/web schedules.

    The web calculation assumes no intermediate web stiffeners because their
    dimensions are not available in the public reference data.  Its result is
    deliberately conservative and identifies whether shell-level detail is
    required before a design can be accepted.
    """
    materials = load_materials()
    glass = get_ply("ELT5500_UD", materials)
    db = get_ply("Saertex_DB", materials)
    x = stations["r_m"].to_numpy()
    chord = stations["chord_m"].to_numpy()
    dcap = stations["cap_separation_m"].to_numpy()
    bcap = stations["cap_width_m"].to_numpy()
    ncap = stations["n_cap_glass"].to_numpy()
    nweb_nominal = stations["n_web_db_per_skin"].to_numpy()
    # Effective DB layers remain an integer count.  The support pitch is the
    # maximum unsupported panel length when ribs or sandwich support are used.
    nweb = np.where(nweb_nominal > 0,
                    np.ceil(nweb_nominal * web_ply_multiplier), 0.0)
    moment = stations["M_park_fault_design_Nm"].to_numpy()
    ei = stations["EI_flap_Nm2"].to_numpy()
    cases, _ = load_envelopes(x)
    shear = cases["U_park_fault"]["Q_design"]
    web_h = 0.65 * stations["tc"].to_numpy() * chord

    cap_t = ncap * glass["t_ply"]
    cap_strain = np.divide(moment * dcap / 2.0, CAP_EI_FRACTION * ei,
                           out=np.zeros_like(moment), where=cap_t > 0)
    cap_sigma = glass["E1"] * cap_strain
    cap_d = glass["E1"] * cap_t**3 / (12.0 * (1.0 - glass["nu12"]**2))
    cap_sigma_cr = np.divide(4.0 * np.pi**2 * cap_d, bcap**2 * cap_t,
                             out=np.full_like(moment, np.inf), where=cap_t > 0)

    web_t = nweb * db["t_ply"]
    web_tau = np.divide(0.5 * shear, 4.0 * web_h * web_t,
                        out=np.zeros_like(shear), where=web_t > 0)
    bay = np.maximum(np.gradient(x), 0.5)
    if web_support_spacing_m is not None:
        if web_support_spacing_m <= 0.0:
            raise ValueError("web_support_spacing_m must be positive")
        bay = np.minimum(bay, web_support_spacing_m)
    aspect = bay / np.maximum(web_h, 1e-6)
    ks = 5.34 + 4.0 / aspect**2
    web_d = db["E1"] * web_t**3 / (12.0 * (1.0 - db["nu12"]**2))
    web_tau_cr = np.divide(ks * np.pi**2 * web_d, web_h**2 * web_t,
                           out=np.full_like(shear, np.inf), where=web_t > 0)

    return pd.DataFrame({
        "r_m": x, "design": design,
        "web_ply_multiplier": web_ply_multiplier,
        "web_support_spacing_m": (np.nan if web_support_spacing_m is None
                                   else web_support_spacing_m),
        "cap_sigma_compression_Pa": cap_sigma,
        "cap_sigma_cr_Pa": cap_sigma_cr,
        "cap_buckling_RF": np.divide(cap_sigma_cr, cap_sigma,
                                      out=np.full_like(cap_sigma, np.inf), where=cap_sigma > 0),
        "web_tau_Pa": web_tau, "web_tau_cr_Pa": web_tau_cr,
        "web_shear_buckling_RF": np.divide(web_tau_cr, web_tau,
                                             out=np.full_like(web_tau, np.inf), where=web_tau > 0),
    })


def build_phase5_model(design_kind: str = "glass") -> tuple[BladeBeamFE, pd.DataFrame, dict]:
    design = build_design(design_kind)
    stiffness = pd.read_csv(AERO_DATA / "reference_blade" / "blade_stiffness.csv", comment="#")
    x = design.stations["r_m"].to_numpy()
    mass = np.interp(x, stiffness["r_m"], stiffness["mass_pm_kg_m"])
    ei_edge = np.interp(x, stiffness["r_m"], stiffness["EI_edge_Nm2"])
    model = BladeBeamFE(x, design.stations["EI_flap_Nm2"].to_numpy(), ei_edge, mass)
    cases, _ = load_envelopes(x)
    return model, design.stations, cases


def write_phase5_outputs() -> dict[str, object]:
    """Run static, modal and local panel-buckling Phase 5 analyses."""
    out = PHASE_ROOT / "results"
    out.mkdir(parents=True, exist_ok=True)
    model, stations, cases = build_phase5_model("glass")
    x = model.x
    static_results: list[StaticResult] = []
    for case_name in ("U_op", "U_park_feather", "U_park_fault"):
        q = distributed_load_from_shear(x, cases[case_name]["Q_design"],
                                        cases[case_name]["M_design"][0])
        static_results.append(model.static("flap", q, case_name))
    mass_df = load_blade_mass(str(AERO_DATA / "reference_blade" / "blade_stiffness.csv"))
    edge_q = GAMMA_F * 9.80665 * np.interp(x, mass_df["r_m"], mass_df["mass_pm_kg_m"])
    static_results.append(model.static("edge", edge_q, "U_gravity"))

    flap_modal = model.modal("flap")
    edge_modal = model.modal("edge")
    buckling = panel_buckling(stations, "glass")

    static_table = pd.DataFrame([{
        "case": r.case, "direction": r.direction, "root_shear_N": r.root_shear_N,
        "root_moment_Nm": r.root_moment_Nm, "tip_displacement_m": r.tip_displacement_m,
    } for r in static_results])
    modal_table = pd.DataFrame([
        {"direction": direction, "mode": i + 1, "frequency_hz": frequency}
        for direction, result in (("flap", flap_modal), ("edge", edge_modal))
        for i, frequency in enumerate(result.frequencies_hz)
    ])
    def write_csv_if_changed(table: pd.DataFrame, path: Path) -> None:
        """Avoid touching an unchanged result that a viewer may hold open."""
        rendered = table.to_csv(index=False)
        if path.exists() and path.read_text(encoding="utf-8") == rendered:
            return
        try:
            path.write_text(rendered, encoding="utf-8", newline="")
        except PermissionError:
            # Windows viewers/sync clients can hold a CSV open. Accept the
            # locked artifact only when it is numerically identical.
            existing = pd.read_csv(path)
            pd.testing.assert_frame_equal(
                existing.reset_index(drop=True),
                table.reset_index(drop=True),
                check_dtype=False,
                rtol=1e-12,
                atol=1e-12,
            )

    write_csv_if_changed(static_table, out / "phase5_static_cases.csv")
    write_csv_if_changed(modal_table, out / "phase5_modal.csv")
    write_csv_if_changed(buckling, out / "phase5_panel_buckling.csv")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from io import BytesIO

    def save_figure_if_changed(fig, path: Path, dpi: int = 170) -> None:
        """Render deterministically and avoid replacing an identical PNG."""
        buffer = BytesIO()
        fig.savefig(buffer, format="png", dpi=dpi)
        rendered = buffer.getvalue()
        if path.exists() and path.read_bytes() == rendered:
            return
        path.write_bytes(rendered)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for result in static_results:
        if result.direction == "flap":
            axes[0].plot(x, result.displacement_m, lw=2, label=result.case)
    axes[0].set(xlabel="Span [m]", ylabel="Flap deflection [m]", title="Phase 5 FE static response")
    axes[0].grid(alpha=.3); axes[0].legend()
    for i in range(3):
        axes[1].plot(x, flap_modal.mode_shapes[:, i], label=f"flap {i + 1}: {flap_modal.frequencies_hz[i]:.2f} Hz")
    axes[1].set(xlabel="Span [m]", ylabel="Normalised modal displacement", title="Flapwise FE modes")
    axes[1].grid(alpha=.3); axes[1].legend()
    fig.tight_layout(); save_figure_if_changed(fig, out / "phase5_static_modal.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    finite_web = np.where(np.isfinite(buckling["web_shear_buckling_RF"]), buckling["web_shear_buckling_RF"], np.nan)
    finite_cap = np.where(np.isfinite(buckling["cap_buckling_RF"]), buckling["cap_buckling_RF"], np.nan)
    ax.semilogy(x, finite_cap, label="spar-cap compression", lw=2)
    ax.semilogy(x, finite_web, label="web shear", lw=2)
    ax.axhline(1.0, color="red", ls=":", label="RF = 1")
    ax.set(xlabel="Span [m]", ylabel="Classical panel buckling RF", title="Local panel buckling screen")
    ax.grid(alpha=.3, which="both"); ax.legend()
    fig.tight_layout(); save_figure_if_changed(fig, out / "phase5_buckling.png"); plt.close(fig)

    # Main visual FEA artefact: factored parked-fault load shape and deformed
    # centreline.  Deflection is amplified only for legibility and labelled.
    parked = next(r for r in static_results if r.case == "U_park_fault")
    q_fault = distributed_load_from_shear(x, cases["U_park_fault"]["Q_design"],
                                          cases["U_park_fault"]["M_design"][0])
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.plot(x, np.zeros_like(x), color="#1c3d34", lw=5, label="undeformed beam axis")
    scale = 2.5
    ax.plot(x, -scale * parked.displacement_m, color="#e76f51", lw=2.8,
            label=f"FE deformed axis ×{scale:g}")
    sample = np.arange(0, len(x), 2)
    ax.quiver(x[sample], np.full_like(sample, 1.15, dtype=float),
              np.zeros_like(sample, dtype=float), -0.85 * q_fault[sample] / q_fault.max(),
              angles="xy", scale_units="xy", scale=1, color="#76c7e5", width=0.004,
              label="factored parked-fault load")
    ax.annotate(f"tip FE deflection = {abs(parked.tip_displacement_m):.2f} m",
                xy=(x[-1], -scale * parked.tip_displacement_m), xytext=(42, -9.5),
                arrowprops={"arrowstyle": "->", "color": "#e76f51"}, color="#9d2e20")
    ax.set(xlabel="Span [m]", ylabel="Illustrative displacement [m]",
           title="Phase 5 FE: parked-fault load and flapwise deformation")
    ax.grid(alpha=.25); ax.legend(loc="lower left")
    fig.tight_layout(); save_figure_if_changed(fig, out / "phase5_parked_fault_fea.png"); plt.close(fig)

    return {
        "model": model, "stations": stations, "static": static_results,
        "flap_modal": flap_modal, "edge_modal": edge_modal, "buckling": buckling,
    }


if __name__ == "__main__":
    result = write_phase5_outputs()
    for static in result["static"]:
        print(f"{static.case:16s} {static.direction:4s} tip={static.tip_displacement_m:.3f} m "
              f"Mroot={static.root_moment_Nm / 1e6:.3f} MN m")
    print("Flap frequencies [Hz]:", np.round(result["flap_modal"].frequencies_hz, 4))
    print("Minimum web buckling RF:", np.nanmin(result["buckling"]["web_shear_buckling_RF"].replace(np.inf, np.nan)))
