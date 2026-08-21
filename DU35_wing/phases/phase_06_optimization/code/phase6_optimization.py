"""Restartable Phase 6 screening DoE for web-buckling closure.

The public data does not include rib/flange mass or adhesive geometry.  The
DoE therefore reports laminate mass separately and treats support density as a
second objective rather than inventing an unsupported support-mass number.
"""

from pathlib import Path
import sys
from io import BytesIO

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures

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

from fea_model import build_phase5_model, panel_buckling
from materials_db import get_ply, load_materials
from structural_model import build_design

BUCKLING_RF_TARGET = 1.5  # Screening design target; must be confirmed with company criteria.
PLY_MULTIPLIERS = (1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)
SUPPORT_SPACINGS_M = (0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)


def _write_csv_if_changed(table: pd.DataFrame, path: Path) -> None:
    """Skip an unchanged checkpoint; verify a locked copy before accepting it."""
    rendered = table.to_csv(index=False)
    if path.exists() and path.read_text(encoding="utf-8") == rendered:
        return
    try:
        path.write_text(rendered, encoding="utf-8", newline="")
    except PermissionError:
        pd.testing.assert_frame_equal(
            pd.read_csv(path).reset_index(drop=True),
            table.reset_index(drop=True),
            check_dtype=False,
            rtol=1e-12,
            atol=1e-12,
        )


def _save_figure_if_changed(fig, path: Path, dpi: int = 170) -> None:
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi)
    rendered = buffer.getvalue()
    if path.exists() and path.read_bytes() == rendered:
        return
    path.write_bytes(rendered)


def web_mass_delta_kg(stations: pd.DataFrame, multiplier: float) -> float:
    """Additional DB laminate mass; ribs/flanges are deliberately excluded."""
    db = get_ply("Saertex_DB", load_materials())
    x = stations.r_m.to_numpy()
    hweb = 0.65 * stations.tc.to_numpy() * stations.chord_m.to_numpy()
    nominal = stations.n_web_db_per_skin.to_numpy()
    resized = np.where(nominal > 0, np.ceil(nominal * multiplier), 0.0)
    return float(np.trapz(4.0 * hweb * (resized - nominal) * db["t_ply"] * db["rho"], x))


def run_web_support_doe() -> pd.DataFrame:
    """Evaluate the full discrete web-support design grid and checkpoint CSV."""
    _, stations, _ = build_phase5_model("glass")
    baseline = build_design("glass")
    rows = []
    for multiplier in PLY_MULTIPLIERS:
        for spacing in SUPPORT_SPACINGS_M:
            buckling = panel_buckling(stations, "glass", multiplier, spacing)
            web_rf = float(buckling.web_shear_buckling_RF.replace(np.inf, np.nan).min())
            added_mass = web_mass_delta_kg(stations, multiplier)
            rows.append({
                "web_ply_multiplier": multiplier,
                "web_support_spacing_m": spacing,
                "support_density_per_m": 1.0 / spacing,
                "min_web_buckling_RF": web_rf,
                "added_web_laminate_mass_kg": added_mass,
                "blade_mass_excluding_supports_kg": baseline.mass_kg + added_mass,
                "feasible_screening": web_rf >= BUCKLING_RF_TARGET,
            })
    result = pd.DataFrame(rows).sort_values(
        ["feasible_screening", "added_web_laminate_mass_kg", "support_density_per_m"],
        ascending=[False, True, True]).reset_index(drop=True)
    out = PHASE_ROOT / "results"
    out.mkdir(parents=True, exist_ok=True)
    _write_csv_if_changed(result, out / "phase6_web_support_doe.csv")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for multiplier, subset in result.groupby("web_ply_multiplier"):
        ax.plot(subset.web_support_spacing_m, subset.min_web_buckling_RF, marker="o",
                label=f"web ply ×{multiplier:g}")
    ax.axhline(BUCKLING_RF_TARGET, color="red", ls=":", label="RF target 1.5")
    ax.set(xlabel="Web support spacing [m]", ylabel="Minimum web buckling RF",
           title="Phase 6 DoE: web thickness and support spacing")
    ax.grid(alpha=.3); ax.legend(ncol=2)
    fig.tight_layout(); _save_figure_if_changed(fig, out / "phase6_web_support_doe.png"); plt.close(fig)
    return result


def validate_surrogate_and_select(doe: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Hold out DoE rows, validate a GP surrogate, and directly recheck choices."""
    features = doe[["web_ply_multiplier", "web_support_spacing_m"]].copy()
    # The plate equation has a strong inverse-square support-spacing term.
    # Exposing it to the surrogate is more data-efficient than a generic GP.
    features["inverse_spacing_squared"] = 1.0 / features.web_support_spacing_m**2
    features = features.to_numpy()
    target = np.log(doe["min_web_buckling_RF"].to_numpy())
    train, test = train_test_split(np.arange(len(doe)), test_size=0.25, random_state=17)
    surrogate = make_pipeline(PolynomialFeatures(degree=2), Ridge(alpha=1e-6))
    surrogate.fit(features[train], target[train])
    predicted = np.exp(surrogate.predict(features[test]))
    validation = doe.iloc[test][["web_ply_multiplier", "web_support_spacing_m", "min_web_buckling_RF"]].copy()
    validation["surrogate_web_buckling_RF"] = predicted
    validation["relative_error_pct"] = 100.0 * (predicted / validation.min_web_buckling_RF - 1.0)

    feasible = doe.loc[doe.feasible_screening].copy()
    mass_min, mass_max = feasible.added_web_laminate_mass_kg.min(), feasible.added_web_laminate_mass_kg.max()
    density_min, density_max = feasible.support_density_per_m.min(), feasible.support_density_per_m.max()
    feasible["balanced_score"] = ((feasible.added_web_laminate_mass_kg - mass_min) / (mass_max - mass_min)
                                  + (feasible.support_density_per_m - density_min) / (density_max - density_min))
    chosen = pd.concat([
        feasible.nsmallest(1, "added_web_laminate_mass_kg").assign(selection="lowest laminate mass"),
        feasible.nsmallest(1, "support_density_per_m").assign(selection="lowest support density"),
        feasible.nsmallest(1, "balanced_score").assign(selection="balanced screening trade"),
    ]).drop_duplicates(subset=["web_ply_multiplier", "web_support_spacing_m"])

    _, stations, _ = build_phase5_model("glass")
    direct = []
    for row in chosen.itertuples(index=False):
        checked = panel_buckling(stations, "glass", row.web_ply_multiplier,
                                 row.web_support_spacing_m)
        direct.append(float(checked.web_shear_buckling_RF.replace(np.inf, np.nan).min()))
    chosen["direct_FE_web_buckling_RF"] = direct
    chosen["direct_FE_pass"] = chosen.direct_FE_web_buckling_RF >= BUCKLING_RF_TARGET

    out = PHASE_ROOT / "results"
    _write_csv_if_changed(validation, out / "phase6_surrogate_validation.csv")
    _write_csv_if_changed(chosen, out / "phase6_selected_candidates.csv")
    return validation, chosen


if __name__ == "__main__":
    doe = run_web_support_doe()
    validation, chosen = validate_surrogate_and_select(doe)
    print(chosen.to_string(index=False))
    print(f"Max held-out surrogate error: {validation.relative_error_pct.abs().max():.1f}%")
