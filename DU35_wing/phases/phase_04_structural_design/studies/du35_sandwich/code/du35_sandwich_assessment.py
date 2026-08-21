"""Reproducible screening assessment for the two upper DU35 sandwich panels.

The unknown transverse-support pitch is retained as an a/b sweep.  Two load
paths are deliberately reported side-by-side: an existing 109 MPa cap-stress
reference and a load mapped from the existing beam-model moment/curvature.
"""

from pathlib import Path
import csv
import sys

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parents[1]
PHASE4_ROOT = Path(__file__).resolve().parents[3]
REPO_ROOT = Path(__file__).resolve().parents[6]
RBS_CODE = PHASE4_ROOT / "studies" / "du35_rbs" / "code"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(RBS_CODE))

from airfoil_panel_laminate import (
    CAP_THICKNESS_M, PANEL_WIDTH_M, UD_GFRP, analyze_du35_upper_panel,
)
from composite_physics.buckling_rbs import critical_uniaxial_compression_rbs
from du35_sandwich_geometry import (
    CORE_THICKNESS_M, FACE_THICKNESS_M, LE_UPPER, TE_UPPER,
    beam_panel_strain, load_inner_beam_snapshot, normalized_upper_position,
)
from composite_physics.reduced_stiffness import compute_reduced_D
from composite_physics.sandwich_panel import SandwichMaterial, SandwichPanel, assess_sandwich_panel

REFERENCE_CAP_STRESS_PA = 109e6
ASPECT_RATIOS = np.linspace(0.5, 4.0, 141)

TRIAX_FACE = {"E1": 27.7e9, "E2": 13.65e9, "G12": 7.20e9, "v12": 0.39}
H100_CORE = {"E1": 0.130e9, "E2": 0.130e9, "G12": 0.049e9, "v12": 0.32}
MATERIAL = SandwichMaterial(
    face=TRIAX_FACE,
    core=H100_CORE,
    face_compression_Pa=350e6,
    core_compression_modulus_Pa=0.130e9,
    core_shear_modulus_Pa=0.049e9,
)


def build_panels():
    return {
        geometry.name: SandwichPanel(
            geometry.name, geometry.midsurface_width_m,
            FACE_THICKNESS_M, CORE_THICKNESS_M, MATERIAL,
        )
        for geometry in (LE_UPPER, TE_UPPER)
    }


def load_cases():
    beam = load_inner_beam_snapshot()
    conservative_strain = REFERENCE_CAP_STRESS_PA / UD_GFRP["E1"]
    cases = {}
    for geometry in (LE_UPPER, TE_UPPER):
        cases[geometry.name] = {
            "109MPa_reference": conservative_strain,
            "beam_mapped": beam_panel_strain(geometry, beam),
        }
    return beam, cases


def _limit_map(result):
    return {limit.mode: limit.critical_Nx_N_per_m for limit in result.limits}


def run_assessment(output_root: Path | None = None):
    output_root = output_root or STUDY_ROOT / "results"
    docs = STUDY_ROOT / "docs"
    output_root.mkdir(parents=True, exist_ok=True)
    docs.mkdir(parents=True, exist_ok=True)

    panels = build_panels()
    beam, strains = load_cases()
    geometry_by_name = {item.name: item for item in (LE_UPPER, TE_UPPER)}
    sweep_rows, comparison_rows = [], []

    for panel_name, panel in panels.items():
        geometry = geometry_by_name[panel_name]
        for case_name, strain in strains[panel_name].items():
            applied = panel.membrane_resultant(strain)
            results = [assess_sandwich_panel(panel, ratio, applied)
                       for ratio in ASPECT_RATIOS]
            for result in results:
                limits = _limit_map(result)
                sweep_rows.append({
                    "panel": panel_name, "load_case": case_name,
                    "a_over_b": result.a_over_b, "a_m": result.a_m,
                    "applied_Nx_N_per_m": applied,
                    "global_CLPT_Nx_N_per_m": result.global_clpt_Nx_N_per_m,
                    "global_shear_corrected_Nx_N_per_m": result.global_shear_corrected_Nx_N_per_m,
                    "face_wrinkling_Nx_N_per_m": limits["face_wrinkling"],
                    "core_shear_crimping_Nx_N_per_m": limits["core_shear_crimping"],
                    "face_compression_Nx_N_per_m": limits["face_compression"],
                    "governing_mode": result.governing_mode,
                    "reserve_factor": result.reserve_factor,
                    "margin_of_safety": result.margin_of_safety,
                })
            worst = min(results, key=lambda item: item.reserve_factor)
            comparison_rows.append({
                "panel": panel_name, "load_case": case_name,
                "width_m": panel.width_m,
                "x_over_c_zone": f"{geometry.x_over_c_start:.2f}-{geometry.x_over_c_end:.2f}",
                "axial_strain": strain, "face_stress_Pa": panel.face_stress(strain),
                "applied_Nx_N_per_m": applied,
                "worst_a_over_b_in_sweep": worst.a_over_b,
                "governing_mode": worst.governing_mode,
                "critical_Nx_N_per_m": worst.governing_critical_Nx_N_per_m,
                "reserve_factor": worst.reserve_factor,
                "margin_of_safety": worst.margin_of_safety,
            })

    # Existing spar-cap RBS result, evaluated on the same a/b sweep.
    cap = analyze_du35_upper_panel()
    d_star = compute_reduced_D(cap.stiffness.A, cap.stiffness.B, cap.stiffness.D)
    cap_strains = {
        "109MPa_reference": REFERENCE_CAP_STRESS_PA / UD_GFRP["E1"],
        "beam_mapped": beam.cap_strain,
    }
    for case_name, strain in cap_strains.items():
        applied = UD_GFRP["E1"] * CAP_THICKNESS_M * strain
        cap_results = [critical_uniaxial_compression_rbs(
            d_star, ratio * PANEL_WIDTH_M, PANEL_WIDTH_M, applied
        ) for ratio in ASPECT_RATIOS]
        worst = min(cap_results, key=lambda item: item.margin_of_safety)
        comparison_rows.append({
            "panel": "upper_spar_cap", "load_case": case_name,
            "width_m": PANEL_WIDTH_M, "x_over_c_zone": "between_webs",
            "axial_strain": strain, "face_stress_Pa": UD_GFRP["E1"] * strain,
            "applied_Nx_N_per_m": applied,
            "worst_a_over_b_in_sweep": worst.a_m / worst.b_m,
            "governing_mode": "global_buckling_RBS",
            "critical_Nx_N_per_m": worst.Nx_cr_N_per_m,
            "reserve_factor": worst.Nx_cr_N_per_m / applied,
            "margin_of_safety": worst.margin_of_safety,
        })

    geometry_rows = [{
        "panel": item.name, "chord_m": 4.491873,
        "x_over_c_start": item.x_over_c_start, "x_over_c_end": item.x_over_c_end,
        "midsurface_width_m": item.midsurface_width_m,
        "centroid_y_m": item.centroid_y_m,
        "normalized_upper_position": normalized_upper_position(item),
        "face_thickness_m": FACE_THICKNESS_M, "core_thickness_m": CORE_THICKNESS_M,
    } for item in (LE_UPPER, TE_UPPER)]
    load_rows = [{
        "beam_station_r_m": beam.r_m, "chord_m": beam.chord_m,
        "moment_Nm": beam.moment_Nm, "EI_Nm2": beam.EI_Nm2,
        "cap_separation_m": beam.cap_separation_m,
        "curvature_per_m": beam.curvature_per_m, "cap_strain": beam.cap_strain,
        "reference_cap_stress_Pa": REFERENCE_CAP_STRESS_PA,
    }]

    _write_csv(output_root / "du35_sandwich_geometry.csv", geometry_rows)
    _write_csv(output_root / "du35_load_cases.csv", load_rows)
    _write_csv(output_root / "du35_sandwich_sweep.csv", sweep_rows)
    _write_csv(output_root / "du35_panel_comparison.csv", comparison_rows)
    _plot_span_sweep(docs / "du35_sandwich_span_sweep.svg", sweep_rows)
    _plot_comparison(docs / "du35_panel_comparison.svg", comparison_rows)
    return beam, comparison_rows


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot_span_sweep(path, rows):
    width, height = 1120, 480
    colors = ("#d97706", "#2563eb")
    panels = ("LE_upper_sandwich", "TE_upper_foam")
    maxima = []
    for panel in panels:
        data = [r for r in rows if r["panel"] == panel]
        maxima.extend(r["global_shear_corrected_Nx_N_per_m"] / 1e6 for r in data)
        maxima.extend(r["applied_Nx_N_per_m"] / 1e6 for r in data)
    ymax = max(maxima) * 1.08
    svg = [_svg_header(width, height),
           '<text x="560" y="28" text-anchor="middle" class="title">DU35 upper sandwich panels — support-spacing sweep</text>']
    for index, (panel, color) in enumerate(zip(panels, colors)):
        left, top, plot_w, plot_h = 70 + index * 535, 65, 455, 330
        base = [r for r in rows if r["panel"] == panel and r["load_case"] == "109MPa_reference"]
        points = []
        for row in base:
            px = left + (row["a_over_b"] - .5) / 3.5 * plot_w
            py = top + plot_h - row["global_shear_corrected_Nx_N_per_m"] / 1e6 / ymax * plot_h
            points.append(f"{px:.1f},{py:.1f}")
        svg += _axes(left, top, plot_w, plot_h, ymax)
        svg.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="3"/>')
        for case, dash in (("109MPa_reference", "9,6"), ("beam_mapped", "2,5")):
            row = next(r for r in rows if r["panel"] == panel and r["load_case"] == case)
            y = top + plot_h - row["applied_Nx_N_per_m"] / 1e6 / ymax * plot_h
            svg.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#374151" stroke-width="2" stroke-dasharray="{dash}"/>')
        svg.append(f'<text x="{left+plot_w/2}" y="{top+plot_h+52}" text-anchor="middle">{panel.replace("_", " ")}</text>')
    svg += ['<text x="560" y="463" text-anchor="middle">Assumed support aspect ratio a/b = 0.5…4.0</text>',
            '<text x="15" y="230" transform="rotate(-90 15 230)" text-anchor="middle">Nx [MN/m]</text>',
            '<line x1="735" y1="425" x2="765" y2="425" stroke="#374151" stroke-width="2" stroke-dasharray="9,6"/><text x="772" y="430" font-size="12">109 MPa reference</text>',
            '<line x1="915" y1="425" x2="945" y2="425" stroke="#374151" stroke-width="2" stroke-dasharray="2,5"/><text x="952" y="430" font-size="12">beam mapped</text>', '</svg>']
    path.write_text("\n".join(svg), encoding="utf-8")


def _plot_comparison(path, rows):
    panels = ("LE_upper_sandwich", "TE_upper_foam", "upper_spar_cap")
    cases = ("109MPa_reference", "beam_mapped")
    values = {(p, c): next(r["reserve_factor"] for r in rows
                           if r["panel"] == p and r["load_case"] == c)
              for p in panels for c in cases}
    ymax = max(values.values()) * 1.18
    width, height = 950, 480
    left, top, plot_w, plot_h = 75, 55, 820, 345
    svg = [_svg_header(width, height),
           '<text x="475" y="28" text-anchor="middle" class="title">DU35 upper panels — governing screening result</text>']
    svg += _axes(left, top, plot_w, plot_h, ymax)
    colors = {cases[0]: "#f59e0b", cases[1]: "#2563eb"}
    group_w, bar_w = plot_w / 3, 72
    for i, panel in enumerate(panels):
        centre = left + (i + .5) * group_w
        for j, case in enumerate(cases):
            value = values[(panel, case)]
            x = centre + (-.55 if j == 0 else .55) * bar_w - bar_w / 2
            y = top + plot_h - value / ymax * plot_h
            h = top + plot_h - y
            svg.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w}" height="{h:.1f}" fill="{colors[case]}"/>')
            svg.append(f'<text x="{x+bar_w/2:.1f}" y="{max(top+14,y-6):.1f}" text-anchor="middle" font-size="12">{value:.2f}</text>')
        svg.append(f'<text x="{centre:.1f}" y="{top+plot_h+28}" text-anchor="middle">{("LE sandwich", "TE foam", "Spar cap")[i]}</text>')
    y_rf = top + plot_h - 1 / ymax * plot_h
    svg += [f'<line x1="{left}" y1="{y_rf:.1f}" x2="{left+plot_w}" y2="{y_rf:.1f}" stroke="#dc2626" stroke-width="2"/>',
            '<rect x="590" y="430" width="14" height="14" fill="#f59e0b"/><text x="612" y="442" font-size="12">109 MPa reference</text>',
            '<rect x="760" y="430" width="14" height="14" fill="#2563eb"/><text x="782" y="442" font-size="12">beam mapped</text>',
            '<text x="18" y="235" transform="rotate(-90 18 235)" text-anchor="middle">Minimum reserve factor</text>', '</svg>']
    path.write_text("\n".join(svg), encoding="utf-8")


def _svg_header(width, height):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
            '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#111827;font-size:14px}.title{font-size:19px;font-weight:600}</style>'
            '<rect width="100%" height="100%" fill="white"/>')


def _axes(left, top, width, height, ymax):
    items = [f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+height}" stroke="#111827"/>',
             f'<line x1="{left}" y1="{top+height}" x2="{left+width}" y2="{top+height}" stroke="#111827"/>']
    for i in range(5):
        y = top + height - i / 4 * height
        value = i / 4 * ymax
        items += [f'<line x1="{left}" y1="{y:.1f}" x2="{left+width}" y2="{y:.1f}" stroke="#d1d5db"/>',
                  f'<text x="{left-7}" y="{y+5:.1f}" text-anchor="end" font-size="11">{value:.1f}</text>']
    return items


if __name__ == "__main__":
    beam, rows = run_assessment()
    print(f"Beam snapshot: r={beam.r_m:.3f} m, M={beam.moment_Nm/1e6:.3f} MNm, cap strain={beam.cap_strain:.6g}")
    for row in rows:
        print(f"{row['panel']:22s} {row['load_case']:18s} "
              f"Nx={row['applied_Nx_N_per_m']/1e6:6.3f} MN/m  "
              f"RF={row['reserve_factor']:6.3f}  {row['governing_mode']}")
