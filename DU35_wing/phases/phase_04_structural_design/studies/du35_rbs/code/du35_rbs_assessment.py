"""Run the traceable RBS stability screening for the measured DU35 panel."""

from airfoil_panel_laminate import CAP_THICKNESS_M, PANEL_WIDTH_M, analyze_du35_upper_panel
from composite_physics.buckling_rbs import critical_uniaxial_compression_rbs
from composite_physics.panel_loads import compressive_resultant_from_stress
from composite_physics.reduced_stiffness import compute_reduced_D

# From the earlier beam-screening calculation: 109 MPa in the cap.  This is
# not a design load case; it is retained solely to keep the same load path
# traceable in the Day-2 presentation.
REFERENCE_CAP_STRESS_PA = 109e6
ASPECT_RATIOS_A_OVER_B = (0.5, 1.0, 2.0, 4.0)


def report_method_applicability(stiffness, d_star, tol=1e-9):
    """Print why the D-only plate formula is rejected and RBS accepted.

    Two separate conditions are checked and must not be confused:

    * ``B == 0``      -- required before a formula written in ``D`` alone may
      be used.  A single-sided skin offsets stiffness from the mid-plane, so
      this fails here.
    * ``D*16 = D*26 = 0`` -- required before the orthotropic simply-supported
      Navier solution may be used on the *reduced* panel.  This one holds.
    """
    b_max = float(abs(stiffness.B).max())
    d16, d26 = float(stiffness.D[0, 2]), float(stiffness.D[1, 2])
    ds16, ds26 = float(d_star[0, 2]), float(d_star[1, 2])
    symmetric = b_max < 1.0
    orthotropic = abs(ds16) <= tol and abs(ds26) <= tol

    print("\n--- method applicability ---")
    print(f"max |B|        = {b_max:,.1f} N")
    print(f"B11            = {stiffness.B[0, 0]:,.1f} N"
          f"   (mid-plane offset B11/A11 = {stiffness.B[0, 0] / stiffness.A[0, 0] * 1e3:.4f} mm)")
    print(f"D16,  D26      = {d16:.4g}, {d26:.4g} N m")
    print(f"D*16, D*26     = {ds16:.4g}, {ds26:.4g} N m")
    print(f"stack mirrored about mid-plane : {symmetric}")
    print(f"reduced panel orthotropic      : {orthotropic}")
    if symmetric:
        print("-> D-only plate formula is admissible for this stack.")
    else:
        print("-> REJECT D-only: the stack is extension-bending coupled, so a")
        print("   formula written in D alone does not apply to this panel.")
    if orthotropic:
        print("-> ACCEPT RBS: D* carries the coupling and leaves no bend-twist")
        print("   terms, so the orthotropic Navier solution applies to D*.")
    else:
        print("-> REJECT RBS: D*16/D*26 are not negligible; a fully coupled")
        print("   eigenvalue solution is required.")
    return symmetric, orthotropic


def run_assessment():
    panel = analyze_du35_upper_panel()
    stiffness = panel.stiffness
    d_star = compute_reduced_D(stiffness.A, stiffness.B, stiffness.D)
    applied = compressive_resultant_from_stress(REFERENCE_CAP_STRESS_PA, CAP_THICKNESS_M)
    results = [critical_uniaxial_compression_rbs(
        d_star, a_m=ratio * PANEL_WIDTH_M, b_m=PANEL_WIDTH_M,
        applied_Nx_N_per_m=applied.Nx_N_per_m,
    ) for ratio in ASPECT_RATIOS_A_OVER_B]
    return panel, d_star, applied, results


if __name__ == "__main__":
    panel, d_star, applied, results = run_assessment()
    print("DU35 upper cap — RBS buckling screening")
    print(f"b = {panel.surface_width_m:.6f} m; h_cap = {CAP_THICKNESS_M * 1e3:.1f} mm")
    print(f"Applied Nx = {applied.Nx_N_per_m / 1e6:.3f} MN/m ({applied.basis})")
    print("D* [kN m]:")
    print(d_star / 1e3)
    report_method_applicability(panel.stiffness, d_star)
    print("\na/b     a [m]    Nx,cr [MN/m]   m,n   MS")
    for result in results:
        print(f"{result.a_m / result.b_m:3.1f}  {result.a_m:8.4f}  "
              f"{result.Nx_cr_N_per_m / 1e6:12.3f}  ({result.m},{result.n})  "
              f"{result.margin_of_safety:+.3f}")
    print("\nScope: RBS approximation; simply supported edges; no certified spanwise support spacing.")
