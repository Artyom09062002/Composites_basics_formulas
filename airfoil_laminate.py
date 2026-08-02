"""ABD stiffness of the DU35 upper panel — the equivalent two-zone stack.

Run from the folder containing lamina_mechanics.py, transformations.py and
laminate.py:

    python du35_panel.py

Geometry is the measured upper panel between the two shear web centrelines.
The material stack is the engineering approximation used in the report:
a 3 mm triaxial skin and a 65 mm unidirectional cap, each represented as one
equivalent 0-degree orthotropic zone.  This is not the manufacturing layup.

>>> REPLACE THE TWO MATERIAL BLOCKS BELOW WITH YOUR OWN NUMBERS <<<
The values here are representative E-glass/epoxy properties, not measured
ones, and they will not reproduce the report's B11 unless they match what
you used.
"""

import numpy as np

from laminate import assemble_laminate_stiffness

np.set_printoptions(precision=4, suppress=False, linewidth=140)

# ---------------------------------------------------------------------------
# Measured panel geometry (FreeCAD, midsurface arc between web centrelines)
# ---------------------------------------------------------------------------
PANEL_WIDTH_B = 0.8108        # m — width b for the plate idealisation
CHORD = 4.492                 # m — section chord, for context only

# ---------------------------------------------------------------------------
# Materials — REPLACE THESE
# index 0 = triaxial skin, index 1 = unidirectional cap
# ---------------------------------------------------------------------------
TRIAX_SKIN = {
    "E1": 27.7e9,    # Pa
    "E2": 13.7e9,    # Pa
    "G12": 7.20e9,   # Pa
    "v12": 0.39,     # -
}

UD_CAP = {
    "E1": 41.8e9,    # Pa
    "E2": 14.0e9,    # Pa
    "G12": 2.63e9,   # Pa
    "v12": 0.28,     # -
}

MATERIALS = [TRIAX_SKIN, UD_CAP]

# ---------------------------------------------------------------------------
# The stack, listed from the -h/2 face to the +h/2 face.
# The skin sits on one side only, which is why B will not vanish.
# ---------------------------------------------------------------------------
SKIN_THICKNESS = 3.0e-3       # m
CAP_THICKNESS = 65.0e-3       # m

layup = [
    {"theta": 0.0, "t": SKIN_THICKNESS, "mat": 0},   # triax skin
    {"theta": 0.0, "t": CAP_THICKNESS, "mat": 1},    # UD cap
]

# ---------------------------------------------------------------------------
result = assemble_laminate_stiffness(layup, MATERIALS)

total_thickness = result.z[-1] - result.z[0]

print("DU35 UPPER PANEL — equivalent two-zone stack")
print(f"  panel width b      : {PANEL_WIDTH_B:.4f} m  (measured midsurface arc)")
print(f"  section chord      : {CHORD:.3f} m")
print(f"  stack thickness h  : {total_thickness * 1e3:.1f} mm"
      f"   ({SKIN_THICKNESS * 1e3:.0f} mm skin + {CAP_THICKNESS * 1e3:.0f} mm cap)")
print(f"  zone interfaces    : {result.z * 1e3} mm\n")

print("A  [N/m]    in-plane stiffness")
print(result.A, "\n")

print("B  [N]      extension-bending coupling")
print(result.B, "\n")

print("D  [N*m]    bending stiffness")
print(result.D, "\n")

print("ABD [6x6]   A block in N/m, B block in N, D block in N*m")
print(result.ABD, "\n")

# ---------------------------------------------------------------------------
# What the numbers mean for the next step
# ---------------------------------------------------------------------------
b_max = abs(result.B).max()
d16, d26 = result.D[0, 2], result.D[1, 2]

print("--- applicability of a D-only plate-buckling formula ---")
print(f"  max |B|      = {b_max:,.1f} N")
print(f"  B11          = {result.B[0, 0]:,.1f} N")
print(f"  D16, D26     = {d16:.4g}, {d26:.4g} N*m")
print()
print(f"  B is zero      : {b_max < 1.0}"
      "   (single-sided skin offsets material from the midplane)")
print(f"  D16=D26 is zero: {abs(d16) < 1e-9 and abs(d26) < 1e-9}")
print()
if b_max >= 1.0:
    print("  -> REJECT: the stack is extension-bending coupled, so a formula")
    print("     written in D alone does not apply to this panel.")
else:
    print("  -> B vanishes; check D16 and D26 before accepting a D-only formula.")
