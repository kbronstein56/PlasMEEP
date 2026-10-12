"""Production feed geometry without Meep/PMMI imports (FEM multimode ports)."""
from __future__ import annotations

C_MPS = 2.99792458e8
A_M = 0.020  # Meep unit = 20 mm pitch

fs_Hz = 3.85e9
a = A_M
wall_thickness = 0.004 / a
width_base = 0.048 / a
clear_width = width_base - 2 * wall_thickness  # 2a in lattice units
feed_length = 0.060 / a

# Nondimensional frequency f*a/c (matches sixport_common.fs_a)
fs_a = fs_Hz * a / C_MPS
