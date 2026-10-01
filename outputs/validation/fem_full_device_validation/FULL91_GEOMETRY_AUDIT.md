# Full 91-bulb geometry audit

Centers and horns are taken from the current production builders:
`Rod_Array_Hexagon_train` and `sixport_common.full_horns`.
The FEM helper `bulb_centers_mesh` is the same lattice translated onto the device domain.

- Bulb count: 91
- Maximum center difference versus the FEM helper: `3.553e-15`
- Nearest-neighbor spacing: `1.000000000000` to `1.000000000000` (pitch `1.0`)
- Lattice basis: `[[0, 1], [sqrt(3)/2, 1/2]]`, pointy-top triangular lattice, side length 6
- Material plasma radius: `0.22999999999999995` a = 4.6000 mm
- Unused metadata `r_plasma`: `0.25` a
- Quartz inner radius: `0.32499999999999996` a
- Quartz outer radius: `0.375` a
- Vacuum gap: `0.095` a between plasma and quartz
- Quartz permittivity: 3.8
- Domain: `30.0` by `28.0` a, PML `2.0` a
- Frequency: 3.85 GHz (`fs_a=0.25684435330257704`)
- Plasma frequency: 8.0 GHz, collision rate 1.0 MHz
- Production B: 0.05 T along +z, `f_c a/c = 0.09337289556608207`
- Ports: outward normals sorted by angle, index 0 nearest +x, then counterclockwise

PEC walls are the four filled horn prisms per port (`left_flare`, `right_flare`, `left_feed`, `right_feed`).
Natural Neumann on Hz is the PEC condition. Metal triangles are removed from the Helmholtz operator by setting rho to 0.

Per-bulb center differences are all at the maximum above. Ordering is the `train_elem_locs` order.

Geometry check: PASS
