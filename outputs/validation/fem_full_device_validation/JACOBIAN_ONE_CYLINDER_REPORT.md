# One-cylinder Jacobian validation

Frozen criteria: `JACOBIAN_PASS_CRITERIA.md` — minimum relative error ≤ 1e-4
with truncation visible (≥10× the minimum at larger `h`).

Script: `scripts/validation/fem_full_device_validation/jacobian_one_cylinder.py`
Data: `jacobian_one_cylinder.json`, `jacobian_fd_comparison.json`

## Constitutive FD (before FEM)

Centered FD of `ε(s)` and `ρ(s)` vs analytic `deps_ds` / `drho_ds`:

| B (T) | min ‖dε_fd − dε‖/‖dε‖ | min ‖dρ_fd − dρ‖/‖dρ‖ | pass |
|-------|------------------------|------------------------|------|
| 0     | 1.2e-15                | 1.5e-10                | yes  |
| +0.05 | 2.7e-15                | 2.1e-10                | yes  |
| −0.05 | 2.7e-15                | 2.1e-10                | yes  |

## FEM forward sensitivity

Direct: `A dx/ds = −(dA/ds) x` with one LU of `A` reused for the sensitivity solve.
Mesh: 4996 nodes, 1098 plasma triangles. Source fixed in air; probe beyond the disk.

| B (T) | min field rel | min sample rel | min power rel | truncation | pass |
|-------|---------------|----------------|---------------|------------|------|
| 0     | 5.6e-10       | 1.8e-10        | 1.1e-09       | yes        | yes  |
| +0.05 | 6.7e-10       | 4.0e-10        | 3.3e-08       | yes        | yes  |
| −0.05 | 7.7e-10       | 5.5e-10        | 3.3e-10       | yes        | yes  |

Sweep shape (all three `B`): large `h` → truncation; intermediate `h` → plateau
descending toward machine precision; no loose threshold.

Real and imaginary parts of the complex probe sample agree separately at the
same order. Magnitude, phase, and `|Hz|²` derivatives also meet ≤ 1e-4 at the
best `h`.

## Verdict

**ONE_CYLINDER_JACOBIAN: PASS** for `B ∈ {0, +0.05, −0.05}`.
