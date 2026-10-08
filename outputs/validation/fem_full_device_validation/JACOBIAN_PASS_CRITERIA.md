# Jacobian pass criteria

Frozen before the finite-difference comparison in `jacobian_check.py`.

The design parameter is the dimensionless density scale `s` on one plasma disk:
`f_p^2(s) = s * f_p,ref^2`, with `f_p,ref = 8 GHz` and `s = 1` at the production density.
`B`, frequency, collision rate, mesh, and the source are fixed.

The direct sensitivity is `A dx/ds = -(dA/ds) x`, with `dε/ds` and `dρ/ds = -ρ (dε) ρ` analytic, and `dA/ds` the stiffness assembly of `dρ` (the mass term does not depend on `s`).

Centered differences use `[x(s+h) - x(s-h)] / (2h)` at

    h = 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5

Acceptance, for both `B = 0` and `B = +0.05 T`:

- the minimum relative error `||dx_fd - dx_direct|| / ||dx_direct||` over that set is `≤ 1e-4`
- some larger `h` has an error at least 10 times that minimum, so the curve is not a single accidental hit
- the same bound on one complex field sample and on `d/ds |Hz|^2` at that sample

A miss stays a miss. The threshold is not loosened after the sweep.
