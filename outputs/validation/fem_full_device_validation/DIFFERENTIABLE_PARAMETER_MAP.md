# Differentiable parameter map

Started only after `FULL91_GYROTROPIC_FORWARD_VALIDATED`.

The first differentiable implementation keeps the mesh and geometry fixed.
Design variables change **material parameters only**. No shape derivatives and
no remeshing.

## Variable

For rod `i`,

    p_i = s_i

where `s_i` is a dimensionless density scale:

    f_p,i^2 = s_i * f_p,ref^2
    f_p,ref = 8.00 GHz   (production)
    s_i = 1                at the production density
    allowed range          s_i ∈ [0, s_max] with s_max set by the optimizer
                           (physical plasma frequency must stay below the
                           solver’s validated regime; s_max = 2 is the
                           initial working bound for tests)

Units: dimensionless. Mapping to electron density:

    n_e,i / n_e,ref = s_i
    because ω_p^2 ∝ n_e.

Collision rate `γ`, lattice geometry, quartz, B, and frequency are held fixed
in the first Jacobian.

## Tensor path

    s_i
      → f_p,i = f_p,ref * sqrt(s_i)
      → ε_i(ω, B) from the validated Lorentz tensor
      → ρ_i = ε_i^{-1}
      → element coefficients on the triangles of rod i
      → A(p)

Analytic derivatives at fixed `ω`, `γ`, `f_c`:

    ∂ε_⊥/∂s = − f_p,ref^2 * u / (f * (u^2 − f_c^2))
    ∂η/∂s   =   f_p,ref^2 * f_c / (f * (u^2 − f_c^2))
    ∂ρ/∂s   = − ρ (∂ε/∂s) ρ

The Hz mass term does not depend on `s`. `dA/ds` is the stiffness assembly of
`dρ` only (`mass_scale = 0` in `assemble_anisotropic`).

## Elements affected

Only plasma triangles whose centroids lie inside rod `i`. Quartz, vacuum,
PEC walls, and other rods are untouched by `∂/∂s_i`.

## Initial test set

1. one plasma disk in a small box
2. two or three disks
3. a reduced multiport only after (1)–(2) pass the frozen Jacobian criteria

The full 91-variable Jacobian is not built until the small-device finite-
difference agreement passes.
