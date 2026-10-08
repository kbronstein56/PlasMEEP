# Differentiability derivation audit

Independent check of the prepared material-only forward sensitivity path.
Mesh and geometry are fixed. No shape derivatives. No remeshing.

Time convention throughout: `e^{-iωt}` (matches `gyrotropic_tensor.py` and
`assemble_anisotropic`).

## 1. Design parameter

For rod `i`,

    p_i ≡ s_i ∈ [0, s_max]

with

    f_{p,i}^2 = s_i · f_{p,ref}^2
    f_{p,ref} = 8.00 GHz   (production)
    s_i = 1                at the production density
    n_{e,i} / n_{e,ref} = s_i

`γ`, lattice geometry, quartz, `B`, and frequency are held fixed for the
first Jacobian. Only plasma triangles whose centroids lie inside rod `i`
feel `∂/∂s_i`.

**Verdict:** definition matches `DIFFERENTIABLE_PARAMETER_MAP.md` and
`plasma_sensitivity.py`.

## 2. Constitutive tensor

Cold-plasma Lorentz model (`tensor_ordinary`), ordinary frequencies:

    u = f + i γ
    D = u² − f_c²
    ε_⊥ = 1 − (f_p²) u / (f D)
    η   = (f_p²) f_c / (f D)
    ε_xy = +i η
    ε_yx = −i η
    ε_yy = ε_⊥

With `f_p² = s f_{p,ref}²`,

    ∂ε_⊥/∂s = − f_{p,ref}² · u / (f D)
    ∂η/∂s   =   f_{p,ref}² · f_c / (f D)
    ∂ε_xy/∂s = +i ∂η/∂s
    ∂ε_yx/∂s = −i ∂η/∂s

These are independent of the instantaneous `f_p` except through the chain
`∂(f_p²)/∂s = f_{p,ref}²`. The “1” in `ε_⊥` drops out under `∂/∂s`.

**Verdict:** `deps_ds` matches the analytic derivative of `tensor_ordinary`.
Off-diagonal signs (`+iη`, `−iη`) match the validated e^{-iωt} convention.

## 3. Resistivity derivative

    ρ = ε^{-1}
    dρ = − ρ (dε) ρ

Component form used in assembly: `dρ_xx, dρ_xy, dρ_yx, dρ_yy`.

**Verdict:** `drho_ds` implements the matrix product correctly. This must
also be checked by centered FD of `ρ(s)` before blaming FEM assembly
(Phase 3).

## 4. Discrete operator

Weak form (`assemble_anisotropic`), stretched anisotropic Helmholtz for `Hz`:

    A = K(ρ) − M(k₀²)

`K` depends on `ρ`; the mass `M` does **not** depend on `s` (plasma
parameter enters only through in-plane `ρ`). Therefore

    dA/ds = K(dρ)     with mass_scale = 0

Empty-row pinning belongs to `A`, not to `dA`. Sensitivity assembly must
use `pin_empty=False` so empty rows of `dA` stay zero; pinned rows of `A`
then give `dx_i = 0` on unused DOFs.

**Verdict:** `jacobian_check.py` assembles `dA` with `mass_scale=0` and
`pin_empty=False`. Correct.

## 5. Forward sensitivity equation

Primal:

    A(p) x(p) = b(p)

Differentiate:

    (dA/dp_k) x + A (dx/dp_k) = db/dp_k

Material-only, fixed nodal source (interpolation weights independent of `s`):

    db/dp_k = 0

hence

    A (dx/dp_k) = − (dA/dp_k) x

**RHS sign:** minus. Implementation: `lu.solve(-(dA @ x))`. Correct.

**Factorization reuse:** one `splu(A)` serves the primal solve and every
sensitivity RHS at fixed `(ω, B, mesh, s_base)`. Do not refactor `A` once
per design variable.

## 6. Checklist answers

| # | Question | Answer |
|---|----------|--------|
| 1 | definition of `s_i` | dimensionless density scale; `f_p² = s f_p,ref²` |
| 2 | `f_p,i² = s_i f_p,ref²` | yes |
| 3 | `dε/ds` | analytic Lorentz derivatives above |
| 4 | `dρ/ds` | `−ρ (dε) ρ` |
| 5 | `dA/ds` | stiffness of `dρ` only |
| 6 | sensitivity RHS sign | `−(dA)x` |
| 7 | time convention | `e^{-iωt}` |
| 8 | tensor off-diagonal signs | `ε_xy = +iη`, `ε_yx = −iη` |
| 9 | source `b` independent of `s` | yes (fixed nodal weights) |
| 10 | geometry/mesh fixed | yes (material-only) |

## 7. Scope of the first implementation

MATERIAL-ONLY differentiation only.

- No remeshing.
- No shape / boundary motion derivatives.
- No differentiation of PML stretch factors (they are geometric).
- No differentiation of quartz or vacuum.

## 8. Audit conclusion

The prepared map, criteria, and `plasma_sensitivity.py` /
`jacobian_check.py` path are mathematically consistent with the discrete
FEM. Numerical pass/fail is decided only by the frozen FD criteria in
`JACOBIAN_PASS_CRITERIA.md`, not by this algebra.
