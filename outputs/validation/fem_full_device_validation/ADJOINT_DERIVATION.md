# Discrete adjoint for a real scalar objective

Derived only after the forward Jacobian path is in place. Material-only
parameters `p_k = s_k`. Fixed mesh and fixed source `b`.

Time convention: `e^{-iωt}`. Complex fields `x ∈ ℂⁿ`.

## Objective

Let `y = L x` be a complex linear monitor (port amplitude), and take the real
objective

    J(x) = |y|² = y* y = (L x)* (L x)

More generally, any real scalar of the form

    J = j(x, x*)

with `J` real-valued. Wirtinger derivatives treat `x` and `x*` as independent:

    dJ = (∂J/∂x) dx + (∂J/∂x*) dx*

For `J = |L x|²`,

    ∂J/∂x  = (L x)* L = y* L
    ∂J/∂x* = (L x) L*   (row-vector sense)

so the differential is

    dJ = 2 Re( y* L dx )

## Forward constraint

    A(p) x(p) = b          (db = 0 for material-only, fixed source)

    A dx = − (dA) x

## Forward gradient (for comparison)

    dJ/dp_k = 2 Re( y* L (dx/dp_k) )
            = 2 Re( y* L A⁻¹ (−∂A/∂p_k x) )

## Adjoint variable

Define the complex adjoint row (or column `λ`) so that

    λ† A = y* L

i.e.

    A† λ = L† y

where `†` is the conjugate transpose (`A^H`). Then

    y* L dx = λ† A dx = λ† (− ∂A/∂p_k x)
            = − λ† (∂A/∂p_k) x

and

    dJ/dp_k = 2 Re( − λ† (∂A/∂p_k) x )
            = −2 Re( λ† (∂A/∂p_k) x )

### Why `A†` and not `Aᵀ`

`J` is a real function of complex `x`. The chain rule that isolates a single
linear solve for all parameters contracts `y* L` against `dx`. Matching

    (y* L) A⁻¹ = λ†

requires `λ† = (y* L) A⁻¹`, hence `A† λ = L† y`.

If the weak form produced a complex-symmetric (not Hermitian) `A`, one could
write a transpose system with a different pairing; our anisotropic plasma
operator with PML stretch is **not** Hermitian in general, and the correct
discrete adjoint of the complex bilinear form for this real objective is the
**conjugate-transpose** system above.

For a holomorphic objective that depends on `x` but not `x*` (rare for power
objectives), a transpose adjoint can appear. `|y|²` depends on both, so the
gradient uses `2 Re(·)` and `A†`.

## Gradient assembly

With `∂A/∂s_k = K(∂ρ/∂s_k)` (stiffness only),

    ∂J/∂s_k = −2 Re( λ† K(∂ρ/∂s_k) x )

Element-wise this is a local contraction of the adjoint and forward fields
against `∂ρ/∂s_k` on rod `k` only — no need to form `K` globally if an
element residual is available.

## Algorithm

1. Solve `A x = b` (forward).
2. Form `y = L x` and RHS `L† y`.
3. Solve `A† λ = L† y` (one adjoint solve per real scalar objective).
4. For each `k`, evaluate `−2 Re(λ† (∂A/∂s_k) x)`.

Cost: one forward factorization (or factor `A` and solve `A†` via the
conjugate-transpose factors), one adjoint RHS, then cheap local reductions
over parameters. Preferred when `N_parameters ≫ N_objectives`.

## Sign checklist

| Term | Sign |
|------|------|
| Forward sensitivity RHS | `−(∂A)x` |
| Adjoint equation | `A† λ = L† y` |
| Gradient | `−2 Re(λ† (∂A) x)` |

## Non-goals of this note

Shape derivatives, mesh motion, and source dependence are excluded. Validation
against FD and against the forward Jacobian is required before the adjoint is
declared correct (`ADJOINT_VALIDATION_REPORT.md`).
