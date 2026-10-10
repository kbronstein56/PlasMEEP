# Port-power adjoint derivation

The Goal-B adjoint for `J = |L x|²` does **not** apply unchanged to
guide-normal Poynting power. This note derives the correct discrete adjoint
for the vacuum-feed monitor used in FULL91.

## Discrete power

At sample `s` on the monitor line (trapezoid weight `w_s`, outward `n̂`):

    H_z = ℓ_s · x          (barycentric, ℓ real)
    ∇H_z · n̂ = g_s · x     (g from P1 gradients, real)
    (S·n̂)_s = ½ Re[ α (g_s·x) (ℓ_s·x)* ],   α = −i/ω

    P = ∑_s w_s (S·n̂)_s

Monitors are in vacuum: `ρ = I` on those elements. **Q does not depend on
the plasma design variables s.** Explicit `∂P/∂s|_x = 0`.

## Differential

    dP = ∑_s (w_s/2) Re[ α (g_s·dx) (H_z)* + α (∇H_z·n̂) (ℓ_s·dx)* ]

This is linear in `dx`. Collecting nodal contributions gives a complex
vector `q` such that

    dP = Re( q† dx )

Explicit assembly (per sample, accumulate onto the three triangle nodes):

    q += (w_s/2) conj(α) H_z  g_s     # from α (g·dx) H_z*
    q += (w_s/2) α (∇H_z·n̂) ℓ_s       # from α (∂n H) (ℓ·dx)*
    # so that dP = Re(q† dx) with NumPy vdot convention q†dx = Σ conj(q_i) dx_i

## Constraint

Material-only, fixed source:

    A dx = − (dA) x

## Adjoint

Choose `λ` by

    A† λ = q

Then

    dP = Re(q† dx) = Re(λ† A dx) = Re(λ† (−dA x)) = − Re(λ† (dA) x)

so

    ∂P/∂s_k = − Re( λ† (∂A/∂s_k) x )

**Sign / dagger:** conjugate transpose `A†`, not plain transpose, because the
pairing is `Re(q† dx)` over complex DOFs (same complex-calculus reason as
in `ADJOINT_DERIVATION.md`).

**Factor relative to `|Lx|²`:** that case used `dJ = 2 Re(y* L dx)` and
`A† λ = L† y` with `∂J/∂s = −2 Re(λ† dA x)`. Here the factor `2` is absorbed
into the definition of `q` via `dP = Re(q† dx)`.

## Full circulator objective

`J = (1/6) ∑_j [ w_des t_des,j − … + w_acc a_j ]` with `t = P/P_ref`,
`a_j = (−P_jj)/P_ref`, and frozen `P_ref`.

For each source field `x_j`:

    q_j = (1/(6 P_ref)) * (
          w_des  q(P_des)
        − w_rev  q(P_rev)
        − w_thr  q(P_thr)
        − w_leak (q(P_leak_a)+q(P_leak_b))
        − w_acc  q(P_jj)     # because a_j = −P_jj/P_ref
    )

Solve `A† λ_j = q_j` (one adjoint RHS per source, same LU as the forward
factorization via `trans='H'`).

    ∂J/∂s_k = ∑_j − Re( λ_j† (∂A/∂s_k) x_j )

For a sixfold orbit variable `q_m` with `s_i = q_m` on the orbit:

    ∂J/∂q_m = ∑_{i ∈ orbit m} ∂J/∂s_i
