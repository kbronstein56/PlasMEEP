# Differentiable FEM — final status

Label: **DIFFERENTIABLE_FEM_VALIDATED**

| Gate | Result |
|------|--------|
| FORWARD_JACOBIAN | PASS |
| FULL91_JACOBIAN | PASS |
| ADJOINT | PASS |
| OPTIMIZATION_PROOF | PASS |

Branch: `agent/eigenmode-ports`

## Exact derivative equations

Material-only density scale `s` with `f_p² = s f_p,ref²`, fixed mesh/source:

    A(s) x(s) = b
    A dx/ds = − (dA/ds) x
    dA/ds = K(dρ/ds),   dρ = −ρ (dε) ρ
    dε_⊥/ds = − f_p,ref² u / (f (u²−f_c²))
    dη/ds   =   f_p,ref² f_c / (f (u²−f_c²))

Real objective `J = |L x|²`:

    A† λ = L† y,   y = L x
    dJ/ds = −2 Re( λ† (dA/ds) x )

## Worst derivative errors (frozen ≤ 1e-4)

| Suite | Worst min rel | Notes |
|-------|---------------|-------|
| One-cylinder field (B=0,±0.05) | 7.7e-10 | ε/ρ FD also pass |
| Multicylinder columns | 5.5e-8 | index OK |
| Multiport amp/power/phase | 5.0e-8 | |
| Full91-F selected rods B=0 | 1.4e-6 (amp) | 6 rods |
| Full91-F selected rods +0.05 T | 1.9e-6 (amp) | 6 rods |
| Adjoint vs forward (1-cyl) | ≤1.4e-13 | |
| Adjoint full91-M spot | ≤1.2e-13 fwd-adj; FD ≤2e-5 | 3 rods |

## Full91 Jacobian timing (grade F, B=0)

- Factorization + primal: 112 s
- 91 sensitivity RHS: 1386 s (15.2 s/column)
- Total 5×91 port Jacobian: 1497 s
- Peak RSS: 42.4 GiB

## Fan 2019

Mathematically equivalent forward-mode framework (`A dx = −(dA)x`, one RHS
per parameter, operator reuse). Discretization is P1 FEM, not FDTD.
See `FAN_FORWARD_MODE_COMPARISON.md`.

## Commits (Goal B sequence)

See `git log` on `agent/eigenmode-ports` from `5480b2e` onward.
