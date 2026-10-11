# Modal-S adjoint gradient (spot check)

Gate: **MODAL_S_ADJOINT_FD** — PASS (FULL91-C)

## Objective

Modal analogue of the Pref FoM, averaged over six sources:

```
J = (1/6) Σ_j [ T_des − 2 T_rev − 1.5 T_thr − T_leak − 0.25 |S_jj|² ]
```

with `T = |S|² = |b/a_j|²` from the validated modal extraction.

## Adjoint

`b_i = α_i · x`, `a_j = α_a · x` with α from the flat-mode (Hz, ∂n) projector.
Quotient rule for `|b/a|²`; solve `Aᴴ λ = q`, then

```
∂J/∂s_k = −(1/6) Σ_j Re(λ_j† (∂A/∂s_k) x_j)
```

Freezing `a_j` (vacuum-feed approximation) is **not** accurate enough here
(~7–37% FD error). Full `da_j` terms are required.

## FD check (uniform s≡1, B=+0.05 T, ε=10⁻³)

| rod | adjoint | FD | rel err |
|----:|--------:|---:|--------:|
| 0 | −0.130197 | −0.130198 | 1.4×10⁻⁵ |
| 7 | −0.034416 | −0.034417 | 8.2×10⁻⁶ |
| 14 | −0.183643 | −0.183644 | 4.3×10⁻⁶ |

No new 91-variable optimization was launched.
