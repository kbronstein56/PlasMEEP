# Circulator objective gradient validation

Three-way check of the **complete** six-port power objective:
centered FD, forward Jacobian chain rule, port-power adjoint.

Predeclared relative tolerance: `0.0001` (with truncation ≥5×).
Weights and normalization frozen in `CIRCULATOR_OBJECTIVE_DERIVATION.md`.

## grade C, B = 0.05 T

- J(s=1) = -6.920525e-01, P_ref = 1.364567e-02

| rod | g_adj | min FD vs adj | fwd vs adj | pass |
|-----|------:|--------------:|-----------:|:----:|
| center(85) | -6.7295e-01 | 4.972e-07 | 1.155e-14 | yes |
| inner(68) | 4.3949e-01 | 1.293e-05 | 1.680e-14 | yes |
| outer(0) | -5.6254e-01 | 1.324e-05 | 1.914e-14 | yes |
| near_horn0(40) | -1.4412e-01 | 5.494e-06 | 5.970e-15 | yes |
| sym_a(35) | 5.0063e-01 | 2.052e-05 | 1.109e-14 | yes |
| sym_b(64) | 5.2609e-01 | 2.075e-05 | 1.477e-15 | yes |

Tied-orbit FD: min rel `8.055e-05` → PASS

**PHYSICAL_PORT_OBJECTIVE_GRADIENT: PASS**
