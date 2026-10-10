# Circulator objective gradient validation

Three-way / adjoint–FD check of the complete six-port power objective.
Predeclared relative tolerance: `1e-4` with truncation ≥5× (same gate as Jacobian).
Weights and `P_ref` frozen in `CIRCULATOR_OBJECTIVE_DERIVATION.md`.

## grade C, B = 0.05 T

FULL91-C: FD + forward + adjoint

- J(s=1) = -6.920525e-01, P_ref = 1.364567e-02

| rod | g_adj | min FD vs adj | pass |
|-----|------:|--------------:|:----:|
| center(85) | -6.7295e-01 | 4.972e-07 | yes |
| inner(68) | 4.3949e-01 | 1.293e-05 | yes |
| outer(0) | -5.6254e-01 | 1.324e-05 | yes |
| near_horn0(40) | -1.4412e-01 | 5.494e-06 | yes |
| sym_a(35) | 5.0063e-01 | 2.052e-05 | yes |
| sym_b(64) | 5.2609e-01 | 2.075e-05 | yes |

Tied-orbit FD min rel `8.055e-05` → PASS

## grade M, B = 0.05 T

FULL91-M: FD + adjoint (forward≡adjoint on C to machine precision)

- J(s=1) = -6.789093e-01, P_ref = 1.369015e-02

| rod | g_adj | min FD vs adj | pass |
|-----|------:|--------------:|:----:|
| center(85) | -6.1604e-01 | 1.074e-07 | yes |
| inner(68) | 4.7649e-01 | 1.398e-05 | yes |
| outer(0) | -5.8520e-01 | 1.374e-05 | yes |
| near_horn0(40) | -1.5022e-01 | 6.269e-06 | yes |
| sym_a(35) | 5.2963e-01 | 2.212e-05 | yes |
| sym_b(64) | 5.5421e-01 | 2.243e-05 | yes |

Tied-orbit FD min rel `9.185e-05` → PASS

**PHYSICAL_PORT_OBJECTIVE_GRADIENT: PASS**
