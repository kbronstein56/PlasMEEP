# Full 91-bulb Jacobian validation

Material-only density scales `s_k` on the validated FULL91-F mesh
(4,137,402 nodes). One factorization of `A` at production `s=1`; each
column is one sensitivity RHS. Frozen gate: relative error ≤ 1e-4 with
truncation visible.

Scripts: `jacobian_full91.py`, `diff_fem_core.py`  
Data: `full91_jacobian_fd_selected.json`, `full91_jacobian_summary.json`,
`full91_jacobian_F_Bp0.0000_src0.json`

## Selected rods (source port 0, B = 0)

| label | rod | min field rel | min amp rel (worst recv) | pass |
|-------|-----|---------------|--------------------------|------|
| center | 85 | 1.03e-6 | 1.03e-6 | yes |
| inner | 68 | 8.00e-7 | 8.68e-7 | yes |
| outer | 0 | 5.17e-7 | 6.52e-7 | yes |
| near_horn0 | 40 | 6.90e-7 | 1.37e-6 | yes |
| sym_a | 23 | 8.34e-7 | 8.57e-7 | yes |
| sym_b | 81 | 8.34e-7 | 8.58e-7 | yes |

All receiving ports (1–5) were checked. Symmetry-related pair `sym_a` /
`sym_b` show matching FD curves at B = 0, as expected.

## Full port Jacobian (B = 0, source 0)

| quantity | value |
|----------|-------|
| shape | 5 recv ports × 91 rods |
| factor + primal solve | 112 s |
| all sensitivity RHS | 1386 s |
| time per column | 15.2 s |
| total Jacobian | 1497 s |
| peak RSS | 42.4 GiB |
| mean \|J\| | 0.158 |
| max \|J\| | 0.526 |

One LU; 91 sensitivity solves; no per-parameter refactorization.

## Magnetized FD

`B = +0.05 T` selected-rod FD is recorded separately when the companion
run finishes (`full91_jacobian_fd_selected` updated / tagged in the
summary). Thresholds are not loosened.

## Verdict

**FULL91_JACOBIAN: PASS** at `B = 0` on grade F for selected rods and the
complete 5×91 complex monitor Jacobian.
