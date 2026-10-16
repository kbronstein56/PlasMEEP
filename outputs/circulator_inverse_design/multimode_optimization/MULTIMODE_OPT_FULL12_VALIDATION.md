# MULTIMODE_OPT_FULL12_VALIDATION

Design: `outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json` label=primary grade=C.

- Onsager r_F = 1.983e-13 (bar 0.005) → **PASS**
- Power balance max residual = 0.1361 (bar 0.05) → **FAIL**

## TE0 metrics

| B | T_des | T_rev | IL | iso | RL | TE1 |
|---|------:|------:|---:|----:|---:|----:|
| +0.05 | 0.8041 | 0.0101 | 0.95 | 19.00 | 16.63 | 0.0187 |
| −0.05 | 0.0101 | 0.8041 | 19.95 | -19.00 | 16.63 | 0.0221 |
| 0 | 0.0510 | 0.0510 | 12.93 | 0.00 | 17.65 | 0.1192 |

B-reversal |T_des(+B) − T_rev(−B)| = 3.4195e-14

## Power balance per channel

| ch | port | mode | P_out | P_mat | residual | PASS |
|---:|-----:|-----:|------:|------:|---------:|:----:|
| 0 | 0 | 0 | 0.9058 | 0.2303 | 0.1361 | False |
| 1 | 0 | 1 | 0.9865 | 0.0335 | 0.0200 | True |
| 2 | 1 | 0 | 0.9058 | 0.2275 | 0.1334 | False |
| 3 | 1 | 1 | 1.0246 | 0.0302 | 0.0548 | False |
| 4 | 2 | 0 | 0.9052 | 0.2295 | 0.1346 | False |
| 5 | 2 | 1 | 0.9946 | 0.0352 | 0.0298 | True |
| 6 | 3 | 0 | 0.9047 | 0.2309 | 0.1356 | False |
| 7 | 3 | 1 | 1.0103 | 0.0320 | 0.0423 | True |
| 8 | 4 | 0 | 0.9060 | 0.2245 | 0.1305 | False |
| 9 | 4 | 1 | 0.9919 | 0.0346 | 0.0265 | True |
| 10 | 5 | 0 | 0.9045 | 0.2314 | 0.1360 | False |
| 11 | 5 | 1 | 0.9985 | 0.0367 | 0.0352 | True |

