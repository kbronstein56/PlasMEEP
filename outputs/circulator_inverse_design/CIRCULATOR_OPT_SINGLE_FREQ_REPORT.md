# Single-frequency sixfold circulator optimization

Operating point: `f = 3.85 GHz`, `B = +0.05 T`, all six sources.
Mesh: FULL91-M. Variables: 16 C6-tied density scales. Optimizer: L-BFGS-B.
Objective / weights / `P_ref` frozen in `CIRCULATOR_OBJECTIVE_DERIVATION.md`.

## Result (uniform start `q≡1`)

| metric | initial | best |
|--------|--------:|-----:|
| J | −0.679 | **5.380** |
| desired avg `t_des` | 0.042 | **4.358** |
| reverse avg | 0.016 | 0.114 |
| through avg | 0.550 | **0.022** |
| leak avg | 0.113 | 0.064 |
| accepted avg | 1.000 | 5.390 |
| isolation avg | 4.11 dB | **15.83 dB** |
| insertion (10 log10 t_des) | −13.8 dB | +6.39 dB |
| worst-port `t_des` | 0.042 | 4.356 |

`t = P / P_ref` with frozen `P_ref = 1.369e-2` from uniform `s≡1`.
Values `t_des > 1` mean the optimized device delivers more desired-port
outward power than the uniform device’s mean accepted power scale.

## Circulation

Best 6×6 outward power matrix is nearly circulant: for every source `j`,
port `(j+1) mod 6` dominates. Parameter range `q ∈ [0.457, 1.404]`.

## Runtime

~3.2 h wall, 35 L-BFGS iterations, 80 function evaluations, peak mesh M.

## Artifact

`opt_single_uniform_M_Bp0.0500.json`
