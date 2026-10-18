# M→F transfer (M-optimized design)

Design: `outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C_warm_M.json`
Criterion: `MESH_TRANSFER_CRITERION.md` (frozen before judgment).

**MESH_TRANSFER: FAIL**
**MESH_CONVERGENCE: FAIL**

| metric | M | F | Δ(F−M) | bar | pass |
|--------|--:|--:|-------:|----:|:----:|
| J | 0.7229 | 0.6698 | -0.0531 | 0.05 | N |
| T_des_avg | 0.8187 | 0.7953 | -0.0234 | 0.04 | Y |
| IL_avg_dB | 0.8690 | 0.9949 | +0.1259 | 0.4 | Y |
| isolation_avg_dB | 20.3038 | 18.4068 | -1.8969 | 1.5 | N |
| RL_avg_dB | 15.9680 | 13.2353 | -2.7327 | 1.5 | N |
| TE1_avg | 0.0265 | 0.0203 | -0.0062 | 0.015 | Y |
| T_thr_avg | 0.0010 | 0.0018 | +0.0008 | 0.01 | Y |
| T_other_avg | 0.0462 | 0.0677 | +0.0215 | 0.03 | Y |

M→F remains large → do **not** call mesh-converged.
Investigate local interface refinement / finer reference rather than
only more optimizer iterations.
