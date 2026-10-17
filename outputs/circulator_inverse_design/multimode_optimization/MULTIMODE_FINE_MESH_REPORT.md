# MULTIMODE_FINE_MESH_REPORT

Primary design: `opt_uniform_C.json` (uniform warm-start C6 multimode optimum).
Mesh transfer only — **no re-optimization** on M or F.

Sources: `fine_mesh_primary_C_to_M.json`, `fine_mesh_primary_M_to_F.json`.
Both transfers `fine_ok: true` (C→M ~2.2 h, M→F ~2.4 h).

## Full C / M / F comparison (TE0-drive means, B=+0.05 T, 3.85 GHz)

| metric | C | M | F | M−C | F−M |
|--------|--:|--:|--:|----:|----:|
| J | 0.7082 | 0.4875 | 0.2661 | −0.2207 | −0.2213 |
| T_des | 0.8041 | 0.6946 | 0.5582 | −0.1095 | −0.1364 |
| IL (dB) | 0.947 | 1.583 | 2.532 | +0.636 | +0.949 |
| isolation (dB) | 19.00 | 15.35 | 16.86 | −3.66 | +1.51 |
| RL (dB) | 16.63 | 11.48 | 12.63 | −5.16 | +1.16 |
| TE1 | 0.0187 | 0.0224 | 0.0439 | +0.0037 | +0.0215 |
| through | 0.0018 | 0.0084 | 0.0267 | +0.0066 | +0.0183 |
| other TE0 | 0.0488 | 0.1137 | 0.1715 | +0.0649 | +0.0578 |
| reflection R | 0.0218 | 0.0712 | 0.0546 | +0.0494 | −0.0166 |
| P_material | 0.2290 | 0.1571 | 0.2958 | −0.0719 | +0.1387 |

## Interpretation

- Evaluation on F **completed successfully**.
- Metrics are **not mesh-converged**: each refinement step drops J by ~0.22 and
  T_des by ~0.11–0.14. Do **not** quote C performance as continuum-accurate.
- Isolation and RL partially recover M→F while transmission continues to fall.
- TE1 rises on finer meshes (0.019 → 0.044) but remains far below Pref-peak ~0.32.

**MULTIMODE_FINE_MESH: EVALUATION_COMPLETE** (C, M, and F evaluations completed).
**MESH_CONVERGENCE: FAIL** (large C→M→F drift; continuum FoMs not achieved).

Do **not** label the fine-mesh result PASS merely because all meshes ran.
`EVALUATION_COMPLETE` ≠ mesh-converged FoMs.
