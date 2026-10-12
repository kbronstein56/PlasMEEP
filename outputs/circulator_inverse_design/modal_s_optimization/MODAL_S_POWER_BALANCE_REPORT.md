# Phase 2 — Modal S power / absorption consistency

Grade: FULL91-M · f=3.85 GHz · B=+0.05 T

Identity (algebraic by construction of `P_PML`):

```
1 = Σ_i |S_ij|² + P_material/P_inc + P_PML_and_higher/P_inc
```

with `P_PML_and_higher = P_inc − Σ|b_fund|² − P_material`.

## Results

| Design | PASS | Σ\|S\|² | P_mat/P_inc | P_PML+higher/P_inc | T_des |
|--------|:----:|-------:|------------:|--------------------:|------:|
| uniform | yes | 0.936 | 0.026 | 0.038 | 0.007 |
| sixfold peak | yes | 0.831 | 0.041 | **0.128** | 0.742 |
| broadband | yes | 0.808 | 0.046 | **0.147** | 0.717 |

## Interpretation

Optimized designs put **~13–15%** of incident power into PML + higher-mode
channels not captured by fundamental `|S|²`, vs ~4% for uniform. That budget
matches the TE1 / thin-PML story in the Onsager audit and is **not** silently
lumped into material absorption.

Material absorption remains ~4–5% on optimized designs.
