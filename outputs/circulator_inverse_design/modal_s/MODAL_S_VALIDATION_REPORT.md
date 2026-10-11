# Modal S-parameter validation report

Branch: `agent/eigenmode-ports`  
Gate: **MODAL_SIXPORT_S_VALIDATED** — PASS

## Results

| Case | Test | Result | Key number |
|------|------|--------|------------|
| 0 | Synthetic a,b + power identity | PASS | err < 1e-15 |
| 1 | Straight matched guide `S21→1`, `S11→0` | PASS | \|S11\|²=0.0000, \|S21\|²=1.0041 |
| 2 | Quarter-wave ε=3.8 slab vs Fresnel | PASS | R=0.3369 (an 0.3403), T=0.6651 (an 0.6597) |
| 3 | Full91-C B=0 reciprocity `S≈Sᵀ` | PASS | ‖S−Sᵀ‖_F/‖S‖_F = 2.5×10⁻⁴ |
| 4 | Onsager `S(+B)≈Sᵀ(−B)` at ±0.05 T | PASS | 1.8×10⁻⁴ |
| 5 | Passive modal power balance | PASS | P_abs_modal/P_inc ≈ 0.062 > 0; \|P_abs_modal−P_abs_mat\|/P_inc ≈ 0.036 |

## Power-normalization check

On the straight guide and synthetic fields, `|a|²` and `|b|²` match ∫S·n̂ modal decomposition (`∫S·n̂ = |b|²−|a|²`) to machine/FEM precision.

## Case 5 note (PML)

`P_abs_modal = |a_j|² − Σ|b_i|²` includes **PML** absorption at the outer boundary.
`P_abs_material` integrates only lossy plasma/quartz. Therefore
`P_abs_material ≤ P_abs_modal` is expected; they need not be equal.
Both are non-negative for the magnetized uniform design, confirming passivity.

## Thresholds

Thresholds were **not** loosened relative to the stated plan (reciprocity/Onsager < 5%;
slab ΔR,ΔT < 0.03; matched |S11|²<0.02, ||S21|−1|<0.05).
