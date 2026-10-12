# Phase 3 — Modal-S adjoint at nonuniform Pref designs

**Verdict: `MODAL_S_NONUNIFORM_GRADIENT = PASS`**

Objective (full quotient rule, `a_j` not frozen):

```
J_modal = (1/6) Σ_j [ T_des − 2 T_rev − 1.5 T_thr − T_leak − 0.25 R ]
```

Grade: FULL91-C · f=3.85 GHz · B=+0.05 T  
Centered FD h-sweep; gate `max_rel ≤ 5×10⁻⁴` at the best h.

## Rod set

| role | rod |
|------|----:|
| center | 85 |
| inner ring | 74 |
| outer ring | 0 |
| near horn 0 | 10 |
| C6 pair | 0, 1 |

## Results

| Design | max best rel err | PASS |
|--------|-----------------:|:----:|
| sixfold peak | **3.4×10⁻⁵** (h=3×10⁻⁶, rod 10) | yes |
| broadband | **2.9×10⁻⁴** (h=10⁻⁴) | yes |

Near-horn rod 10 on the peak design is nonlinear in `s`: coarse h (10⁻³–10⁻⁴)
gives O(10%) FD bias; refining to h∈{10⁻⁵, 3×10⁻⁶} recovers adjoint agreement
to 3×10⁻⁵. Functional linearity `α·x` vs `extract_ab` remains ~10⁻¹⁴.

C6 pair (rods 0,1) on the peak design: gradients agree to ~0.5% (`s` tied).

## Artifacts

- `grad_valid_summary_C.json` (broadband + initial peak)
- `grad_valid_opt_single_uniform_M_Bp0.0500_C.json` (refined peak)
- `grad_valid_C_peak_refined.log`
- `grad_rod10_debug.json`
