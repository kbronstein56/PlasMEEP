# Stage-1 vs Stage-2 comparison (same 28-frequency set)

**Designs**
- Stage-1: `opt_broadband_production_C.json` (Objective A, 21-freq train; rescored here)
- Stage-2: `opt_broadband_stage2_A2_C.json` (Objective A2_hybrid, 28-freq train)

**Frequency set (frozen Stage-2):**  
`4.0, 4.1, … + {4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95}` (28 pts)

A2 scores for Stage-1 use offline `J_A2_hybrid` on Stage-1 spectra (obj-grid + dense audit).  
Stage-2 A2 scores are the true production per-source-mean `objective_A2_from_S` from the final rec.

Soft-min aggregation: τ=0.05, α=0.25 (identical).

---

## Broadband A2 objective

| Metric | Stage-1 design @ 28 | Stage-2 best @ 28 |
|--------|-------------------:|------------------:|
| \(J_{BB}\) (A2 soft-min) | −9.557 | **−6.229** |
| mean A2 | −4.054 | **−3.130** |
| worst A2 | −11.558 @ **4.1** | **−7.397 @ 5.65** |

Δ\(J_{BB}\) ≈ **+3.33** (higher is better).

---

## Engineering at 4.1 GHz (on-grid for both)

| Metric | Stage-1 | Stage-2 |
|--------|--------:|--------:|
| \(T_\mathrm{des}\) | 3.36×10⁻⁴ | **0.141** |
| IL (dB) | 39.94 | **8.51** |
| isolation (dB) | −18.80 | **+7.07** |
| RL (dB) | 1.60 | 4.23 |
| TE1 | 0.220 | 0.326 |
| A2 | −11.6 / −13.9 | **−0.84** |

Stage-2 **fixed the Objective-A pathology** at 4.1 GHz that Stage-1 left dead.

---

## Stage-1 engineering on the 28-set (from Stage-1 dense audit / obj metrics)

| Metric | mean | worst |
|--------|-----:|------:|
| \(T_\mathrm{des}\) | 0.077 | 3.36×10⁻⁴ @ 4.1 |
| IL (dB) | 14.95 | 39.94 @ 4.1 |
| isolation (dB) | −1.41 | −18.80 @ 4.1 |
| RL (dB) | 5.48 | (min RL ~0.6 class) |
| TE1 | 0.248 | — |

## Stage-2 engineering

### On final objective rec

| Location | \(T_\mathrm{des}\) | IL | iso | RL | TE1 | A2 |
|----------|-------------------:|---:|----:|---:|----:|---:|
| 4.1 GHz | 0.141 | 8.51 | +7.07 | 4.23 | 0.326 | −0.84 |
| worst A2 @ 5.65 | 0.018 | 18.50 | −11.89 | 6.79 | 0.132 | −7.40 |

### Full band 50 MHz audit (Stage-2 design)

| Metric | mean | worst | f (GHz) |
|--------|-----:|------:|--------:|
| \(T_\mathrm{des}\) | 0.102 | 0.0018 | 4.3 (on-grid) |
| IL (dB) | 13.01 | 27.47 | 4.3 (on-grid) |
| isolation (dB) | 2.58 | **−20.14** | **5.75 (off-grid)** |
| RL (dB) | 6.47 | — | — |
| TE1 | 0.286 | — | — |
| A2 | −2.85 | **−10.98** | **5.75 (off-grid)** |

---

## Scientific reading

1. Objective A2 + enrichment **worked as intended**: catastrophic 4.1 GHz transmission failure is gone.
2. Soft-min parked on **5.65**, but dense audit found a **worse between-grid hole at 5.75 GHz** (A2 −10.98, iso −20 dB) — same sampling lesson as Stage-1’s 5.625 GHz.
3. **4.3 GHz** remains a severe on-grid IL failure (27.5 dB).
4. Stage-3 is scientifically motivated to add **≥5.75 GHz**, but must **not** auto-start; confirm peak with 25 MHz local scan first (`BB_STAGE2_AUDIT_REPORT.md`).
