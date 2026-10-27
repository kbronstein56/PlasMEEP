# BROADBAND_PRODUCTION_SPEC_STAGE2_PROPOSED

**Status:** SUPERSEDED by `BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md` after
`A2_GRADIENT_FD_VALIDATION.json` overall PASS (2026-10-27).  
Kept for history. Still does **not** authorize auto-launch.

Evidence: `BB_OBJECTIVE_ALIGNMENT_AUDIT.md`, `BB_OBJECTIVE_A_VS_A2.md`, `BB_ENRICHMENT_RECOMMENDATION.md`.

---

## Proposed choices (delta vs stage-1 freeze)

| Item | Stage-1 (frozen) | Stage-2 (proposed) |
|------|------------------|--------------------|
| Symmetry | C6, 16 vars | **unchanged** |
| Incident / desired | TE0 / TE0 mode-preserving | **unchanged** |
| Per-frequency FoM | Objective A (linear \|S\|²) | **Objective A2_hybrid** |
| Broadband combiner | soft-min τ=0.05, α=0.25 | **unchanged** |
| Frequency set | 21 × 100 MHz | **28** = 21 + enrichment |
| Mesh | FULL91-C | **unchanged** |
| Gradient | discrete \(A^\top\) adjoint | **unchanged** (A2 channel weights) |
| Warm start | uniform | **`opt_broadband_production_C.json` q_best** |

### Objective A2_hybrid

\[
J_{A2} = 2\log(T_\mathrm{des}+\varepsilon) - \log(T_\mathrm{rev}+\varepsilon)
 - 1.5\,T_\mathrm{thr} - T_\mathrm{other} - 0.25\,R - T_\mathrm{TE1}
\]

with \(\varepsilon=10^{-6}\). No hard IL/isolation dB thresholds. Implementation: `objective_A2.py`.

### Frequency set (28)

```
4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.8, 4.85, 4.9,
5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.6, 5.625, 5.65, 5.7,
5.8, 5.85, 5.9, 5.95, 6.0
```

Added vs stage-1: `4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95`.

---

## Proposed run parameters (manual approval only)

| Parameter | Value | Rationale |
|-----------|------:|-----------|
| maxiter | **30** | match stage-1 L-BFGS budget |
| bounds | [0.05, 1.8] | match stage-1 |
| grade | C | match stage-1 |
| start design | `opt_broadband_production_C.json` | warm-start |

### Expected wall-time

Stage-1: **~44 h** wall for 21 freqs × 30 iters (observed `elapsed_s` ≈ 1.58×10⁵).

Scale ≈ \(28/21 \approx 1.33\):

| Estimate | Hours |
|----------|------:|
| Nominal (30 iter, 28 freq) | **~55–60 h** |
| Conservative (+line-search / LU variance) | **~65–75 h** |

Do **not** auto-launch. After a completed stage-2 run: dense ≤50 MHz audit; do not auto-chain further opts.

---

## Freeze gate

Promote this file to `BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md` only when:

1. `A2_GRADIENT_FD_VALIDATION.json` → `pass: true` (SF and 3-freq soft-min, rel ≤ 1e-3)
2. Human approval to start the enrichment+A2 optimization
