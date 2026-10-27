# Objective A vs Engineering-Aligned A2 (decision)

**Date:** 2026-10-27  
**Spectra:** same stage-1 design (`opt_broadband_production_C.json`) — **no new FEM optimization**.  
**Evidence:** `BB_OBJECTIVE_ALIGNMENT_AUDIT.md`, `BB_OBJECTIVE_ALIGNMENT_OFFLINE.json`, `objective_A2.py`.  
**FD:** `A2_GRADIENT_FD_VALIDATION.json` — **OVERALL PASS**
(SF rel \(1.75\times10^{-5}\); 3-freq soft-min rel \(4.07\times10^{-5}\); bar \(10^{-3}\)).

---

## Why did stage-1 improve J_BB while mean T_des≈0.094, IL≈14.3 dB, isolation≈1.4 dB?

Stage-1:

| Metric | Start | End |
|--------|------:|----:|
| \(J_{BB}\) (soft-min A) | −0.801 | **−0.529** |
| mean \(T_\mathrm{des}\) (final) | — | **0.094** |
| mean IL (final) | — | **14.3 dB** |
| mean isolation (final) | — | **1.4 dB** |
| worst IL on grid | — | **39.9 dB @ 4.1 GHz** |
| worst iso on grid | — | **−18.8 dB @ 4.1 GHz** |

Optimizer did its job on the **scalar it was given**. Soft-min A pushed up the worst *absolute-power* \(J_A\) frequencies (large reverse / TE1). It did **not** force useful transmission, because:

1. Absolute \(+T_\mathrm{des}\) rewards are linear and saturate at O(1); a dead channel (\(T_\mathrm{des}\sim10^{-4}\)) costs almost nothing relative to a −0.5 TE1/reverse event.
2. Absolute \(−2\,T_\mathrm{rev}\) only hurts when reverse power is *large*. Collapsing both desired and reverse toward zero (reflector / absorber / mode-scrambler) escapes the reverse penalty.
3. Soft-min therefore “improved” by trading one absolute-power failure mode for another that is still an engineering failure (4.1 GHz).

Re-scoring the **same** 21 final spectra under A2_hybrid:

| Objective | mean \(J_f\) | \(\min J_f\) | \(J_{BB}\) (τ=0.05, α=0.25) | Rank of 4.1 GHz |
|-----------|-------------:|-------------:|----------------------------:|----------------:|
| A (current) | −0.442 | −0.657 @ 5.7 | −0.529 | 11 / 21 |
| **A2_hybrid** | −3.28 | **−11.56 @ 4.1** | **−9.37** | **1 / 21** |

Under A2, stage-1’s “success” is revealed as a large residual soft-min hole at 4.1 GHz — consistent with 40 dB IL.

---

## Would enrichment alone plausibly fix this?

**No.**

- 4.1 GHz was **already** an objective frequency.
- Dense audit added 5.625 GHz as a between-grid soft-min hole for *A*; enrichment fixes sampling, not the FoM.
- Offline: A_current still ranks 4.1 as **34th-worst** of 52 dense points (not top-5). A2 ranks it **#1**.

Enrichment remains necessary for off-grid spikes (5.625, 5.65, …) but must be paired with A2.

---

## Decision

| Question | Answer |
|----------|--------|
| Is Objective A misaligned with engineering quality? | **YES** |
| Preferred replacement | **A2_hybrid** = \(2\log(T_\mathrm{des}+\varepsilon)-\log(T_\mathrm{rev}+\varepsilon)\) − linear thr/other/R/TE1 |
| Hard IL/iso dB thresholds? | **Not required**; \(\varepsilon=10^{-6}\) only |
| Keep C6 / TE0 mode-preserving / soft-min τ,α? | **YES** |
| Overwrite `BROADBAND_PRODUCTION_SPEC_FROZEN.md`? | **NO** — historical stage-1 freeze |
| Stage-2 authorization | Proposed spec only after FD PASS; **do not auto-launch** |

Code freeze of A2 (pending FD): `scripts/.../objective_A2.py` (`A2_FREEZE`).

---

## Gradient path (Phase 4)

Channel weights for column \(n=2j_\mathrm{src}\):

\[
\frac{\partial J_j}{\partial |S_{m,n}|^2}
=
\begin{cases}
2/(T_\mathrm{des}+\varepsilon) & m=\mathrm{des\,TE0}\\
-1/(T_\mathrm{rev}+\varepsilon) & m=\mathrm{rev\,TE0}\\
-W_\mathrm{THR} & m=\mathrm{thr\,TE0}\\
-W_\mathrm{LEAK} & m=\mathrm{other\,TE0}\\
-W_\mathrm{REFL} & m=\mathrm{src\,TE0}\\
-W_\mathrm{TE1} & m=\mathrm{any\,TE1}
\end{cases}
\]

Same bilinear discrete adjoint as production (`A^\top`, 1 LU + 6+6 RHS). No JAX.

FD bars: relative error ≤ 1e-3 on (1) single-frequency orbit FD and (2) 3-frequency soft-min E2E FD — **both PASS** (`A2_GRADIENT_FD_VALIDATION.json`).

Stage-2 production freeze: `BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md` (does not launch optimization).
