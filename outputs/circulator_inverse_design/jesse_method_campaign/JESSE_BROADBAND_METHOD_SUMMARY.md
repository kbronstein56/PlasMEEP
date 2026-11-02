# Jesse broadband method summary

**Stage-2 COMPLETE** (A2_hybrid + 28-freq enrichment).  
Stage-1 historical freeze remains `BROADBAND_PRODUCTION_SPEC_FROZEN.md`.  
Stage-2 freeze: `BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md`.  
**Do not auto-start Stage-3.**

Continuum SF quote (not M-mesh): \(J\approx0.626\pm0.007\).

---

### 1. How many frequency points?

Stage-1 trained on **21×100 MHz**. Dense audit found soft-min hole at **5.625 GHz**.  
Stage-2 used **28** fixed frequencies (21 + enrichment). Post-Stage-2 **50 MHz** audit
is the next sampling check (in progress / see `bb_stage2_audit_50MHz.json`).

### 2. How are narrow between-grid failures prevented?

Fixed set per L-BFGS → dense audit → add worst freqs → restart.  
Stage-1 → Stage-2 already applied this once. Further enrichment only after Stage-2 audit.

### 3–6. Differentiable BB / adjoint / efficiency

Soft-min τ=0.05 α=0.25 retained after A2 scale sensitivity
(`A2_SOFTMIN_SCALE_SENSITIVITY.md`): at Stage-2 start τ∈[0.05,0.50] was a weighting
no-op; α=0.25 supplies multi-freq floor. FD discrete adjoint; no JAX production.

### 7–9. Mode / experiment / C6

**Objective A2_hybrid** (mode-preserving TE0; log transmission + log isolation).  
C6 unchanged. Objective A retained only as Stage-1 historical FoM.

### 10. Time-domain

Still planned; not blocking post-Stage-2 audit.

### Stage-1 production result (historical)

| | Value |
|--|------:|
| FoM | Objective A |
| J_BB | −0.801 → **−0.529** |
| mean T_des / IL / iso | 0.094 / 14.3 dB / 1.4 dB |
| pathology | 4.1 GHz IL≈40 dB under-penalized by A |

### Stage-2 production result

| | Value |
|--|------:|
| FoM | **A2_hybrid** |
| freqs | **28** |
| J_BB | −11.33 → **−6.229** |
| 4.1 GHz | T_des 3.4e-4→**0.141**, IL 40→**8.5 dB**, iso −19→**+7.1 dB** |
| final worst | A2 −7.40 @ **5.65 GHz** (IL≈18.5, iso≈−11.9) |
| early flags | **none** |
| wall | ~78 h / 30 L-BFGS iters |

Details: `BB_STAGE2_RESULTS.md`, `BB_STAGE1_VS_STAGE2.md`.

### Enrichment / Stage-3?

**50 + 25 MHz audits COMPLETE — enrichment needed.**  
Off-grid: worst A2/IL @ **5.75** (A2≈−11.0, IL≈26 dB, iso≈−20 dB); worst iso @
**5.825** (≈−21.3 dB). On-grid residual: **4.3 GHz** IL≈27.5 dB.

If Stage-3 approved (`BB_STAGE3_PROPOSAL.md`):
- why: close 5.75/5.825 between-grid hole
- add: **5.75, 5.825** → 30 freqs
- objective: A2_hybrid; τ=0.05; α=0.25; C6
- start: Stage-2 `q_best`
- question: close hole without undoing 4.1 recovery?

**Do not auto-start Stage-3.**
