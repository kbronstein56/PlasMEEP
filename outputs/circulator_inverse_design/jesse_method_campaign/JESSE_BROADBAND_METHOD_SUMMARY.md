# Jesse broadband method summary

**Stage-2 COMPLETE** (A2_hybrid + 28-freq enrichment).  
**Stage-3 COMPLETE** (A2_hybrid + 30-freq enrichment; holes at 5.75/5.825 repaired).  
Stage-1 historical freeze remains `BROADBAND_PRODUCTION_SPEC_FROZEN.md`.  
Stage-2 freeze: `BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md`.

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

### Stage-3 production result

| | Value |
|--|------:|
| FoM | **A2_hybrid** |
| freqs | **30** (Stage-2 28 + 5.75 + 5.825) |
| J_BB | −9.02 → **−3.659** |
| 5.75 | A2 −10.98→**−2.95**, iso −20→**+3.1 dB** |
| 5.825 | A2 −10.71→**+1.38**, iso −21→**+17.9 dB** |
| 4.1 GHz | T_des 0.141→**0.289**, IL 8.5→**5.4 dB** (preserved + improved) |
| new 50 MHz worst | A2 **−4.57 @ 4.3 GHz** |
| wall | ~11.0 h / 30 L-BFGS / `n_workers=12` |

Details: `BB_STAGE3_RESULTS.md`, `BB_STAGE3_AUDIT_REPORT.json`.

### Post-Stage-3 dense audit (COMPLETE)

Full-band **25 MHz** + local **12.5 MHz** refine (`BB_STAGE3_DENSE_AUDIT.md`):

- Exact enrichment nodes **5.75 / 5.825 remain repaired**.
- **Hidden continuum failures found:** worst at **5.8125 GHz** (A2≈−12.0, iso≈−18.7 dB, IL≈32 dB) between objective nodes 5.800 and 5.825; also **4.775 GHz** (A2≈−8.6, IL≈26 dB).
- Recommendation: **C** — substantial hidden broadband failure remains.
- Stage-4: **PROPOSAL ONLY** (`BB_STAGE4_PROPOSAL.md`) — **NOT STARTED**.
