# Jesse broadband method summary

Production stage-1 BB opt **completed** (nit=30, J_BB −0.801 → **−0.529**).
Dense audit shows **enrichment needed** before any restart.
Design: `opt_broadband_production_C.json` (`q_best`).

Continuum SF quote (not M-mesh): \(J\approx0.626\pm0.007\).

---

### 1. How many frequency points?

Initial production used **21×100 MHz**. Post-opt 50 MHz: worst iso/IL already
on-grid @ 4.1. Soft-min fails between grids: worst J moves **5.7 → 5.65 → 5.625**
(J_A to **−1.15**). Add priority freqs before restart (see § enrichment).

### 2. How are narrow between-grid failures prevented?

Fixed set per L-BFGS → dense audit → add worst freqs → restart. Stage-1 audit
identified the list below; **do not auto-start** another 30-iter run.

### 3–6. Differentiable BB / adjoint / efficiency

Unchanged: soft-min τ=0.05 α=0.25; FD adjoint; latest cost assembly+dA dominant
(LU ~9 s historical-84%-LU obsolete).

### 7–9. Mode / experiment / C6

Frozen: Objective A, C6.

### 10. Time-domain

Still planned; not blocking enrichment decision.

### Stage-1 production result

| | Value |
|--|------:|
| J_BB | −0.801 → **−0.529** |
| mean / worst J (obj) | −0.442 / −0.657 @ 5.7 |
| mean T_des / IL / iso / TE1 | 0.094 / 14.3 dB / 1.4 dB / 0.256 |
| stop | iteration limit (30), 44 h |

### Enrichment needed?

**YES.** Priority frequencies to add to the objective set:

**`4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95`** (→ 28 freqs total)

Full spike list and rationale: `BB_ENRICHMENT_RECOMMENDATION.md`.
