# Stage-5 dense spectral audit

**Design:** Stage-5 `q_best` (`opt_broadband_stage5_A2_C.json`)  
**Forward-only** (no optimization). Stage-6 **NOT STARTED**.

Companion: `BB_STAGE5_DENSE_AUDIT_SUMMARY.json`

## 1. Baseline — full-band 6.25 MHz (n = 321)

Same grid as the Stage-5 objective.

| Metric | f (GHz) | Value |
|--------|--------:|------:|
| Worst A2 | **5.1375** | **−7.585** |
| Worst IL | 5.09375 | 27.57 dB |
| Worst isolation | 4.575 | −10.33 dB |

Objective-grid softmin worst matches the audit baseline worst (5.1375).

### Historical watch (on 6.25 MHz nodes)

| f (GHz) | Stage-4 dense A2 (approx) | Stage-5 A2 | Stage-5 IL | Stage-5 iso |
|--------:|--------------------------:|-----------:|-----------:|------------:|
| 4.73125 | ≈ −10.15 | **−5.88** | 19.7 | −3.4 |
| 4.775 | repaired on Stage-4 | **−2.45** | 16.7 | **+8.4** |
| 5.08125 | **−11.04** | **−6.17** | 22.7 | −2.0 |
| 5.5875 | bad in Stage-3 | **−6.54** | 19.6 | −6.0 |
| 5.8125 | −12.02 → Stage-4 ≈−4.5 | **−5.19** | 16.2 | −4.9 |

Historical Stage-4 corpses are **improved on-node**, but not eliminated.

## 2. Refine — 3.125 MHz local neighborhoods

Centers: all historical failure regions + every 6.25 MHz local A2 minimum with A2 ≤ −3 + true baseline worst.  
n ≈ 192 new solves (seeded with 6.25 where overlapping).

### True worst (authoritative for this audit resolution)

| Quantity | Value |
|----------|------:|
| **f** | **4.734375 GHz** |
| **A2** | **−8.741** |
| **IL** | **25.04 dB** |
| **isolation** | **−11.03 dB** |
| **T_des** | 3.32×10⁻³ |

This frequency is **between** objective nodes 4.73125 (A2≈−5.88) and 4.7375 (A2≈−7.32). Soft-min therefore under-weighted the true floor during optimization.

### Other serious 3.125 MHz features (A2 ≲ −8)

| f (GHz) | A2 | IL (dB) | iso (dB) |
|--------:|---:|--------:|---------:|
| 4.734375 | −8.74 | 25.0 | −11.0 |
| 4.228125 | −8.28 | 22.4 | −11.3 |
| 5.084375 | −8.24 | 26.6 | −7.2 |
| 5.571875 | −8.05 | 21.9 | −10.5 |
| 5.559375 | −7.97 | 20.6 | −10.8 |

Multiple sharp local minima with A2 ≤ −5 remain across 4–6 GHz at 3.125 MHz resolution.

## 3. Comparison vs Stages 3–4

| Stage | True dense-audit worst | A2 | Notes |
|------:|------------------------|---:|-------|
| 3 | 5.8125 (12.5 MHz local) | ≈ −12.02 | between sparse nodes |
| 4 | 5.08125 (6.25 local) | ≈ −11.04 | new hole after enrichment |
| **5** | **4.734375 (3.125 local)** | **≈ −8.74** | between **6.25 MHz objective** nodes |

Worst absolute A2 improved (less catastrophic), but the **between-grid notch mechanism is not gone**.

## 4. Maximum gap / sampling implication

- Objective max gap: **6.25 MHz** (by construction).  
- Observed residual notch depth at **3.125 MHz** offset from a sampled node proves that **6.25 MHz spacing is still insufficient** for a fully defensible continuum claim over 4–6 GHz.
- A scientifically safer Stage-6 (if ever authorized) would need **full-band ≤ 3.125 MHz** or an equivalent max-gap rule — **not** isolated enrichment of 4.734375 alone.

## 5. Conclusion

### **B — materially improved but narrower residual(s) remain**

Not **A**: continuum not resolved at the tested audit resolution.  
Not full **C**: the full-band 6.25 MHz grid did stop the Stage-2→4 pattern of *completely silent* multi-25 MHz gaps and improved many historical corpses; soft-min now sees a dense set of competing on-grid failures. But features **narrower than 6.25 MHz** still deepen between nodes (true worst A2≈−8.7 off-grid).

```
STAGE-6: NOT STARTED
```
