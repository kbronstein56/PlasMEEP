# Broadband frequency sampling study (4–6 GHz)

## Physics setup

- Target band: **4.0–6.0 GHz**
- Throughout: TE0 and TE1 propagate; TE2 remains cutoff (~7.49 GHz)
- **No** propagating-channel-count transition in-band

## Frozen grid-sufficiency bars (before seeing results)

File: `GRID_SUFFICIENCY_BARS_FROZEN.json`

| Metric | Bar when refining grid |
|--------|------------------------|
| Worst isolation change | **< 0.5 dB** |
| Worst IL change | **< 0.2 dB** |
| Worst RL change | **< 0.2 dB** |
| Worst-iso frequency location jump | **≤ 100 MHz** |

## Staged reconnaissance — RESULTS (SF design `opt_uniform_C_warm_M`, FULL91-C)

### 200 MHz grid (11 points) — complete

| Worst metric | Value | Frequency |
|--------------|------:|----------:|
| Isolation | **−24.2 dB** | 5.8 GHz |
| IL | **40.8 dB** | 6.0 GHz |
| RL | **0.85 dB** | 4.4 GHz |
| J_A | **−1.03** | 5.8 GHz |

At the SF design frequency (~3.85 GHz) the device is strong; across 4–6 GHz it is
**not** a broadband isolator. This is precisely Jesse’s spike/sampling concern:
a design good at one frequency can fail badly between / outside sample points.

### 100 MHz grid (21 points) — complete

| Worst metric | Value | Frequency |
|--------------|------:|----------:|
| Isolation | **−33.1 dB** | 5.9 GHz |
| IL | **43.5 dB** | 5.9 GHz |
| RL | **0.80 dB** | 4.3 GHz |

### 200 vs 100 comparison — **FAIL frozen bars**

From `spectral_grid_comparison.json`:

| Check | Δ | Pass? |
|-------|--:|:-----:|
| Worst isolation | **8.90 dB** | FAIL (<0.5) |
| Worst IL | **2.69 dB** | FAIL (<0.2) |
| Worst RL | 0.05 dB | PASS |
| Worst-iso location | 0.10 GHz | FAIL (≤0.1 strict) |

**Conclusion:** 200 MHz objective sampling is **not** grid-sufficient for this
spectrum. Need ≤100 MHz with adaptive enrichment toward 50/25 MHz near 5.7–6.0 GHz
and other features. Soft-min on 11 points would miss the 5.9 GHz catastrophe.

Optimization enrichment (when BB opt starts):

1. Fixed objective set per L-BFGS run  
2. Dense audit  
3. Add worst between-grid frequencies  
4. Restart from previous optimum  
5. Never change frequencies mid line-search  

## Answer template for Jesse (updated with data)

> The SF-optimized design already fails broadband isolation on a 200 MHz audit
> of 4–6 GHz (worst iso ≈ −24 dB near 5.8 GHz). So the question is not only
> “how many points for a good device,” but “how we sample while optimizing.”
> We freeze grid-sufficiency bars, use soft-min + mean on a fixed set, then
> adaptively enrich from dense audits. Soft-min alone does not prevent
> between-grid failures.

**Do not** present “11 points” as sufficient for a finished BB device.
