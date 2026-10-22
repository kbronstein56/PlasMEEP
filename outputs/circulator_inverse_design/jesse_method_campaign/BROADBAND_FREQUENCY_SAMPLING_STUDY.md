# Broadband frequency sampling study (4–6 GHz)

## Frozen grid-sufficiency bars (before results)

`GRID_SUFFICIENCY_BARS_FROZEN.json`: Δiso < 0.5 dB, ΔIL/RL < 0.2 dB, location ≤100 MHz.

## Full-band reconnaissance (SF design, FULL91-C)

| Grid | Worst iso | @ f | Worst IL | Pass vs finer? |
|------|----------:|-----:|---------:|:--------------:|
| 200 MHz (11) | −24.2 dB | 5.8 | 40.8 dB | **FAIL** vs 100 |
| 100 MHz (21) | −33.1 dB | 5.9 | 43.5 dB | **FAIL** vs 50 |
| 50 MHz (41) | **−52.8 dB** | **5.75** | **58.1 dB** | see 25 MHz band |

200→100: Δiso **8.9 dB**. 100→50: Δiso **19.8 dB**.

## Targeted 25 MHz refinement — 5.50–6.00 GHz

Artifact: `spectral_25MHz_5p50_6p00.json` (21 points, 25 MHz spacing; midpoints
computed fresh, 50 MHz nodes reused).

| Source | Worst iso | @ f | Worst IL | @ f |
|--------|----------:|-----:|---------:|-----:|
| 50 MHz restricted to band | −52.8287 dB | **5.75** | 58.1160 dB | **5.75** |
| 25 MHz band scan | −52.8287 dB | **5.75** | 58.1160 dB | **5.75** |

Nearby 25 MHz midpoints (e.g. 5.725, 5.775) are **less** extreme than 5.75
(next-worst iso ≈ −31 dB). Feature peak is resolved on the 50 MHz node.

### 50 → 25 band comparison vs frozen bars

From `spectral_grid_comparison.json` → `50_vs_25_band_5p50_6p00`:

| Check | Δ | Pass? |
|-------|--:|:-----:|
| Worst isolation | **0.0 dB** | **PASS** (<0.5) |
| Worst IL | **0.0 dB** | **PASS** (<0.2) |
| Worst-iso location | **0.0 GHz** | **PASS** (≤0.1) |
| Worst-IL location | **0.0 GHz** | **PASS** (≤0.1) |
| **pass_all** | — | **PASS** |

**Verdict:** In the critical upper band, worst-isolation and worst-IL frequencies
have **converged** under the frozen sampling bars at 50 MHz relative to 25 MHz.
No concrete reason to change the frozen production objective, C6, Objective A,
τ/α, or the initial **100 MHz** optimization grid. Keep adaptive enrichment
(≤50 MHz audit + 25 MHz locally when features appear) as already frozen.

## Answer for Jesse

> Full-band coarse grids (200/100 MHz) are not validation-sufficient. The
> catastrophic hole near 5.75 GHz is captured at 50 MHz and unchanged at 25 MHz
> (Δiso = ΔIL = 0). Production starts with a **100 MHz objective set** and
> adaptively enriches from dense audits — never mid line-search.

Freeze unchanged: `BROADBAND_PRODUCTION_SPEC_FROZEN.md`.
