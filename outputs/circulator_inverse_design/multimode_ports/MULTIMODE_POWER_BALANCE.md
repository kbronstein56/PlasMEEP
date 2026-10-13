# MULTIMODE_POWER_BALANCE

**Status: PASS** (projected flat-port measure; bar 5%)

Design `opt_single_uniform_M_Bp0.0500`, grade C, B = +0.05 T, f = 3.85 GHz.

## Gate

For each of 12 unit-power incident channels:

\[
|1 - \sum_{mn}|S_{mn}|^2 - P_{\mathrm{material}}| \le 0.05
\]

with `P_material ≥ −0.02`.

## Result (frozen measure = projected)

| Drive | mean Σ\|S\|² | mean TE0 | mean TE1 | mean P_mat | mean residual proxy |
|-------|-------------|----------|----------|------------|---------------------|
| TE0   | 0.9835 | 0.6659 | 0.3176 | 0.0347 | ≤ 0.0216 |
| TE1   | 0.9911 | 0.3195 | 0.6716 | 0.0324 | ≤ 0.0377 |

Worst channel residual: **0.0377** (ch3, port1 TE1).
Former arclength failure ch7 (port3 TE1): **0.0673 → 0.0352**.

`MULTIMODE_POWER_BALANCE = PASS`

See `TE1_POWER_BALANCE_DIAGNOSTIC.md` and `PROJECTED_MEASURE_EXPERIMENT.md`.
