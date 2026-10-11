# Fine-mesh validation and B-reversal

Design: best sixfold single-frequency result
(`opt_single_uniform_M_Bp0.0500.json`), re-evaluated on FULL91-F.

## Mesh transfer (B = +0.05 T, 3.85 GHz)

| metric | M (opt) | F (eval) | |Δ| |
|--------|--------:|---------:|----:|
| J | 5.380 | 5.373 | 0.007 |
| t_des avg | 4.358 | 4.355 | 0.003 |
| t_rev avg | 0.114 | 0.116 | 0.002 |
| t_thr avg | 0.022 | 0.022 | <1e-5 |
| isolation avg | 15.83 dB | 15.74 dB | 0.09 dB |

Well inside the frozen 0.10 dB port gate for the dominant observables.

## B-reversal (same material, FULL91-F)

| B | t_des (CCW) | t_rev (CW) | through | note |
|---|------------:|-----------:|--------:|------|
| +0.05 T | 4.355 | 0.116 | 0.022 | CCW wins |
| −0.05 T | 0.127 | 4.072 | 0.022 | CW wins (direction reverses) |
| 0 | 0.042 | 0.101 | 0.122 | no circulation |

**B_REVERSAL: PASS** (power-matrix test: at +B, P_des>P_rev for all sources;
at −B, P_rev>P_des for all sources).

## Local frequency sweep on F (+0.05 T)

Single-frequency design is narrowband: strong at 3.85 GHz, degraded by
3.90 GHz. Motivates the broadband stage.

| f (GHz) | t_des | t_rev | t_thr | iso (dB) |
|--------:|------:|------:|------:|---------:|
| 3.70 | 1.20 | 0.022 | 0.155 | 17.4 |
| 3.75 | 0.84 | 0.068 | 0.009 | 10.9 |
| 3.80 | 1.71 | 0.026 | 0.237 | 18.1 |
| 3.85 | 4.36 | 0.116 | 0.022 | 15.7 |
| 3.90 | 0.19 | 0.269 | 0.520 | −1.5 |
| 3.95 | 0.94 | 0.222 | 0.644 | 6.3 |
| 4.00 | 2.15 | 0.036 | 0.009 | 17.8 |
