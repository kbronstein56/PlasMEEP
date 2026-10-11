# Broadband sixfold optimization

Objective frozen in `CIRCULATOR_BROADBAND_OBJECTIVE.md`:

    J_bb = mean_n J(f_n) + 0.5 min_n J(f_n)

on `f ∈ {3.75…3.95} GHz`, warm-started from the single-frequency design.
Search mesh: FULL91-C (speed). Validation: FULL91-M band eval + F spot.

## Optimization (C)

| | start (warm) | best |
|--|-------------:|-----:|
| J_bb | 2.463 | **3.811** |
| J_mean | — | 3.214 |
| J_min (worst f) | — | 1.193 |

Wall ~2.0 h, 12 L-BFGS iterations.

## Band comparison on FULL91-M (same frozen P_ref)

| f (GHz) | single-freq t_des / iso_dB | broadband t_des / iso_dB |
|--------:|---------------------------:|-------------------------:|
| 3.75 | 1.77 / 14.1 | 2.89 / 12.8 |
| 3.80 | 1.60 / 17.8 | 2.73 / 24.0 |
| 3.85 | **4.36 / 15.8** | 4.17 / 17.2 |
| 3.90 | 0.23 / −0.4 | **0.65 / 8.3** |
| 3.95 | 1.12 / 8.0 | **1.64 / 13.4** |

| score on M | single-freq | broadband |
|------------|------------:|----------:|
| J_bb | 1.701 | **3.602** |
| J_min | −0.438 | **1.185** |

Broadband trades a little peak `t_des` at 3.85 GHz for a large rescue of the
3.90 GHz band edge.

FULL91-F at 3.85 GHz for the broadband design: `J=5.188`, `t_des=4.163`
(close to M).

## Artifact

`opt_broadband_C_Bp0.0500.json`, `broadband_vs_single_band_eval.json`
