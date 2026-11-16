# Stage-5 results — full-band 6.25 MHz A2 broadband

**Status:** COMPLETE (optimizer iteration limit)  
**Stage-6:** **NOT STARTED**  
**Branch:** `agent/eigenmode-ports`  
**Warm start:** Stage-4 `opt_broadband_stage4_A2_C.json` `q_best`  
**Kernel:** optimized A2 path (`PLASMEEP_A2_FAST_KERNEL=1`)

## Configuration (frozen)

| Item | Value |
|------|------:|
| Objective | A2_hybrid |
| τ / α | 0.05 / 0.25 |
| Bounds | [0.05, 1.8] |
| Design | C6 / 16 vars, FULL91-C |
| Frequency rule | `[round(4.0+i*0.00625,6) for i in range(321)]` |
| n_freq | **321** |
| max gap | **6.25 MHz** |
| n_workers / OMP | 12 / 1 |
| maxiter | 30 |

## Optimization outcome

| Quantity | Value |
|----------|------:|
| Final J_BB | **−6.055965** |
| Accepted iters | **30 / 30** |
| BB evaluations | **77** |
| Wall time | **76571 s ≈ 21.3 h** |
| Termination | L-BFGS-B iteration limit |
| Objective-grid worst | **5.1375 GHz**, A2 ≈ **−7.585** |
| IL / iso at obj worst | 25.78 dB / −4.71 dB |
| max softmin weight | 0.279 |
| N_eff | 7.43 |

Start (Stage-4 q on 321-grid): J_BB ≈ −8.609, worst A2 ≈ −11.04 @ 5.08125 GHz.

## Watch frequencies (final objective grid)

| f (GHz) | A2 | IL (dB) | iso (dB) |
|--------:|---:|--------:|---------:|
| 4.10000 | +0.145 | 5.83 | +9.29 |
| 4.73125 | −5.875 | 19.71 | −3.38 |
| 4.77500 | −2.446 | 16.70 | +8.35 |
| 5.08125 | −6.171 | 22.74 | −1.98 |
| 5.13750 | −7.585 | 25.78 | −4.71 |
| 5.81250 | −5.192 | 16.22 | −4.92 |

## Post-audit (mandatory)

See `BB_STAGE5_DENSE_AUDIT.md`. True worst under 3.125 MHz refine:

**4.734375 GHz, A2 ≈ −8.741, IL ≈ 25.04 dB, iso ≈ −11.03 dB**

(between 6.25 MHz nodes 4.73125 / 4.7375).

## Conclusion (preview)

**B — materially improved, but narrower residuals remain.**  
The 6.25 MHz full-band objective removed the Stage-2→4 “silent gap” class for many historical corpses and improved J_BB / worst A2, but 3.125 MHz auditing still finds A2 ≲ −8 features between objective nodes. Stage-6 is **not** authorized.
