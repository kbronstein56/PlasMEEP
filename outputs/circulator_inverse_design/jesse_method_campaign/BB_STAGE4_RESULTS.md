# Stage-4 broadband optimization — results

**Status:** complete (resumed after external interrupt at EVAL #67)  
**Do not launch Stage-5 automatically.**

## Optimization outcome

| Quantity | Value |
|----------|------:|
| Final `J_BB` (47-freq soft-min) | **−4.131096** |
| Start (Stage-3 `q_best` on Stage-4 grid) | −9.514531 |
| Accepted L-BFGS iterations | **30 / 30** |
| Completed BB objective evaluations | **74** |
| Wall time | **≈22.6 h** (incl. pre-interrupt) |
| Optimizer message | `STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT` |
| Warm start | Stage-3 `q_best` |
| Package | A2_hybrid, τ=0.05, α=0.25, C6, FULL91-C, n_workers=12, omp=1 |

Resume verification: log shows `RESUME Stage-4 from nit=29` and `remaining L-BFGS iters=1 (done=29)` before completion.

## Dense post-audit (forward-only)

| Scan | Worst f (GHz) | Worst A2 | IL (dB) | isolation (dB) |
|------|--------------:|---------:|--------:|---------------:|
| Full 4–6 @ 25 MHz | **4.725** | **−9.16** | 26.24 | −11.09 |
| Local 12.5 MHz | 4.725 | −9.16 | 26.24 | −11.09 |
| Local 6.25 MHz | **5.08125** | **−11.04** | 32.61 | −13.59 |

**True dense-audit worst (across available Stage-4 scans):**  
**5.08125 GHz**, A2 = **−11.04**, IL = **32.61 dB**, isolation = **−13.59 dB**.

Secondary serious hole: **4.725–4.73125 GHz** (A2 ≈ −9.2 to −10.2), adjacent to the repaired 4.775 node.

## Targeted Stage-3 → Stage-4 comparison

| f (GHz) | Stage-3 A2 / IL / iso | Stage-4 A2 / IL / iso | Assessment |
|--------:|-----------------------|-----------------------|------------|
| **5.8125** | −12.02 / 31.96 / −18.67 | −4.48 / 15.56 / −2.75 | **Improved, not fully repaired** (iso still negative) |
| **4.775** | −8.59 / 25.90 / −8.35 | −2.09 / 16.20 / **+9.37** | **Largely repaired** at the exact node |
| **5.5875** | −2.72 / 11.42 / +1.78 | −5.09 / 15.54 / −3.83 | **Regressed** (became objective-grid worst) |
| **4.725** | −2.44 / 11.43 / +3.21 | −9.16 / 26.24 / −11.09 | **New serious off-grid hole** |
| **5.08125** | (not a Stage-3 corpse) | −11.04 / 32.61 / −13.59 | **New sharp off-grid hole** |

4.1 GHz remains healthy (A2≈+0.03, IL≈5.93 dB, iso≈+8.9 dB, T_des≈0.255).

## Repair verdicts

- **5.8125:** improved substantially (ΔA2 ≈ +7.5) but **not fully repaired** — still A2≈−4.5, iso≈−2.7 dB.
- **4.775:** **repaired at the sampled node** (iso −8.4 → +9.4 dB), but the nearby **4.725** pedestal became a deep corpse.
- **Serious new off-grid holes:** **yes** — especially **5.08125 GHz** (true dense worst) and **4.725 GHz**. Soft-min again parked failures between objective nodes.

## Recommendation

**Broadband optimization is not yet sufficient.** Stage-4 fixed the authorized notch samples partially, but dense audit reveals new continuum holes off the 47-frequency grid. Do **not** start Stage-5 without a new localization + frequency-set proposal analogous to the Stage-3→Stage-4 process.

Artifacts:
- `opt_broadband_stage4_A2_C.json`
- `bb_stage4_checkpoint.json` (`status=complete`)
- `BB_STAGE4_DENSE_AUDIT_SUMMARY.json` / `.md`
- `bb_stage4_dense_{25MHz,12p5MHz_local,6p25MHz_local}.json`
