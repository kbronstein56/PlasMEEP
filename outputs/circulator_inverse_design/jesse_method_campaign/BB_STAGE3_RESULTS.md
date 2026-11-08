# Stage-3 results (AUTHORIZED run complete)

**Branch:** `agent/eigenmode-ports`  
**Runner:** `production_broadband_opt_stage3.py` (`n_workers=12`, `omp=1`)  
**Warm start:** `opt_broadband_stage2_A2_C.json` `q_best`  
**Grid:** 30 freqs = Stage-2 28 + **5.75** + **5.825**  
**Objective:** A2_hybrid, τ=0.05, α=0.25, C6, FULL91-C, maxiter=30

## Optimization summary

| Quantity | Value |
|----------|------:|
| Start \(J_{BB}\) (30-freq softmin) | −9.019 |
| Final \(J_{BB}\) | **−3.659** |
| Accepted L-BFGS iters | 30 / 30 |
| BB evals | 54 |
| Wall time | **39 553 s ≈ 11.0 h** |
| Final on-objective worst A2 | **−4.572 @ 4.3 GHz** |

## Did the 5.75 / 5.825 holes get repaired?

**Yes — both enrichment frequencies improved strongly** (exact 25 MHz audit points, same as Stage-2 hole report):

| f (GHz) | Stage-2 A2 | Stage-3 A2 | Stage-2 iso (dB) | Stage-3 iso (dB) | Stage-2 IL (dB) | Stage-3 IL (dB) |
|--------:|----------:|----------:|-----------------:|-----------------:|----------------:|----------------:|
| **5.75** | −10.98 | **−2.95** | −20.1 | **+3.1** | 26.0 | **14.5** |
| **5.825** | −10.71 | **+1.38** | −21.3 | **+17.9** | 24.4 | **10.1** |

ΔA2 ≈ **+8.0** at 5.75 and **+12.1** at 5.825; isolation flipped from deep reverse to forward sense at both nodes.

## 4.1 GHz recovery — preserved and improved

| Metric | Stage-2 | Stage-3 |
|--------|--------:|--------:|
| \(T_{des}\) | 0.141 | **0.289** |
| IL (dB) | 8.51 | **5.39** |
| iso (dB) | +7.07 | **+9.97** |
| A2 | −0.84 | **+0.41** |

## New broadband worst (50 MHz dense audit)

| Stage | Worst A2 | f (GHz) | IL (dB) | iso (dB) |
|-------|----------|--------:|--------:|---------:|
| Stage-2 | −10.98 | **5.75** | 26.0 | −20.1 |
| Stage-3 | **−4.57** | **4.3** | 15.8 | −0.66 |

Stage-3 moves the continuum bottleneck from the mid-band isolation corpse back to the known **on-grid 4.3 GHz** pressure (already noted in the Stage-3 proposal as FoM-visible). Local 25 MHz audit also shows a mild off-node dip at **5.775 GHz** (A2≈−4.73, iso≈−2.8 dB) — much milder than Stage-2’s −20 dB holes.

Other Stage-2 pain points on the 50 MHz grid also improved, e.g. **5.65 GHz** A2 −7.40→−0.82, iso −11.9→+10.1 dB.

## Artifacts

| File | Role |
|------|------|
| `opt_broadband_stage3_A2_C.json` | Final Stage-3 design |
| `bb_stage3_checkpoint.json` / `bb_stage3_history.json` | Recoverable history |
| `bb_stage3_audit_50MHz.json` | Dense 50 MHz audit |
| `bb_stage3_audit_25MHz_local.json` | 5.50–6.00 GHz @ 25 MHz |
| `BB_STAGE3_AUDIT_REPORT.json` | Machine-readable S2 vs S3 + holes |
| `bb_stage3_nohup.out` | Full run log (untracked) |

## Scientific takeaway

Soft-min A2 on the enriched 30-frequency set **closed the Stage-2 between-grid isolation holes at 5.75 and 5.825 GHz without undoing (and while improving) the 4.1 GHz recovery**. Residual risk concentrates at **4.3 GHz** (on-grid) and a mild local notch near **5.775 GHz**.
