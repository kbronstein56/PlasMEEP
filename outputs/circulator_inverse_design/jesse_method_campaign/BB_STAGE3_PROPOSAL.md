# Stage-3 proposal — **AUTHORIZED / LAUNCHED**

Generated after Stage-2 dense 50 MHz + targeted 25 MHz audits.  
**Authorization:** 2026-11-07 — launch with `n_workers=12`, `omp_threads_per_worker=1`.

## Why needed

Stage-2 fixed the Stage-1 4.1 GHz transmission corpse and improved \(J_{BB}\)
(−11.33 → −6.23). Dense audit then found a **new between-grid hole**:

| f (GHz) | A2 | IL (dB) | iso (dB) | in Stage-2 obj? |
|--------:|---:|--------:|---------:|:---------------:|
| 5.65 | −7.40 | 18.5 | −11.9 | yes (soft-min focus) |
| **5.75** | **−10.98** | **26.0** | **−20.1** | **no** |
| **5.825** | −10.71 | 24.4 | **−21.3** | **no** |

Same sampling failure mode as Stage-1’s 5.625 GHz discovery.

## Exact package

| Item | Value |
|------|-------|
| Objective | A2_hybrid |
| τ | 0.05 |
| α | 0.25 |
| Symmetry | C6 (16 vars) |
| Mesh | FULL91-C |
| Bounds | [0.05, 1.8] |
| Start design | `opt_broadband_stage2_A2_C.json` `q_best` |
| Frequencies | Stage-2 28 + **5.75 + 5.825** = **30** |
| maxiter | 30 L-BFGS |
| Parallelism | `--n-workers 12 --omp-threads-per-worker 1` |
| Expected wall | ~16–17 h (measured 5.06× vs serial ~84 h) |

### Exact frequency set (GHz)

```
4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.8, 4.85, 4.9,
5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.6, 5.625, 5.65, 5.7,
5.75, 5.8, 5.825, 5.85, 5.9, 5.95, 6.0
```

## Scientific question

Can soft-min A2 on the 30-frequency set close the 5.75/5.825 isolation/A2
hole without undoing the 4.1 GHz recovery (T_des≈0.14, IL≈8.5 dB, iso≈+7 dB)?

## Explicit non-goals

- Does not change A2 / τ / α / C6 mid-run.
- Does not loosen criteria in response to poor performance.
- Does not claim Stage-2 continuum physics; F/Q/X2 still required later.
- Notes residual **on-grid** 4.3 GHz IL≈27.5 dB as a separate pressure target.

## Runner / artifacts

- Runner: `production_broadband_opt_stage3.py`
- Pipeline: `run_stage3_authorized.sh`
- Checkpoint: `bb_stage3_checkpoint.json`
- History: `bb_stage3_history.json`
- Final: `opt_broadband_stage3_A2_C.json`
- Post-audit: `audit_bb_stage3_design.py` → `bb_stage3_audit_50MHz.json`, `BB_STAGE3_AUDIT_REPORT.json`
