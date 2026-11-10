# Stage-4 — **AUTHORIZED** (2026-11-10)

Authorized after Stage-3 notch localization (`BB_STAGE3_NOTCH_LOCALIZATION.json`).
Use the FINAL 47-frequency set below. Do **not** add further frequencies unless a
concrete numerical issue appears during setup. Do **not** launch Stage-5 automatically.

## Localization summary (fixed Stage-3 `q_best`, no optimization)

### 5.8 GHz feature (5.775–5.850 @ 6.25 MHz)

| Quantity | Value |
|----------|------:|
| Notch center | **5.8125 GHz** |
| A2 / IL / iso | **−12.02 / 31.96 dB / −18.67 dB** |
| Deep core (A2 < −7) | 5.8125–5.81875 GHz (**≈6–12.5 MHz**) |
| Bad core (A2 < −5) | 5.80625–5.81875 GHz (**12.5 MHz**) |
| Local minima | **1 deep** @ 5.8125; 1 weak @ 5.8375 (A2≈−2.5, not a corpse) |
| Shoulders | 5.800 A2=+5.53; 5.825 A2=+1.38 |

3.125 MHz refine **not required**: 6.25 MHz already resolves a single deep minimum with depressed neighbors on both sides.

### 4.775 GHz feature (4.725–4.825 @ 6.25 MHz)

| Quantity | Value |
|----------|------:|
| Notch center | **4.775 GHz** |
| A2 / IL / iso | **−8.59 / 25.90 dB / −8.35 dB** |
| Deep core (A2 < −7) | single sample @ 4.775 |
| Pedestal (A2 < −5) | ~4.74–4.825 (**≈80–90 MHz** wide) |
| Structure | primary min @ 4.775 + weaker local mins @ 4.744 / 4.756 |

**12.5 MHz** sampling with flanks at 4.75 / 4.775 / 4.8 / 4.825 is appropriate; full-band 6.25 MHz is not justified.

## Why not 12.5 MHz alone through 5.8 GHz

A 12.5 MHz ladder with only `{5.800, 5.8125, 5.825}` leaves the measured bad samples **5.80625** and **5.81875** off-grid. Soft-min could again park a deep floor between objective nodes. Therefore Stage-4 uses **6.25 MHz only inside the measured 5.800–5.825 core**, and **12.5 MHz** across the broader 5.70–5.90 guard band.

## FINAL Stage-4 frequency set (GHz) — AUTHORIZED

```
4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.75, 4.775, 4.8, 4.825, 4.85, 4.9,
5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.5875, 5.6, 5.625, 5.65,
5.7, 5.7125, 5.725, 5.7375, 5.75, 5.7625, 5.775, 5.7875,
5.8, 5.80625, 5.8125, 5.81875, 5.825,
5.8375, 5.85, 5.8625, 5.875, 5.8875, 5.9, 5.95, 6.0
```

| Item | Value |
|------|------:|
| **n_freq** | **47** |
| 12-worker batches | **4** (`ceil(47/12)`) |
| Est. wall / BB eval | **~970 s ≈ 16 min** |
| Est. Stage-4 total | **~14–16 h** (~50–55 BB evals) |

### Unchanged package

| Item | Value |
|------|-------|
| Objective | A2_hybrid |
| τ | 0.05 |
| α | 0.25 |
| Bounds | [0.05, 1.8] |
| Symmetry | C6 |
| Mesh | FULL91-C |
| FEM / adjoint | unchanged |
| Optimizer | L-BFGS, maxiter=30 |
| Warm start | Stage-3 `q_best` |
| Parallelism | `n_workers=12`, `omp_threads_per_worker=1` |

## Status

**AUTHORIZED 2026-11-10. Launch via `run_stage4_authorized.sh`.**
After optimization: dense post-audit + Stage-3 comparison + viz. **No Stage-5.**
