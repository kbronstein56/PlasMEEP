# Stage-4 proposal — **NOT AUTHORIZED / DO NOT LAUNCH**

Generated from Stage-3 full-band 25 MHz + 12.5 MHz refine audits  
(`BB_STAGE3_DENSE_AUDIT.md`).

## Why needed

Stage-3 repaired the *exact* enrichment nodes 5.75 and 5.825 GHz and improved 4.1 GHz.  
Dense continuum sampling then found **new between-grid corpses**, worst at:

| f (GHz) | A2 | IL (dB) | iso (dB) | relative to Stage-3 obj grid |
|--------:|---:|--------:|---------:|------------------------------|
| **5.8125** | **−12.02** | **32.0** | **−18.7** | between 5.800 and 5.825 |
| **4.775** | −8.59 | 25.9 | −8.4 | off 50/100 MHz historical grids |
| 5.5875 | −8.53 | 24.0 | −10.5 | mid-band |
| 5.7125 | −8.10 | 27.7 | −6.2 | near 5.75 neighborhood |

Same sampling failure mode as Stage-1→2 and Stage-2→3.

## Proposed package (proposal only)

| Item | Value |
|------|-------|
| Objective | A2_hybrid |
| τ | 0.05 |
| α | 0.25 |
| Symmetry | C6 (16 vars) |
| Mesh | FULL91-C |
| Bounds | [0.05, 1.8] |
| Start design | `opt_broadband_stage3_A2_C.json` `q_best` |
| Frequencies | Stage-3 30 + **5.8125, 4.775, 5.5875, 5.7125** (= **34**) |
| maxiter | 30 L-BFGS |
| Parallelism | `--n-workers 12 --omp-threads-per-worker 1` (verified Stage-3 config) |

### Exact proposed frequency set (GHz)

```
4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.775, 4.8, 4.85, 4.9,
5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.5875, 5.6, 5.625, 5.65, 5.7,
5.7125, 5.75, 5.8, 5.8125, 5.825, 5.85, 5.9, 5.95, 6.0
```

(If a leaner Stage-4 is preferred: add only **5.8125** and **4.775** first.)

## Scientific question

Can soft-min A2 on the enriched set close the **5.8125** and **4.775** continuum corpses without undoing Stage-3’s 5.75/5.825/4.1 gains?

## Explicit non-goals / constraints

- Does **not** change A2 / τ / α / C6 mid-run.
- Does **not** loosen criteria in response to poor performance.
- Does **not** claim continuum physics; F/Q/X2 still required later.
- **Does not authorize launch.** Wait for explicit user authorization after reviewing this proposal and `BB_STAGE3_DENSE_AUDIT.md`.

## Expected cost

Roughly Stage-3 wall × (34/30) ≈ **~12–13 h** at `n_workers=12`, assuming similar L-BFGS eval count.

## Status

**PROPOSAL ONLY. Stage-4 NOT STARTED.**
