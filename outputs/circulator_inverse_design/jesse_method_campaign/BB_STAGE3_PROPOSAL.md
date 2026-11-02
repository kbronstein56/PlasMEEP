# Stage-3 proposal (NOT authorized / NOT launched)

Generated after Stage-2 dense 50 MHz + targeted 25 MHz audits.

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
| Expected wall | ~80–85 h |

## Scientific question

Can soft-min A2 on the 30-frequency set close the 5.75/5.825 isolation/A2
hole without undoing the 4.1 GHz recovery (T_des≈0.14, IL≈8.5 dB, iso≈+7 dB)?

## Explicit non-goals of this proposal

- Does **not** authorize auto-start.
- Does not change A2 / τ / α / C6 mid-run.
- Does not claim Stage-2 continuum physics; F/Q/X2 still required later.
- Notes residual **on-grid** 4.3 GHz IL≈27.5 dB as a separate pressure target
  the FoM already sees — enrichment alone will not invent a cure if A2 tradeoffs
  prefer other frequencies.

## Go / no-go

Await human approval before launching `production_broadband_opt_stage2.py`
variant with the 30-freq set (or a dedicated Stage-3 runner).
