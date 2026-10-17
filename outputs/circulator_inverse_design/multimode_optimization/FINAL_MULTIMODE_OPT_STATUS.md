# FINAL_MULTIMODE_OPT_STATUS

Primary opt JSON: `outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json`

Campaign finish log: `campaign_finish.log` ended
`SF CAMPAIGN COMPLETE 2026-10-16T01:23:44-07:00`.

## Gates

| Gate | Status | Notes |
|------|--------|-------|
| MULTIMODE_C6_OPT | **PASS** | Best = uniform warm start, J=0.708 on C |
| OPTIMIZED_FULL12_ONSAGER | **PASS** | r_F = 1.98×10⁻¹³ on C |
| OPTIMIZED_POWER_BALANCE | **PASS** (post-fix) | Was FAIL 0.136 from buggy `Im(ε)` absorption; fixed to `(ε−ε†)/(2i)` → residual 0.0046 on C. See `OPTIMIZED_POWER_BALANCE_DIAGNOSTIC.md` |
| MULTIMODE_FINE_MESH | **EVALUATION_COMPLETE** | C, M, F evaluated successfully (≠ mesh-converged) |
| MESH_CONVERGENCE | **FAIL** | FoMs drift strongly C→M→F; continuum not achieved |
| FREE91_MULTIMODE | **FAILED** | `IndexError` in singleton-orbit masks; deprioritized until power/M/F done |
| C6_TIED_FD | **PASS** | |
| BROADBAND_READY | **CHANNEL_OK_4_6GHZ** | Advisor 4–6 GHz: TE0+TE1 only; no TE1-cutoff derivation. Prior cutoff warning was for 3.70–4.00 GHz char only. |

## Why this file looked “missing” earlier

`finish_sf_campaign.sh` **did** write `FINAL_MULTIMODE_OPT_STATUS.md` at campaign end
(01:23). During earlier monitoring, M→F was still running and the file did not
exist yet; afterward it sat **uncommitted/unpushed** on disk until this cleanup.

## Primary design (optimized / reported on C)

- Start: `uniform` → best J **0.7082** (from −0.4876)
- T_des 0.804 / IL 0.95 dB / iso 19.0 dB / RL 16.6 dB / TE1 0.0187
- q ∈ [0.804, 1.207]; runtime ~5.2 h on C

## Mesh transfer (existing results only)

| | C | M | F |
|--|--:|--:|--:|
| J | 0.708 | 0.487 | 0.266 |
| T_des | 0.804 | 0.695 | 0.558 |
| TE1 | 0.019 | 0.022 | 0.044 |

See `MULTIMODE_FINE_MESH_REPORT.md` for the full table.

## Power balance by mesh

| Mesh | Status | max residual |
|------|--------|-------------:|
| C | **FAIL** | 0.136 |
| M | **NOT_RUN** | — |
| F | **NOT_RUN** | — |

## Free-91

Attempted after freq char; crashed:

```
IndexError: list index out of range
  in gradient_q_adjoint ... o_masks[m]
```

No `opt_free91_*.json` was written → earlier auto-status said `NOT_ATTEMPTED`.

## SINGLE_FREQ_VALIDATION summary

**PARTIAL**: C6 opt + Onsager + fine-mesh *evaluations* succeeded; power balance
failed on C; FoMs drift strongly C→M→F; free-91 failed; no M/F power-balance.

## Recommended next action

Re-optimize (or warm-start refine) the multimode C6 design on **FULL91-M** under
the same frozen objective, then re-check power balance on M — do **not** treat
the C optimum as mesh-converged.
