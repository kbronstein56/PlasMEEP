# FINAL_MULTIMODE_OPT_STATUS

Primary SF opt (C): `opt_uniform_C.json` (J=0.708 on C).
Best mesh-refined design: `opt_uniform_C_warm_M.json` (J=0.723 on M).

Campaign finish log: `campaign_finish.log` ended
`SF CAMPAIGN COMPLETE 2026-10-16T01:23:44-07:00`. Follow-on: power-balance fix + M refine + M→F.

## Gates

| Gate | Status | Notes |
|------|--------|-------|
| MULTIMODE_C6_OPT | **PASS** | C: J=0.708; M refine recovers to J=0.723 |
| OPTIMIZED_FULL12_ONSAGER | **PASS** | r_F = 1.98×10⁻¹³ on C |
| OPTIMIZED_POWER_BALANCE | **PASS** (post-fix) | Fixed `(ε−ε†)/(2i)` → residual 0.0046 on C |
| MULTIMODE_FINE_MESH | **EVALUATION_COMPLETE** | C/M/F evals + M refine + M→F done |
| MESH_CONVERGENCE | **FAIL** | M→F of M-opt: ΔJ=-0.053, Δiso=-1.9 dB, ΔRL=-2.7 dB (frozen bars). Better than C→F but not converged |
| FREE91_MULTIMODE | **FAILED** | IndexError cause recorded; deprioritized |
| C6_TIED_FD | **PASS** | |
| BROADBAND_READY | **CHANNEL_OK_4_6GHZ** | Advisor 4–6 GHz: TE0+TE1 only; no TE1-cutoff derivation |

## Why this file looked “missing” earlier

`finish_sf_campaign.sh` **did** write `FINAL_MULTIMODE_OPT_STATUS.md` at campaign end
(01:23). During earlier monitoring, M→F was still running and the file did not
exist yet; afterward it sat **uncommitted/unpushed** on disk until this cleanup.

## Primary design (optimized / reported on C)

- Start: `uniform` → best J **0.7082** (from −0.4876)
- T_des 0.804 / IL 0.95 dB / iso 19.0 dB / RL 16.6 dB / TE1 0.0187
- q ∈ [0.804, 1.207]; runtime ~5.2 h on C

## Mesh transfer

### C optimum transferred (no re-opt)

| | C | M | F |
|--|--:|--:|--:|
| J | 0.708 | 0.487 | 0.266 |
| T_des | 0.804 | 0.695 | 0.558 |
| TE1 | 0.019 | 0.022 | 0.044 |

### M-refined design (`opt_uniform_C_warm_M`)

| | M opt | same on F | Δ(F−M) |
|--|------:|----------:|-------:|
| J | 0.723 | 0.670 | −0.053 |
| T_des | 0.819 | 0.795 | −0.023 |
| iso dB | 20.30 | 18.41 | −1.90 |
| RL dB | 15.97 | 13.24 | −2.73 |
| TE1 | 0.0265 | 0.0203 | −0.006 |

See `M_REFINE_REPORT.md`, `M_TO_F_TRANSFER_REPORT.md`, `MESH_TRANSFER_CRITERION.md`.

## Power balance by mesh

| Mesh | Status | notes |
|------|--------|-------|
| C | **PASS** | max residual 0.0046 (fixed abs.); was 0.136 with Im(ε) bug |
| M (M-opt) | **PASS** | P_out+P_mat≈0.997 (rem≈0.0028) |
| F (M-opt) | **PASS** | P_out+P_mat≈0.997 (rem≈0.0027) |

## Free-91

Attempted after freq char; crashed:

```
IndexError: list index out of range
  in gradient_q_adjoint ... o_masks[m]
```

No `opt_free91_*.json` was written → earlier auto-status said `NOT_ATTEMPTED`.

## SINGLE_FREQ_VALIDATION summary

C6 opt + Onsager + absorption fix + M refine: **done**. Mesh continuum:
**not achieved** (M→F fails frozen bars on J/iso/RL). Free-91: failed (cause
recorded). Do **not** start long broadband optimization until mesh strategy
(interface refinement / finer reference) is decided.

## Recommended next action

Investigate local interface / plasma–quartz refinement or an even finer
reference mesh — not more optimizer iterations alone. Broadband channel count
is already OK for 4–6 GHz.
