# BROADBAND_PRODUCTION_SPEC_FROZEN

Frozen **2026-10-22** after Jesse method campaign studies.
Long 4–6 GHz optimization may start **only after** this file exists.
This freeze does **not** itself launch that optimization.

**Historical stage-1 freeze — do not overwrite.** Post-stage-1 objective
realignment lives in `BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md`
(Objective A2_hybrid + 28-freq enrichment). Stage-2 is also not auto-launched.

## Frozen choices

| Item | Freeze |
|------|--------|
| Symmetry | **C6** (16 vars). C3 orbits validated but short-opt lost (0.444 vs 0.727). |
| Incident mode | **TE0** drive on each of 6 ports (keep all six sources on FULL91-C). |
| Desired-output definition | **Mode-preserving Objective A**: reward desired TE0; penalize all TE1 + reverse/through/other/reflection TE0. Weights: `W_DES=1, W_REV=2, W_THR=1.5, W_LEAK=1, W_REFL=0.25, W_TE1=1`. |
| Broadband objective | \(J_{BB}=\alpha\,\mathrm{mean}(J_f)+(1-\alpha)\,\mathrm{softmin}_\tau(J_f)\) |
| τ | **0.05** |
| α | **0.25** |
| Initial frequency set | **100 MHz** spacing on [4.0, 6.0] GHz → **21** points. (200 MHz fails frozen bars by 8.9 dB iso.) |
| Adaptive enrichment | After each opt stage: dense ≤50 MHz audit; add worst-iso / worst-IL frequencies (and neighbors at 25 MHz if feature width demands); **restart** L-BFGS; never change set mid line-search. |
| Optimization mesh | FULL91-**C** |
| Continuum validation | Freeze design; evaluate worst + representative freqs on F/Q/X2 (or current continuum ladder); quote continuum brackets, not M-only J. |
| Gradient method | Efficient frequency-domain discrete adjoint (`A^\top`, 1 LU + 6+6 RHS). No JAX production port. |
| Time-domain | Meep cross-check per `TIME_DOMAIN_BROADBAND_VALIDATION.md` after gyrotropic gates — independent, not ground truth. |

## Evidence pointers

| Study | Artifact |
|-------|----------|
| Forward vs adjoint | `forward_mode_bench.json` (agree 2e-14; prefer adjoint) |
| Sampling | `spectral_*MHz.json`, `spectral_grid_comparison.json` (200→100→50 all FAIL bars) |
| Soft-min E2E | `softmin_e2e_fd.json` PASS |
| Mode A vs B | `opt_mode_A.json` / `opt_mode_B.json` → freeze A |
| C3 vs C6 | `c3_c6_short_opt.json` → freeze C6 |
| Literature | `docs/circulator_literature/` |

## Post-freeze 25 MHz audit (2026-10-22)

Targeted 25 MHz scan on 5.50–6.00 GHz vs band-restricted 50 MHz:
worst iso and worst IL both remain at **5.75 GHz** with **Δ = 0** under frozen
bars (`50_vs_25_band_5p50_6p00` PASS). **No change** to C6, Objective A, τ, α,
or the initial 100 MHz objective grid.

## Explicit non-goals of this freeze

- Does not reopen DtN / gyrotropic absorption validation.
- Does not authorize quoting M-mesh J=0.723 as physics.
- Does not claim 21 frequencies are final validation density.
