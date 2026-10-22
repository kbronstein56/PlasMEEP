# Mode-preserving vs port-power objectives

## Frozen weights (before seeing results)

See `MODE_OBJECTIVE_WEIGHTS_FROZEN.json`.

**Objective A — mode-preserving (production-equivalent)**

- Good: desired-port **TE0**
- Bad: reverse/through/other TE0, reflection TE0, **all outgoing TE1**
- Weights: `W_DES=1, W_REV=2, W_THR=1.5, W_LEAK=1, W_REFL=0.25, W_TE1=1`

**Objective B — port-power circulator**

- Good: **total** desired power \(T_{\mathrm{des,TE0}}+T_{\mathrm{des,TE1}}\)
- Bad: total reverse / through / other / reflected (TE0+TE1)
- **No** TE1-as-TE1 penalty at the desired port
- Same port-role weights; `W_TE1_B=0`

## Status

| Step | Status |
|------|--------|
| Weights frozen | Done |
| Baseline eval on existing designs | Running / see `mode_baselines.json` |
| Short equal-budget opts A vs B | Running / see `opt_mode_A.json`, `opt_mode_B.json` |
| Spectra (des0/des1/rev0/rev1) | From baselines + post-opt |
| B-reversal check | After short opts |

## Questions to answer (fill from JSON)

1. How much performance is sacrificed by TE0 preservation?
2. Does port-power naturally exploit TE1?
3. Broadband implications (needs spectral audit)?
4. Experimental meaning (see experimental mode memo)?
5. Does B reversal still reverse circulation?
6. Is output mode stable vs frequency?

## Interim recommendation

Treat TE0 preservation as a **design choice**, not physics law (Zhao 2025 /
Wang 2021). Freeze production BB objective only after A/B short-opt numbers exist.
