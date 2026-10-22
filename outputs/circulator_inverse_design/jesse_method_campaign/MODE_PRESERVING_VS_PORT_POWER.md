# Mode-preserving vs port-power objectives

## Frozen weights (before seeing results)

See `MODE_OBJECTIVE_WEIGHTS_FROZEN.json`.

**Objective A — mode-preserving:** reward desired TE0; penalize all outgoing TE1.  
**Objective B — port-power:** reward desired TE0+TE1; no TE1-as-TE1 penalty.

## Short controlled opts (FULL91-C, maxiter 5, start `opt_uniform_C_warm_M`)

| | Objective A | Objective B |
|--|------------:|------------:|
| J_A | **0.727** | 0.699 |
| J_B | 0.755 | **0.755** |
| des0 | 0.825 | 0.830 |
| des1 | 0.0079 | 0.0177 |
| TE1 / des_tot | 0.9% | 2.1% |
| iso TE0 (dB) | **20.2** | 18.0 |
| iso tot (dB) | 19.5 | 16.9 |

Artifacts: `opt_mode_A.json`, `opt_mode_B.json`, `mode_ab_short_opt.json`.

## Answers (SF, short budget)

1. **Sacrifice from TE0 preservation:** small at this SF point — A wins J_A by ~0.03; des1 stays <1%.
2. **Does B exploit TE1?** Mildly (2% vs 1%); not a mode-converter redesign in 5 iters.
3. **Broadband?** Neither SF short-opt addresses 4–6 GHz (spectral study shows SF design fails BB).
4. **Experimentally meaningful?** B aligns with total VNA power; A needs mode resolution to verify TE0.
5. **B reversal:** pending dedicated ±B check on both finals.
6. **Mode stability vs f:** needs spectral of A vs B designs (des0/des1 plots).

## Recommendation

For production BB freeze: **default Objective A** (mode-preserving) unless a longer B
campaign shows clear broadband / experimental benefit. Do not claim TE0 preservation
is free forever — only that on this SF warm start it is cheap.
