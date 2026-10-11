# In-silico circulator inverse design — status

Branch: `agent/eigenmode-ports`

## Gates

| Gate | Result |
|------|--------|
| PHYSICAL_PORT_OBJECTIVE_GRADIENT | **PASS** |
| SINGLE_FREQ_CIRCULATOR_OPT | **PASS** |
| FINE_MESH_DESIGN_VALIDATION | **PASS** |
| B_REVERSAL | **PASS** |
| BROADBAND_OPT | **PASS** |
| 91_VARIABLE_REFINEMENT | **PASS** (small gain; near-C6 retained) |
| MODAL_SIXPORT_S_VALIDATED | **PASS** (see `modal_s/`) |

**Terminology:** Pref-campaign `insertion_avg_dB = 10 log10(P/P_ref)` is **not**
conventional insertion loss. True modal IL for the sixfold peak is ≈1.3 dB
(`T_des≈0.74`); see `modal_s/MODAL_S_DESIGN_EVAL_REPORT.md`.

## Primary peak design (sixfold, 3.85 GHz, +0.05 T)

Source: `opt_single_uniform_M_Bp0.0500.json` (FULL91-M opt → FULL91-F validated).

| metric | initial `s≡1` | optimized |
|--------|--------------:|----------:|
| J | −0.679 | **5.380** |
| t_des avg | 0.042 | **4.358** (F: 4.355) |
| t_rev avg | 0.016 | 0.114 |
| t_thr avg | 0.550 | **0.022** |
| isolation avg | 4.1 dB | **15.8 dB** |
| worst-port t_des | 0.042 | 4.356 |
| q range | 1 | [0.457, 1.404] |

`t = P/P_ref` with frozen `P_ref = 1.369×10⁻²` (uniform accepted-power scale).
Optimizer uses linear powers; dB are report-only.

## Broadband candidate

`opt_broadband_C_Bp0.0500.json`, validated on M:

| score | single-freq design | broadband design |
|-------|-------------------:|-----------------:|
| J_bb on M | 1.70 | **3.60** |
| worst-f J | −0.44 (3.90 GHz) | **1.19** |

## Free-91 refinement

From the sixfold peak design: J 5.380 → **5.476**, desired_std across
sources still ~0.003 (rotational equivalence essentially kept). Prefer the
simpler sixfold design unless the tiny gain is needed.

## Multi-start

Perturbed C6 start on C reached J≈5.04 (C mesh), below the uniform-start
M design. Uniform start retained.

## Robustness (peak design, M)

Sensitive to global density scale (±5%) and to frequency/B shifts (narrowband).
2% random per-rod noise: isolation stays ~14–15 dB.

## Runtime / memory (order of magnitude)

- Gradient validation M: ~52 min
- Single-freq opt M: ~3.2 h
- F eval + sweep: ~23 min
- Broadband opt C: ~2.0 h
- Free-91 refine M: ~1.9 h

## Commits (this campaign)

See `git log` from `4ffc851` onward on `agent/eigenmode-ports`.
