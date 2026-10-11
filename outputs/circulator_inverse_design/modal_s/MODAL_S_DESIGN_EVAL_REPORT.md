# Existing designs under true modal S-parameters

Mesh: FULL91-M · f = 3.85 GHz · no re-optimization

## Terminology correction

| Pref-campaign label | True meaning | Conventional name |
|---------------------|--------------|-------------------|
| `t_des = P_des/P_ref` | Outward net port power over frozen uniform accepted-power scale | **not** \|S\|² |
| `insertion_avg_dB = 10 log10(t_des)` | Pref-relative gain | **not** insertion loss |
| — | `IL = −10 log10(\|S_des\|²)` | conventional insertion loss |

For the sixfold peak design, Pref reporting gave `t_des≈4.36` / “+6.4 dB insertion”
while the structure accepts ≈5.39 `P_ref`. Under modal S the same design has
**T_des≈0.74** and **IL≈1.3 dB**.

## +B = +0.05 T summary (FULL91-M)

| Design | T_des | IL (dB) | Iso (dB) | RL (dB) | T_jj | T_thr | T_leak | P_abs/P_inc | Pref t_des | Pref “ins.” (dB) | Strong? |
|--------|------:|--------:|---------:|--------:|-----:|------:|-------:|------------:|-----------:|-----------------:|:-------:|
| Sixfold peak | 0.742 | 1.30 | 16.0 | 12.4 | 0.058 | 0.0037 | 0.009 | 0.169 | 4.358 | +6.39 | **yes** |
| Broadband | 0.717 | 1.44 | 17.5 | 12.6 | 0.056 | 0.0022 | 0.020 | 0.192 | 4.174 | +6.21 | **yes** |
| Uniform s≡1 | 0.007 | 21.4 | 6.9 | 1.1 | 0.783 | 0.121 | 0.023 | 0.064 | 0.042 | −13.8 | no |

Strong-circulator criteria: `T_des>0.5`, iso>10 dB, `T_thr<0.1`, `T_jj<0.3`.

## Bias dependence (sixfold peak)

| B | T_des | Iso (dB) | T_jj | Notes |
|---|------:|---------:|-----:|-------|
| +0.05 T | 0.742 | +16.0 | 0.058 | CCW circulation |
| −0.05 T | 0.020 | −15.5 | 0.048 | Circulation reverses (desired↔reverse) |
| 0 | 0.005 | +1.6 | 0.803 | Reflective; no circulation |

Onsager on peak design: ‖S(+B)−Sᵀ(−B)‖_F/‖S‖_F ≈ 0.040 (formulation gate on uniform was 1.8×10⁻⁴).

## Judgment

The Pref inverse-design campaign remains a **successful proof**. Under true modal
normalization both the sixfold peak and broadband designs are **strong circulators**
(~1.3–1.4 dB IL, ~16–17.5 dB isolation, ~12.5 dB return loss, through <0.4%).
Uniform s≡1 is not.

Full complex 6×6 S matrices (+B, −B, B=0) are in
`eval_*_M.json` under this directory.

## Prior optimization preserved

No rods were re-optimized. Design vectors remain those in
`opt_single_uniform_M_Bp0.0500.json` and `opt_broadband_C_Bp0.0500.json`.
