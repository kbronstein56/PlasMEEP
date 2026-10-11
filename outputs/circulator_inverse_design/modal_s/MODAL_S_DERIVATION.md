# Modal six-port scattering formulation (FEM)

Branch: `agent/eigenmode-ports`

## Audit of existing machinery

| Piece | Location | Finding |
|-------|----------|---------|
| Horn geometry | `sixport_common` | Source more outward than monitor: `src·n̂ − mon·n̂ = +1.2a`. Monitor sits between source and device. |
| Feed mode | `full91_solve.port_lines` | Flat Hz fundamental (`weights` = trapezoid ds). Matches “Hz PEC mode is flat”. |
| Injection | `full91_solve.inject` / circulator `rhs` | Consistent nodal load of flat profile on `source_xy`. |
| Monitor power | `circulator_core.port_power_one` / `guide_normal_flux` | Net ∫ S·n̂ only. At the driven port this mixes incident and reflected. |
| Pref campaign | `metrics_from_power` | `t = P/P_ref` with frozen uniform accepted-power `P_ref`. `insertion_avg_dB = 10 log10(t_des)` is **not** conventional insertion loss. |

## Definitions

Port outward normal `n̂` points away from the device.

At a vacuum-feed monitor:

- **Incident** amplitude `a`: wave traveling `−n̂` (into the device)
- **Outgoing** amplitude `b`: wave traveling `+n̂` (out of the device)

Flat L2-normalized mode `φ = 1/√W` on the monitor span `W`:

```
Hz     = (a_f + b_f) φ
∂_n Hz = i β (b_f − a_f) φ ,   β = k0
```

Power-normalized amplitudes (`e^{-iωt}`, `c = 1`):

```
â = a_f √(β/(2ω)) ,   ˆb = b_f √(β/(2ω))
```

so that `|â|²` = incident modal power, `|ˆb|²` = outgoing modal power, and

```
∫ S·n̂ ds = |ˆb|² − |â|²
```

for a pure bidirectional fundamental (verified synthetically and on the straight guide).

Scattering matrix for a fixed unit incident at source port `j` (other external incidents ≈ 0):

```
S_ij = ˆb_i / â_j ,   â_j ≡ A[j,j]
```

Reflection at the driven port uses `b_j/a_j`, **not** `−P_jj/P_ref`.

## Conventional circulator metrics (CCW desired at +B)

| Quantity | Definition |
|----------|------------|
| Desired transmission | `T_des = \|S_{j+1,j}\|²` |
| Insertion loss | `IL = −10 log10(T_des)` |
| Reverse isolation | `10 log10(T_des / T_rev)` |
| Return loss | `RL = −10 log10(\|S_jj\|²)` |
| Through leakage | `T_thr = \|S_{j+3,j}\|²` |
| Other-port leakage | `T_leak = \|S_{j+2,j}\|² + \|S_{j+4,j}\|²` |
| Absorbed | `P_abs = \|a_j\|² − Σ_i \|b_i\|²` |

## Correction of prior wording

In the Pref campaign, `insertion_avg_dB = 10 log10(P_des/P_ref)` measured **gain relative to a frozen accepted-power scale**, not insertion loss. Optimized `t_des ≈ 4.36` means the structure accepts/transmits ~4.36 `P_ref`, not `\|S21\|² ≈ 4.36`. That label is retired for modal-S reporting.

## Validation gates

See `validate_modal_s.py` and `VALIDATION_SUMMARY.json`.

## Code

- `scripts/circulator_inverse_design/modal_s_core.py`
- `scripts/circulator_inverse_design/straight_guide_fem.py`
- `scripts/circulator_inverse_design/validate_modal_s.py`
- `scripts/circulator_inverse_design/evaluate_modal_s_designs.py`
