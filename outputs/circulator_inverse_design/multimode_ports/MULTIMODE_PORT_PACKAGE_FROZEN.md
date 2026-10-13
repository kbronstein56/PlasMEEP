# MULTIMODE_PORT_PACKAGE_FROZEN

**Frozen after:** DIRECT_DTN_FULL12_ONSAGER = PASS and MULTIMODE_POWER_BALANCE = PASS.

Do not make silent changes to the items below without a new validation campaign.

## Boundary geometry

- Outer DtN planes from `port_outer_center` + exterior excision beyond ports.
- True physical-boundary edges via `build_port_from_boundary` (`plane_tol=0.12`).
- Clear width `W = clear_width = 2a = 2.0` (unchanged physical device).
- Staircase mesh cut accepted; **not** remeshed (Phase 4 skipped).

## Modal basis

- Parallel-plate Neumann modes: `φ_m(y) = cos(m π (y+W/2)/W)`, m = 0,1,…
- Propagating TE0/TE1 at 3.85 GHz; two evanescent Robin modes in `A`.

## Quadrature / measure (frozen)

- **`measure = projected`**: edge weight `|e · t̂|` (flat-port `dy`).
- Arclength `|e|` retained only for diagnostics (`--measure arclength`).
- Same `c_m` for Robin, incident RHS, and `Q` columns (reciprocal).

## Normalization

- Discrete L2: `∫ φ_m² dℓ = 1` with projected measure on port edges.
- Power waves: `â = a_f √(β/(2ω))`, unit-power drives in S columns.

## S definition

\[
S = I + 2i\,Q,\qquad
Q_{mn}=\sqrt{\beta_m\beta_n}\,c_m^\top A^{-1} c_n.
\]

Frozen FEM pair: Robin `(-iβ) c cᵀ`, RHS `+2iβ a_f c`.

## Power balance

- Identity: `|1 − Σ|S|² − P_material| ≤ 0.05` per channel.
- `P_material` from `material_absorption` (field integral); discrete `Ax=f` closes.

## Validated gates (this freeze)

| Gate | Result |
|------|--------|
| Straight-guide TE0/TE1 | PASS |
| DIRECT_DTN_FULL12_ONSAGER | PASS ~1e-14 |
| MULTIMODE_POWER_BALANCE | PASS (max residual 0.0377 on C) |
| TE1 mesh C→M ch7 | 0.0352 → 0.0253 (converging) |
| Absorption cross-check | discrete identity OK; field residual ≤ 5% |
| DIRECT_DTN_ADJOINT_FD | PASS (forward≡adjoint ~1e-13; FD rods 0,30,45) |

## Cause of original TE1 failure

Arclength staircase measure inflated TE1 modal loads on port 3 (ch7);
`Σ|S|²≈1.036` with residual 0.067. Projected flat-port measure fixes it
without breaking Onsager. Not a geometry or material-asymmetry bug.

## Implementation entry points

- `boundary_ports.py` — edges, projected measure, `modal_c_vector`
- `direct_dtn_s.py` — `S = I + 2i Q`
- `full91_direct_dtn.py` — FULL91 assembly / Onsager
- `power_balance_direct.py` — 12-channel balance
