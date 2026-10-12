# Multimode port derivation (parallel-plate Hz feeds)

Time convention: **e^{-iωt}**. The validated FULL91 FEM solves the scalar Hz equation with natural Neumann conditions on ρ→0 PEC metal (no Hz Dirichlet on air–metal interfaces).

## Geometry

Production clear aperture width:

    W = clear_width = width_base − 2·wall_thickness = 2a  (40 mm clear)

Transverse coordinate on the port line: y ∈ [−W/2, W/2] with PEC at y = ±W/2.

## Transverse eigenmodes (Neumann Hz on PEC sidewalls)

    φ_m(y) = cos(m π (y + W/2) / W) ,   m = 0, 1, 2, …

    ∂φ_m/∂y|_{y=±W/2} = 0

Cutoff wavenumbers and longitudinal constants at angular frequency ω = k₀ (c=1):

    k_{c,m} = m π / W
    β_m = +√(k₀² − k_{c,m}²)     (propagating, real β_m ≥ 0)
    α_m = +√(k_{c,m}² − k₀²)     (evanescent; use exp(−α_m |Δx|) along guide)

### Cutoffs at 3.85 GHz (numerical)

With k₀ = 2π f a / c and W = 2:

| m | f_cut (GHz) | at 3.85 GHz |
|---|-------------|-------------|
| 0 | —           | propagating (β₀ = k₀) |
| 1 | ≈ 3.747     | propagating |
| 2 | ≈ 7.494     | evanescent  |

TE0 and TE1 are **open propagating channels**; TE2+ are evanescent at the design frequency.

## Normalization and orthogonality

L² normalization on the clear aperture (trapezoid quadrature on port samples):

    ∫ φ_m φ_n dy ≈ Σ_q w_q φ_m(y_q) φ_n(y_q) = δ_{mn}

Mode profiles used in code: `modes.normalize_on_line`.

## Bidirectional expansion on a port plane

At a port with outward normal **n̂** (pointing out of the device):

    H_z = Σ_m (a_{m,f} + b_{m,f}) φ_m

    ∂_n H_z = Σ_m i β_m (b_{m,f} − a_{m,f}) φ_m     (propagating; β_m real > 0)

**Sign convention (aligned with `modal_s_core.extract_ab`):**

- **Incident** modal amplitude a travels **−n̂** (into the device).
- **Outgoing** modal amplitude b travels **+n̂** (out of the device).

Power-normalized amplitudes (vacuum feed):

    â_m = a_{m,f} √(β_m / (2ω)),   ω = k₀

    P_inc,m = |â_m|² ,   P_out,m = |b̂_m|²

For propagating modes, net guide-normal Poynting flux through the port is Σ_m (|b̂_m|² − |â_m|²) when cross-mode coupling is negligible.

## Evanescent branch

For β_m = i α_m (α_m > 0), use the same (a,b) decomposition with β_m in the projection:

    a_f = ½(c_H − c_d/(i β_m)),   b_f = ½(c_H + c_d/(i β_m))

Evanescent modes improve local DtN accuracy but are **not** counted as propagating S-matrix ports unless explicitly listed in the 12×12 channel map.

## 12×12 channel ordering (Phase 4)

    index = 2·port + mode_index
    mode_index = 0 → TE0 (m=0), 1 → TE1 (m=1)

Columns of S_full are **unit power** on one incoming propagating mode; rows are outgoing power-normalized modes.
