# Phase 3 — Projected-measure experiment

## Question

Does replacing staircase arclength `ds = |e|` with the flat-port projected
measure `dℓ = |e · t̂|` fix TE1 power balance without breaking Onsager?

## Mathematical consistency with flat-port DtN

The continuous DtN / modal theory is posed on a **flat** waveguide cross-section
with transverse coordinate `y` and measure `dy`:

\[
(c_m)_p = \int_{-W/2}^{W/2} \phi_m(y)\, N_p(y)\, dy,
\qquad
\int_{-W/2}^{W/2} \phi_m\phi_n\, dy = \delta_{mn}.
\]

The FULL91 outer port cut is a **staircase** approximation to that plane.
Using Euclidean edge length `|e|` integrates along the jagged polyline, not
along the intended aperture. For TE0 (`φ₀` constant) the bias is mild; for TE1
(`φ₁(y)` varies across the aperture) it systematically mis-normalizes `c₁`
and over-counts `|S|²` on some mesh-registration draws.

The projected weight `|e · t̂|` is the absolute contribution of each edge to
the **nominal aperture coordinate** `y`. It is the discrete analogue of `dy`
on the same node set and the same `φ_m(y)` samples. Excitation and observation
use the **same** `c_m`, so the Robin / RHS / `Q` construction remains reciprocal
and Onsager-compatible.

**Conclusion:** adopting projected measure is consistent with the flat-port
derivation; arclength was a mesh-artifact measure, not the theory measure.

## Results (grade C, design `opt_single_uniform_M_Bp0.0500`)

| Gate | arclength | projected |
|------|-----------|-----------|
| Straight-guide TE0/TE1 | PASS | PASS |
| FULL12 Onsager `r_F` | ~9e-14 PASS | ~8e-14 PASS |
| Power ch7 (port3 TE1) residual | **0.0673 FAIL** | **0.0352 PASS** |
| Max 12-channel residual | 0.0673 | 0.0377 |
| MULTIMODE_POWER_BALANCE | FAIL | **PASS** |

All 12 channels ≤ 5% with projected measure. Onsager remains at numerical precision.

## Decision

**Adopt `measure=projected` as the frozen production default.**
Arclength retained only for A/B diagnostics (`--measure arclength`).
Phase 4 (straight-port remesh) is **not required** for the 5% gate.
