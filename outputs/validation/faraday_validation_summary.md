# Faraday Benchmark Validation Record

## Status
**PASS** at resolution 32 and 64.

## Fixes applied (validation script only)
- π-periodic Stokes unwrap: `ψ = ½ unwrap(atan2(S2, S1))`
- Pre-flight propagation check for B=0, +B, and −B
- PASS gates use |κ| and B-reversal only; absolute sign reported separately
- Tightened gates: B=0 ≤ max(0.002, 3% scale); ±B asymmetry ≤ 10%; theory |κ| error ≤ 15%; R² / RMSE fit gates
- DFT z monotonicity + length match checks
- `sim.run(until=...)` retained; DFT via `add_dft_fields([Ex,Ey], [fs], ...)`

## Defaults
- a = 0.028 m
- fs = 5.0 GHz, fp = 2.0 GHz, γ = 1 MHz, |B| = 0.05 T
- Theory |κ| = 0.078460 rad/a

## Results

| res | run_time | κ(B=0) | κ(+B) | κ(−B) | \|κ\| err vs theory | asymmetry | status |
|-----|----------|--------|-------|-------|---------------------|-----------|--------|
| 32  | 80       | 0      | −0.078528 | +0.078528 | 0.087% | 0 | PASS |
| 64  | 100      | 0      | −0.078479 | +0.078479 | 0.024% | 0 | PASS |

Global sign mapping (not a gate): **measured_opposite_theory_kappa_sign**
(κ_meas(+B) ≈ −κ_theory(+B)). B-reversal holds.

## Convergence
|κ| moves from 0.078528 (res 32) → 0.078479 (res 64) toward theory 0.078460.
RMSE improves slightly; R² = 1.000 for ±B at both resolutions.

## Artifacts
- `outputs/validation/faraday_res32.json`
- `outputs/validation/faraday_res64.json`
- corresponding `.log` files

## Remaining uncertainties
- Absolute Faraday sign convention (Stokes ψ vs (k₊−k₋)/2 labeling) is opposite here; deferred as a convention issue, not a magnitude/reversal failure.
- ContinuousSource DFT includes startup transient; sufficient here but not a settled-time DFT protocol.
- Homogeneous Faraday PASS does not by itself clear six-port reciprocity residuals.
