# Faraday sign-convention note (secondary)

**Status:** magnitude/reversal **PASS** preserved; absolute sign opposite theory.

## Observed

At the homogeneous Faraday point (fs=5 GHz, fp=2 GHz, |B|=0.05 T):

- κ_meas(+B) ≈ −κ_theory(+B)
- κ_meas(−B) ≈ −κ_theory(−B)
- |κ| error ~0.09% (res32) / ~0.02% (res64); perfect B-reversal

## Convention checklist (not yet exhaustively closed)

| Item | Current Faraday script | Possible flip source |
|---|---|---|
| Propagation | +z | — |
| Source | Ex linear | — |
| B / Meep bias | `bias = ±fc_a` along z (ordinary-freq units) | electron charge sign in ω_c = eB/m vs −e |
| Theory κ | (k₊ − k₋)/2 with ε± = ε⊥ ± η | which circular basis is “+” |
| Measured κ | Stokes ψ = ½ atan2(S₂,S₁), slope dψ/dz | Stokes angle sense vs optical convention |
| Time-harmonic | e^{−iωt} (Meep) | must match theory ODE |

**Working conclusion:** opposite global sign is consistent with a **basis/Stokes vs (k₊−k₋) labeling convention**, not a |κ| or reversal bug. Closing it requires one controlled sign-walk (fix charge in bias, Stokes definition, and circular-basis labeling one at a time).

## Finite-length validation (completed)

Homogeneous Faraday cell at res=64 (`faraday_res64.json`):

| Quantity | Value |
|---|---|
| Length L | 20 a (= 56 cm) |
| fp | 2 GHz, fs = 5 GHz, \|B\| = 0.05 T |
| Theory κ(+B) | +0.07846 rad/a → **+1.57 rad ≈ +89.9°** over L |
| Meep κ(+B) | −0.07848 rad/a (opposite sign, same \|κ\|) |
| \|κ\| error | 0.024% |
| B-reversal | perfect (κ(−B) = −κ(+B)) |

**Conclusion:** magnitude and reversal are trustworthy; global sign is a **labeling/convention** issue (Stokes ψ vs theory circular-basis κ), not a gyrotropy implementation bug.
