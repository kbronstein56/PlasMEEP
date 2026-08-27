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

## Finite-length reporting plan (later)

For length L (a-units) at fixed (fp, B, fs):

- theory_deg = degrees(κ_theory × L)
- meep_deg = degrees(ψ(L) − ψ(0)) with same Stokes unwrap

Report: “For this fp, B, and L, theory predicts X° and Meep gives Y° (and −B → −Y°).”

Do **not** block the high-resolution reciprocity study on this.
