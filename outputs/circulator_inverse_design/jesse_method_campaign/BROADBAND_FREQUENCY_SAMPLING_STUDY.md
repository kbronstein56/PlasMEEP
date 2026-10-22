# Broadband frequency sampling study (4–6 GHz)

## Physics setup

- Target band: **4.0–6.0 GHz**
- Throughout: TE0 and TE1 propagate; TE2 remains cutoff (~7.49 GHz)
- **No** propagating-channel-count transition in-band

## Frozen grid-sufficiency bars (before seeing results)

File: `GRID_SUFFICIENCY_BARS_FROZEN.json`

| Metric | Bar when refining grid |
|--------|------------------------|
| Worst isolation change | **< 0.5 dB** |
| Worst IL change | **< 0.2 dB** |
| Worst RL change | **< 0.2 dB** |
| Worst-iso frequency location jump | **≤ 100 MHz** |

Justification: continuum SF isolation uncertainty ~0.13 dB; experimental VNA
uncertainty makes sub-0.5 dB grid jitter sub-dominant to mesh/continuum error.

## Staged reconnaissance

Design: `opt_uniform_C_warm_M` on practical mesh (FULL91-**C**).

| Stage | Spacing | Points | Artifact |
|-------|---------|-------:|----------|
| A1 | 200 MHz | 11 | `spectral_200MHz.json` |
| A2 | 100 MHz | 21 | `spectral_100MHz.json` |
| A3 | 50 MHz (+25 MHz local if needed) | ≤41 | `spectral_50MHz.json` |
| Compare | — | — | `spectral_grid_comparison.json` |

Forward-only eval (no adjoint) to locate features; verify worst points later on
continuum-quality meshes.

## Adaptive enrichment rule (for optimization — not during line search)

1. Choose fixed objective frequency set
2. Optimize (L-BFGS) with that fixed set
3. Dense audit
4. Add worst between-grid frequencies
5. Restart from previous optimum
6. **Never** change the frequency set during a line search

## Answer template for Jesse

> “\(N\) objective frequencies at spacing \(\Delta f\) are initially sufficient
> when they pass the frozen grid-sufficiency bars against a denser audit; we
> still use adaptive enrichment, and final validation uses much denser sampling.
> Soft-min alone does not prevent between-grid spikes.”

**Do not** present “11 points” as universal truth until comparison JSON passes.
