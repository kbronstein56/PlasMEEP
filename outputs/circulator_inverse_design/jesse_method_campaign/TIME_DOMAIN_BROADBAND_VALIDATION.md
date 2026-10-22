# Time-domain broadband validation

Jesse asked to validate broadband with time domain. This is an **independent
cross-check**, not the optimizer ground truth.

Related prior plan: `../broadband_method/TIME_DOMAIN_BROADBAND_VALIDATION_PLAN.md`.

## Principles

1. Optimize in frequency-domain FEM (validated DtN / adjoint).
2. Validate with Meep broadband pulse + DFT / modal extraction.
3. **Do not tune FEM to match Meep.** Report disagreement honestly.
4. Known context: Meep↔FEM plasma-interface disagreements exist at B=0 in prior
   campaigns — interpret full-device TD residuals in that light.

## Pre-device gates (must pass first)

| Gate | Test |
|------|------|
| 1. Magnetized plasma tensor in Meep | Small Voigt/Faraday benchmark vs analytic / FEM |
| 2. B-sign / nonreciprocity | Flip \(B\); circulation direction must reverse consistently |
| 3. Modal / port power extraction | Compare TD modal powers to FD on a simple guide |
| 4. TE0/TE1 decomposition (if needed) | Project DFT aperture fields onto same basis as FEM |

## Full-device protocol

- Broadband pulse spanning **4–6 GHz**.
- For each necessary incident source class (C6: theoretically 1 after rotation
  check; C3: 2 classes — verify first):
  - Accumulate DFT fields / modal amplitudes.
  - Extract dense spectrum.
- Compare to FD FEM:
  - desired total / TE0 / TE1
  - reverse total
  - through / other / reflection
  - absorption where robust
- Convergence:
  - spatial resolution (ppc)
  - run duration / field decay (resonance-limited; use decay, not only \(1/\Delta f\))

## Status

| Item | Status |
|------|--------|
| Plan frozen | **This document + prior plan** |
| Small Meep gyrotropic gate | Pending execution |
| Full-device 4–6 GHz TD | Pending (after gates) |
| Long BB optimization | **Blocked** until production spec frozen |

## Deliverables when executed

- `meep_vs_fem_spectrum_*.json/png`
- Residual tables per metric
- Honest FAIL/PASS vs stated bars (bars to be set after small-gate residuals known)
