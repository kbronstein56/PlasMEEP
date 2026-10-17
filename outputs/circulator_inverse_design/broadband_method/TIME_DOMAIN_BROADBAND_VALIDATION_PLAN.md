# TIME_DOMAIN_BROADBAND_VALIDATION_PLAN

Independent check of a broadband FEM candidate using **Meep / PlasMEEP** FDTD.
Meep is **not** treated as unquestionable ground truth (known B=0 plasma-interface
disagreement with FEM). Report agreement, disagreement, and convergence of each side.

## Goal

Independent cross-check of the frequency-domain FEM broadband spectrum — **not**
assumed ground truth (known Meep↔FEM plasma-interface disagreements).

For each TE0 source horn, recover vs frequency over **4–6 GHz**:

- desired TE0 transmission
- reverse TE0
- through TE0
- other TE0 leakage
- TE1 conversion (all ports)
- reflection (where DFT/modal extraction is reliable)

Throughout 4–6 GHz both TE0 and TE1 propagate (TE2 cut off), so Meep modal
extraction should use the same two propagating channels as the FEM 12×12 package.

## Setup (proposed)

1. **Geometry / materials:** same FULL91 device, same \(s\) vector as the FEM candidate,
   \(B=+0.05\) T (and optionally \(B=-0.05\) for Onsager-style swap check in time domain
   is harder; prefer FEM for B-reversal).
2. **Source:** broadband pulse with spectrum covering 4–6 GHz (e.g. Gaussian centered
   near 5 GHz with adequate bandwidth).
3. **Monitors:** modal (or flux + mode decomposition) DFT at each port plane for TE0 and TE1.
4. **DFT samples:** initially **25–50 MHz** spacing (41–81 frequencies). Tighten to 25 MHz
   near FEM-identified sharp features.
5. **Runtime:** continue until port transient energy decays enough for the required
   spectral accuracy (document residual energy / DFT convergence vs runtime).
6. **Resolution:** start from prior quartz-agreed ppc; treat plasma interfaces with
   known caution — if plasma B≠0 Meep is still untrusted, label the comparison as
   **partial / diagnostic** and prioritize FEM dense sweep + power balance.

## Comparison metrics

At each DFT frequency \(f\):

| Quantity | FEM | Meep | Report |
|----------|-----|------|--------|
| \(T_{\mathrm{des}}\) | ✓ | ✓ | Δ and ΔdB |
| isolation | ✓ | ✓ | ΔdB |
| TE1 total | ✓ | ✓ | absolute Δ |
| \(R\) | ✓ | if reliable | flag DFT noise |

Also: mesh/ppc convergence on Meep; FEM M↔F transfer of the same design.

## Staging

1. **Empty horns / quartz-only** broadband pulse (expect good Meep↔FEM from prior B=0 work).
2. **Uniform plasma B=0** (expect possible disagreement — document).
3. **Optimized magnetized design** only after (1) is clean; interpret (2)–(3) in light of
   known plasma FDTD issues.

## Deliverables

- Script sketch under `scripts/circulator_inverse_design/broadband_method/`
  (implementation optional in this campaign if RAM/time limited).
- JSON of DFT curves + markdown comparison table.
- Explicit statement: validation = dense FEM audit **plus** time-domain cross-check;
  neither alone is sufficient given prior cross-code history.
