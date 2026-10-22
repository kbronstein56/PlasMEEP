# Literature action items

Only ideas sufficiently relevant to potentially **test** in PlasMEEP.
Updated: 2026-10-21 (seed pass).

| ID | Source | Action | Priority | Status |
|----|--------|--------|----------|--------|
| A1 | LIT-001, LIT-003 | Implement/validate **C3** tied orbits (~31 vars) and controlled C6 vs C3 opt | High | In campaign |
| A2 | LIT-004 | Benchmark **forward-mode** sensitivity vs adjoint for 16-param C6 after shared LU | High | In campaign |
| A3 | LIT-005 | Optional JAX `custom_vjp` wrapping validated adjoint (API only); no FULL91 rewrite | Medium | Prototype if cheap |
| A4 | LIT-006 | Time-domain Meep validation plan for 4–6 GHz; Meep ≠ ground truth | High | Docs + staged TD |
| A5 | LIT-007 | After soft-min sampling study, consider autocorrelation / smoothness penalty if dense audits show spike gaming | Medium | After Phase 7 |
| A6 | LIT-008, LIT-009 | Short **mode-preserving vs port-power** objective comparison | High | In campaign |
| A7 | LIT-002, LIT-003 | Keep continuum SF quotes (not M-mesh J) when comparing to experimental plasma isolation | Ongoing | Docs |

**Not actioned (NOT_APPLICABLE or premature):**

- Replacing SuperLU with an unbenchmarked solver based on paper claims alone.
- Porting production FEM to JAX because ceviche uses autograd.
- Claiming experimental TE0/TE1 resolution without a measurement plan.
