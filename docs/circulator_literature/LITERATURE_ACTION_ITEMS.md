# Literature action items

Only ideas sufficiently relevant to potentially **test** in PlasMEEP.
Updated: 2026-11-01 (post Stage-2).

| ID | Source | Action | Priority | Status |
|----|--------|--------|----------|--------|
| A1 | LIT-001, LIT-003 | Implement/validate **C3** tied orbits (~31 vars) and controlled C6 vs C3 opt | High | Done (C6 frozen) |
| A2 | LIT-004 | Benchmark **forward-mode** sensitivity vs adjoint for 16-param C6 after shared LU | High | Done (prefer adjoint) |
| A3 | LIT-005 | Optional JAX `custom_vjp` wrapping validated adjoint (API only); no FULL91 rewrite | Medium | Assessed; no production JAX |
| A4 | LIT-006 | Time-domain Meep validation plan for 4–6 GHz; Meep ≠ ground truth | High | Docs; not blocking |
| A5 | LIT-007 / LIT-010 | Dense-audit enrichment loop; minimax/soft-min sampling hygiene | High | Stage-2 done; Stage-2 50 MHz audit in flight |
| A6 | LIT-008, LIT-009 | Short **mode-preserving vs port-power** objective comparison | High | Done → A then A2_hybrid |
| A7 | LIT-002, LIT-003 | Keep continuum SF quotes (not M-mesh J) when comparing to experimental plasma isolation | Ongoing | Docs |
| A8 | LIT-011 | After Stage-2 audit, continuum F/Q/X2 at worst+representative freqs before quoting physics | High | Pending audit |
| A9 | LIT-012 | Do not conflate frequency-broadband robustness with fabrication-stochastic TO in talks | Ongoing | Docs |

**Not actioned (NOT_APPLICABLE or premature):**

- Replacing SuperLU with an unbenchmarked solver based on paper claims alone.
- Porting production FEM to JAX because ceviche uses autograd.
- Claiming experimental TE0/TE1 resolution without a measurement plan.
