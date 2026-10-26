# Jesse method / science questions — status

Branch `agent/eigenmode-ports`. Hub:
`jesse_method_campaign/JESSE_BROADBAND_METHOD_SUMMARY.md`.

| Item | Status |
|------|--------|
| Method questions 1–6 | Answered; production spec frozen |
| Stage-1 BB opt (21×100 MHz, maxiter 30) | **COMPLETE** — J_BB −0.801 → −0.529; `opt_broadband_production_C.json` |
| Dense 50/25 MHz audit of final design | **COMPLETE** |
| Adaptive enrichment | **NEEDED** — add `4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95` |
| Next L-BFGS restart | **Not started** (await explicit approval) |

Do not rerun `production_broadband_opt.py` with the old 21-point grid.
