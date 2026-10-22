# Jesse method / science questions — status

**Campaign complete** on `agent/eigenmode-ports` (checkpoint `9304ea3`).
Hub: `jesse_method_campaign/JESSE_BROADBAND_METHOD_SUMMARY.md`.
Production freeze: `jesse_method_campaign/BROADBAND_PRODUCTION_SPEC_FROZEN.md`.
Long 4–6 GHz opt **not started**.

| # | Question | Status |
|---|----------|--------|
| 1 | Broadband sampling | **Answered** — 200/100 MHz fail bars; −53 dB hole @ 5.75 GHz on 50 MHz; freeze 100 MHz + enrichment |
| 2 | Time-domain validation | **Plan frozen** — Meep cross-check; gates pending execution |
| 3 | Differentiable BB / JAX | **Answered** — soft-min+mean; E2E FD PASS; JAX tiny PASS; no FULL91 port |
| 4 | Gradient efficiency | **Answered** — already FD adjoint; forward-mode agree 2e-14, slower for scalar J |
| 5 | Mode preservation | **Answered** — freeze A; B mild TE1 use in short opt |
| 6 | C6 vs C3 | **Answered** — C6 choice not physics; freeze C6 (short-opt win) |

Literature: `docs/circulator_literature/`.
