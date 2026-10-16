# JESSE_BROADBAND_METHOD_SUMMARY

Plain answers to the advisor questions. Details live in the sibling docs in this folder.
Single-frequency multimode freeze is **untouched**.

## 1. How many objective frequencies for 4–6 GHz isolation?

**Adaptive — not a single fixed number.**

- Initial optimization grid: **11** frequencies (4.0, 4.2, …, 6.0 GHz).
- Dense audit: **41** frequencies (0.05 GHz).
- Add worst / unresolved audit points to the objective and re-optimize until the
  dense response no longer hides failures between samples.
- Final FD validation: ≤50 MHz spacing (25 MHz near sharp features).

See `ADAPTIVE_BROADBAND_SAMPLING_REPORT.md` and `BROADBAND_DIFFERENTIABLE_OBJECTIVE.md`.

TE1 is propagating throughout 4–6 GHz (cutoff ~3.747 GHz); TE2 remains cut off
(~7.49 GHz) → same 12-channel TE0+TE1 DtN definition across the band.

## 2. Is the broadband objective differentiable?

**Yes.** Per-frequency \(J(f)\) is the validated multimode score. Broadband
aggregation is a **soft-min (log-sum-exp)** plus a small mean term (\(\tau=0.05\),
\(\alpha=0.25\), frozen). Gradient:

\[
\partial J_{\mathrm{BB}}/\partial s
= \sum_i w_i\,\partial J(f_i)/\partial s
\]

with analytic softmax weights \(w_i\). No JAX required for this chain rule.

## 3. JAX / autograd?

**Assessment only — not ported.** Reverse-mode AD through a linear solve **is** an
adjoint; it does **not** remove Maxwell/adjoint RHS work. Porting FULL91 sparse
complex FEM + DtN to JAX would be a large rewrite with unclear wall-clock gain.
**Retain the validated analytic frequency-domain adjoint** for production
(`JAX_AUTOGRAD_FEASIBILITY.md`).

## 4. Is there a frequency-domain adjoint?

**Yes — that is the current method.**

\[
A(f,s)x=b,\qquad
\partial S/\partial s
\propto
(A^{-\top}c)^\top(\partial A/\partial s)(A^{-1}c).
\]

Production already reuses one SuperLU factor via `lu.solve(..., trans="T")`
(`FREQUENCY_DOMAIN_ADJOINT_METHOD.md`).

## 5. Most efficient implementation found?

**Measured on FULL91-C** (`opt_uniform_C`, `ADJOINT_EFFICIENCY_BENCHMARK.md`):

| Step | Count | Time (s) |
|------|------:|---------:|
| Assembly | 1 | ~115 |
| LU factorization | **1** | **~1008** (dominates) |
| TE0 forward RHS | **6** | ~2.3 |
| Scalar-objective transpose RHS | **6** | ~2.9 |
| Orbit \(\partial A\) contractions | 16 | ~65 |
| **Total efficient eval+grad** | | **~1195 s** |
| Baseline (12+12 + duplicate S) | | ~1397 s |
| Speedup | | **1.17×** |
| Grad agreement vs frozen | | \(\|\Delta g\|/\|g\|\sim 9\times10^{-14}\) |
| FD orbit spot check | | rel \(6.9\times10^{-5}\) → **PASS** |

Eliminates TE1-incident solves during opt (TE1 **outputs** remain in \(J\)) and
removes duplicate `direct_S_from_factor` solves. Further gains require faster
factorization or multi-frequency parallelism (LU still ~84% of runtime).

Parallel-frequency strategy: see `PARALLEL_FREQUENCY_STRATEGY.md`.

## 6. How will broadband be validated?

1. Dense independent **frequency-domain** sweep (41+ points, then ≤50/25 MHz).
2. Adaptive enrichment until hidden dips are gone.
3. Independent **time-domain Meep** broadband pulse + DFT/modal monitors
   (`TIME_DOMAIN_BROADBAND_VALIDATION_PLAN.md`), interpreted in light of known
   Meep↔FEM plasma-interface issues — Meep is not sole ground truth.

## Campaign hygiene

- Branch: `agent/eigenmode-ports` only.
- Frozen SF objective / port package: **not modified**.
- New work lives under `broadband_method/` (scripts + outputs).
