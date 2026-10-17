# JESSE_BROADBAND_METHOD_SUMMARY

Plain answers to the advisor questions. Details live in the sibling docs in this folder.
Single-frequency multimode freeze is **untouched**.

## 0. Modal channels on the advisor target band (4–6 GHz)

| mode | cutoff | status in 4–6 GHz |
|------|--------|-------------------|
| TE0 | — | propagates |
| TE1 | ≈ 3.747 GHz | propagates across the **entire** band |
| TE2 | ≈ 7.49 GHz | cutoff (does not propagate) |

**Therefore the validated 12-channel TE0+TE1 port formulation applies unchanged
over the whole 4–6 GHz target.** There is **no** modal-channel discontinuity and
**no** TE1-cutoff derivation required for broadband on this band.

The earlier `NEEDS_CUTOFF_DERIVATION` / TE1-cutoff warning applied only to the
**3.70–4.00 GHz** single-frequency characterization window (which straddles TE1
cutoff). Do not carry that warning into 4–6 GHz planning.

Remaining broadband questions are frequency sampling density, robust aggregation
(τ/α), and compute — not channel count.

## 1. How many objective frequencies for 4–6 GHz isolation?

**Adaptive — not a single fixed number. Do not pick an arbitrary N.**

Procedure (once a trustworthy single-frequency design exists):

1. Densely sample 4.00–6.00 GHz at ≤50 MHz.
2. Measure actual spectral feature widths (isolation / TE1 / IL dips).
3. Compare visibility on 500 / 250 / 200 / 100 MHz grids.
4. Choose an **initial** optimization grid from the narrowest relevant features.
5. Adaptive enrichment: sparse opt → dense audit → add worst between-grid
   frequencies → re-optimize until the dense audit no longer hides failures.

Placeholder starting guess (until feature-width study): 11 pts @ 200 MHz for
opt, 41 pts @ 50 MHz for audit — then refine from measured widths.

See `ADAPTIVE_BROADBAND_SAMPLING_REPORT.md` and `BROADBAND_DIFFERENTIABLE_OBJECTIVE.md`.

## 2. Is the broadband objective differentiable?

**Yes.** Per-frequency \(J(f)\) is the validated multimode score. Broadband
aggregation is a **soft-min (log-sum-exp)** plus a mean term. Prototype
hyperparameters (FD-checked, **not** chosen only because FD passed):

- **τ = 0.05** — how closely soft-min approximates worst-frequency \(J\)
- **α = 0.25** — mean-performance vs worst-frequency weighting

Fixed-spectra sensitivity (`SOFTMIN_TAU_ALPHA_SENSITIVITY.md`) shows τ≈0.02
tracks sharp dips more tightly; keep 0.05/0.25 until real SF spectra are used
to retune. Soft-min weight FD: rel ≈ 5.7×10⁻¹¹ PASS.

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
