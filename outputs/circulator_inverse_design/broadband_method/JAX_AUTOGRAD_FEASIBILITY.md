# JAX / Autograd feasibility (Phase 4)

**Conclusion up front:** retain the validated analytic frequency-domain adjoint for
production unless a small-problem benchmark shows a clear wall-clock win. AD does
**not** eliminate Maxwell or adjoint linear solves.

## Options

### A. HIPS autograd / Ceviche-style

- Reverse-mode AD with custom VJPs for sparse complex solves.
- Proven in photonics (Ceviche) for frequency-domain Maxwell.
- Would require wrapping our SciPy/SuperLU (or PETSc) solve with an autograd
  primitive whose VJP is exactly \(A^\top\lambda = g\).

### B. JAX

- Excellent for dense GPU kernels and custom `linear_solve` VJPs.
- Our bottleneck is **large sparse complex LU** (~0.8M–4M DOF), not dense matmul.
- JAX’s sparse story for complex FEM + SuperLU/MUMPS is immature relative to SciPy.
- Porting cost: rewrite assembly, gyrotropic tensor, DtN ports, and mesh plumbing.

### C. Current analytic direct-DtN adjoint (recommended production)

- Already FD-validated; SuperLU `trans="T"` reuse; C6 orbit contraction.
- Phase 2 reduces RHS from ~24–36 to **6+6** per frequency without changing math.

## What would have to be rewritten for JAX/autograd

| Component | Effort | Need custom VJP? |
|-----------|--------|------------------|
| Sparse FEM assembly (complex anisotropic) | High | No (or rule for `coo` build) |
| Gyrotropic \(\rho(s,B,f)\) | Medium | No |
| Sparse factorization / solve | High | **Yes** — VJP = adjoint solve |
| Direct-DtN ports / projected measure | Medium | No |
| Objective \(J(S)\) | Low | No |
| Broadband soft-min aggregation | Low | No |

## Small-problem test policy

If useful later: a 2-port straight guide or single plasma disk with dense JAX
`jnp.linalg.solve` can verify that reverse-mode AD matches our analytic \(\partial S/\partial s\).
That does **not** demonstrate wall-clock gains on FULL91.

**No production FEM port in this campaign.**

## Recommendation

| Criterion | Winner |
|-----------|--------|
| Correctness / already validated | Analytic adjoint (C) |
| Wall-clock on FULL91 today | Analytic + Phase-2 efficiency |
| Maintenance of derivative code | AD could help *if* assembly were JAX-native — it is not |
| Broadband multi-frequency | Analytic adjoint × \(N_f\) (+ parallel) |

**Decision:** keep analytic frequency-domain adjoint. Revisit JAX only if (i) a GPU
sparse complex direct solver path exists in-stack and (ii) a measured FULL91-scale
prototype beats the efficient SciPy path by a clear margin.
