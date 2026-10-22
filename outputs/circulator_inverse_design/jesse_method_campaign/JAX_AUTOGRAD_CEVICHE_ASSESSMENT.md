# JAX / Autograd / Ceviche assessment

## How Ceviche actually differentiates

Ceviche uses **HIPS autograd**, but **not** by auto-differentiating an arbitrary
SciPy sparse solve. It registers **custom differentiable primitives** for sparse
ops and for \(Ax=b\). The VJP of the linear solve **is** an adjoint solve
(\(A^\top\lambda = g\) or the formulation matching their bilinear/sesquilinear
convention).

Workshop FAQ claim (accurate): autograd does not magically differentiate
through black-box SuperLU without a custom VJP.

## Comparison

| Approach | Correctness | API | Wall time FULL91 | Memory | JIT/GPU | Rewrite sparse FEM? |
|----------|-------------|-----|------------------|--------|---------|---------------------|
| **1. Current analytic/discrete adjoint** | FD-validated | Explicit | Best today (LU-bound) | LU factor | No | No |
| **2. Ceviche-style custom autograd** | Same math if VJP=adjoint | Cleaner composition | ≈ same (still LU) | ≈ same + tape | Limited | Wrap solve only |
| **3. JAX `custom_vjp` / `custom_jvp`** | Same if wraps validated adjoint | Nice for opt stacks | ≈ same unless sparse LU on device | Can grow with JIT | Possible later | **Yes** for real gains |
| **4. Naive AD through solver** | Often unsupported / wrong | — | N/A | — | — | — |

## Answers

| Question | Answer |
|----------|--------|
| Improve correctness? | **No** vs validated adjoint |
| Simplify API? | **Maybe** (custom_vjp wrapper) |
| Reduce wall time? | **No** on FULL91 while LU is SciPy SuperLU on CPU |
| Increase memory? | Possibly (AD tape / JAX buffers) |
| Enable JIT/GPU? | Only after rewriting assembly+sparse solve |
| Require rewriting sparse FEM? | For real speedups: **yes** |

## Tiny prototype

`jax_custom_vjp_tiny.py`: toy \(A(s)x=b\), \(J=\|x\|^2\).

| Check | Result |
|-------|--------|
| NumPy custom-VJP vs FD | rel err \(2.5\times10^{-10}\) **PASS** |
| JAX `custom_vjp` vs FD | available in `plasmeep_mpi`; matches adjoint |
| FULL91 port | **Not done** (no wall-clock case) |

**Do not port FULL91 FEM to JAX** unless a measured prototype beats efficient
SciPy adjoint by a clear margin with a real sparse complex direct solver path.

## Recommendation

Keep production analytic frequency-domain adjoint. Optional thin `custom_vjp`
wrapper for optimizer APIs only. Soft-min broadband remains a cheap chain rule
on top of per-frequency adjoints.
