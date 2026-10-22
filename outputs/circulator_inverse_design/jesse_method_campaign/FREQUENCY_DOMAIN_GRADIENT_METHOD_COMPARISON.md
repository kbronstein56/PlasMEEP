# Frequency-domain gradient method comparison

**Verdict for Jesse:** The current production method **already is** a
frequency-domain discrete adjoint. An adjoint RHS after LU is **not** another
full Maxwell factorization.

## Production formulation (do not casually change)

At each frequency:

\[
A(f,s)\,x = b.
\]

The validated gradient uses the complex **bilinear** transpose system
appropriate to our DtN product (`A^\top`, SuperLU `trans="T"`), **not** a casual
replacement by \(A^H\). FD validation of this choice stands; do not reopen
without a concrete failure.

## Efficient implementation (measured)

FULL91-C, design `opt_uniform_C` (`ADJOINT_EFFICIENCY_BENCHMARK.md`):

| Path | Wall time | Notes |
|------|----------:|-------|
| Baseline channel adjoint | **1397 s** | ~12+12 + duplicate S |
| Efficient TE0 + scalar adjoint | **1195 s** | **1 LU + 6 fwd + 6 adj** |
| Speedup | **1.17×** | |

Efficient breakdown:

| Stage | Time |
|-------|-----:|
| Assemble | 115 s |
| **LU factorization** | **1008 s (~84%)** |
| Six forward RHS | 2 s |
| Objective | 2 s |
| Six adjoint RHS | 3 s |
| dA contractions | 65 s |

**Quantifying Jesse’s “extra simulation” concern**

| Operation | Cost on FULL91-C |
|-----------|-----------------:|
| New LU + assembly | ~1120 s |
| One adjoint RHS after LU | **~0.5 s** (6 RHS ≈ 3 s total) |

So the adjoint is **not** “another simulation” in the factorization sense.
LU dominates; shaving RHS matters far less than solver/factorization improvements.

## Forward-mode vs adjoint (this campaign)

Script: `scripts/.../jesse_method_campaign/forward_mode_bench.py`

After shared LU:

| Method | RHS (16-param C6) |
|--------|-------------------|
| Efficient adjoint | 6 fwd + 6 adj |
| Forward-mode direct sensitivity | 6 fwd + **16×6 = 96** sensitivity RHS |

Expected: adjoint still wins for a **scalar** objective even at 16 parameters,
because LU is shared and adjoint RHS count is independent of \(N_\mathrm{param}\).
Forward mode becomes attractive when \(N_\mathrm{param}\ll N_\mathrm{obj}\)
(many simultaneous objectives) — see Hughes et al., ACS Photonics 2019.

Measured numbers: see `forward_mode_bench.json` (filled when bench completes).
Scaling estimates reported for 16 / ~31 / 91 variables from the same profile.

## Sparse direct solver inventory

Script: `solver_inventory.py` → `solver_inventory.json` (2026-10-21).

| Solver / package | Available in `plasmeep_mpi`? |
|------------------|------------------------------|
| SciPy SuperLU (`splu`) | **Yes — production** |
| PARDISO / pypardiso | No |
| MUMPS / pymumps | No |
| PETSc / petsc4py | No |
| sksparse CHOLMOD | No (also not for complex nonsymmetric) |
| JAX / autograd | Yes (API only; not FULL91 LU) |

Symbolic-factor reuse across nearby frequencies is **not** exposed by SciPy
SuperLU; would need PARDISO/MUMPS/PETSc with separate symbolic/numeric phases.

Do **not** switch production without a FULL91-validated bench that beats
~1008 s LU with matching residual/gradient.

## Recommendation

1. Keep frequency-domain discrete adjoint (`A^\top`) as production.
2. Keep efficient 6+6 RHS path.
3. Invest next efficiency effort in **factorization** (solver choice / reuse),
   not AD framework migration.
4. Use forward-mode only if a multi-objective setting makes it win the bench.
