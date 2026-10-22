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

## Efficient implementation

### Historical FULL91-C profile (prior campaign only)

From `ADJOINT_EFFICIENCY_BENCHMARK.md`. **Do not treat “LU ≈84%” as the current
machine profile** — that was a load/era-specific observation.

| Path | Wall time | Notes |
|------|----------:|-------|
| Baseline channel adjoint | **1397 s** | ~12+12 + duplicate S |
| Efficient TE0 + scalar adjoint | **1195 s** | **1 LU + 6 fwd + 6 adj** |
| Speedup | **1.17×** | |

| Stage | Time (historical) |
|-------|------------------:|
| Assemble | 115 s |
| LU factorization | **1008 s (~84% then)** |
| Six forward RHS | 2 s |
| Six adjoint RHS | 3 s |
| dA contractions | 65 s |

### Latest benchmark (authoritative for current cost structure)

FULL91-C, `opt_uniform_C_warm_M`, `forward_mode_bench.json`:

| Stage | Time |
|-------|-----:|
| Assembly | **77.8 s** |
| LU | **9.0 s** |
| Six forward RHS | 0.6 s |
| Six adjoint RHS | **1.8 s** |
| dA contraction | **36.8 s** |
| **Total** | **~126.0 s** |

**Current dominant costs:** assembly + derivative (`dA`) construction — **not**
adjoint RHS solves and **not** LU on this machine. Adjoint RHS remains ≪ 1% of
wall time.

| Forward-mode (after shared LU) | Value |
|--------------------------------|------:|
| Sensitivity RHS | **96** (16×6) |
| Sensitivity wall | 62 s (RHS 24 s + dA 37 s) |
| Full forward est. | **150 s** |
| ‖g_fwd−g_adj‖/‖g‖ | **2.0×10⁻¹⁴ PASS** |

**Quantifying Jesse’s “extra simulation” concern (either machine):**

| Operation | Cost |
|-----------|------|
| New LU + assembly | dominant factorization/assembly |
| One adjoint RHS after LU | **≪ 1 s** on this machine; ~0.5 s in historical profile |

An adjoint RHS is **not** another full simulation.

## Scaling estimate (from this campaign’s profile)

| Case | n_param | Est. adjoint | Est. forward | Prefer |
|------|--------:|-------------:|-------------:|--------|
| C6 | 16 | ~48 s* | ~70 s* | **adjoint** (1.45×) |
| C3 | 31 | ~83 s* | ~126 s* | **adjoint** |
| free91 | 91 | ~221 s* | ~352 s* | **adjoint** |

\*Post-LU + dA model from `forward_mode_bench.json` (excludes full assemble in the
scaling helper’s LU-only slice; use full-wall 126 vs 150 for absolute compare).

Hughes et al. (2019): forward mode wins when \(N_\mathrm{param}\ll N_\mathrm{obj}\).
For our **scalar** FoM, adjoint remains preferred at 16–91 parameters.

## Sparse direct solver inventory

`solver_inventory.json` (2026-10-21):

| Solver / package | Available in `plasmeep_mpi`? |
|------------------|------------------------------|
| SciPy SuperLU (`splu`) | **Yes — production** |
| PARDISO / pypardiso | No |
| MUMPS / pymumps | No |
| PETSc / petsc4py | No |
| JAX / autograd | Yes (API only) |

SciPy SuperLU has **no** public symbolic-factor reuse across nearby frequencies.
Do not switch production without a FULL91-validated win over SuperLU.

## Recommendation

1. Keep frequency-domain discrete adjoint (`A^\top`) as production.
2. Keep efficient 6+6 RHS path.
3. Forward-mode agrees to machine precision but is slower for scalar J.
4. Next efficiency lever: factorization/solver and/or dA contraction cost — not AD frameworks.
