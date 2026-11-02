# Broadband performance audit (pre–Stage-3)

**Status:** Stage-3 optimization **NOT STARTED / AWAITING AUTHORIZATION**  
**Branch:** `agent/eigenmode-ports`  
**Date:** 2026-11-02  
**Machine:** ~252 GiB RAM; host reports 128 logical CPUs (`nproc --all`).  
**Companion JSON:** `BROADBAND_PERFORMANCE_AUDIT.json`

Physics, mesh grade, A2 objective, softmin (τ,α), ports, B sign, and SuperLU path were **not** changed. Sparse-solver migration was **not** performed.

---

## 1. What Stage-2 actually spent time on

| Quantity | Value |
|----------|------:|
| Wall time | 281 238 s ≈ **78.1 h** |
| Accepted L-BFGS iters | 30 |
| Broadband evals (`EVAL #`) | **83** |
| Frequencies / eval | **28** (serial) |
| Implied time / freq | **≈121 s** |
| Implied time / BB eval | ≈3388 s |

Frequencies were evaluated **strictly sequentially** in `eval_bb_A2` (for-loop). `OMP_NUM_THREADS=8` (etc.) did **not** make frequencies concurrent — measured below.

---

## 2. Single-frequency timing breakdown (measured)

Design: Stage-2 `q_best`, grade **C**, **5.65 GHz**, `OMP/OPENBLAS/MKL=8`.

| Phase | Wall (s) | Share |
|-------|--------:|------:|
| `build_system` assemble | **94.5** | **62.5%** |
| `dA` orbit contraction (16×) | **43.2** | **28.5%** |
| SuperLU `splu` factor | **9.9** | **6.5%** |
| Adjoint 6 RHS | 2.55 | 1.7% |
| Forward 6 RHS | 0.92 | 0.6% |
| S + objective | 0.07 | — |
| **Total** | **≈151** | 100% |

| Resource | Measured |
|----------|----------|
| DOFs | 823 166 |
| LU factorizations / freq | **1** |
| Forward RHS / freq | **6** |
| Adjoint RHS / freq | **6** |
| RSS after factor | ≈2.25 GiB |
| RSS end (with dA path) | ≈3.07 GiB |
| CPU during heavy phase | **≈1.00 equivalent core** (despite OMP=8) |

**Bottleneck:** serial per-frequency work dominated by **mesh/system assembly** and **orbit `dA` contractions**, not multicore factorization. SuperLU factor is only ~10 s of ~151 s.

---

## 3. What OMP / OPENBLAS / MKL actually do

Thread sweep at **5.0 GHz** (warmup + timed), same design:

| Threads | Wall (s) | Factor (s) | ≈cores (sample) |
|--------:|--------:|-----------:|----------------:|
| 1 | 152.8 | 10.2 | 1.00 |
| 2 | 153.8 | 10.4 | 1.00 |
| 4 | 154.2 | 10.1 | 1.00 |
| 8 | 153.6 | 9.9 | 1.00 |

**Conclusion:** Raising OMP/OPENBLAS/MKL from 1→8 does **not** accelerate the production A2 path on this mesh. Those variables do **not** imply frequency parallelism. For process-parallel runs, set them to **1 per worker** to avoid oversubscription.

Production path: `scipy.sparse.linalg.splu` → **SuperLU**. No UMFPACK / PyPardiso in the env. Factorization is effectively **serial**.

---

## 4. Work classification

| Work | Class |
|------|--------|
| Frequencies | **Independent** — safe process-parallel target |
| Assemble / mesh materials | Serial within a process |
| SuperLU factor | Largely serial (1× / frequency) |
| Forward 6 RHS | Serial backsolves after factor (reuse LU) |
| Adjoint 6 RHS | Serial backsolves after factor (reuse LU) |
| `dA` over 16 orbits | Serial loop in-process |
| OMP/BLAS threads | Negligible benefit here; risk of nested oversubscription |

---

## 5. Process-parallel frequency benchmark (A2)

Probe frequencies: **4.1, 5.0, 5.65, 5.75 GHz** (Stage-2 `q_best`).  
Workers use `spawn` + `OMP/OPENBLAS/MKL=1`.

| Workers | Wall (s) | Speedup | Rel. err J_BB | Rel. err ‖g‖ | Equivalence |
|--------:|--------:|--------:|--------------:|-------------:|:-----------:|
| 1 (serial) | 626.5 | 1.00× | — | — | — |
| 2 | 341.3 | **1.84×** | 8.9e−15 | 3.3e−13 | **PASS** |
| 4 | 213.5 | **2.93×** | 8.9e−15 | 3.3e−13 | **PASS** |
| 6 | 215.1 | 2.91× | (match) | (match) | **PASS** |

Notes:
- On a **4-frequency** probe, 6 workers cannot beat 4 (only 4 independent tasks).
- Peak child RSS sum with 4 workers observed ~**15 GiB** transient; steady ~3 GiB / worker.
- **8 workers** would be RAM-safe (~25 GiB) on this host for a **30-frequency** Stage-3 eval; not separately timed on the 4-freq probe.

Implementation (default **off**):

- Core: `scripts/.../a2_efficient_core.py` → `eval_bb_A2(..., n_workers=1)`
- Opt runner: `production_broadband_opt_stage2.py --n-workers N --omp-threads-per-worker 1`

---

## 6. Safest fastest verified configuration

**Verified (numerical equivalence):** `--n-workers 4 --omp-threads-per-worker 1`  
**Parent env:** `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`

Recommended launch sketch (when Stage-3 is authorized — **not now**):

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u production_broadband_opt_stage3.py \  # when authorized
  --n-workers 4 --omp-threads-per-worker 1
```

Conservative next step after authorization: try **`--n-workers 6`** (projected better on 30 freqs; RAM still ≪ 252 GiB). Do **not** set OMP=64 or n_workers≈64.

Default remains **serial** (`n_workers=1`) until explicitly enabled.

---

## 7. Stage-3 runtime estimate (proposal frequency set)

Stage-3 proposal: **30** frequencies (Stage-2 28 + **5.75** + **5.825**).  
Assume same L-BFGS effort as Stage-2 (**83** BB evals) and historical **121 s/freq**.

| Config | Est. wall | Active cores | Est. peak RAM |
|--------|----------:|-------------:|--------------:|
| Current serial (as Stage-2) | **≈84 h** | ~1 | ~3 GiB |
| Verified `--n-workers 4`, OMP=1 | **≈29 h** | ~4 | ~12 GiB |
| Projected `--n-workers 6` | **≈19 h** | ~6 | ~19 GiB |
| Projected `--n-workers 8` | **≈15 h** | ~8 | ~25 GiB |

Line-search eval count can vary; treat ±20% as plausible.  
**STAGE-3 IS NOT STARTED.**

---

## 8. Future multicore sparse-solver plan (do not migrate yet)

SuperLU is only ~6–7% of SF time today; **assembly + dA** dominate. A faster sparse direct solver is still worthwhile later (and for leaner forward-only audits), but it is **not** the first Stage-3 lever.

When authorized to benchmark (not migrate production):

1. Candidates that support **complex nonsymmetric** CSC + multi-RHS: **oneAPI MKL PARDISO** (via `pypardiso` or `PyPardiso`-style wrappers), optionally **CHOLMOD** only if structure becomes SPD (it is **not** for this Helmholtz/DtN system).
2. Acceptance gates (all required):
   - Same S / A2 / softmin J to tight tolerance vs SuperLU on grade-C golden freqs
   - Gradient FD or adjoint cross-check unchanged
   - Multi-RHS forward + adjoint reuse after one factor
   - Peak RAM per process documented
   - Threaded factor does not regress when `n_workers>1` (set solver threads = 1 under process parallel, or serial freqs + threaded factor — pick one level of parallelism)
3. Keep SuperLU behind a switch until golden diffs pass on the Stage-2 `q_best` + uniform designs.

**Higher priority for wall-time than solver swap:** frequency process-parallel (done, opt-in) and, later, speeding **`build_system`** / **`assemble_dA_rod`** (profile-guided, still physics-preserving).

---

## 9. Artifacts

| File | Role |
|------|------|
| `profile_broadband_performance.py` | Audit driver |
| `a2_efficient_core.py` | Shared timed A2 + optional process-parallel |
| `production_broadband_opt_stage2.py` | `--n-workers` opt-in (default 1) |
| `BROADBAND_PERFORMANCE_AUDIT.json` | Raw measurements |
| `perf_audit_nohup.out` | Run log (untracked) |
