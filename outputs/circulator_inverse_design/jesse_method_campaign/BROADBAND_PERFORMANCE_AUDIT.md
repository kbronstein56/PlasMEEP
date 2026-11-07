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

### 5a. Initial sweep (earlier audit)

| Workers | Wall (s) | Speedup | Rel. err J_BB | Rel. err ‖g‖ | Equivalence |
|--------:|--------:|--------:|--------------:|-------------:|:-----------:|
| 1 (serial) | 626.5 | 1.00× | — | — | — |
| 2 | 341.3 | **1.84×** | 8.9e−15 | 3.3e−13 | **PASS** |
| 4 | 213.5 | **2.93×** | 8.9e−15 | 3.3e−13 | **PASS** |
| 6 | 215.1 | 2.91× | (match) | (match) | **PASS** |

### 5b. Dedicated 6/8 scaling re-measure (2026-11-02 evening)

Same probe, instrumented peak process-tree RSS + integrated CPU. Companion: `BROADBAND_PARALLEL_SCALING_6_8.json`.

| Workers | Wall (s) | Speedup | ≈cores | Peak RAM (GiB) | ΔJ vs serial | Δ‖g‖ vs serial | Equivalence |
|--------:|--------:|--------:|-------:|---------------:|-------------:|---------------:|:-----------:|
| 1 | 613.1 | 1.00× | 1.00 | 4.13 | — | — | — |
| **4** | **203.7** | **3.01×** | **3.99** | **16.14** | **0** | **0** | **PASS** |
| 6 | 206.0 | 2.98× | 3.93 | 15.96 | 0 | 0 | **PASS** |
| 8 | 204.0 | 3.00× | 3.98 | 15.19 | 0 | 0 | **PASS** |

**12 workers:** skipped on the 4-freq probe — 8 did **not** beat 4. That probe **cannot** decide Stage-3 worker count (only 4 concurrent jobs).

### 5c. Stage-3-representative 12-frequency scaling (authoritative for Stage-3)

**Date:** 2026-11-07  
**Companion:** `BROADBAND_PARALLEL_SCALING_STAGE3_SUBSET.json`  
**Driver:** `benchmark_parallel_scaling_stage3_subset.py`  
**Probe (12 GHz ⊆ Stage-3 30-set):**  
`[4.1, 4.3, 4.5, 4.7, 5.0, 5.2, 5.4, 5.6, 5.65, 5.75, 5.825, 6.0]`  
Includes required problem freqs **4.1, 4.3, 5.65, 5.75, 5.825**.  
`OMP=OPENBLAS=MKL=1` per worker. Design: Stage-2 `q_best`, grade C. **No Stage-3 optimization.**

| Workers | Wall (s) | Speedup | ≈cores | Peak RAM (GiB) | ΔJ | Δ‖g‖ | Equivalence |
|--------:|--------:|--------:|-------:|---------------:|---:|-----:|:-----------:|
| 1 | 1435.0 | 1.00× | 1.00 | 4.19 | — | — | — |
| 4 | 458.1 | 3.13× | 3.97 | 16.46 | 0 | 0 | **PASS** |
| 6 | 376.7 | 3.81× | 5.95 | 23.89 | 0 | 0 | **PASS** |
| 8 | 377.3 | 3.80× | 6.41 | 30.19 | 0 | 0 | **PASS** |
| **12** | **283.8** | **5.06×** | **11.59** | **43.34** | **0** | **0** | **PASS** |

Notes:
- With **12** independent frequencies, 6 and 8 workers both need **2 batches**, so their walls match (~377 s); this does **not** mean 8 is useless for the real 30-freq Stage-3 (there 8→4 batches vs 6→5).
- **12 workers** clearly wins: **5.06×**, ~11.6 effective cores, **≈43 GiB** peak (≪ 252 GiB).
- All configs numerically equivalent to serial (J and ‖g‖ exact in this run).

Implementation (default **off**):

- Core: `scripts/.../a2_efficient_core.py` → `eval_bb_A2(..., n_workers=1)`
- Opt runner: `production_broadband_opt_stage2.py --n-workers N --omp-threads-per-worker 1`
- Scaling drivers: `benchmark_parallel_scaling.py`, `benchmark_parallel_scaling_stage3_subset.py`

---

## 6. Safest fastest **measured** configuration (Stage-3 planning)

**Winner (12-freq Stage-3 subset — authoritative):**  
`--n-workers 12 --omp-threads-per-worker 1`  
**Parent env:** `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`

| Metric | Value |
|--------|------:|
| Speedup vs serial | **5.06×** |
| Effective cores | **≈11.6** |
| Peak RAM | **≈43 GiB** |
| J / ‖g‖ vs serial | exact match |

The earlier 4-freq “recommend 4 workers” conclusion is **superseded** for Stage-3; that probe was concurrency-limited.

Recommended Stage-3 launch (when authorized — **not now**):

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u <stage3_runner>.py \
  --n-workers 12 --omp-threads-per-worker 1
```

Do **not** set OMP=64 or n_workers≈64. Default remains **serial** (`n_workers=1`) until explicitly enabled.

---

## 7. Stage-3 runtime estimate (proposal frequency set)

Stage-3 proposal: **30** frequencies (Stage-2 28 + **5.75** + **5.825**).  
Assume same L-BFGS effort as Stage-2 (**83** BB evals) and historical **121 s/freq** → serial ≈ **83.7 h**.

Estimates apply **measured** speedups from the **12-freq Stage-3 subset**:

| Config | Est. wall | Basis | ≈cores (measured) | Peak RAM (measured) |
|--------|----------:|-------|------------------:|--------------------:|
| Serial | **≈84 h** | Stage-2 history | ~1 | ~4 GiB |
| `--n-workers 4` | ≈27 h | measured 3.13× | ~4.0 | ~16 GiB |
| `--n-workers 6` | ≈22 h | measured 3.81× | ~6.0 | ~24 GiB |
| `--n-workers 8` | ≈22 h | measured 3.80× (≈6 on 12-freq) | ~6.4 | ~30 GiB |
| **`--n-workers 12` (recommended)** | **≈16.6 h** | **measured 5.06×** | **~11.6** | **~43 GiB** |

±20% plausible from line-search eval count. On the full 30-freq Stage-3 set, 8 workers may pull slightly ahead of 6 (4 vs 5 batches); **12 remains the measured fastest safe choice**.

**STAGE-3 IS NOT STARTED / AWAITING AUTHORIZATION.**

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
| `profile_broadband_performance.py` | Initial audit driver |
| `benchmark_parallel_scaling.py` | 4-freq 4/6/8 scaling (concurrency-limited) |
| `benchmark_parallel_scaling_stage3_subset.py` | **Authoritative** 12-freq Stage-3-subset scaling |
| `a2_efficient_core.py` | Shared timed A2 + optional process-parallel |
| `production_broadband_opt_stage2.py` | `--n-workers` opt-in (default 1) |
| `BROADBAND_PERFORMANCE_AUDIT.json` | Aggregated measurements |
| `BROADBAND_PARALLEL_SCALING_6_8.json` | 4-freq 6/8 results |
| `BROADBAND_PARALLEL_SCALING_STAGE3_SUBSET.json` | 12-freq scaling + Stage-3 recommendation |
| `*_nohup.out` | Run logs (untracked) |
