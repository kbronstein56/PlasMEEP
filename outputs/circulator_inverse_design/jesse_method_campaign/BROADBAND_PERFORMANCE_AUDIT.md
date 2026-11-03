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

**12 workers:** skipped — 8 workers did **not** beat 4 (204.0 s vs 203.7 s). Expected: a 4-frequency probe has at most **4** concurrent independent solves, so n_workers∈{6,8} cannot improve wall time here.

Notes:
- 6 and 8 are **numerically exact** vs serial on this probe (J and ‖g‖ bit-identical in this run) and memory-safe (~16 GiB ≪ 252 GiB).
- They are **not faster** than 4 workers on this benchmark.
- Peak RAM for 6/8 stayed ~16 GiB (only 4 workers busy at once on a 4-freq job).

Implementation (default **off**):

- Core: `scripts/.../a2_efficient_core.py` → `eval_bb_A2(..., n_workers=1)`
- Opt runner: `production_broadband_opt_stage2.py --n-workers N --omp-threads-per-worker 1`
- Scaling driver: `benchmark_parallel_scaling.py`

---

## 6. Safest fastest **measured** configuration

**Winner (measured wall time + numerical equivalence + RAM-safe):**  
`--n-workers 4 --omp-threads-per-worker 1`  
**Parent env:** `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`

- Speedup vs serial on the 4-freq probe: **3.01×**
- Peak RAM: **≈16 GiB**
- Effective cores: **≈4**
- Objective / gradient vs serial: **exact match** in the scaling re-measure

6 and 8 workers pass equivalence and are RAM-safe but are **not** faster in measurement on this probe. Do **not** use them as the authorized Stage-3 default without a separate multi-frequency (>4) scaling study.

Recommended Stage-3 launch (when authorized — **not now**):

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u <stage3_runner>.py \
  --n-workers 4 --omp-threads-per-worker 1
```

Do **not** set OMP=64 or n_workers≈64. Default remains **serial** (`n_workers=1`) until explicitly enabled.

---

## 7. Stage-3 runtime estimate (proposal frequency set)

Stage-3 proposal: **30** frequencies (Stage-2 28 + **5.75** + **5.825**).  
Assume same L-BFGS effort as Stage-2 (**83** BB evals) and historical **121 s/freq**.

| Config | Est. wall | Basis | Active cores | Peak RAM |
|--------|----------:|-------|-------------:|---------:|
| Serial | **≈84 h** | Stage-2 history | ~1 | ~4 GiB |
| **`--n-workers 4` (recommended)** | **≈28 h** | **measured 3.01×** on 4-freq probe | ~4 | ~16 GiB |
| `--n-workers 6` | ≈19 h | *projected only* (not faster on 4-freq) | ~6* | ~19 GiB† |
| `--n-workers 8` | ≈15 h | *projected only* (not faster on 4-freq) | ~8* | ~26 GiB† |

\*Core counts for 6/8 assume Stage-3’s 30 frequencies keep that many workers busy; the 4-freq probe only exercised ~4 cores.  
†Peak RAM for 6/8 on Stage-3 would rise vs the 4-freq probe once all workers are concurrently loaded.

**Authorized recommendation until further measurement:** Stage-3 with **`--n-workers 4`**, estimated **≈28 h** (±20% for line-search variability).

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
| `benchmark_parallel_scaling.py` | Dedicated 4/6/8(+12) scaling + resource monitor |
| `a2_efficient_core.py` | Shared timed A2 + optional process-parallel |
| `production_broadband_opt_stage2.py` | `--n-workers` opt-in (default 1) |
| `BROADBAND_PERFORMANCE_AUDIT.json` | Initial audit measurements |
| `BROADBAND_PARALLEL_SCALING_6_8.json` | Measured 6/8 scaling + Stage-3 recommendation |
| `perf_audit_nohup.out` / `parallel_scaling_6_8_nohup.out` | Run logs (untracked) |
