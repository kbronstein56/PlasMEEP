# Broadband frequency-parallel scaling (production notes)

## Machine

- AMD Ryzen Threadripper PRO 3995WX: **64 cores / 128 threads**
- WSL available RAM ≈ **440 GiB**

## Established facts

1. Independent frequency solves parallelize correctly and **numerically reproduce** serial J and ∇J (within tight relative tolerances).
2. Parent and per-worker `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` are required to avoid oversubscription.
3. Useful wall-time scaling for the FULL91-C A2 adjoint stack **saturates near 12 process workers**. Beyond that, per-frequency FEM time inflates under memory/bandwidth contention and wall time stays flat while RAM grows linearly.

## Benchmark history

| Study | Probe | Workers tested | Recommendation |
|-------|-------|----------------|----------------|
| `BROADBAND_PARALLEL_SCALING_6_8.json` | small | 6, 8 | early Stage-3 planning |
| `BROADBAND_PARALLEL_SCALING_STAGE3_SUBSET.json` | 12 freqs | 1–12 | **12 workers** for Stage-3/4 |
| **`BROADBAND_PARALLEL_SCALING_48FREQ.json`** | **48 freqs** | **12, 24, 32, 48** | **still 12 workers** (24 ≈ +5% wall, 2× RAM) |

### 48-frequency high-worker result (Stage-4 `q_best`)

| n_workers | wall (s) | vs 12w | peak RAM (GiB) | equiv cores | J/grad OK |
|----------:|---------:|-------:|---------------:|------------:|:---------:|
| 12 | 1094.7 | 1.00× | 44.7 | 11.6 | yes |
| 24 | 1038.8 | 1.05× | 86.3 | 23.6 | yes |
| 32 | 1063.1 | 1.03× | 111.7 | 26.6 | yes |
| 48 | 1083.9 | 1.01× | 163.5 | 47.3 | yes |

**Do not assume 48 workers will turn a 321-frequency Stage-5 eval into a ~15-minute job.**  
Expect ≈ `ceil(n_freq/12) × ~274 s` wall per BB eval at 12 workers.

## Production default

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u ... --n-workers 12 --omp-threads-per-worker 1
```

See also: `BB_STAGE5_SAMPLING_STRATEGY.md` (Stage-5 **not** authorized).
