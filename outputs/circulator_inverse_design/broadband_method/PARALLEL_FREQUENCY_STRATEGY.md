# PARALLEL_FREQUENCY_STRATEGY

Each frequency: 1 LU + 6 TE0 fwd + 6 adj (efficient path). OMP_NUM_THREADS=1 per worker to avoid BLAS oversubscription. Peak aggregate RAM for multi-process is ~workers × single-process footprint; parent VmRSS does not capture children.

| workers | N_f | wall (s) | mean s/freq |
|--------:|----:|---------:|------------:|
| 1 | 2 | 3038.7 | 1518.0 |
| 2 | 2 | 4010.3 | 4009.0 |

**Recommended concurrency (this machine / grade C):** 1 worker(s) for the tested subset.

For full 11-frequency objective+grad, expect roughly `11 × (efficient single-freq time) / workers` wall seconds per L-BFGS step, plus orbit dA cost already inside each worker.

