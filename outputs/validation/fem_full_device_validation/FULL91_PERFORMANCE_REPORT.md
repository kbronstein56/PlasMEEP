# Full 91-bulb solver performance

Timings are wall-clock seconds from the solve records. Factorizations used
8 OpenMP / OpenBLAS threads. One factorization is reused for all six source
ports. A new factorization is required when B, frequency, or the plasma
tensor changes.

Observed free memory before these solves was about 244 GiB. The code’s
LU estimate (`8e-6` GiB per node) is conservative: FULL91-F was estimated
at 33 GiB and the resident drop after factorization was about 12 GiB
(free memory near 232 GiB).

## Meshes

| Grade | h plasma/quartz | h horn | h air | h PML | Nodes | Triangles | nnz | Matrix | Est. LU | Factor | Assemble | Six RHS |
|-------|-----------------|--------|-------|-------|-------|-----------|-----|--------|---------|--------|----------|---------|
| C | 0.012 a | 0.040 a | 0.060 a | 0.12 a | 823,166 | 1,645,362 | 5,648,690 | 0.13 GiB | 6.6 GiB | 9.7 s | 4.6 s | 1.5 s |
| M | 0.008 a | 0.025 a | 0.040 a | 0.08 a | 1,846,108 | 3,690,764 | 12,630,911 | 0.28 GiB | 14.8 GiB | 28–30 s | 10 s | 4.0 s |
| F | 0.005 a | 0.016 a | 0.030 a | 0.06 a | 4,137,402 | 8,272,868 | 28,243,948 | 0.63 GiB | 33.1 GiB | 82–83 s | 24 s | 9.4 s |
| A | 0.005 a | 0.016 a | 0.022 a | 0.06 a | 5,061,363 | 10,120,790 | 34,712,645 | 0.78 GiB | 40.5 GiB | 128 s | — | ~10 s |
| Q | 0.005 a | 0.016 a | 0.018 a | 0.06 a | 6,052,164 | 12,102,392 | 41,646,904 | 0.93 GiB | 48.4 GiB | 190 s | — | ~12 s |
| X2 | 0.0042 a | 0.016 a | 0.022 a | 0.06 a | 6,233,120 | 12,464,304 | — | — | ~50 GiB | 165 s | — | ~12 s |
| R | 0.004 a | 0.016 a | 0.055 a | 0.11 a | 4,902,903 | 9,804,748 | 33,599,913 | 0.75 GiB | 39.2 GiB | 82 s | — | ~10 s |
| X | 0.004 a | 0.016 a | 0.022 a | 0.06 a | 6,641,826 | 13,281,716 | — | — | — | aborted | — | — |
| S | 0.0035 a | 0.012 a | 0.030 a | 0.06 a | 7,154,769 | 14,307,602 | — | — | — | aborted | — | — |

Quartz wall thickness is `0.050 a`. Requested elements across the wall:
C 4.17, M 6.25, F 10.0. Longest edge after the air-seed fix: C `0.181 a`,
M `0.120 a`, F `0.091 a`. Median edge: C `0.0124 a`, M `0.0083 a`, F `0.0049 a`.

One RHS is about one sixth of the six-RHS time: C `0.25 s`, M `0.67 s`, F `1.6 s`.
A FULL91-M case (assemble + factor + six ports) is about 59 s. A FULL91-F
case is about 120 s. Twelve FULL91-M cases (B samples and the frequency
samples) finished in 707 s.

## What can be reused

| Change | Factorization |
|--------|----------------|
| Another source port, same geometry, frequency, and B | Reuse the LU |
| Plasma density, collision rate, or radius | New factorization |
| Applied B, including the sign | New factorization |
| Frequency | New factorization |

## Consequence for a later optimization

These numbers are feasibility data only. No adjoint and no optimizer were built.

- A serial sweep of N candidates on FULL91-M is about `N` minutes if each candidate is one frequency and one B.
- The same sweep on FULL91-F is about `2N` minutes.
- Six-port characterization after the factor is a few seconds.
- FULL91-S is not a usable production mesh on this SuperLU build. The abort was an index-workspace failure with 242 GiB still available, the same class as the earlier uniform mesh near `6.7e6` degrees of freedom.
