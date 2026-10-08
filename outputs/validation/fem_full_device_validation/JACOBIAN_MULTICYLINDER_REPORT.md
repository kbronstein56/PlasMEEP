# Multi-cylinder Jacobian validation

Three independent plasma disks, each with its own density scale `s_k`.
Direct columns `dx/ds_k` from one factorization of `A(s=1)` versus centered FD
that perturbs only rod `k`. Frozen relative tolerance ≤ 1e-4 with truncation.

Script: `jacobian_multicylinder.py`  
Data: `jacobian_multicylinder.json`

## Indexing check

At the best `h`, the FD field for column `k` must match direct column `k` and
disagree with every other direct column by ≥10× in relative residual. This
guards against swapped rod masks.

## Results

| B (T) | rod | min field rel | index OK | pass |
|-------|-----|---------------|----------|------|
| 0     | 0   | 2.9e-8        | yes      | yes  |
| 0     | 1   | 2.2e-8        | yes      | yes  |
| 0     | 2   | 3.0e-8        | yes      | yes  |
| +0.05 | 0   | 5.3e-8        | yes      | yes  |
| +0.05 | 1   | 4.0e-8        | yes      | yes  |
| +0.05 | 2   | 5.5e-8        | yes      | yes  |

Probe-sample relative errors track the field errors. Cross-column residuals
stay O(1) while own-column residuals fall below 1e-7.

## Verdict

**MULTICYLINDER_JACOBIAN: PASS**
