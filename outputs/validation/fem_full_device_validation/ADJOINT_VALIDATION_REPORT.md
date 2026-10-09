# Adjoint validation (three-way)

Objective: `J = |H_z(probe)|²` on the one-cylinder mesh.

Gradient routes:

1. centered finite differences of `J(s)`
2. forward sensitivity `dx/ds` then `dJ = 2 Re(y* L dx/ds)`
3. discrete adjoint `A† λ = L† y`, `dJ/ds = −2 Re(λ† (dA/ds) x)`

Frozen relative tolerance ≤ 1e-4 with truncation visibility.
Script: `adjoint_check.py`  
Data: `adjoint_validation.json`

## Forward vs adjoint

| B (T) | ‖dJ_fwd − dJ_adj‖ / ‖dJ_fwd‖ |
|-------|------------------------------|
| 0     | 1.2e-15                      |
| +0.05 | 1.4e-13                      |
| −0.05 | 9.2e-16                      |

## FD agreement (minimum over `h` sweep)

| B (T) | FD vs forward | FD vs adjoint | pass |
|-------|---------------|---------------|------|
| 0     | 1.1e-09       | 1.1e-09       | yes  |
| +0.05 | 3.3e-08       | 3.3e-08       | yes  |
| −0.05 | 3.3e-10       | 3.3e-10       | yes  |

Large `h` shows truncation; intermediate `h` forms the agreement plateau.

## Verdict

**ADJOINT: PASS** on the one-cylinder problem for `B ∈ {0, +0.05, −0.05}`.
The conjugate-transpose adjoint matches the forward Jacobian to machine
precision and both match centered FD inside the frozen gate.
