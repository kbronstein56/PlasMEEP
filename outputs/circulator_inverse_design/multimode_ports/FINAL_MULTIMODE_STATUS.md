# Final multimode port status (checkpoint)

## Gates

| Gate | Status |
|------|--------|
| MULTIMODE_PORT_VALIDATED (straight empty guide) | **PASS** |
| FULL12_ONSAGER (FULL91 + DtN, 12×12) | **UNRESOLVED** (probe solve runs; full 12×12×B sweep in progress) |
| MULTIMODE_POWER_BALANCE | **NOT_REACHED** |
| MULTIMODE_ADJOINT_FD | **NOT_REACHED** |
| MULTIMODE_CIRCULATOR_OPT | **NOT_REACHED** |

## Straight-guide DtN (frozen)

- Robin stiffness: **−i β_m** Galerkin blocks on port DOFs
- Incident RHS: **+2 i β a φ** (field amplitude for unit power)
- Validation mesh: L=24, h=0.025, 2 propagating + 2 evanescent modes
- See `MULTIMODE_PORT_VALIDATION.md`, `straight_dtn_validation.json`

## FULL91 DtN probe (grade C, uniform s=1, B=+0.05 T)

- Port boundary nodes per horn: ≈112–116 on outer feed plane
- Unit TE0 drive at port 0: |a_in| ≈ 1.21 at monitor extractors (legacy monitor planes)
- Full 12-column scattering + Onsager sweep: run `full91_dtn_solve.py`

## Root cause (unchanged)

Production clear width 2a supports **TE0 + TE1** at 3.85 GHz. Fundamental-only 6×6 S is an open subblock; 12×12 propagating S is required for Onsager with gyrotropic interior.

## Commits

Pending checkpoint commit on `agent/eigenmode-ports` (multimode DtN implementation + straight-guide PASS).
