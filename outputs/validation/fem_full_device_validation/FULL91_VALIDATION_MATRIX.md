# Full 91-bulb forward validation matrix

Branch: `agent/eigenmode-ports`

Production reference mesh for the magnetized solves: **FULL91-F**
(4,137,402 nodes, 8,272,868 triangles).

Overall label:

    FULL91_GYROTROPIC_FORWARD_PARTIALLY_VALIDATED

The magnetized operator, Onsager–Casimir identity, absorption sign, and
B-reversal all pass on the factorable meshes. The model is not cleared
for adjoint development. The B = 0 port powers are still moving between
FULL91-M and FULL91-F by more than 0.10 dB and 1 degree, and a tighter
interface mesh (FULL91-S, 7,154,769 nodes) cannot be factored with this
32-bit SuperLU build.

| # | Gate | Status | Evidence |
|---|------|--------|----------|
| 1 | Exact production geometry | PASS | Max bulb-center difference `3.553e-15`. Plasma radius `0.230 a`. See `FULL91_GEOMETRY_AUDIT.md`. |
| 2 | Gyrotropic path at B = 0 equals the scalar operator | PASS | On FULL91-C, `\|\|A - A_scalar\|\| / \|\|A\|\| = 0`. Reaction symmetry on C/M/F is `3.3e-14`, `2.1e-14`, `2.5e-14`. |
| 3 | B = 0 mesh convergence | FAIL | FULL91-M to FULL91-F, source P1: through `+0.226 dB` and `+1.03 deg`; adjacent ports about `-0.14 dB` and `+1.87 deg`; next ports about `-0.31 dB` and `+1.24 deg`. Trend C to M to F is monotone. FULL91-S aborted in SuperLU (`malloc fails for local dworkptr`) with 242 GiB still free. |
| 4 | B ≠ 0 mesh convergence | PASS | FULL91-M to FULL91-F at `±0.05 T`. Every receiving port is within `0.083 dB` and `0.87 deg`. Weak-port absolute power differences are below `1.3e-5`. |
| 5 | B reversal | PASS | Adjacent-port circulation reverses between `+0.05 T` and `-0.05 T`. Through power is unchanged at the `1e-6` level. Absorption is equal for both signs to about `3e-9`. |
| 6 | Onsager–Casimir `S_ij(+B) = S_ji(-B)` | PASS | Reaction residual on FULL91-F at `0.05 T`: max `1.173e-13`, RMS `2.71e-14` (matrix scale `0.294`). Same residual on FULL91-M is below `8e-14` at `0.0125`, `0.025`, `0.0375`, and `0.05 T`. Same-B reaction nonsymmetry is `0.108` of the max entry, so the device is nonreciprocal. |
| 7 | Passivity / power balance | PASS for absorption; UNRESOLVED for a closed port sum | Plasma absorption is positive for every source, both B signs, every B sample, and every frequency sample. Horn-monitor sum plus absorption leaves `11%` of the driven-port net flux at B = 0 and `23%` at `±0.05 T`. The six horn lines are not a closed surface. |
| 8 | B → 0 continuity | PASS | P1 minus P5 power contrast on FULL91-M is `-2.0e-6`, `+8.77e-3`, `+3.75e-3`, `+7.42e-4`, `+3.65e-4` at `0`, `0.0125`, `0.025`, `0.0375`, `0.05 T`, and the opposite sign at `-B`. No absorption sign change. |
| 9 | Frequency response | PASS as a successful sweep; the spectrum is steep | Five frequencies at `+0.05 T` on FULL91-M, each with a new factorization and a recomputed tensor. Absorption stays positive. Through power changes by a factor of about five between `3.754 GHz` and `3.850 GHz`. That is recorded as device structure, not a solver failure. |
| 10 | Solver stability on the production meshes | PASS | C, M, and F factorized repeatedly. Peak observed memory drop on F is about 12 GiB from 244 GiB free. FULL91-S is outside this SuperLU index range and was not repeated. |

Adjoint, gradient, and Bayesian work were not started.
