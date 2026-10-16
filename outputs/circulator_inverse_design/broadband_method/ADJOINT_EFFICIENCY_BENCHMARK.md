# ADJOINT_EFFICIENCY_BENCHMARK

Grade `C`, design `opt_uniform_C`.

## Wall time (one eval+grad)

| path | seconds | RHS / LU |
|------|--------:|----------|
| baseline (`_channel_basis` + `direct_S_from_factor`) | 1397.0 | ~12+12 fwd/adj + duplicate S solves |
| efficient (TE0-only + scalar adjoint) | 1195.0 | **6 fwd + 6 adj, 1 LU** |
| speedup | **1.17×** | |

## Efficient path breakdown

| stage | s |
|-------|--:|
| assemble_s | 114.763 |
| factor_lu | 1007.795 |
| forward_te0_rhs | 2.252 |
| S_and_objective | 1.820 |
| adjoint_te0_rhs | 2.857 |
| dA_contraction | 65.436 |
| total | 1194.924 |

## Accuracy vs frozen baseline

| check | value | bar |
|-------|------:|----:|
| ΔJ | 2.442e-15 | ~0 |
| max ‖ΔS[:,TE0]‖ | 6.351e-15 | 1e-8 |
| ‖Δg‖/‖g‖ | 8.697e-14 | 0.001 |
| FD orbit 1 rel | 6.882e-05 | 0.001 |

**EFFICIENT_ADJOINT: PASS**

Redundant solves removed: TE1 incident columns not factored into opt;
`direct_S_from_factor` not called (S from U); channel-wise V replaced by
six scalar-objective adjoint RHS.

