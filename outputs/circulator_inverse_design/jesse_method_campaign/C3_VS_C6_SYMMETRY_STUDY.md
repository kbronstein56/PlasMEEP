# C3 vs C6 symmetry study

## Answer to Jesse (explicit)

**C6 is not required by circulator physics.** It follows from imposing full
60° rotational symmetry on our **six-port** cyclic target. Fan/Wang is **C3**
because that device has **three** ports.

| Symmetry | Variables | Orbit construction |
|----------|----------:|--------------------|
| C6 | 16 | production |
| C3 | **31** | **PASS** (`THREEFOLD_DESIGN_ORBITS.md`) |

## Rotation / source reduction

`rotation_source_check.json`: predicting other TE0 sources by 60° port permutation
gives ~**5%** relative column error on FULL91-C → **do not** drop to one source
on this mesh. LU/dA still dominate anyway.

## Controlled short opt (FULL91-C, Objective A, maxiter 5)

Same physical warm start: C6 `opt_uniform_C_warm_M` best → C3 via orbit average
(map |Δs|=0 to machine precision for a C6-symmetric vector).

| | C6 | C3 |
|--|---:|---:|
| n_var | 16 | 31 |
| J_best | **0.727** | 0.444 |
| T_des | 0.825 | 0.708 |
| isolation (dB) | **20.2** | 14.6 |
| RL (dB) | 16.3 | 15.9 |
| TE1 | 0.030 | 0.045 |
| IL (dB) | 0.84 | 1.50 |

Artifacts: `opt_c6_short.json`, `opt_c3_short.json`, `c3_c6_short_opt.json`.

**Reading:** Under this short equal budget, **C6 wins**. Extra C3 freedom did not
help from a C6-shaped start in 5 L-BFGS iterations. A longer C3 campaign from a
C3-native start could still win — not ruled out.

## Recommendation

- Tell Jesse: C6 is a symmetry **choice** matching six-port 60°; Fan is C3 for
  three ports.
- Empirically freeze **C6** for the first BB production run unless a longer C3
  study from a C3 start beats C6 on continuum metrics.
- Keep six incident sources on FULL91-C.
