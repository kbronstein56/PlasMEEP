# C3 vs C6 symmetry study

## Answer to Jesse (explicit)

**C6 is not required by circulator physics.** It follows from imposing full
60° rotational symmetry on our **six-port** cyclic target
\(1\to2\to3\to4\to5\to6\to1\), which is C6-invariant.

Fan/Wang’s classic magneto-optic circulator is **naturally C3** because it has
**three** waveguides/ports — not because threefold is universally “more physical.”

| Symmetry | Variables | Status |
|----------|----------:|--------|
| C6 | **16** (1+15×6) | Production |
| C3 | **31** (1+30×3) | Orbits **PASS** (`THREEFOLD_DESIGN_ORBITS.md`) |

Under exact C3, ports form `{1,3,5}` and `{2,4,6}`; C6 makes all six equivalent.

## Rotation / source reduction (numerical)

`rotation_source_check.json` on `opt_uniform_C_warm_M` / FULL91-C:

| Check | Result |
|-------|--------|
| Predict sources 1…5 from source 0 by 60° port permutation | max rel column err **≈5.0%** |
| Pass bar 1e−6 | **FAIL** |

**Interpretation:** the triangular FULL91-C mesh is not exactly C6-symmetric, so
fields do not rotate to machine precision. **Do not** drop to one representative
source in production on this mesh. Theoretical RHS savings would be tiny versus
LU/dA anyway.

C3 still needs **two** inequivalent source classes when/if a mesh is exactly C3.

## Controlled short opt

Same Objective A, mesh C, mapped warm start, maxiter 5 — running after spectral
(`c3_c6_short_opt.json`). Neutral map C6→C3: average rods within each C3 orbit.

## Recommendation (so far)

- Answer Jesse: C6 because we imposed 60° on six ports; Fan is C3 because of three ports.
- Empirically test C3 with short opts before freezing BB symmetry.
- Keep six sources on FULL91-C despite C6 design tying.
