# C3 vs C6 symmetry study

## Answer to Jesse (explicit)

**C6 is not required by circulator physics.** It follows from imposing full
60° rotational symmetry on our **six-port** cyclic target
\(1\to2\to3\to4\to5\to6\to1\), which is C6-invariant.

Fan/Wang’s classic magneto-optic circulator is **naturally C3** because it has
**three** waveguides/ports — not because threefold is universally “more physical.”

Sixfold is a **strong constraint**: 91 rods → **16** variables (1 center + 15×6).

Threefold: 91 rods → **~31** variables (1 center + 30×3). Under exact C3, ports
form two rotational classes `{1,3,5}` and `{2,4,6}`; C6 makes all six equivalent.

## Orbit construction

Script: `threefold_orbits.py` → `THREEFOLD_DESIGN_ORBITS.md` / `threefold_orbits.json`.

## Controlled comparison

Same objective (mode-preserving A), mesh, mapped warm start, iteration budget:

| Case | Variables | Artifact |
|------|----------:|----------|
| C6 | 16 | `opt_c6_short.json` |
| C3 | ~31 | `opt_c3_short.json` |
| free91 | 91 | optional / skip if too expensive |

C6 design mapped into C3 as neutral warm start (`map_c6_q_to_c3`).

## Source reduction

| Symmetry | Theoretical sources | Verify first |
|----------|---------------------|--------------|
| C6 | 1 representative TE0 | `rotation_source_check.json` |
| C3 | 2 inequivalent classes | same |

**Do not overstate speedup:** LU dominates (~84%); cutting 5 of 6 RHS saves seconds,
not the factorization.

## Status / results

Fill from `c3_c6_short_opt.json` and rotation check when complete.
