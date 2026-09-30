# Adversarial audit of the B = 0 FEM validation

Branch `agent/eigenmode-ports`. This audit looks for a reason the label could still be wrong. It does not add B = 0 tests for reassurance, and it does not start B ≠ 0.

Decision: **B0_FEM_VERIFIED_AND_VALIDATED**.

**B0 VALIDATION FROZEN**

No FEM physics error was found. No validation row was downgraded. The count stays 130 PASS, 13 MEEP_CROSS_CODE_FAIL, 1 SUPERSEDED_TEST_HARNESS, 0 UNRESOLVED.

## 1. Clean-checkout reproduction

The solver, the analytic references, and the stored JSON are committed. The sentinel calls those functions and writes `sentinel_repro.json`. It does not overwrite the stored validation files. Untracked field dumps were not inputs.

Commands, from the repository root, with `OMP_NUM_THREADS=8`:

```
python scripts/validation/fem_physics_validation/sentinel_repro.py mms
python scripts/validation/fem_physics_validation/sentinel_repro.py planar
python scripts/validation/fem_physics_validation/sentinel_repro.py oblique
python scripts/validation/fem_physics_validation/sentinel_repro.py recip
python scripts/validation/fem_physics_validation/sentinel_repro.py power
python scripts/validation/fem_physics_validation/sentinel_repro.py bare
python scripts/validation/fem_physics_validation/sentinel_repro.py three
python scripts/validation/fem_physics_validation/sentinel_repro.py seven
python scripts/validation/fem_physics_validation/sentinel_repro.py guide
python scripts/validation/fem_physics_validation/sentinel_repro.py coated
python scripts/validation/fem_physics_validation/sentinel_compare.py
```

`coated` is the graded 5.30 GHz mesh, 5,940,235 DOFs, about 2266 s. `guide` is the acceptance mesh h = 0.00064, 3,909,063 DOFs, about 385 s. The other groups are the stored configurations named below.

Exact differences, new minus stored:

| Sentinel | Stored configuration | Quantity | Delta |
| --- | --- | --- | --- |
| Full-complex off-diagonal manufactured solution, n = 64 | `mms.json` | L2 | +2.44e-18 |
| | | H1 | −1.39e-17 |
| | | max nodal | 0 |
| Analytic PEC guide, h = 0.025 | `planar_fem.json` | L2 | 0 |
| Propagating PEC guide, h = 0.00064 | `guide_eh_close.json` | Hz L2, gradient L2, element E/H max, power relative | 0 |
| Plasma slab 0.40 at 3.85 GHz, h = 0.02 | `planar_fem.json` | T dB | 0 |
| | | T phase, degrees | −2.60e-14 |
| | | R dB | +3.86e-15 |
| Bare plasma cylinder, h = 0.04 | `scatter_partial.json` `bare_prod` | forward dB and phase | 0 |
| Three coated cylinders, h = 0.04 | `scatter_partial.json` `coated_3` | forward dB and phase | 0 |
| Seven coated bulbs, h = 0.02, h_edge = 0.0045 | `coated7_local.json` | forward dB | −2.40e-11 |
| | | forward phase, degrees | −4.01e-10 |
| | | near-quartz dB | +3.74e-11 |
| | | DOFs | 0 |
| Coated cylinder, 5.30 GHz, graded | `coated_530_graded.json` | forward dB | −1.43e-12 |
| | | forward phase, degrees | −9.92e-12 |
| | | near-quartz dB | +1.67e-12 |
| | | ring L2 | +9.26e-14 |
| | | DOFs | 0 |
| Oblique PML, dpml = 1.2, h = 0.02 | `oblique_pml_fit.json` | \|R\| at 0°, 25°, 60° | +4.5e-17, +4.2e-16, −2.9e-15 |
| Quartz-disk Green reciprocity | `conservation.json` | relative difference | 0 |
| Plasma-slab power | `conservation.json` | \|r\|² + \|t\|² | 0 |

The largest phase difference is 4.0e-10 degree on the seven-bulb forward probe. That is roundoff, not a changed observable.

## 2. Weak-form derivation

`B0_WEAK_FORM_INDEPENDENT_AUDIT.md` re-derives the formulation from Maxwell's equations with the phasor `e^{-iωt}` and then compares it with the assembler.

The scalar equation is `div(ρ ∇ Hz) + k0² Hz = 0` with `ρ = ε^{-1}` and `k0 = ω`. The weak form is `∫ (ρ ∇ Hz) · ∇ v − k0² ∫ Hz v`. The code assembles stiffness minus `k0²` times the consistent mass (`area/6` and `area/12`). Off-diagonal `ρ` is contracted as `(ρ ∇ φ_j) · ∇ φ_i`.

`E = (i/ω) ρ (∂y Hz, −∂x Hz)`. The Poynting components are `Sx = ½ Re(Ey Hz*)` and `Sy = −½ Re(Ex Hz*)`. Guide power for the cosine mode is `½ Re(β) (width/2) / k0`, which is the expression in `guide_eh_close`. Natural Neumann data are PEC for this polarization. Interface continuity of `Hz` and of `ρ ∂n Hz` is the natural transmission condition. Passive media have `Im(ε) > 0`, and the Drude value used here does.

The scatterer load is a unit consistent delta. It is not the factor `i/4` in front of `H_0^{(1)}`. Ratio tests cancel that scale. No sign discrepancy was found.

## 3. PML derivation

The same note derives the complex-coordinate weights. With `s = 1 + i σ/ω` and `σ ≥ 0`, an outgoing wave `exp(i kx x̃)` decays. The physical-coordinate weights are `sy/sx` on `ρ_xx`, `sx/sy` on `ρ_yy`, `1` on both off-diagonal entries, and `sx sy` on the mass. The assembler matches that table, including the unstretched off-diagonal terms. `σ` is quadratic, zero on the inner face, and `σ_max = −ln(1e-15) · 3 / (2 dpml)`. The same stretch is used at normal and oblique incidence. The Bloch test restricts the mesh so only the `+x` layer is entered.

## 4. Analytic-reference independence

`analytic_maxwell.py` imports `json`, `pathlib`, `numpy`, and `scipy.special` only. It does not import the FEM assembler, a mesh, or Meep. `analytic_sweeps.py` and `cluster_reference.py` call that module. They do not read FEM matrices or fit a scale to a FEM field.

An alternate calculation, written without copying those coefficient formulas, was then compared with the library:

| Check | Absolute difference |
| --- | --- |
| Graf addition of `H_0` at nmax = 12 | 6.73e-14 |
| Bare-cylinder field at (2.8, 0) from an independent Mie coefficient | 5.55e-17 |
| Coated n = 0 coefficient from an independent cylindrical transfer | 6.94e-18 |
| Plasma slab `\|r\|²` and `\|t\|²` after correcting the audit transfer | 2.22e-16 and 5.55e-17 |

The first draft of the slab transfer used `cos + sin` where `exp(i kd) = cos + i sin`. That draft was wrong. The corrected transfer matches `stack_response`. The library was not changed. The Mie and coated comparisons agreed before that slab correction, so they do not depend on it.

Production `ε(3.85 GHz)` recomputed from `f = 3.85 GHz`, `fp = 8.00 GHz`, `γ = 1.00e6 Hz`, and `a = 0.020 m` differs from the FEM value by `1.78e-15`. At 5.30 GHz the same formula gives `−1.2783908053260666 + 0.00042988505760869176i`.

## 5. Invariance

`invariance_audit.json`, PEC guide, 525 DOFs, h = 0.05, except the hex test which is the analytic coated seven-bulb field.

| Change | Difference |
| --- | --- |
| COLAMD versus NATURAL ordering | 5.27e-15 relative |
| Reversed local node order | 1.23e-15 relative |
| Source phase `exp(i 0.7)` | 1.17e-15 relative, after undoing the phase |
| Whole-geometry translation by (0.37, −0.21) | 0 |
| Mirror through the guide centerline | 1.74e-15 relative to the expected sign flip of the half-wave cosine |
| Opposite quad diagonal, error versus the exact mode | both triangulations 0.022682237, norms equal to 9e-15 |
| 60° rotation of the hex, source, and sample | 1.17e-15 |
| Reordering the hex centers by 60° | 1.93e-15 |

On this coarse guide, two homogeneous monitor nodes differ from the exact mode by 1.86% and 3.30%. That is the h = 0.05 discretization error. The acceptance guide at h = 0.00064 has element E/H max 0.0897% and reproduced exactly in the sentinel.

## 6. Units

`a = 0.020 m`. `f_a = f_Hz · a / c`, so `fs_a = 0.2568443533025771` at 3.85 GHz. `k0 = 2π f_a = ω` in units `c = 1`. Plasma frequency and collision rate use the same ordinary-frequency conversion. `ε = 1 − fp² / (f² + i f γ)`. Monitor coordinates in the scatterer suite are in units of `a`. The Poynting factor `1/2` is the phasor convention, and the guide cosine contributes another `1/2` through `∫ cos² = width/2`. Both are in `p_exact`.

`sixport_common.r_plasma = 0.005 / a = 0.250`. The validated bulbs use radius `0.230 a = 4.60 mm`. Quartz radii `0.325` and `0.375` are 13 mm and 15 mm. The validation rows name `0.230` explicitly. They certify that geometry. They do not by themselves certify a 5.00 mm core.

## 7. Semantic review

Every PASS row was read against its evidence file. No row was downgraded.

Rows that share one solve say so (`6f` with `6a`, `7e` with `7a`). The bare frequency row `6g` matches `freq_fem.json`: worst probes `−0.0029 dB` at 3.20 GHz through `−0.0248 dB / −0.147°` at 5.30 GHz on the finest bare mesh. The coated 5.30 GHz entry in that same file, uniform `h = 0.02`, is `−0.293 dB` and `−12.3°`. That mesh is not a pass. Row `7f` is the graded 5,940,235-DOF solve, and the sentinel reproduced it.

Four markdown rows (`6g`, `12g`, `14i`, `15a`) contained raw `|` characters, so extra cells displayed a leftover `UNRESOLVED` or `FAIL` beside the real status. The status column and `audit_row_status.json` already had the correct label. Those lines were rewritten so each row has one status. The count is unchanged.

`14i` reports a FEM volume about 95% below a free-space Hankel volume. The row already excludes that comparison: the FEM source is a unit delta. The FEM evidence is the 0.28% change from `h = 0.02` to `h = 0.01` and the 1.5% same-mesh contour balance. The analytic balance on the Hankel field is a 0.24% residual.

Row `12h` says 11 cells across the quartz wall. That figure is `(0.375 − 0.325) / 0.0045`. On the default mesher a numeric area constraint can ignore per-region areas, so this is the requested boundary spacing. The pass itself is the probe errors, which the sentinel reproduced.

## 8. Frozen thresholds

`PASS_CRITERIA.md` was added in `d7b1aee` on 2026-09-29, while rows `3d` and `7f` were still failing. `git diff d7b1aee HEAD` for that file is empty. It is an ancestor of the guide closeout `9921fa6`, the coated closeout `405a234`, and the PML closeout `f3eb58a`.

The file still says: smooth cases `≤ 0.01 dB`, `≤ 0.1°`, `≤ 0.1%`; one cylinder `≤ 0.05 dB`, `≤ 0.5°`; several cylinders `≤ 0.10 dB`, `≤ 1°`. Those bars were not loosened. The guide element max `0.0897%` is under `0.1%`. The graded coated probes are under `0.05 dB` and `0.5°`. The seven-bulb probes are under `0.10 dB` and `1°`.

## 9. Meep separation

The 13 `MEEP_CROSS_CODE_FAIL` rows are `6m`, `6n`, `7h`, `7i`, `9c`, `10c`, `11c`, `12g`, `20g`, `20h`, `20i`, `20j`, `20k`. For each, the analytic reference has its own truncation or identity check, and the FEM probe or slab result is inside the frozen tolerance. `meep_closeout.json` keeps registrations separate. Examples, Meep only:

| Case | Meep result |
| --- | --- |
| Coated, res 40, ox = 0 | −0.133 dB, −2.40° |
| Coated, res 40, ox = 0.25 | −0.233 dB, −3.57° |
| Pair, 0°, res 30, ox = 0 | +1.07 dB, +16.1° |
| Three cylinders, res 40, ox = 0 | −0.29 dB, −5.96° |
| Seven bare, res 36, ox = 0 | −6.32 dB, −41.4° |
| Seven bare, res 36, ox = 0.25 | −0.53 dB, −5.20° |
| Seven coated, prior grid | about 0.6 dB and 15° |

The straight guide at Meep resolution 80 is `+0.070°` and `−3.6e-6 dB`. That row is PASS and is not one of the 13. Meep was not rerun and registrations were not averaged.

## 10. Evidence chain

Maxwell, `e^{-iωt}` → the weak form in the independent audit, matched term by term to `assemble_anisotropic`.

That continuum operator → independent solutions: Graf and Mie field difference `5.6e-17`, coated `n = 0` difference `6.9e-18`, slab power difference below `1e-15`, manufactured-solution L2 order `1.998`.

Those solutions → FEM convergence: guide Hz L2 `7.78e-5` and element E/H max `0.0897%` at 3,909,063 DOFs; coated forward `−0.0118 dB / −0.461°` at 5,940,235 DOFs; seven-bulb near-quartz `−0.020 dB / −0.335°` at 618,177 DOFs. The sentinel reproduced these values.

FEM solutions → power, reciprocity, and invariance: Green relative difference `0` on the quartz disk and `2.86e-14` on the guide block; lossless contour `1e-12`; factorization reorder `1e-15`; hex rotation `1e-15`.

Those checks → the production-like cluster: seven coated bulbs at the stated `0.230 / 0.325 / 0.375` radii, several probes, last local step `0.08°`.

The analytic field is not taken from the FEM solution, and the FEM error is not used to set the analytic coefficients. Agreement is two calculations of the same Maxwell problem.

## 11. Decision and limitations

**B0_FEM_VERIFIED_AND_VALIDATED**

**B0 VALIDATION FROZEN**

Known limitations, none of which failed a frozen FEM test:

- Thirteen Meep rows remain outside the analytic tolerance.
- The uniform 6.71e6-DOF coated mesh aborted in 32-bit SuperLU. The graded 5.94e6-DOF solve is the 5.30 GHz evidence.
- `sixport_common.r_plasma` is `0.250 a`. The validated core is `0.230 a`.
- The coated 5.30 GHz phase slope is about `−865 deg/GHz`. The `0.5°` bar was not relaxed.
- The old three-wall oblique box stays `SUPERSEDED_TEST_HARNESS`.
- Absolute FEM absorption is not compared with a free-space Hankel source.
- On the default circle mesher, “cells across the quartz” computed from `h_edge` is a boundary-spacing request.
