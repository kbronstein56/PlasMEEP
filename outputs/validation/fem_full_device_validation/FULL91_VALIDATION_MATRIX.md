# Full 91-bulb forward validation matrix

Branch: `agent/eigenmode-ports`

Production reference meshes:

- B = 0 continuum estimate: **FULL91-Q** (6,052,164 nodes)
- Magnetized production reference: **FULL91-Q**, with **FULL91-F** retained as the prior magnetized six-port matrix

Overall label:

    FULL91_GYROTROPIC_FORWARD_VALIDATED

Mandatory FAIL = 0. Mandatory UNRESOLVED = 0.
The Hamid & Cooray ellipse echo-width attempt is a separate literature FAIL
and does not contradict the already-validated Lorentz cylinder suite.

| # | Gate | Status | Evidence |
|---|------|--------|----------|
| 1 | Exact production geometry | PASS | Max bulb-center difference `3.553e-15`. Plasma radius `0.230 a`. See `FULL91_GEOMETRY_AUDIT.md`. |
| 2 | Gyrotropic path at B = 0 equals the scalar operator | PASS | On FULL91-C, `\|\|A - A_scalar\|\| / \|\|A\|\| = 0`. Reaction symmetry on C/M/F is `3.3e-14`, `2.1e-14`, `2.5e-14`. |
| 3 | Full-device B = 0 mesh convergence | PASS | After air-seed and air-refinement grades A then Q (plasma/quartz held at F size): A→Q source P1 max `|ΔdB| = 0.060`, max `|Δphase| = 0.35°`. Frozen target `≤ 0.10 dB`, `≤ 1°`. Earlier M→F exceeded the target; that failure is superseded by the A/Q pair, not erased. Interface-only mesh R moved ports toward M and is diagnostic evidence that empty-cavity air mattered. FULL91-S (7.15e6) still cannot be factored and was not repeated. |
| 4 | Full-device B ≠ 0 mesh convergence | PASS | F→Q at `+0.05 T`: every receiving port `≤ 0.047 dB` and `≤ 0.53°`. Prior M→F at `±0.05 T` was already `≤ 0.083 dB` and `≤ 0.87°`. |
| 5 | +B / −B physics | PASS | Adjacent-port circulation reverses. Through power unchanged at the `1e-6` level. Absorption equal for both signs to about `3e-9`. |
| 6 | Onsager–Casimir | PASS | Reaction residual on FULL91-F: max `1.173e-13`. Same-B nonsymmetry `0.108` of the max entry. |
| 7 | Closed energy / passivity | PASS | Control volume = PML-inner rectangle with a thin source slot excised. Frozen relative residual `≤ 0.02`. Mesh F, all six sources, B = 0 / `+0.05` / `−0.05`: max relative residual `8.36e-3`. Absorption ≥ 0 everywhere. See `full91_balance_F.json`. |
| 8 | Port extraction uncertainty | PASS | Quantified on mesh F for P2/P4/P5 at B = 0 and `+0.05 T`. Fixed production span `0.92`: station shifts of `±0.30 a` change flux by `≤ 0.023 dB`; sample count 21–81 changes flux by `≤ 0.003 dB`. Changing the span to 0.80 or 0.98 moves flux by about `0.6–1.1 dB` because the extractor is not mode-normalized; that is recorded as the cost of changing the definition, not as mesh noise. See `full91_port_uncertainty_F_*.json`. |
| 9 | B continuation | PASS | Dense samples include `0.045` and `0.055 T` around production. Contrast and absorption stay smooth and positive. Nonreciprocity → 0 as B → 0. |
| 10 | Frequency continuation | PASS | Dense samples at `3.80`, `3.82`, `3.84`, `3.86`, `3.88 GHz` around `3.85 GHz` at `+0.05 T`. `Re ε_xx` monotone. Through power is steep but mesh-factorable and absorption-positive. |
| 11 | Solver stability | PASS | Meshes through Q (6.05e6 nodes) factorize. Observed free-memory drop on Q about 21 GiB from 245 GiB. S remains outside this SuperLU. |
| 12 | Literature benchmarks / physics consistency | PASS for inconsistency; paper reproduction separate | Fan 2019/2018/2020: `NOT_REPRODUCIBLE_FROM_PAPER` (no digit tables). Hamid & Cooray 2016 ellipse: Mie calibration max relative error `0.021` against frozen `0.005`, so the ellipse comparison is not interpretable under the frozen criteria and is **FAIL** as a paper reproduction. Independent gyrotropic-cylinder `T_n` agreement (`1.0e-4` dB) is unchanged. No contradiction with the validated Lorentz operator. |

Adjoint / Jacobian work may begin only after this label. Differentiability is not claimed by this matrix alone.
