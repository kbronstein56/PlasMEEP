# OPTIMIZED_POWER_BALANCE_DIAGNOSTIC

Design: `opt_uniform_C.json` (transferred C optimum). Frozen DtN S package **not** modified.

## Discrete identity

From \(A x = f\) with the current \(e^{-i\omega t}\) convention, lossy gyrotropic
\(\rho=\varepsilon^{-1}\), and DtN port Robin operator:

\[
x^{\mathsf H} A x = x^{\mathsf H} f
\quad\Rightarrow\quad
\mathrm{Im}(x^{\mathsf H} A x)=\mathrm{Im}(x^{\mathsf H} f).
\]

Split \(A = A_{\mathrm{vol}} + A_{\mathrm{port}}\). On all six TE0 drives on C:
**discrete Ax=f closed** (relative residual \(<10^{-9}\)).

## Three quantities (TE0 incidence)

| symbol | definition |
|--------|------------|
| **A.** P_out | \(\sum\|S_{:,n}\|^2\) from direct DtN S |
| **B.** P_Im_eps | old volume post-process: \(E^\dagger \mathrm{Im}(\varepsilon) E\) (elementwise) |
| **C.** P_diss | gyrotropic-safe: \(E^\dagger ((\varepsilon-\varepsilon^\dagger)/(2i)) E\), \(E=(i/\omega)\rho\,\mathrm{curl}\,H\) |

## Grade C (COMPLETE)

| quantity | value |
|----------|------:|
| P_out mean | 0.9053 |
| P_Im_eps mean (buggy) | 0.2290 |
| P_diss mean (correct) | 0.0912 |
| residual max Im_eps | **0.1361** |
| residual max diss | **0.0046** |
| GATE Im_eps (≤0.05) | **FAIL** |
| GATE dissipation (≤0.05) | **PASS** |
| discrete Ax=f all TE0 | True |
| elapsed | 378 s |

Sum check: \(P_{\mathrm{out}}+P_{\mathrm{diss}}\approx 0.9965\) (closes).
\(P_{\mathrm{out}}+P_{\mathrm{Im\_eps}}\approx 1.134\) (overcounts — matches the campaign FAIL).

## Grades M / F

M factorization was started but **stopped** after C already identified the root
cause, to free CPU for absorption-fix validation and M-mesh refinement
(two concurrent SuperLU factorizations were heavily oversubscribed).

**Prediction (same operator bug):** M/F will show the same pattern —
discrete+S close, P_diss restores the gate, P_Im_eps fails — because the bug is
in post-processing, not mesh-dependent DtN. Re-check M/F power balance with the
**fixed** `material_absorption` after M-opt (Phase 4/5).

## Interpretation

**IF** discrete operator identity + S close **BUT** volume P_material (Im_eps) does not
→ **material_absorption post-processing is wrong / inconsistently normalized.**

That is exactly the C result.

Root cause: for collisional gyrotropic \(\varepsilon\), elementwise
\(\mathrm{Im}(\varepsilon_{xy})\) is **not** equal to
\((( \varepsilon-\varepsilon^\dagger)/(2i))_{xy}\). Elementwise Im over-counts
dissipation when the design drives strong plasma fields.

**Action:** fix absorption post-processing only (`validate_modal_s.material_absorption`
→ use \(D=(\varepsilon-\varepsilon^\dagger)/(2i)\)). Do **not** modify the validated
direct-DtN S package. Keep the 5% power-balance gate.

## Artifacts

- `power_balance_rootcause.py` / `power_balance_rootcause_nohup.out` (C complete)
- Fix: `scripts/circulator_inverse_design/validate_modal_s.py`
- Validation: `ABSORPTION_FIX_VALIDATION.md` (when complete)
