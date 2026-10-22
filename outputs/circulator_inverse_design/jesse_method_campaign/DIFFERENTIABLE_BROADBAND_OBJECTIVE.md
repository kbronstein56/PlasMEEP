# Differentiable broadband objective

Extends `../broadband_method/BROADBAND_DIFFERENTIABLE_OBJECTIVE.md`.

## Definition

Per-frequency scalar \(J(f)\) from the chosen single-frequency FoM (A or B).

Normalized soft-min + mean:

\[
J_{\mathrm{softmin}}
= -\tau\log\Bigl(\tfrac{1}{N}\sum_i e^{-J_i/\tau}\Bigr)
\quad\text{(max-shift stable log-sum-exp)},
\]

\[
J_{\mathrm{BB}}
= \alpha\,\mathrm{mean}(J_f) + (1-\alpha)\,J_{\mathrm{softmin}}.
\]

Frozen defaults (until production spec): \(\tau=0.05\), \(\alpha=0.25\).

| Symbol | Meaning |
|--------|---------|
| \(\tau\) | Softness; smaller → closer to worst-\(J\) |
| \(\alpha\) | Weight on mean vs soft-min (higher α → more average performance) |
| Sign | We **maximize** \(J\) / \(J_{\mathrm{BB}}\) (same as SF) |

Weights \(w_i=\partial J_{\mathrm{BB}}/\partial J_i\) → chain rule
\(\nabla_q J_{\mathrm{BB}}=\sum_i w_i\nabla_q J(f_i)\).

## Validations

| Check | Result |
|-------|--------|
| Soft-min weight FD | ~5.7e-11 PASS (prior) |
| τ/α sensitivity on fixed spectra | `../broadband_method/SOFTMIN_TAU_ALPHA_SENSITIVITY.md` |
| **End-to-end multi-freq design FD** | `softmin_e2e_fd.json` (3-freq set; this campaign) |

Do **not** choose τ/α from which setting gives the nicest production J.

## Status

E2E FD script: `softmin_e2e_fd.py`. Require relative error < 1e-3 on tested orbits.
