# Differentiable broadband objective

## Definition

\[
J_{\mathrm{softmin}}
= -\tau\log\Bigl(\tfrac{1}{N}\sum_i e^{-J_i/\tau}\Bigr)
\quad\text{(max-shift stable)},
\]

\[
J_{\mathrm{BB}}
= \alpha\,\mathrm{mean}(J_f) + (1-\alpha)\,J_{\mathrm{softmin}}.
\]

Defaults: \(\tau=0.05\), \(\alpha=0.25\). Maximize \(J_{BB}\).

\(\nabla_q J_{BB}=\sum_i w_i\nabla_q J(f_i)\) with
\(w_i=\partial J_{BB}/\partial J_i\).

## Validations

| Check | Result |
|-------|--------|
| Soft-min weight FD | ~5.7e-11 PASS (prior) |
| τ/α sensitivity (fixed spectra) | `../broadband_method/SOFTMIN_TAU_ALPHA_SENSITIVITY.md` |
| **E2E multi-freq design FD** | **PASS** (`softmin_e2e_fd.json`) |

E2E details (FULL91-C, 3 freqs 4/5/6 GHz, orbit 1, h=3e-4):

| Quantity | Value |
|----------|------:|
| J_BB | −0.643 |
| Soft-min weight FD | PASS |
| Orbit-1 rel err analytic vs FD | **3.08×10⁻⁴** (< 1e-3 bar) |

Do **not** retune τ/α from production niceness.
