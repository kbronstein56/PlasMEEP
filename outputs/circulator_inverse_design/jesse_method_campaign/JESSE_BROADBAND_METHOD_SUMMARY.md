# Jesse broadband method summary

Plain-language answers with measured numbers. **Long 4–6 GHz optimization has
not started** (`BROADBAND_PRODUCTION_SPEC_STATUS.md`).

Continuum SF quote (not M-mesh): \(J\approx0.626\pm0.007\), isolation
\(\approx16.9\pm0.13\) dB.

---

### 1. How many frequency points?

No universal \(N\). Frozen bars: Δiso < 0.5 dB, ΔIL/RL < 0.2 dB vs denser audit.

**Measured on SF design:** 200 MHz worst iso **−24 dB @ 5.8 GHz**; 100 MHz finds
**−33 dB @ 5.9 GHz**. Refining 200→100 changes worst iso by **8.9 dB** →
**FAIL** sufficiency bars. So 11 points @ 200 MHz are **not** enough even as a
stable objective grid without enrichment; use ≤100 MHz + adaptive 50/25 MHz near
features.

### 2. How are narrow between-grid failures prevented?

Fixed FoM set per L-BFGS → dense audit → add worst freqs → restart. Soft-min
helps but does not replace audit (Gedeon 2026).

### 3. How is broadband differentiable?

\(J_{BB}=\alpha\,\mathrm{mean}(J_f)+(1-\alpha)\,\mathrm{softmin}_\tau(J_f)\),
\(\tau=0.05\), \(\alpha=0.25\). Chain rule through FD adjoints. Soft-min weight
FD ~5.7e-11; E2E design FD queued (`softmin_e2e_fd.py`).

### 4. JAX/autograd vs current adjoint?

Ceviche = custom VJP = adjoint. Tiny JAX `custom_vjp` **PASS** (rel
\(2.5\times10^{-10}\) vs FD). Does **not** beat FULL91 LU/dA cost. Keep analytic
adjoint.

### 5. Is the current method a frequency-domain adjoint?

**Yes** — discrete bilinear \(A^\top\) (validated).

### 6. What dominates compute time?

| Benchmark | LU | Adjoint RHS | Total efficient |
|-----------|---:|------------:|----------------:|
| Historical FULL91-C | **1008 s (~84%)** | ~3 s | 1195 s (1.17× vs baseline 1397 s) |
| This campaign (same grade) | 9 s | 1.8 s | **126 s** (dA 37 s) |

Forward-mode (16×6=96 RHS): **150 s** est., grad agree **2×10⁻¹⁴**. Prefer
**adjoint** at 16/31/91 params for scalar J.

Adjoint RHS ≠ another factorization.

### 7. Does mode preservation matter?

Short equal-budget opts (5 iters, FULL91-C): A reaches J_A=**0.727** with des1≈0.9%;
B reaches J_B=**0.755** with des1≈2.1% and worse TE0 isolation (18 vs 20 dB).
**Default A** for production unless longer B proves broadband value.

### 8. Can mode content be measured experimentally?

Not with ordinary VNA/horn alone. Near-field scan / dual-probe / measurement-aware
receive FoM — see `EXPERIMENTAL_MODE_RESOLUTION_FEASIBILITY.md`.

### 9. Why C6 versus C3?

C6 = choice for six-port 60° tying (16 vars). Fan is C3 because **three** ports.
C3 orbits **31 vars PASS**. Mesh is not exactly C6: rotating one source predicts
others only to **~5%** — keep six sources on FULL91-C.

### 10. How will time-domain validation be done?

Meep pulse 4–6 GHz after gyrotropic gates; compare to FEM; do not tune FEM to
Meep. `TIME_DOMAIN_BROADBAND_VALIDATION.md`.

### 11. Which papers most affected choices?

Wang & Fan 2005 (C3); Hughes 2019 (forward vs adjoint); Hammond 2022 (TD/FD);
Gedeon 2026 (spike hiding); Zhao 2025 / Wang PRL 2021 (mode conversion);
Mehrpour Bernety 2022 / Li PRApplied 2022 (plasma microwave circulators).

Log: `docs/circulator_literature/LITERATURE_LOG.md`.
