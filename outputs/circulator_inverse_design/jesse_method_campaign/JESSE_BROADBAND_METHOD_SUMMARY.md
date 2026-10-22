# Jesse broadband method summary

All six advisor questions answered with measured numbers.
**Production BB spec is frozen** (`BROADBAND_PRODUCTION_SPEC_FROZEN.md`).
25 MHz upper-band audit **PASS** vs 50 MHz (worst iso/IL stay at 5.75 GHz, Δ=0).

Continuum SF quote (not M-mesh): \(J\approx0.626\pm0.007\), iso \(\approx16.9\pm0.13\) dB.

---

### 1. How many frequency points?

No universal \(N\). Full-band: 200→100 and 100→50 fail frozen bars. Targeted
**25 MHz on 5.50–6.00 GHz** vs band-restricted 50 MHz: worst iso/IL both at
**5.75 GHz** with **Δ = 0** → **PASS**. Production: **100 MHz (21 pts)** initial
objective + adaptive ≤50/25 MHz enrichment.

### 2. How are narrow between-grid failures prevented?

Fixed FoM set per L-BFGS → dense audit → add worst freqs → restart. Soft-min
alone is insufficient.

### 3. How is broadband differentiable?

\(J_{BB}=\alpha\mathrm{mean}+(1-\alpha)\mathrm{softmin}_\tau\), τ=0.05, α=0.25.
E2E design FD **PASS** (rel \(3.08\times10^{-4}\)).

### 4. JAX/autograd vs current adjoint?

Ceviche = custom VJP = adjoint. Tiny JAX `custom_vjp` PASS. No FULL91 JAX port.

### 5. Is the current method a frequency-domain adjoint?

**Yes** — discrete bilinear \(A^\top\).

### 6. What dominates compute time?

**Latest** FULL91-C efficient adjoint (~126 s): assembly **77.8 s**, LU **9.0 s**,
adjoint RHS **1.8 s**, dA **36.8 s**. Dominant costs are **assembly + dA**, not
adjoint RHS. Older “LU ≈84%” is **historical only**. Forward-mode agrees to
**2×10⁻¹⁴** but is slower for scalar J.

### 7. Does mode preservation matter?

Short opts: freeze **Objective A** (mode-preserving).

### 8. Can mode content be measured experimentally?

Not with VNA/horn alone. Near-field / dual-probe / measurement-aware receive FoM.

### 9. Why C6 versus C3?

C6 = six-port 60° choice; Fan is C3 for three ports. Freeze **C6**.

### 10. How will time-domain validation be done?

Meep 4–6 GHz pulse after gyrotropic gates; compare to FEM; do not tune FEM to Meep.

### 11. Which papers most affected choices?

Wang & Fan 2005; Hughes 2019; Hammond 2022; Gedeon 2026; Zhao 2025 / Wang PRL 2021;
Mehrpour Bernety 2022; Li PRApplied 2022.

Log: `docs/circulator_literature/LITERATURE_LOG.md`.
