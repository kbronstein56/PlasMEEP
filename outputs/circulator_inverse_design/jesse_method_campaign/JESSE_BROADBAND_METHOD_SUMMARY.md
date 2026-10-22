# Jesse broadband method summary

All six advisor questions answered with measured numbers.
**Production BB spec is frozen** (`BROADBAND_PRODUCTION_SPEC_FROZEN.md`).
Long 4–6 GHz optimization is **gated** by that freeze (not auto-started here).

Continuum SF quote (not M-mesh): \(J\approx0.626\pm0.007\), iso \(\approx16.9\pm0.13\) dB.

---

### 1. How many frequency points?

No universal \(N\). Frozen bars failed for 200→100 (Δiso 8.9 dB) and 100→50
(Δiso 19.8 dB). Finest audit: worst iso **−52.8 dB @ 5.75 GHz**. Production
starts at **100 MHz (21 pts)** + adaptive enrichment from ≤50/25 MHz audits.

### 2. How are narrow between-grid failures prevented?

Fixed FoM set per L-BFGS → dense audit → add worst freqs → restart. Soft-min
alone is insufficient (Gedeon 2026 + our −53 dB hole).

### 3. How is broadband differentiable?

\(J_{BB}=\alpha\mathrm{mean}+(1-\alpha)\mathrm{softmin}_\tau\), τ=0.05, α=0.25.
E2E design FD **PASS** (rel \(3.08\times10^{-4}\)).

### 4. JAX/autograd vs current adjoint?

Ceviche = custom VJP = adjoint. Tiny JAX `custom_vjp` PASS. No FULL91 JAX port.

### 5. Is the current method a frequency-domain adjoint?

**Yes** — discrete bilinear \(A^\top\).

### 6. What dominates compute time?

Historical: LU **~84%** (1008 s of 1195 s). This campaign: adjoint 126 s
(dA 37 s, LU 9 s). Forward-mode 150 s est., agree **2×10⁻¹⁴**; prefer adjoint.
Adjoint RHS ≪ factorization.

### 7. Does mode preservation matter?

Short opts: A J_A=0.727 (des1≈0.9%); B J_B=0.755 (des1≈2.1%, worse TE0 iso).
**Freeze Objective A.**

### 8. Can mode content be measured experimentally?

Not with VNA/horn alone. Near-field / dual-probe / measurement-aware receive FoM.

### 9. Why C6 versus C3?

C6 = six-port 60° choice; Fan is C3 for three ports. Short opt: C6 **0.727** vs
C3 **0.444**. **Freeze C6.** Keep six sources (~5% rotation residual on mesh).

### 10. How will time-domain validation be done?

Meep 4–6 GHz pulse after gyrotropic gates; compare to FEM; do not tune FEM to Meep.

### 11. Which papers most affected choices?

Wang & Fan 2005; Hughes 2019; Hammond 2022; Gedeon 2026; Zhao 2025 / Wang PRL 2021;
Mehrpour Bernety 2022; Li PRApplied 2022.

Log: `docs/circulator_literature/LITERATURE_LOG.md`.
Snapshot: `CAMPAIGN_RESULTS_SNAPSHOT.md`.
