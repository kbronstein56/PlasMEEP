# Jesse broadband method summary

Plain-language answers. Numbers from validated benches; numerical campaign
artifacts fill remaining gaps. **Long 4–6 GHz optimization has not started.**

Continuum SF quote (not M-mesh): \(J\approx0.626\pm0.007\), isolation
\(\approx16.9\pm0.13\) dB (`MESH_CONVERGENCE_MREFINED_DESIGN.md`).

---

### 1. How many frequency points?

Not a fixed universal \(N\). Staged recon on 4–6 GHz (200→100→50 MHz) against
**frozen** grid bars (Δiso < 0.5 dB, ΔIL/RL < 0.2 dB). Soft-min + adaptive
enrichment; dense audit for validation. See `BROADBAND_FREQUENCY_SAMPLING_STUDY.md`.

### 2. How are narrow between-grid failures prevented?

Fixed objective set per L-BFGS run → dense audit → add worst frequencies →
restart. Soft-min helps but does **not** replace dense audit (Gedeon 2026).

### 3. How is broadband differentiable?

\(J_{BB}=\alpha\,\mathrm{mean}(J_f)+(1-\alpha)\,\mathrm{softmin}_\tau(J_f)\) with
chain rule through per-frequency frequency-domain adjoints. Defaults
\(\tau=0.05\), \(\alpha=0.25\). Soft-min weight FD ~5.7e-11; E2E design FD in
`softmin_e2e_fd.json`.

### 4. JAX/autograd vs current adjoint?

Ceviche = custom VJP around sparse solve (= adjoint). JAX/`custom_vjp` would not
remove LU cost on FULL91. Keep analytic adjoint; optional API wrapper only.
See `JAX_AUTOGRAD_CEVICHE_ASSESSMENT.md`.

### 5. Is the current method a frequency-domain adjoint?

**Yes.** Discrete FD adjoint with bilinear \(A^\top\) (validated; do not casually
switch to \(A^H\)).

### 6. What dominates compute time?

**LU factorization ~84%** (~1008 s of ~1195 s efficient path on FULL91-C).
Six adjoint RHS ≈ 3 s after LU — **not** another full simulation.
Baseline 1397 s → efficient 1195 s (1.17×).

### 7. Does mode preservation matter?

Open scientifically. Objective A (TE0-preserving) vs B (port-power total) with
**frozen** weights; short opts underway. Literature: mode-converting circulators
exist (Zhao 2025, Wang 2021).

### 8. Can mode content be measured experimentally?

Not with ordinary VNA/horn alone. Feasible via near-field scan or dual-probe;
measurement-aware receive-mode FoM is attractive. See
`EXPERIMENTAL_MODE_RESOLUTION_FEASIBILITY.md`.

### 9. Why C6 versus C3?

C6 matches six-port 60° symmetry (16 vars) — a **choice**, not a law. Fan’s
device is C3 because it has **three** ports. C3 (~31 vars) under test.

### 10. How will time-domain validation be done?

Meep broadband pulse 4–6 GHz after small gyrotropic gates; compare spectra to
FEM; do **not** tune FEM to Meep. `TIME_DOMAIN_BROADBAND_VALIDATION.md`.

### 11. Which papers most affected choices?

| Paper | Impact |
|-------|--------|
| Wang & Fan OL 2005 | C3 three-port classic → C6≠required |
| Hughes et al. ACS Photonics 2019 | Forward vs reverse mode scaling |
| Hammond et al. OE 2022 | TD/FD hybrid broadband TO (Meep) |
| Gedeon et al. 2026 preprint | Spike-hiding → dense audit + enrichment |
| Zhao et al. ACS Photonics 2025 | Mode-converting circulation → Objective B |
| Wang et al. PRL 2021 | Modal circulation channels |
| Mehrpour Bernety et al. 2022 | Plasma microwave circulator context |

Literature log: `docs/circulator_literature/LITERATURE_LOG.md`.

---

**Gate:** `BROADBAND_PRODUCTION_SPEC_FROZEN.md` must exist before long BB opt.
