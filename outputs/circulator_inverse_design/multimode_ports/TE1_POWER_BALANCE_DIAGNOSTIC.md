# TE1 power balance diagnostic

Design `opt_single_uniform_M_Bp0.0500`, grade C, bar 0.05.
Boundary measure in production S: **arclength** (`ds`).

## Verdict (Phase 1–2)

| Item | Result |
|------|--------|
| Failing channel | **ch=7, port=3, TE1** |
| Residual | **0.0673** (bar 0.05) |
| Failure signature | `Σ\|S\|² = 1.0357` (over-unity outgoing) + `P_mat = 0.0315` → identity miss **−6.7%** |
| TE0 all ports | PASS (max residual 0.0234) |
| Other TE1 ports | PASS (max residual 0.0364 on port 4) |
| Tracks port orientation alone? | **No** — opposite port 0 TE1 residual 0.0286 PASS |
| Tracks max staircase length? | **No** — port 4 has largest `Σds/W=1.246` but TE1 residual 0.0364 PASS |
| Tracks projected aperture miss? | Weak — all ports have `Σ\|e·t\|/W ∈ [0.944, 0.981]` |
| FEM energy identity | `Im(xᴴ A x) = Im(xᴴ f)` closes for ch7 (−11.027); discrete operator self-consistent |
| Likely cause | Staircase `ds` measure inflates TE1 modal loads / S power relative to flat-port `dy`; port 3 is a bad mesh-registration draw (not uniquely worst geometry, but TE1-sensitive combination of staircase + field) |

## Per-channel residuals

| ch | port | mode | Σ\|S\|² | TE0 | TE1 | P_mat | residual | PASS |
|----|------|------|--------|-----|-----|-------|----------|------|
| 0 | 0 | 0 | 0.9859 | 0.6010 | 0.3849 | 0.0363 | 0.0222 | True |
| 1 | 0 | 1 | 1.0000 | 0.2275 | 0.7725 | 0.0286 | 0.0286 | True |
| 2 | 1 | 0 | 0.9824 | 0.7367 | 0.2457 | 0.0374 | 0.0198 | True |
| 3 | 1 | 1 | 0.9974 | 0.3749 | 0.6224 | 0.0322 | 0.0296 | True |
| 4 | 2 | 0 | 0.9845 | 0.6508 | 0.3337 | 0.0376 | 0.0222 | True |
| 5 | 2 | 1 | 0.9838 | 0.2485 | 0.7353 | 0.0296 | 0.0134 | True |
| 6 | 3 | 0 | 0.9889 | 0.5981 | 0.3909 | 0.0344 | 0.0234 | True |
| 7 | 3 | 1 | 1.0357 | 0.3149 | 0.7208 | 0.0315 | **0.0673** | **False** |
| 8 | 4 | 0 | 0.9813 | 0.6570 | 0.3244 | 0.0342 | 0.0155 | True |
| 9 | 4 | 1 | 1.0034 | 0.3966 | 0.6069 | 0.0329 | 0.0364 | True |
| 10 | 5 | 0 | 0.9812 | 0.7587 | 0.2225 | 0.0398 | 0.0210 | True |
| 11 | 5 | 1 | 0.9808 | 0.3440 | 0.6368 | 0.0314 | 0.0122 | True |

## Port geometry (staircase audit)

Nominal `W = clear_width = 2a = 2.0`. Modal theory assumes flat section with measure `dy`.

| port | n_edges | Σds | Σ\|e·t\| | Σds/W | Σproj/W | y_span | along_std | along [min,max] | TE1 nrm²_ds | TE1 nrm²_proj | \|c1\| | c1ᵀc1 |
|------|---------|-----|---------|-------|---------|--------|-----------|-----------------|-------------|---------------|--------|-------|
| 0 | 50 | 2.4532 | 1.9524 | 1.2266 | 0.9762 | 2.0000 | 0.00885 | [−0.0215, 0.0169] | 1.2104 | 0.9794 | 0.21585 | 0.04659 |
| 1 | 48 | 2.3251 | 1.9252 | 1.1625 | 0.9626 | 2.0000 | 0.00711 | [−0.0133, 0.0192] | 1.1873 | 0.9486 | 0.21251 | 0.04516 |
| 2 | 43 | 2.3400 | 1.9196 | 1.1700 | 0.9598 | 2.0000 | 0.00718 | [−0.0123, 0.0141] | 1.1620 | 0.9510 | 0.22298 | 0.04972 |
| **3** | 47 | 2.4087 | 1.9259 | 1.2043 | 0.9630 | 2.0000 | 0.00792 | [−0.0149, 0.0164] | **1.2215** | 0.9618 | 0.22274 | 0.04961 |
| 4 | 47 | 2.4928 | 1.9612 | **1.2464** | 0.9806 | 2.0000 | 0.00772 | [−0.0166, 0.0155] | 1.2419 | 0.9677 | 0.21894 | 0.04793 |
| 5 | 46 | 2.3347 | 1.8884 | 1.1674 | 0.9442 | 2.0000 | 0.00762 | [−0.0122, 0.0162] | 1.1894 | 0.9249 | 0.21241 | 0.04512 |

Observations:

1. **Every port** has staircase inflation: `Σds/W ≈ 1.16–1.25` while projected measure is slightly **short**: `Σ\|e·t\|/W ≈ 0.94–0.98`.
2. Midpoint deviation from the nominal plane is small (`along_std ~ 0.007–0.009`, range ~±0.02) — not a gross plane-miss.
3. `y_span = W` exactly on all ports — aperture span is correct; the error is in **edge measure**, not missing ends.
4. TE1 raw norm² under `ds` is ~1.16–1.24 (should be 1 for continuous ∫ φ₁² dy); under projected measure it drops to ~0.92–0.98 — closer to theory.
5. Opposite pairs `(0,3)`, `(1,4)`, `(2,5)` have similar staircase stats; **only port-3 TE1 fails**. Failure is therefore **mesh-registration / staircase sampling of φ₁(y)** on that cut, not a global orientation or material asymmetry (TE0 on port 3 still PASS).

## Rotational / causal triage

| Hypothesis | Evidence |
|------------|----------|
| Port orientation | FAIL tracks port 3 only; opposite port 0 TE1 OK → not orientation alone |
| Mesh registration | PASS — port 3 TE1 nrm²_ds is 2nd-highest; staircase edges sample φ₁ unevenly |
| Edge staircase geometry | Global: all ports inflated; local: port 3 not unique worst, but TE1 residual peaks there |
| Mode normalization | TE1 more sensitive than TE0 (constant φ₀); over-unity `Σ\|S\|²` on ch7 |
| Material field asymmetry | Unlikely primary: P_mat~0.03 similar across TE1 channels; FEM `xᴴAx` identity closes |

## Discrete vs field power (ch7)

- `Im(xᴴ A x) = Im(xᴴ f) = −11.027` → discrete anti-Hermitian identity closes.
- Engineering residual uses `1 − Σ\|S\|² − P_material` with S from the same modal `c` as the drive.
- Over-unity `Σ\|S\|²` implies **modal power extraction / normalization** (staircase `ds`) rather than a broken FEM energy identity.

## Phase 3 outcome

Projected measure `dℓ = |e · t̂|` adopted: Onsager ~1e-14, all 12 residuals ≤ 5%
(ch7: 0.0673 → 0.0352). See `PROJECTED_MEASURE_EXPERIMENT.md`.
Package frozen in `MULTIMODE_PORT_PACKAGE_FROZEN.md`.

Artifacts: `te1_power_balance_diagnostic.json`, this file.
