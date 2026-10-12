# Modal-S Onsager audit on the Pref sixfold peak design

**Verdict: `MODAL_S_OPTIMIZED_ONSAGER = FAIL`**

Frozen criterion: `ONSAGER_CRITERION_FROZEN.md` (written **before** these experiments; not modified).

## Frozen thresholds (recap)

| # | Requirement | Result on peak @ 3.85 GHz |
|---|-------------|---------------------------|
| 1 | `r_F ≤ 5×10⁻³` on finest mesh | **FAIL** — C/M/F all ≈ 0.039–0.041 |
| 2 | Significant-channel \|ΔS\|, Δ\|S\| ≤ 0.02 | **FAIL** — worst \|ΔS\|≈0.028, Δ\|S\|≈0.022 on desired channels |
| 3 | Monitor-shift \|S\| invariance ≤ 0.015 | **PASS** — Δ≤0.008 (C), 0.006 (M) |
| 4 | Mesh trend C→M decreasing or both ≤5e-3 | **FAIL** — C 0.0390 → M 0.0400 → F 0.0407 (no improvement) |

## What is *not* the cause

| Hypothesis | Evidence |
|------------|----------|
| Spatial mesh error | C/M/F residuals identical to ~0.001 |
| Monitor station | ±0.3a and +0.7a shifts leave `r_F`≈0.039; \|S\| stable |
| Sample count | 41→81: Δ\|S\|≤0.0003 |
| ∂n noise / one-plane vs two-plane | Two-plane (d=0.25, 0.40) unchanged |
| Port phase gauge | Optimizing 12 port phases does not reduce `r_F` |
| Reaction / FEM operator | Reaction Onsager `r_F ≈ 9×10⁻¹⁴` on the same design |
| Weak PML σ | `PML_SIGMA_SCALE` 1→30: no improvement |

## Root cause (established)

### 1. Propagating TE1 in the production feed (dominant)

At 3.85 GHz: `k0≈1.614`, `clear_width=2` ⇒ `kc₁=π/W≈1.571 < k0` ⇒ **TE1 propagates** (`β₁≈0.37`, `λ₁≈17a`).

PML thickness ≪ `λ₁`, so TE1 is poorly absorbed at the outer boundary.

| Frequency | TE1 | Peak `r_F` (complex) | Peak \|S\|-only `r_F` | Uniform `r_F` |
|-----------|-----|---------------------:|----------------------:|--------------:|
| 3.70 GHz | **cutoff** | 0.0074 | **0.0043** | 1.98×10⁻⁴ |
| 3.75 GHz | just on | 0.0067 | — | 2.05×10⁻⁴ |
| 3.80 GHz | on | 0.0088 | — | 4.05×10⁻⁴ |
| 3.85 GHz | on | **0.0390** | **0.0328** | 1.84×10⁻⁴ |

Removing TE1 (frequency below cutoff) collapses the peak residual by ~5×. Uniform stays excellent at all frequencies because TE0↔TE1 conversion is tiny.

Desired-channel check at 3.85 GHz: `|S_{j+1,j}(+B)|≈0.859` vs `|S_{j,j+1}(−B)|≈0.839` (Δ≈0.021) — same order as TE0→TE1 conversion power (~0.018 mean `|S|²`).

**Physics:** For a multimode junction, only the **full** mode-port S satisfies Onsager. The fundamental block alone need not when TE1 is open. Reaction Onsager remains exact; fundamental modal S is a truncation.

### 2. Residual ~7×10⁻³ floor with TE1 cut off

Even at 3.70 GHz, peak complex `r_F≈0.0074` (above 5×10⁻³). Magnitude Onsager is 0.0043 with worst significant Δ\|S\|=0.001 — **power metrics are fine**; leftover error is mostly phase/near-field. F mesh does not reduce it. Likely soft-source + flare reactive fields on the optimized (nonuniform) plasma.

### 3. TE1 feed constrictions

Symmetric feed constrictions (`W_filter≈1.45–1.6`) reduce TE1 at the monitor and can cut `r_F` toward ~0.014–0.015, but:
- they perturb the device (`T_des` changes);
- they still do **not** meet 5×10⁻³;
- aggressive filters destroy circulation.

Not accepted as a production fix without remeshing/re-optimization under the filter.

### 4. Multimode 12×12 S

Building TE0+TE1 drives and `S=B A⁻¹` does **not** yet yield `r_F≤5×10⁻³` (uniform ~0.046, peak ~0.12). TE1 soft-source / thin-PML columns are not trustworthy. Fixing that needs thicker PML and/or true waveguide port BCs (out of scope for a quick patch).

## Significant-channel residuals (FULL91-M, 3.85 GHz)

Worst offenders are the **desired** transmissions `|S|≈0.86` and reflections `|S|≈0.24`, not tiny leak channels — so this is not a noise-floor artifact.

## Implications for optimization gate

Frozen criterion **FAIL**. Per campaign rules: **do not start true modal-S inverse design** until extraction meets the frozen bar (or a *new* criterion is frozen *before* further experiments for a revised extractor, without relaxing 5×10⁻³ on the same definition).

Engineering note (not a pass): monitor-invariant `|S|²` metrics on the Pref peak remain useful as Pref-campaign characterization; they are not Onsager-clean at 3.85 GHz.

## Artifacts

- `onsager_audit_summary.json`, `onsager_audit_*_{C,M}.json`, `onsager_audit_F.json`
- `extractor_variants_onsager.json`, `pml_scale_onsager_C.json`
- `multimode_onsager_*.json`, `te1_filter_*.json`
- `ONSAGER_CRITERION_FROZEN.md`
