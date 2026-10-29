# A2 Soft-min Scale Sensitivity (offline)

**Date:** 2026-10-28  
**Scope:** Offline weighting only. No FEM. Objective A2_hybrid unchanged. Stage-2 **not** launched.  
**Spectra:** Stage-1 design rescored with A2_hybrid (`opt_broadband_production_C` 21-grid + dense audit + proposed 28-freq set).  
**Machine-readable:** `A2_SOFTMIN_SCALE_SENSITIVITY.json`

Aggregation (unchanged formula):

\[
J_{BB}=(1-\alpha)\,J_{\mathrm{softmin}}+\alpha\,\mathrm{mean}(J_f),\quad
J_{\mathrm{softmin}}=-\tau\log\!\Big(\mathrm{mean}_i e^{-J_i/\tau}\Big)
\]

Combined gradient weights: \(w_i=(1-\alpha)\,p_i+\alpha/n\) with softmax \(p_i=\partial J_{\mathrm{softmin}}/\partial J_i\).

---

## 1. Score-scale fact: A2 is not on Objective A’s scale

| Spectrum | Objective | min \(J\) | range | gap worst→2nd | gap worst→median |
|----------|-----------|----------:|------:|--------------:|-----------------:|
| Stage-1 21-grid | A | −0.657 | 0.52 | **0.029** | 0.23 |
| Stage-1 21-grid | A2 | −11.56 | 13.9 | **4.43** | 8.68 |
| Proposed 28-freq | A2 | −11.56 | 13.9 | **3.25** | 7.39 |

Range ratio A2/A ≈ **27×**. Worst-to-2nd gap ratio ≈ **150×** on the 21-grid.

At frozen \(\tau=0.05\), Objective A still had multi-frequency soft-min:

| FoM @ τ=0.05, α=0.25 | \(\max p\) | \(N_{\mathrm{eff}}(p)=1/\sum p_i^2\) | \(\max w\) | \(N_{\mathrm{eff}}(w)\) |
|----------------------|----------:|-------------------------------------:|----------:|----------------------:|
| Objective A (21) | 0.35 | **4.96** | 0.27 | 7.46 |
| A2 (21) | **1.00** | **1.00** | 0.76 | 1.71 |
| A2 (28) | **1.00** | **1.00** | 0.76 | 1.73 |

So the FD 3-frequency example with combined weights ≈ `[0.083, 0.083, 0.833]` is **not** a soft soft-min: pure softmax is one-hot on the worst frequency; the 0.083 floor is exactly \(\alpha/n\).

**Answer to the scale concern:** \(\tau=0.05\) *is* much sharper for A2 than it was for A, **because logarithmic A2 gaps are O(1–10) while A gaps were O(0.03)**. That is a real scale mismatch, not an illusion.

---

## 2. τ sweep at α = 0.25 (proposed 28-freq A2 scores)

Worst frequency is always **4.1 GHz** (\(J_{A2}\approx-11.56\)).

| τ | \(J_{BB}\) | worst A2 | \(\max p\) (soft-min) | \(\max w\) (combined) | \(N_{\mathrm{eff}}(p)\) | \(N_{\mathrm{eff}}(w)\) | \(p>0.99\)? | \(w>0.75\)? | \(w>0.90\)? |
|--:|----------:|---------:|----------------------:|----------------------:|------------------------:|------------------------:|:-----------:|:-----------:|:-----------:|
| 0.05 | −9.557 | −11.558 | 1.000 | 0.759 | 1.00 | 1.73 | yes | yes | no |
| 0.10 | −9.432 | −11.558 | 1.000 | 0.759 | 1.00 | 1.73 | yes | yes | no |
| 0.25 | −9.057 | −11.558 | 1.000 | 0.759 | 1.00 | 1.73 | yes | yes | no |
| 0.50 | −8.433 | −11.558 | 0.997 | 0.757 | 1.01 | 1.74 | yes | yes | no |

Soft-min top mass at τ=0.05 (28-freq):

| f (GHz) | softmax \(p\) |
|--------:|--------------:|
| **4.1** | **1.000** |
| 5.85 | ~10⁻²⁸ |
| 5.65 | ~10⁻³¹ |

At τ=0.50, runner-up 5.85 GHz still has \(p\approx0.0015\) only.

**Within {0.05, 0.10, 0.25, 0.50}, τ does not change optimization weights in any material way** while the Stage-1 A2 outlier gap (~3.25 to 2nd-worst) persists. Only \(J_{BB}\)’s numerical value softens (less negative) as τ grows — that is *not* a reason to pick τ.

### Top-5 combined weights (α=0.25)

When soft-min is one-hot, every non-worst frequency gets the **same** floor \(\alpha/n\approx0.00893\). Apparent “top-5” after 4.1 at τ=0.05 is therefore a **tie artifact** (sorted frequency order), not soft-min sharing.

| τ | Top frequencies by combined \(w\) | Notes |
|--:|-----------------------------------|-------|
| 0.05 | 4.1 (0.759), then tied floor on others | soft-min one-hot |
| 0.10 | 4.1 (0.759); tiny soft-min bleed → 5.85, 5.65 appear in sort | still \(\max p=1\) |
| 0.25 | same structure | still one-hot for practical purposes |
| 0.50 | 4.1 (0.757); 5.85 gets \(p\approx0.0015\) | still \(\max p>0.99\) |

Aggregate **ranking of the worst frequency never changes** (always 4.1). Soft-min internal ranking among #2…#5 only becomes nonzero for τ ≳ 0.5, and remains negligible through τ=0.5.

Same pattern on the 21-grid and dense 52-point set: τ∈[0.05,0.50] → \(\max p\approx1\), \(\max w\approx0.76\), \(N_{\mathrm{eff}}(w)\approx1.7\).

---

## 3. α cross-sweep (28-freq A2; representative τ)

Because soft-min is one-hot for all tested τ, **α alone** sets the multi-frequency floor:

\[
\max w \approx (1-\alpha)+\alpha/n,\qquad N_{\mathrm{eff}}(w)\approx 1\big/\big([(1-\alpha)+\alpha/n]^2+(n-1)(\alpha/n)^2\big)
\]

| τ | α | \(J_{BB}\) | \(\max w\) | \(N_{\mathrm{eff}}(w)\) | \(w>0.75\) | \(w>0.90\) | Soft-min still one-hot? |
|--:|--:|----------:|----------:|------------------------:|:----------:|:----------:|:-----------------------:|
| 0.05 | 0.10 | −10.71 | **0.904** | 1.22 | yes | **yes** | yes |
| 0.05 | **0.25** | −9.56 | **0.759** | **1.73** | yes | no | yes |
| 0.05 | 0.50 | −7.80 | 0.518 | 3.73 | no | no | yes |
| 0.50 | 0.10 | −9.44 | 0.901 | 1.23 | yes | yes | nearly |
| 0.50 | 0.25 | −8.43 | 0.757 | 1.74 | yes | no | nearly |
| 0.50 | 0.50 | −6.74 | 0.517 | 3.74 | no | no | nearly |

α=0.10 → single frequency owns >90% of the combined gradient (true one-freq chase).  
α=0.25 → worst owns ~76%; every other frequency still receives a uniform α/n share.  
α=0.50 → softer mean-heavy blend; weakens worst-case pressure.

---

## 4. What τ would restore Objective-A-like soft-min Neff on *this* A2 spectrum?

On the proposed 28-freq A2 scores (gap worst→2nd ≈ 3.25):

| Target \(N_{\mathrm{eff}}(p)\) | Approximate τ needed | \(\max p\) at that τ |
|-------------------------------:|---------------------:|---------------------:|
| ≥ 2 | ~1.5 | ~0.68 |
| ≥ 3 | ~1.75 | ~0.56 |
| ≥ 5 (≈ Objective A’s Neff) | ~2.3 | ~0.38 |

Matching Objective A’s *worst→2nd gap sharpness* formally would need \(\tau_{A2}\approx\tau_A\times(\Delta_{A2}/\Delta_A)\approx 7.5\) — far outside a soft-min that still prefers the catastrophic outlier.

**Important dynamics note:** if Stage-2 lifts 4.1 until it ties the former 2nd-worst, then at τ=0.05 soft-min immediately becomes two-way (\(N_{\mathrm{eff}}(p)\approx2\)). Soft-min sharpness is gap-driven; it is not permanently locked to one frequency after the pathology is fixed.

---

## 5. Decision

### Is τ=0.05 intentionally appropriate for A2 (approximate minimax)?

**Yes — as the soft-min temperature for Stage-2 start — with the following precise reading:**

1. **Scale mismatch is real.** τ=0.05 is “too sharp” *relative to how it behaved under Objective A*, solely because A2’s log terms create O(1–10) score gaps. The FD example’s ~83% combined weight on the worst frequency is the expected one-hot soft-min + α/n floor, not a mild soft-min.

2. **Within the tested set {0.05, 0.10, 0.25, 0.50}, changing τ does not buy multi-frequency soft-min sharing** on present Stage-1 A2 spectra. Amending the freeze to τ=0.10 or 0.25 would be a **no-op for L-BFGS weights** while 4.1 remains ~3 units worse than #2. That is not a scientific improvement.

3. **Desired optimization behavior** is robust worst-case broadband performance with enough smoothness that L-BFGS does not *only* chase one frequency:
   - Soft-min one-hot on a **catastrophic** outlier (IL~40 dB at 4.1 GHz) is the correct worst-case pressure for A2.
   - Multi-frequency smoothness today comes from **α=0.25**, which caps combined worst weight at \(\approx(1-\alpha)+\alpha/n\sim76\%\) and keeps \(N_{\mathrm{eff}}(w)\approx1.7\) with a uniform α/n floor on all other frequencies — not from τ.
   - Restoring soft-min-internal Neff≳3 *while the outlier gap remains* would require τ~O(1–2), which deliberately dilutes pressure on the dead channel before it is fixed. That fights the Stage-2 engineering goal.

4. **Do not retune τ to make \(J_{BB}\) look prettier.** Larger τ makes \(J_{BB}\) less negative without changing who is optimized when gaps are large.

### Recommendation

| Parameter | Stage-2 freeze | Action |
|-----------|----------------|--------|
| **τ = 0.05** | keep | **No amendment.** Documented intentional approximate-minimax soft-min under A2’s log scale; α provides the multi-freq floor. |
| **α = 0.25** | keep | Still the smoothness knob; α=0.10 would be a true >90% single-freq chase and is rejected. |
| Objective A2_hybrid | keep | unchanged |
| 28-freq set | keep | unchanged |

`BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md` is **left unchanged** (no versioned τ amendment).

---

## 6. Practical takeaway for Stage-2 launch (when approved)

Expect early L-BFGS steps to look nearly minimax on 4.1 GHz in the soft-min channel, with a steady α/n pull from all 28 frequencies. That is consistent with A2’s purpose (do not hide dead transmission). After the worst A2 gap falls to ≲τ, soft-min automatically redistributes — same formula, healthier spectrum.

Do **not** launch Stage-2 from this document alone; this only clears the τ/α scale concern for the existing freeze.
