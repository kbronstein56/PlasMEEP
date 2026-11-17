# Circulator feasibility baseline — Stage-5 design space

**Status:** diagnostic only  
**Stage-6:** **NOT STARTED / AWAITING AUTHORIZATION**  
**Design:** Stage-5 `q_best` (`opt_broadband_stage5_A2_C.json`, commit `7d517cb`)  
**Physics package:** FULL91-C, B = **+0.05 T**, C6 / 16 vars, bounds **[0.05, 1.8]**, A2_hybrid τ=0.05 α=0.25

Companion: `CIRCULATOR_FEASIBILITY_BASELINE_QSTATS.json`

## 1. Bound saturation (16 C6 variables)

| Stage | near lower (≤0.07) | near upper (≥1.78) | interior | q_min | q_max | mean |
|------:|-------------------:|-------------------:|---------:|------:|------:|-----:|
| 2 | 0 | 0 | 16 | 0.460 | 1.330 | 0.841 |
| 3 | 1 | 0 | 15 | 0.058 | 1.658 | 1.197 |
| 4 | 1 | 0 | 15 | 0.059 | 1.664 | 1.175 |
| **5** | **1** | **0** | **15** | **0.0585** | **1.662** | **1.176** |

Stage-5 has **one** variable near the lower floor (~0.0585) and **none** at the upper wall. The design is **not** strongly bound-saturated. Widening `[0.05, 1.8]` is unlikely to be the primary fix based on this evidence alone.

## 2. Optimizer termination state (Stage-5)

| Quantity | Value |
|----------|------:|
| Accepted iters | 30 / 30 (iteration limit) |
| BB evaluations | 77 |
| Wall | 76571 s ≈ 21.3 h |
| Final grad_norm ‖∇_q J_BB‖ | **≈ 352** |
| Start (nit=1) grad_norm | ≈ 2325 |
| Mid (nit=16) grad_norm | ≈ 160 |

Termination is the **L-BFGS iteration limit**, not a small-gradient stationary point. Soft-min / multi-frequency landscape remains nonsmooth; a large residual gradient does **not** by itself prove physical infeasibility, but it also means Stage-5 did not converge to a calm local minimum.

## 3. q_best similarity across Stages 2–5

| Pair | ‖Δq‖₂ | cosine similarity |
|------|------:|------------------:|
| S2–S3 | 2.63 | 0.889 |
| S3–S4 | 0.140 | 0.9998 |
| S4–S5 | **0.026** | **1.0000** |

Stage-2 → Stage-3 moved the design substantially (new enrichment / different soft-min pressures). Stages **3→4→5 kept essentially the same density pattern** while the frequency grid densified. Stage-5 largely **re-tuned a Stage-3/4-shaped device** rather than discovering a new operating regime.

## 4. Failure-mode composition (Stage-5 dense audit)

Deep failures (A2 ≤ −7) simultaneously show:

- high IL (IL ≥ 15 dB for **all** A2 ≤ −7 samples)
- non-positive isolation (iso ≤ 0 dB for **all** A2 ≤ −7 samples)
- very small desired transmission (true worst T_des ≈ 3.3×10⁻³)

So notches are **not** isolation-only or IL-only — they are **joint transmission zeros / circulation collapses**.

## 5. Spectral distribution of failures

On the Stage-5 3.125 MHz local audit (A2 ≤ −5):

| Band | count |
|------|------:|
| 4.0–4.5 GHz | 50 |
| 4.5–5.0 | 17 |
| 5.0–5.5 | 18 |
| 5.5–6.0 | 15 |

Bad features are **distributed across the full 4–6 GHz window**, denser at the low end on this audit set, not confined to one modal pocket. The four A2 ≲ −8 residuals sit near **4.23, 4.73, 5.08, 5.57 GHz** — roughly a ~0.5 GHz spacing family, suggestive of recurring resonant structure rather than a single accidental hole.

## 6. Resonant character (qualitative)

Prior Stage-3/4/5 localization shows:

- deep cores only a few MHz wide
- rapid recovery of A2 / isolation on either shoulder
- soft-min weight migrating between competing corpses as the grid densifies

This is consistent with **narrow resonant transmission zeros** in a multiport gyrotropic scatterer, not with a smooth wideband mismatch. Field movies at failure frequencies (PlasMEEP-viz Stage-4/5) show disrupted circulation patterns vs repaired nodes — qualitative modal pathology, not a soft global penalty.

## 7. Documented B / hardware context (for later phases)

| Source | B |
|--------|--:|
| Production `B_PLUS` | **+0.05 T** |
| Sketchbook / KNOWN_TRUTHS sign tests | ±0.05 T |
| Related literature already logged (Li et al., PRA 2022; Y-junction plasma circulator) | ~**0.10 T** (104 mT) center field |

Production optimization has only been authorized at **0.05 T**. A higher-B probe is scientifically motivated by logged literature, but must be labeled as **extrapolation beyond the production bias** unless/until hardware authorization expands the range.

## 8. Preliminary feasibility reading (baseline only)

Evidence **against** “Stage-5 failed only because of bound walls”:
- 15/16 variables interior; one soft lower-bound touch.

Evidence **for** “4–6 GHz remains hard for this fixed Stage-3/4/5 density family”:
- S3–S5 q nearly frozen while denser sampling kept finding new A2 ≲ −8 notches;
- failures are joint IL+isolation collapses distributed over the band.

Evidence **that optimizer limits are still real**:
- large final gradient; soft-min / off-grid migration history across Stages 2–5.

**Does not prove impossibility.** Motivates the bandwidth ladder, B probe, and DOF probe in the feasibility study.

## 9. Stage-5 notch-width characterization (added)

Forward-only local scans on Stage-5 `q_best` (`BB_STAGE5_NOTCH_WIDTHS.json`):

| Region | min f (GHz) | min A2 | A2&lt;−5 width (MHz) | A2&lt;−7 width (MHz) | notes |
|--------|------------:|-------:|--------------------:|--------------------:|-------|
| 4p734 | 4.734375 | −8.74 | 18.8 | **4.7** (n=4) | matches Stage-5 true-worst seed |
| 4p228 | **4.189063** | **−9.22** | 4.7 | **2.3** (n=4 @ 0.781 MHz) | deeper neighbor of ~4.228 family |
| 5p084 | 5.084375 | −8.24 | 23.4 | **1.56** (n=2) | deep core barely resolved |
| 5p572 | **5.576563** | **−11.91** | 3.9 | **3.1** (n=5 @ 0.781 MHz) | **missed by 3.125 MHz audit** |

**Narrowest reliably resolved important (A2&lt;−7, ≥2 samples) width ≈ 1.56 MHz** (5.084 family).  
A still-deeper floor at **5.576563 GHz (A2≈−11.9)** shows that even the Stage-5 3.125 MHz dense audit was **not** a true continuum bound.

This strengthens the feasibility question: failures are **few-MHz resonant collapses**, distributed in frequency, while `q` remains nearly frozen from Stage-3 onward.

```
STAGE-6: NOT STARTED / AWAITING AUTHORIZATION
```
