# BB Objective Alignment Audit (Phase 1–3)

**Date:** 2026-10-27  
**Design:** `opt_broadband_production_C.json` (stage-1)  
**Spectra:** `bb_production` rec_best (21×100 MHz) + `bb_final_audit_50MHz.json` + `bb_final_audit_25MHz_local.json`  
**Machine-readable:** `BB_OBJECTIVE_ALIGNMENT_OFFLINE.json`  
**No FEM in this audit.**

---

## 1. Why does 4.1 GHz (IL≈40 dB, isolation≈−19 dB) get only J_A≈−0.43?

Production Objective A (maximize):

\[
J_A = T_\mathrm{des} - 2\,T_\mathrm{rev} - 1.5\,T_\mathrm{thr} - T_\mathrm{other} - 0.25\,R - T_\mathrm{TE1}
\]

At **4.1 GHz** (already on the objective grid):

| Quantity | Value |
|----------|------:|
| \(T_\mathrm{des}\) | 3.36×10⁻⁴ |
| \(T_\mathrm{rev}\) | 7.86×10⁻³ |
| IL | 39.94 dB |
| isolation | −18.80 dB |
| \(R\) | 0.692 |
| TE1 | 0.220 |
| **J_A** | **−0.428** |

Term decomposition:

| Term | Contribution |
|------|-------------:|
| \(+T_\mathrm{des}\) | +0.0003 |
| \(−2\,T_\mathrm{rev}\) | −0.0157 |
| \(−1.5\,T_\mathrm{thr}\) | −0.0033 |
| \(−T_\mathrm{other}\) | −0.0164 |
| \(−0.25\,R\) | −0.173 |
| \(−\mathrm{TE1}\) | −0.220 |
| **sum** | **−0.428** |

**Mechanism:** IL ~40 dB means desired power is ~10⁻⁴. Isolation −19 dB means reverse is only ~20× larger than desired — still **tiny in absolute power**. Absolute reverse/through penalties therefore nearly vanish. The score is dominated by reflection and TE1 (~−0.4), which look “merely bad,” not catastrophic.

Contrast **5.625 GHz** (off-grid dense audit, worst J_A = −1.153): \(T_\mathrm{des}\approx0.037\), \(T_\mathrm{rev}\approx0.490\). Reverse absolute power is large → \(−2\,T_\mathrm{rev}\approx−0.98\) dominates. Engineering quality is poor but **not** as dead as 4.1 GHz; Absolute-A ranks it worse because reverse power is large, not because transmission is useful.

**Conclusion:** Absolute-power Objective A allows both desired and reverse powers to collapse toward zero without a sufficiently severe penalty. Frequency enrichment alone cannot fix that pathology.

---

## 2. Other frequencies where J_A ranking disagrees with engineering quality

On the 21-point objective grid, largest |rank(IL) − rank(J_A)| gaps (rank 1 = worst):

| f (GHz) | IL (dB) | iso (dB) | J_A | rank IL | rank J_A | gap |
|--------:|--------:|---------:|----:|--------:|---------:|----:|
| 5.9 | 18.0 | −1.4 | −0.281 | 5 | 19 | −14 |
| 4.5 | 17.8 | +8.9 | −0.319 | 6 | 18 | −12 |
| **4.1** | **39.9** | **−18.8** | **−0.428** | **1** | **11** | **−10** |
| 5.2 | 11.4 | +2.2 | −0.599 | 14 | 4 | +10 |
| 4.7 | 6.6 | +3.9 | −0.480 | 18 | 9 | +9 |

Pattern:

- High-IL / low-\(T_\mathrm{des}\) points are **under-penalized** by J_A (4.1, 5.9, 4.5, 5.8).
- Moderate-transmission points with large TE1 / leakage can rank as “worst J_A” despite better IL (5.2, 4.7).

Obj-grid correlations (higher score should track better engineering):

| Pair | Pearson | Spearman |
|------|--------:|---------:|
| J_A vs \(T_\mathrm{des}\) | 0.53 | — |
| J_A vs (−IL) | 0.29 | 0.37 |
| J_A vs isolation | 0.44 | 0.55 |

Weak coupling of J_A to insertion loss is the smoking gun.

---

## 3. Does absolute-power A hide tiny-desired / tiny-reverse pathology?

**Yes.**

Synthetic point \(T_\mathrm{des}=10^{-4}\), \(T_\mathrm{rev}=10^{-3}\), \(R=0.9\), TE1=0.05:

- \(J_A \approx -0.277\) (better than many stage-1 frequencies with usable \(T_\mathrm{des}\)).
- Same point under A2_hybrid: \(J \approx -11.8\) (among the worst conceivable scores).

So both numerator and denominator of isolation can approach zero while J_A remains O(−0.3). A ratio/isolation term alone would still be insufficient without a separate floor on useful transmission; A2 uses **both** \(\log(T_\mathrm{des}+\varepsilon)\) and \(\log((T_\mathrm{des}+\varepsilon)/(T_\mathrm{rev}+\varepsilon))\).

---

## Phase 2 — Candidate engineering-aligned objectives

All remain **mode-preserving** (TE0 desired, TE1 undesirable). C6 unchanged. Not port-power Objective B.

| Candidate | Formula (maximize) | Thresholds? |
|-----------|-------------------|-------------|
| **A_current** | linear A as frozen | none |
| **A2_log_ratio / A2_hybrid** | \(2\log(T_\mathrm{des}+\varepsilon)-\log(T_\mathrm{rev}+\varepsilon) - 1.5\,T_\mathrm{thr}-T_\mathrm{other}-0.25\,R-\mathrm{TE1}\) | none (\(\varepsilon=10^{-6}\)) |
| **A2_softplus_floor** | \(-\mathrm{softplus}(T_\mathrm{floor}-T_\mathrm{des})-\mathrm{softplus}(\log\rho_\mathrm{floor}-\log(T_\mathrm{des}/T_\mathrm{rev}))-\ldots\) | needs \(T_\mathrm{floor}\), \(\rho_\mathrm{floor}\) |

### On thresholds (not frozen here)

If softplus floors were used:

| Symbol | Example | Physical meaning |
|--------|---------|------------------|
| \(T_\mathrm{floor}=0.1\) | IL ≈ 10 dB | minimum useful TE0 through power |
| \(\rho_\mathrm{floor}=10\) | isolation ≈ 10 dB | minimum desired/reverse ratio |

These are **illustrative engineering targets**, not derived from a validation bar. Preferring A2_hybrid avoids freezing arbitrary dB cutoffs; \(\log(T_\mathrm{des})\) already diverges to \(-\infty\) as transmission dies.

---

## Phase 3 — Offline ranking on dense spectra (n=52)

| Metric | A_current | A2_hybrid | A2_softplus_floor |
|--------|----------:|----------:|------------------:|
| Rank of **4.1 GHz** (1=worst) | **34** | **1** | **1** |
| Rank of **5.625 GHz** | 1 | 11 | 6 |
| Spearman(score, −IL) | 0.27 | **0.89** | 0.70 |
| Spearman(score, isolation) | 0.60 | **0.92** | 0.99 |
| 4.1 in worst-5? | no | **yes** | **yes** |
| Pathology synth worse than “good”? | **no** | **yes** | **yes** |
| **Decision** | **REJECT** | **ACCEPT (recommended)** | ACCEPT (secondary) |

A2_hybrid worst-10 head: 4.1, 4.075, 5.85, 5.65, 5.725, 5.7, 5.35, 4.85, 5.15, 4.15 — dominated by high-IL / dead-transmission points. 5.625 remains badly ranked (11) because reverse is large and \(T_\mathrm{des}\) is only moderate — still correctly treated as a serious failure, but no longer allowed to overshadow a 40 dB IL corpse.

**Recommended production candidate:** `A2_hybrid` (`objective_A2.py`), pending Phase-4 FD validation.

---

## Implications for enrichment

Adaptive enrichment (add 4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95) remains **necessary** for between-grid soft-min holes, but is **not sufficient**: 4.1 GHz was already an objective frequency with catastrophic IL and was under-weighted. Stage-2 must combine **A2 + enrichment**.
