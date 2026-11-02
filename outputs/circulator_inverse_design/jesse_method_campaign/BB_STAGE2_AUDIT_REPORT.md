# Stage-2 dense audit report

**Design:** `opt_broadband_stage2_A2_C.json`  
**50 MHz:** `bb_stage2_audit_50MHz.json` (41 pts, 4.0–6.0 GHz) — COMPLETE  
**25 MHz local:** `bb_stage2_audit_25MHz_local.json` (5.50–5.90 GHz) — COMPLETE  
**Between-grid:** `bb_stage2_between_grid_report.json`  
**No Stage-3 launched.**

---

## Full-band engineering (50 MHz)

| Metric | mean | worst | f (GHz) | on Stage-2 obj grid? |
|--------|-----:|------:|--------:|:--------------------:|
| \(T_\mathrm{des}\) | 0.102 | 0.00183 | **4.3** | yes |
| IL (dB) | 13.01 | **27.47** | **4.3** | yes |
| isolation (dB) | 2.58 | **−20.14** | **5.75** | **no** |
| RL (dB) | 6.47 | — | — | — |
| TE1 | 0.286 | — | — | — |
| A2 | −2.85 | **−10.98** | **5.75** | **no** |

On-grid A2 worst: **5.65 GHz** (A2=−7.40) — matches optimizer soft-min focus.

---

## 25 MHz refinement (5.50–5.90)

| Metric | worst | f (GHz) |
|--------|------:|--------:|
| A2 | **−10.98** | **5.75** |
| IL | **25.97** | **5.75** |
| isolation | **−21.28** | **5.825** |
| \(T_\mathrm{des}\) | 0.0027 | 5.75 |

Top-5 worst A2 in the local window: **5.75, 5.825, 5.725, 5.775, 5.875**.

50 vs 25 (band-restricted):  
`d_worst_iso≈1.14 dB`, location shift 5.75→5.825 (0.075 GHz) → **iso bar FAIL**; IL location unchanged.

---

## Enrichment decision

| Question | Answer |
|----------|--------|
| Enrichment needed? | **YES** |
| Primary add | **5.75 GHz** (worst A2 + worst IL in window; off-grid vs Stage-2 28-set) |
| Secondary add | **5.825 GHz** (worst isolation on 25 MHz; adjacent spike) |
| On-grid residual | **4.3 GHz** IL≈27.5 dB (already in objective — FoM/pressure issue, not sampling alone) |

If Stage-3 approved → objective grid = Stage-2 28 + `{5.75, 5.825}` → **30 frequencies**.

---

## Proposed Stage-3 package (NOT started)

| Item | Value |
|------|-------|
| Why | Close between-grid A2/iso hole ~5.75–5.825 that Stage-2 soft-min never saw |
| Scientific question | Can A2 soft-min on 30 freqs suppress the 5.75/5.825 failure without undoing 4.1 GHz recovery? |
| Objective | A2_hybrid (unchanged) |
| τ / α | 0.05 / 0.25 (unchanged) |
| C6 / mesh / bounds | unchanged |
| Start | Stage-2 `q_best` |
| Add | **5.75, 5.825** |
| maxiter | 30 (suggested) |
| Expected wall | ~Stage-2 × 30/28 ≈ **80–85 h** |

**Do not auto-launch.** Human approval required.

Continuum F/Q/X2 and Meep remain after any Stage-3 decision.
