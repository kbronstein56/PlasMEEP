# Stage-2 Broadband A2 Optimization — Results

**Status:** COMPLETE (2026-11-01)  
**Spec:** `BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN.md`  
**Runner:** `production_broadband_opt_stage2.py`  
**Design:** `opt_broadband_stage2_A2_C.json` (`q_best`)  
**Do not auto-start Stage-3.**

---

## Run summary

| Item | Value |
|------|------:|
| Objective | A2_hybrid |
| τ / α | 0.05 / 0.25 |
| Frequencies | 28 (enriched) |
| Symmetry / mesh | C6 / FULL91-C |
| Warm start | Stage-1 `opt_broadband_production_C.json` |
| Accepted L-BFGS iters | **30 / 30** |
| Wall time | **~78.1 h** (281238 s) |
| Stop | iteration limit |
| Early flags (iter 5) | **none** (`bb_stage2_early_flags.json`) |

---

## Objective trajectory

| Milestone | \(J_{BB}\) | worst A2 | worst f | \(T_\mathrm{des}\)@4.1 | IL@4.1 | iso@4.1 | max \(w\) | \(N_{\mathrm{eff}}\) |
|-----------|----------:|---------:|--------:|-----------------------:|-------:|--------:|----------:|---------------------:|
| Start | −11.331 | −13.881 | **4.1** | 3.36×10⁻⁴ | 39.94 | −18.80 | 0.759 | 1.73 |
| Accepted 1 | −9.520 | −11.470 | 5.8 | 0.131 | 8.83 | 7.13 | 0.759 | 1.73 |
| Accepted 5 | −7.044 | −8.337 | 5.4 | 0.136 | 8.68 | 7.16 | 0.734 | 1.84 |
| Best / final | **−6.229** | **−7.397** | **5.65** | **0.141** | **8.51** | **7.07** | 0.407 | 3.75 |

Early diagnostic (first 5 accepted): 4.1 GHz \(T_\mathrm{des}\) improved **~400×** (3.4×10⁻⁴ → 0.136); A2@4.1 −13.88 → −0.85. No stall / mismatch / oscillation flags.

---

## Expected minimax behavior — observed

1. **Initial dominance:** start soft-min one-hot on 4.1 GHz (max \(w\approx0.76\)).
2. **4.1 recovery:** by accepted iter 1, worst frequency already left 4.1; engineering at 4.1 became usable (IL~9 dB, iso~7 dB).
3. **Migration:** worst-frequency pressure moved among 5.8 → 5.7 → 5.4 → 4.3 → **5.65** (enrichment frequency). This is healthy minimax, not a failure.
4. **Late multi-frequency:** final \(N_{\mathrm{eff}}\approx3.75\), max \(w\approx0.41\) — soft-min no longer single-frequency once the catastrophic outlier closed.

---

## Best-design engineering (objective grid)

### At 4.1 GHz (former corpse)

| Metric | Stage-1 | Stage-2 best |
|--------|--------:|-------------:|
| \(T_\mathrm{des}\) | 3.36×10⁻⁴ | **0.141** |
| IL (dB) | 39.94 | **8.51** |
| isolation (dB) | −18.80 | **+7.07** |
| A2 | ≈ −11.6 (offline avg) / −13.9 (true) | **−0.84** |

### Current worst (5.65 GHz)

| Metric | Value |
|--------|------:|
| A2 | −7.40 |
| \(T_\mathrm{des}\) | 0.018 |
| IL | 18.5 dB |
| isolation | −11.9 dB |
| TE1 | 0.132 |

Mean A2 on 28-grid (final): **−3.13**.

---

## Artifacts

| File | Role |
|------|------|
| `opt_broadband_stage2_A2_C.json` | final / best design |
| `bb_stage2_checkpoint.json` | restartable checkpoint |
| `bb_stage2_history.json` | full history |
| `bb_stage2_early_diagnostic.json` | first 5 accepted |
| `bb_stage2_early_flags.json` | empty flags |
| `BB_STAGE1_VS_STAGE2.md` | same-grid comparison |
| `bb_stage2_audit_50MHz.json` | dense audit (running / post) |

---

## Next (not auto-started)

1. Dense **50 MHz** audit of Stage-2 design (in progress / complete — see audit artifacts).
2. Targeted **25 MHz** only if between-grid hole found.
3. Decide enrichment / Stage-3 only after documenting scientific need.
4. Continuum F/Q/X2 and Meep remain later independent checks.
