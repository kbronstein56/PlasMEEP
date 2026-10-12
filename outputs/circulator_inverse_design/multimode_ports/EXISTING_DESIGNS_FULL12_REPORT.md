# Existing designs under FULL12 DtN S (Phase 7)

Grade C, f=3.85 GHz, B=+0.05 T. Extract: inset `(Hz,∂n)` 0.08a inside DtN cut; `S=B A⁻¹`.

## FULL12 Onsager (frozen bar 5×10⁻³)

| Design | r_F (+B/−Bᵀ) | Status |
|--------|--------------|--------|
| uniform | 0.157 | **FAIL** |
| Pref sixfold peak | 0.183 | **FAIL** |
| Pref broadband | 0.186 | **FAIL** |

Prior “PASS” at r_F=0.004 was retracted (S≈I from ∂n≈0 on the cut).

## TE0 figures (DtN 12×12, mean over 6 sources)

| Design | T_des | IL (dB) | isolation (dB) | RL (dB) | R_TE0 | TE1 total | rem |
|--------|-------|---------|----------------|---------|-------|-----------|-----|
| uniform | 0.011 | 19.6 | 7.0 | ~1.3 | 0.749 | 0.084 | −0.014 |
| Pref sixfold peak | **0.548** | **2.62** | **12.8** | ~11 | 0.078 | **0.319** | 0.011 |
| Pref broadband | **0.548** | **2.62** | **13.7** | ~12 | 0.070 | **0.316** | 0.018 |

`rem = 1 − Σ|S|²` over 12 channels (material + non-port lump). TE1 is **not** absorption.

## Soft-source comparison (fundamental 6×6)

| Design | soft T_des | soft r_F | soft TE1/P_inc |
|--------|------------|----------|----------------|
| uniform | 0.0074 | 1.8×10⁻⁴ | small |
| Pref peak | 0.739 | 0.039 FAIL | **0.106** |

DtN Pref T_des (0.55) is in the same ballpark as soft (0.74) but lower; DtN TE1 (0.32) is higher than soft TE1 (0.11) — extractor still not frozen for optimization.

## Gate

**FULL12_ONSAGER: FAIL** — do not open multimode optimization.
**Design FoMs:** informative only until Onsager passes.
