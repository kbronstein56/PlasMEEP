# Multimode power balance (direct DtN)

## Pref sixfold peak (grade C, +0.05 T)

| Drive | Σ\|S\|² | TE0 out | TE1 out | material | rem (1−out−mat) |
|-------|--------|---------|---------|----------|-----------------|
| TE0 (mean) | 0.984 | 0.667 | **0.317** | 0.037 | −0.021 |
| TE1 (mean) | 1.000 | 0.318 | 0.682 | 0.031 | −0.031 |

Frozen bar: engineering residual ≤ 0.05.

- All **6 TE0** drive channels: **PASS**
- **5/6 TE1** channels: PASS; channel 7 (port 3, TE1): residual **0.067** (column Σ\|S\|²=1.036)

Empty-guide direct S is unitary to \(10^{-10}\). Full91 staircase ports show ~few‑percent column-norm overshoot on some TE1 drives — numerical, not “PML hiding TE1”.

**TE1 is explicit in S** (≈32% of TE0-incident power on Pref peak), matching the soft-source TE1≈11–32% story; it is **not** lumped into absorption.

## Soft-source Pref TE1 cross-check

Earlier soft-source extract: TE1 out / P_inc ≈ **0.106**. Direct DtN Pref TE1 ≈ **0.317**. Different port reference / normalization; both show TE1 is a major open channel.

**MULTIMODE_POWER_BALANCE: FAIL** (strict: one TE1 channel over bar). TE0-driven balance **PASS**. Adjoint/opt gated off until TE1 column norms are tightened.
