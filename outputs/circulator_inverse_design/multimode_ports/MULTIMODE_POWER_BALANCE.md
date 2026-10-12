# Multimode power balance (Phase 6)

## Soft-source Pref sixfold peak (grade C, +0.05 T)

From `modal_s_core` TE1 report on soft-source solves:

| Term | Fraction of P_inc |
|------|-------------------|
| Fundamental modal outgoing (Σ\|b₀\|²) | ≈ 0.74 (des) + R + leak … |
| **TE1 outgoing (all ports)** | **≈ 0.106** |
| Prior “PML/higher-mode” bucket | ≈ 0.13–0.15 |

**Conclusion:** most of the previously reported ~13% PML/higher-mode term on the Pref peak is **outgoing TE1**, not material loss and not pure PML. Do not call TE1 absorption.

## DtN power identity

With inset extract, smoke test on uniform: `P_out/P_inc ≈ 0.99` (propagating TE0+TE1). Remaining ≤ few % is material + discretization.

Frozen engineering bar (when extract trusted): `|1 − Σ|S|² − P_mat/P_inc − P_nonport/P_inc| ≤ 0.05`.

**MULTIMODE_POWER_BALANCE: PARTIAL** (soft-source TE1 split quantified; full DtN+material ledger pending trusted FULL12).
