# Existing designs under FULL12 DtN S (Phase 7)

## Critical caveat on the first Onsager PASS

The first exterior-excision run reported `r_F=0.00419` **PASS**, but raw `|a|≈|b|` on every channel (∂n≈0 on the cut). That forces **S ≈ I**, which is trivially Onsager-symmetric and **not** a validated multimode S.

Inset `(Hz, ∂n)` extract (0.08a into the feed) restores diagonal dominance (~0.98) and `P_out/P_inc≈0.99`, and matches soft-source **uniform** modal metrics (R≈0.75–0.78, T_des≈0.007–0.01).

**FULL12_ONSAGER status with inset extract: FAIL / under re-measurement** (uniform smoke `r_F≈0.16`). Do not optimize.

## Soft-source reference (grade C, fundamental 6×6, for context)

| Design | mean T_des | mean R | r_F Onsager | mean TE1 out/P_inc |
|--------|------------|--------|-------------|--------------------|
| uniform | 0.0074 | 0.783 | 1.8×10⁻⁴ | (small) |
| Pref sixfold peak | **0.739** | 0.065 | **0.039** FAIL | **0.106** |
| Pref broadband | (see JSON) | | | |

≈10.6% of Pref-peak incident power leaves as TE1 — most of the old ~13% “PML/higher-mode” bucket.

## DtN re-eval (first report JSON — **invalid S≈I**)

| Design | claimed r_F | T_des (bogus) | note |
|--------|-------------|---------------|------|
| uniform | 0.00419 PASS | ~0 | S≈I artifact |
| Pref peak | 0.0138 FAIL | 0.011 | invalid extract |
| Pref broadband | 0.0097 FAIL | 0.006 | invalid extract |

Re-run with inset extract → `eval_designs_full12.py` / updated JSON.
