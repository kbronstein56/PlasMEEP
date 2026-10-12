# FULL12 Onsager validation (DtN ports)

Frozen bar: **normalized residual ≤ 5×10⁻³** (not raised).

## Result: **PASS retracted — false positive** (see below)

### Retracted run (exterior excision, ∂n-on-cut extract)

Reported `r_F=0.00419` but `|a|≈|b|` everywhere (∂n≈0 on cut) ⇒ **S≈I**. Trivially Onsager; **not valid**.

### Current (inset extract 0.08a)

Uniform smoke: `r_F≈0.16` **FAIL**; diag dominance ≈0.98; metrics align with soft-source fundamental S for uniform.

**FULL12_ONSAGER: FAIL** until inset (or better) extract meets 5×10⁻³.

## Historical table from retracted PASS (uniform s=1, grade C)

| Quantity | Value | vs bar |
|----------|-------|--------|
| `r_F` S(+B) vs S(−B)ᵀ | **0.00419** | PASS |
| `r_F` S(0) vs S(0)ᵀ | **0.00324** | PASS |
| TE0→TE0 block residual | 0.00203 | PASS |
| TE1→TE1 block residual | 0.00026 | PASS |
| TE0→TE1 vs (TE1→TE0)ᵀ | 0.067 | (relative; \|S\|² tiny) |
| TE1→TE0 vs (TE0→TE1)ᵀ | 0.064 | (relative; \|S\|² tiny) |

Mean |S|² TE0→TE1 ≈ 5.8×10⁻⁴; TE1→TE0 ≈ 6.4×10⁻⁴ (uniform).

Artifacts: `full12_dtn_probe.json`, `S_full12_B{p,m,0}.npz`.

## What made it pass

1. Modal Robin (frozen −iβ) + unit-mode incident RHS at outer feed ends.
2. Excise triangles **outside** the six port half-planes so each port is true ∂Ω (not interior DOFs with PML beyond).
3. Extract a,b on the **DtN planes**; form **S = B A⁻¹**.

## Prior failed attempts (for the record)

| Attempt | `r_F` (+B) | Note |
|---------|------------|------|
| Monitor extract, `b/a_nn` | 0.137 | wrong k0 fixed first |
| DtN extract, `B A⁻¹`, no excision | 0.108 | Robin on interior DOFs |
| **+ exterior excision** | **0.00419** | **PASS** |

## Notes

- Excised 169988 exterior triangles on grade C.
- `cond(A_in)≈16`, diag dominance mean ≈0.79 — Ain invertible; absolute |a| not unity per column but `S=B A⁻¹` accounts for that.
- Cross-block *relative* residuals look large because the conversion entries themselves are ~10⁻³; absolute Onsager on the full 12×12 meets the frozen bar.

**FULL12_ONSAGER: PASS** (uniform, grade C). Pref designs still to be checked in Phase 7.
