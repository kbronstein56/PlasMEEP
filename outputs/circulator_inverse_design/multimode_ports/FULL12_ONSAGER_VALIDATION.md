# FULL12 Onsager validation (DtN ports)

Frozen bar: **normalized residual ≤ 5×10⁻³** (not raised).

## Run (grade C, uniform s=1, f=3.85 GHz)

| Quantity | Value | vs bar |
|----------|-------|--------|
| `r_F` S(+B) vs S(−B)ᵀ | **0.137** | FAIL |
| `r_F` S(0) vs S(0)ᵀ | **0.230** | FAIL |
| TE0→TE0 block residual | **0.00310** | PASS |
| TE1→TE1 block residual | **0.00826** | FAIL (close) |
| TE0→TE1 vs (TE1→TE0)ᵀ | **0.528** | FAIL |
| TE1→TE0 vs (TE0→TE1)ᵀ | **0.932** | FAIL |

Mean |S|² TE0→TE1 ≈ 0.0096; TE1→TE0 ≈ 0.0035.

Artifact: `full12_dtn_probe.json`, `S_full12_B{p,m,0}.npz`.

## Interpretation

1. **Same-mode TE0 block already satisfies Onsager** under outer-feed DtN — consistent with reaction Onsager and the hypothesis that the old 6×6 failure was an open multimode subblock.
2. **B=0 reciprocity also fails** (`r≈0.23`), so the current TE1 drive/extract path is not yet validated (not a gyrotropic-only issue).
3. Likely causes under investigation:
   - Extraction at **legacy monitor** planes while DtN lives at the **outer feed end**
   - Column-wise `b / a_nn` instead of `S = B A⁻¹`
   - TE1 phase / power normalization between DtN and monitor

**FULL12_ONSAGER: FAIL** (extractor not frozen; diagnose before optimization).
