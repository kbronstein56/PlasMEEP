# Frozen modal-S Onsager acceptance criterion

**Frozen before Phase-1 experiments. Do not edit after seeing results.**

Date frozen: 2026-10-12 (campaign start).

## Context

- Reaction / FEM operator Onsager is already validated to ~roundoff.
- Uniform `s≡1` modal-S Onsager residual was ~1.8×10⁻⁴ (FULL91-C).
- Sixfold Pref-optimized design showed ~0.040 on FULL91-M — treated as
  an **extraction / discretization** issue until proven otherwise.
- Feed at 3.85 GHz: `k0 ≈ 1.614`, `clear_width = 2`. The n=1 transverse
  mode has `kc = π/W ≈ 1.571 < k0` and **propagates** (`β₁ ≈ 0.37`).
  One-mode S is still the production target, but TE1 contamination of the
  `(Hz,∂n)` projector must be quantified.

## Primary gate (PASS requires all)

Evaluated on the **sixfold peak** design at `f=3.85 GHz`, `B=±0.05 T`,
production extractor (flat fundamental, default monitor station), on the
finest mesh completed in the audit (at least FULL91-M; FULL91-F if run):

1. **Frobenius Onsager residual**
   ```
   r_F = ||S(+B) − Sᵀ(−B)||_F / ||S(+B)||_F  ≤  5×10⁻³
   ```

2. **Significant-channel elementwise residual**
   For every entry with `|S_ij(+B)| ≥ 0.10`:
   ```
   |S_ij(+B) − S_ji(−B)|  ≤  0.02
   ```
   and the same bound on `| |S_ij(+B)| − |S_ji(−B)| |` (magnitude).

3. **Monitor-station magnitude invariance** (same mesh as (1), peak design)
   Shift each monitor by `Δζ ∈ {−0.3, 0, +0.3}` along the local outward
   normal within the uniform feed (still between source and flare):
   ```
   max_{significant ij}  | |S_ij(Δζ)| − |S_ij(0)| |  ≤  0.015
   ```

4. **Mesh trend**
   `r_F` on FULL91-M must be **strictly smaller** than on FULL91-C for the
   same design/extractor, or both already ≤ 5×10⁻³.
   If FULL91-F is available, `r_F(F) ≤ r_F(M) + 1×10⁻³` (no large regression).

## Secondary diagnostics (inform, do not alone decide PASS)

- TE1 power fraction on monitors (should be reported; large TE1 is allowed
  only if fundamental S still meets (1)–(3)).
- Sample-count doubling (`n_points` 41→81) must change significant `|S|`
  by ≤ 0.005.
- Driven-port `a_j` phase/magnitude consistency across the six ports under C6.

## Fail / UNRESOLVED

- If (1)–(4) are not met after reasonable extractor fixes (monitor shift,
  sample count, optional TE1-orthogonalized projection), declare
  `MODAL_S_OPTIMIZED_ONSAGER = FAIL` or `UNRESOLVED` and **do not** start
  true modal-S optimization.
- Do **not** raise the 5×10⁻³ / 0.02 thresholds after seeing numbers.
