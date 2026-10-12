# FULL12 Onsager validation — **direct DtN S**

Frozen bar: **\(r_F\le 5\times10^{-3}\)** (unchanged).

## Result: **PASS** (grade C)

| Design | \(r_F\) (+B/−Bᵀ) | worst \|Δ\| | \(r_F\) (B=0) | T_des | TE1 | rem |
|--------|------------------|-------------|---------------|-------|-----|-----|
| uniform | **3.9×10⁻¹⁴** | 3.0×10⁻¹⁴ | 2.2×10⁻¹⁴ | 0.0085 | 0.061 | 0.007 |
| Pref sixfold peak | **8.7×10⁻¹⁴** | 8.9×10⁻¹⁴ | 9.2×10⁻¹⁵ | 0.556 | 0.317 | 0.016 |
| Pref broadband | **6.2×10⁻¹⁴** | 5.9×10⁻¹⁴ | 6.5×10⁻¹⁵ | 0.550 | 0.315 | 0.021 |

Method: true boundary-edge \(c_m\), \(S=I+2iQ\), \(Q_{mn}=\sqrt{\beta_m\beta_n}\,c_m^\top A^{-1}c_n\).

Artifact: `direct_dtn_full12_report.json`.

## Prior invalid results (do not use)

| Method | \(r_F\) | Note |
|--------|---------|------|
| Node-cloud + cut \(\partial_n\) | ~0 / false PASS | \(S\approx I\) |
| Node-cloud + inset Ain/Bout | 0.16–0.19 FAIL | See `OLD_RF_FAILURE_EXPLANATION.md` |

**DIRECT_DTN_FULL12_ONSAGER: PASS**
