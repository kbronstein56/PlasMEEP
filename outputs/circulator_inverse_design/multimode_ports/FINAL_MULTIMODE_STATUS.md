# FINAL_MULTIMODE_STATUS

| Gate | Status |
|------|--------|
| TE1_BOUNDARY_DIAGNOSTIC | PASS |
| DIRECT_DTN_FULL12_ONSAGER | PASS (~1e-14, projected) |
| MULTIMODE_POWER_BALANCE | PASS (max 0.0377 ≤ 5%) |
| DIRECT_DTN_ADJOINT_FD | PASS (fwd≡adj ~1e-13; FD ≤1e-3) |
| MULTIMODE_OPT | NOT_REACHED |

**Original failing channel:** ch7 = port 3 TE1, residual 0.0673 under
arclength measure. **Cause:** staircase `ds` vs flat-port `dy` TE1
normalization / power extraction. **Fix:** frozen `measure=projected`
(`|e·t̂|`). Onsager preserved at numerical precision. Optimization not started.
