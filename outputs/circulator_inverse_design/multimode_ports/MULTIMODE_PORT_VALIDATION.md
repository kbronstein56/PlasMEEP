# Multimode port validation (Phase 3)

## Frozen straight-guide setup

| Parameter | Value |
|-----------|-------|
| L | 24 (code units) |
| h | 0.025 |
| n_modes Robin | 2 propagating |
| n_evanescent | 2 |
| f | 3.85 GHz |

## Frozen thresholds

| Check | Threshold |
|-------|-----------|
| Transmission \|S\| | ≥ 0.995 |
| Reflection \|S\| | ≤ 0.015 |
| Cross-mode \|S\| | ≤ 0.02 |

## Results (`straight_dtn_validation.json`)

| Test | Status |
|------|--------|
| TE0 incidence | PASS (T≈1.000, R≈0.004, cross-mode ≲0.012) |
| TE1 incidence | PASS (T≈0.999, R≈0.014, cross-mode ≲0.013) |
| Mode orthogonality proxy | PASS |
| Lossless partial | PASS |

**MULTIMODE_PORT_VALIDATED (straight empty guide): PASS**

Full91 DtN placement and 12×12 Onsager are tracked separately in `FULL12_ONSAGER_VALIDATION.md`.
