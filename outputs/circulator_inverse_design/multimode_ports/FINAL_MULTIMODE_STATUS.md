# Final multimode / DtN status

## Gates

| Gate | Status |
|------|--------|
| DTN_MATRIX_ONSAGER | **PASS** (\(\sim10^{-17}\)) |
| TRUE_BOUNDARY_PORT_ASSEMBLY | **PASS** (y_span=W; staircase edges; empty-guide TE0/TE1 PASS) |
| DIRECT_DTN_FULL12_ONSAGER | **PASS** (\(r_F\sim10^{-14}\) uniform + Pref peak + broadband) |
| MULTIMODE_POWER_BALANCE | **FAIL** (TE0 PASS; 1/6 TE1 channels at 6.7% > 5%) |
| DIRECT_DTN_ADJOINT_FD | **NOT_REACHED** |
| MULTIMODE_OPT | **NOT_REACHED** |

## Cause of prior \(r_F\approx0.16\)–0.19

See `OLD_RF_FAILURE_EXPLANATION.md`: node-cloud Robin + inset Ain/Bout extract, **not** \(A(+B)\neq A(-B)^\top\).

## Pref peak FoMs (direct DtN, TE0 mean)

- T_des ≈ 0.556, IL ≈ 2.6 dB, TE1 ≈ 0.317, material ≈ 0.037

## Commits

`32da092` and follow-ups on `agent/eigenmode-ports`.
