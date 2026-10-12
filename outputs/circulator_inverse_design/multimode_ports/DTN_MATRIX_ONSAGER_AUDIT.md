# DtN matrix Onsager audit (Phase 1)

Grade **C**, uniform s=1, f=3.85 GHz. Artifact: `dtn_matrix_onsager_audit.json`.

## Complete operator (bulk + exterior excision + all modal Robin)

| Test | Residual |
|------|----------|
| \(\|A(+B)-A(-B)^T\|/\|A(+B)\|\) | **1.47×10⁻¹⁷** |
| \(\|A(0)-A(0)^T\|/\|A(0)\|\) | **1.35×10⁻¹⁷** |
| Bulk alone (+B/−B) | 1.47×10⁻¹⁷ |
| Port TE0 / TE1 / evanescent \(\|A-A^T\|/\|A\|\) | **0** |

**DTN_MATRIX_ONSAGER: PASS**

The gyrotropic bulk + modal Robin (rank-one \(c c^T\) with real \(c\)) preserve complex transpose Onsager at assembly level. The previous full12 \(r_F\approx0.16\)–0.19 is **not** a broken \(A(+B)\neq A(-B)^T\).

## Port geometry finding (root cause for S failure)

`mesh_port_nodes(..., plane_tol=0.06)` on each port:

| Port | n_nodes | distinct along-layers | y_span | Σw |
|------|---------|----------------------|--------|-----|
| 0–5 | 112–116 | **106–111** | ≈2.0 (=W) | ≈2.0 |

Nodes span `along ∈ [-0.06,+0.06]` with **~110 distinct layers** — a **2-D node cloud**, not a 1-D boundary.

Modal Robin was assembled with 1-D trapezoid weights along \(y\) on that cloud. That is **not** a consistent P1 boundary-edge integral on \(\partial\Omega\). Combined with inset \((H_z,\partial_n)\) re-extraction, this yields a non-reciprocal discrete \(S\) even though \(A\) itself is Onsager.

## Next

Replace node-cloud quadrature with **true exposed boundary edges** after the port cut (Phase 2), then define **direct DtN \(S\)** without inset Ain/Bout (Phase 3).
