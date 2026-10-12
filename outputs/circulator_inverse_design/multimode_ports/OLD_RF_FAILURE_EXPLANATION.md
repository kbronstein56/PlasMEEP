# Why the old FULL12 \(r_F\approx0.16\)–0.19 failed

## Not a bulk Onsager failure

Phase 1 audit (`DTN_MATRIX_ONSAGER_AUDIT.md`):

\[
\|A(+B)-A(-B)^\top\|/\|A(+B)\| \sim 10^{-17}
\]

including gyrotropic bulk, exterior excision, and **all** modal Robin terms.
Port blocks \(c c^\top\) are exactly complex-symmetric. The FEM operator was fine.

## Cause 1 — node-cloud “boundary”

`mesh_port_nodes(plane_tol=0.06)` collected **~110 distinct along-layers**
(~112 nodes) per port — a slab of interior+boundary nodes, not \(\partial\Omega\).

Modal Robin used 1-D trapezoid weights along \(y\) on that cloud. That is not a
consistent P1 edge integral on the cut. The discrete channel coupling \(Q\) used
for \(S\) was therefore not the same reciprocal boundary functional that entered
\(A\).

## Cause 2 — inset \((H_z,\partial_n)\) → Ain/Bout → \(S=B A^{-1}\)

Primary \(S\) was rebuilt by:

1. solving with DtN RHS,
2. moving **0.08a inward**,
3. estimating \(\partial_n\) from one-sided / inset elements,
4. reconstructing \(a,b\), forming Ain/Bout,
5. \(S=\mathrm{Bout}\,\mathrm{Ain}^{-1}\).

On the cut, \(\partial_n\) was previously ~0 (forcing \(a\approx b\), \(S\approx I\),
false Onsager PASS). After inset, \(\partial_n\) worked enough for FoMs to look
plausible, but:

- reconstructed \(a\) ≠ prescribed DtN incident,
- Ain mixing / conditioning polluted reciprocity,
- observation functional ≠ excitation functional \(c\) in \(A\).

Hence \(S(+B)\neq S(-B)^\top\) at \(r_F\sim0.16\) even though \(A\) obeyed Onsager.

## Fix

1. True exposed **boundary edges** after excision; P1 edge loads \(c_m\).
2. Direct \(S=I+2iQ\) with \(Q_{mn}=\sqrt{\beta_m\beta_n}\,c_m^\top A^{-1}c_n\)
   (same \(c\) as Robin/RHS).

This forces \(A(+B)=A(-B)^\top\Rightarrow S(+B)=S(-B)^\top\) up to solve error.
