# Modal DtN / Robin port weak form

Time convention: **e^{-iωt}**.

## Independent IBP (validated volume operator)

From \(-\nabla\cdot(\rho\nabla H)-k_0^2 H=0\) (isotropic sketch \(\rho=1\)):

\[
\int_\Omega \nabla H\cdot\nabla v - k_0^2 H v
\;-\;
\int_{\partial\Omega}(\partial_n H)\,v
\;=\;0.
\]

Natural Neumann omits the boundary integral. Robin \(\partial_n H=\lambda H\)
contributes **\(-\lambda\)** to the Galerkin matrix.

## Outgoing first-order ABC

Outgoing w.r.t. outward \(\hat n\): \(\partial_n H = i\beta H\) (\(\lambda=i\beta\)).

Matrix contribution:

\[
\boxed{\text{add }(-i\beta)\text{ on the port (modal rank-one }c c^\top\text{)}}
\]

This is the **frozen** implementation sign (straight-guide validated). Comments
that once wrote \(+i\beta\) as the matrix term contradicted IBP; the code sign
is correct.

## Incident forcing

\(\partial_n H = i\beta H - 2i\beta a_f\phi\) ⇒ with the same IBP, the known term
appears as RHS \(+\int 2i\beta a_f\phi\,v\) when the Robin \((-i\beta)\) is on \(A\)
**or** equivalently RHS \(-\int 2i\beta\ldots\) if the overall residual sign is
flipped. The frozen FEM pair that matches empty-guide transmission is:

- Robin: \((-i\beta)\,c c^\top\)
- RHS: \(+2i\beta\, a_f\, c\)

## Discrete direct \(S\)

With true P1 boundary edges and that frozen pair:

\[
S = I + 2i\,Q,
\qquad
Q_{mn}=\sqrt{\beta_m\beta_n}\,c_m^\top A^{-1} c_n.
\]

See `DISCRETE_DTN_SCATTERING_DERIVATION.md`.

## Geometry

Assemble \(c\) on **exposed boundary edges** after exterior excision — not on a
`plane_tol` node cloud.
