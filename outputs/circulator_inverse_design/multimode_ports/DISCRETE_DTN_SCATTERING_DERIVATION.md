# Discrete DtN scattering derivation

Time: **e^{-iωt}**. Validated volume weak form (after IBP):

\[
\int_\Omega \rho\nabla H\cdot\nabla v - k_0^2 H v
\;-\;
\int_{\partial\Omega}(\partial_n H)\,v
\;=\;0
\]

Natural Neumann \(\partial_n H=0\) omits the boundary integral. For Robin
\(\partial_n H=\lambda H\) the boundary term contributes **\(-\lambda\)** to the
Galerkin matrix.

## Outgoing first-order ABC

Outgoing (away from the device) with outward \(\hat n\):

\[
\partial_n H = i\beta\, H
\qquad\Rightarrow\qquad
\lambda=i\beta
\qquad\Rightarrow\qquad
\text{add }(-i\beta)\text{ to }A.
\]

This is the **frozen** `dtn_assembly` / straight-guide sign. (Comments that
wrote \(+i\beta\) as the matrix contribution were inconsistent with IBP; the
code sign is the one required by this derivation.)

## Incident mode

With \(H=(a_f+b_f)\phi\) and \(\partial_n H=i\beta(b_f-a_f)\phi\):

\[
\partial_n H = i\beta H - 2i\beta a_f\phi.
\]

Weak form ⇒ Robin \((-i\beta)\) on \(A\) and RHS

\[
b_p += \int_\Gamma 2i\beta\, a_f\,\phi\, N_p\,ds.
\]

## Discrete channels

P1 boundary loads (true edges):

\[
(c_m)_p = \int_\Gamma \phi_m N_p\,ds,
\qquad
\int\phi_m\phi_n=\delta_{mn}.
\]

Drive channel \(n\) with field amplitude \(a_{f,n}\):

\[
A x = (2i\beta_n a_{f,n})\, c_n,
\qquad
h_m := c_m^\top x = a_{f,m}+b_{f,m}.
\]

Power normalization \(\hat a = a_f\sqrt{\beta/(2\omega)}\). For unit-power drives
\(a_{f,n}=1/\mathrm{scale}_n\), \(\hat a_n=1\):

Frozen code pair: Robin \((-i\beta)\) and incident RHS \((+2i\beta a_f c)\).
On a matched empty-guide port this pair yields

\[
2i\beta\,(c^\top A^{-1} c) = -1,
\]

so the power-wave map consistent with \(R=0\), \(|T|=1\) is

\[
S_{mn}
=
\delta_{mn}
+
2i\sqrt{\beta_m\beta_n}\,(c_m^\top A^{-1} c_n),
\]

i.e.

\[
S = I + 2i\, Q,
\qquad
Q_{mn}=\sqrt{\beta_m\beta_n}\, c_m^\top A^{-1} c_n.
\]

(The continuum bookkeeping with the opposite RHS sign produces \(S=-I+2iQ\);
we keep the FEM pair that already passes straight-guide field tests and use
\(S=I+2iQ\).)

## Automatic Onsager

Port vectors \(c_m\) are **real** and **B-independent**. Modal Robin blocks are
complex-symmetric. Validated bulk satisfies \(A(+B)=A(-B)^\top\). Therefore

\[
Q(+B)=Q(-B)^\top
\quad\Rightarrow\quad
S(+B)=S(-B)^\top
\]

up to linear-solve / roundoff — **without** reconstructing \(a\) from fields.

## What not to do

Do **not** define primary \(S\) by inset \((H_z,\partial_n)\) → Ain/Bout →
`Bout @ inv(Ain)`. That path caused the \(r_F\approx0.16\) failure despite
matrix Onsager holding at \(10^{-17}\).
