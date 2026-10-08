# Fan 2019 forward-mode comparison

Reference: Hughes, Williamson, Minkov, Fan,
“Forward-Mode Differentiation of Maxwell’s Equations,”
ACS Photonics 6, 3010–3016 (2019), DOI 10.1021/acsphotonics.9b01238,
arXiv:1908.10507.

This note is equation-level. The paper does not publish digit tables that
can be reproduced as a numerical pass row.

## What the paper does

The paper differentiates an FDTD update. For a design variable `p` and a
time-domain state, one additional forward-mode field is marched with the
same linearized update. Cost scales with the number of design variables,
not with the number of outputs. That is attractive when many field samples
or spectral channels are differentiated with respect to a few parameters.

The two demonstrations are:

1. near-field intensity of a dielectric scatterer versus permittivity
2. grating-coupler spectral power / efficiency versus fill factor

Neither demonstration supplies a mesh size, fill-factor table, or efficiency
digit that can be matched here.

## What our discrete FEM does

Frequency-domain P1 Helmholtz for `Hz`:

    A(p) x(p) = b(p)

With a fixed source and material-only `p = s` (density scale),

    A dx/ds = − (dA/ds) x
    db/ds = 0

`dε/ds` and `dρ/ds = −ρ (dε) ρ` are analytic. `dA/ds` is the stiffness of
`dρ` only. One factorization of `A` is reused for the primal solve and for
every sensitivity RHS at that `(ω, B, mesh)`.

## Method-independent agreement

| Item | Fan 2019 (FDTD) | This FEM |
|------|-----------------|----------|
| Stationarity | linearized Maxwell update | `A dx = −(dA) x` |
| Solves per design variable | one extra forward-mode field | one extra solve with the same `A` |
| Scaling with N_outputs | one forward-mode field supplies all outputs | one `dx/dp` supplies all linear outputs of `x` |
| Scaling with N_parameters | N_param forward-mode runs | N_param sensitivity RHS (same LU) |
| Exact vs finite difference | exact discrete derivative | exact discrete derivative of this FEM |

## What differs

- Discretization: Yee FDTD versus P1 FEM on triangles.
- Time convention and constitutive model: their examples are reciprocal
  dielectrics; ours is a cold-plasma Lorentz tensor with optional `B`.
- Memory: FDTD stores time history or checkpointed fields; our direct
  SuperLU stores one complex LU per `(mesh, f, B)`.
- Suitability for ~91 plasma controls and six ports: forward mode needs
  about 91 sensitivity solves after one factorization. An adjoint needs
  one extra solve per real scalar objective. With many ports and one
  scalar figure of merit, adjoint wins; with many field samples per rod
  density, forward mode can still be competitive after the LU exists.

## Claim

The discrete sensitivity equation for our linear FEM is the
frequency-domain counterpart of the paper’s forward-mode idea. It is not
a numerical reproduction of their FDTD examples.
