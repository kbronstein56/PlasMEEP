# Fan 2019 forward-mode comparison

Reference: Hughes, Williamson, Minkov, Fan,
“Forward-Mode Differentiation of Maxwell’s Equations,”
ACS Photonics 6, 3010–3016 (2019), DOI 10.1021/acsphotonics.9b01238,
arXiv:1908.10507.

Equation-level comparison. The paper does not publish digit tables that can
be reproduced as a numerical pass row. Our spatial discretization is P1 FEM,
not Yee FDTD; the claim is mathematical equivalence of the *forward-mode
differentiation framework*, not identical numerics.

## Discretized Maxwell equation

| | Fan 2019 | This work |
|-|----------|-----------|
| Form | Time-domain Yee update (linearized FDTD) | Frequency-domain weak form for `Hz` |
| Discrete equation | Linear update / effective `A x = b` in the paper’s notation | `A(p) x(p) = b` with anisotropic `ρ(ε(ω,B))` and PML stretch |
| Time convention | Paper’s FDTD convention | `e^{-iωt}` Lorentz plasma |

## Differentiation of `A x = b`

Both differentiate the discrete linear system. With design parameter `p`:

    (∂A/∂p) x + A (∂x/∂p) = ∂b/∂p

For material-only parameters and a fixed source, `∂b/∂p = 0`, so

    A (∂x/∂p) = − (∂A/∂p) x

That is exactly the forward-mode sensitivity equation used here and the
frequency-domain counterpart of the paper’s extra forward-mode field.

## Parameter derivative

| | Fan 2019 | This work |
|-|----------|-----------|
| Typical `p` | permittivity or geometric fill | density scale `s` with `f_p² = s f_p,ref²` |
| `∂A/∂p` | derivative of the FDTD update operators | stiffness assembly of `∂ρ/∂s = −ρ (∂ε/∂s) ρ` (`mass_scale=0`) |

## Reuse of the forward operator

Both reuse the same linearized forward operator for the sensitivity field:
one extra forward-mode solve per parameter, not a new nonlinear model.
Here that is literal factorization reuse: one `splu(A)` serves the primal
solve and every sensitivity RHS at fixed `(ω, B, mesh, s_base)`.

## One sensitivity RHS per parameter

Yes, in both frameworks. Cost scales as `N_parameters` forward-mode solves
after the base factorization / base forward run.

## Scaling

| Quantity | Fan FDTD forward-mode | This FEM forward-mode |
|----------|----------------------|------------------------|
| vs `N_parameters` | linear (one forward-mode field each) | linear (one RHS each, same LU) |
| vs `N_outputs` | one forward-mode field yields all outputs of that field | one `dx/dp` yields all linear functionals of `x` |
| vs adjoint | forward better when few params, many outputs | same tradeoff; ~91 rods and one scalar FoM favor adjoint |

## What is independent of FDTD vs FEM

- Stationarity of the discrete residual and its total derivative.
- The identity `A dx = −(dA)x` when `db = 0`.
- Exactness of the discrete derivative relative to finite differences on
  *the same* discrete system.
- The parameter/output scaling argument above.

## What differs (and is not claimed equal)

- Yee staggered grids vs P1 triangles.
- Time-domain broadband updates vs single-frequency complex Helmholtz.
- Reciprocal dielectric examples vs gyrotropic cold-plasma `ε(ω,B)`.
- Memory: FDTD checkpoints vs one complex SuperLU factor.

## Claim

Our implementation is **mathematically equivalent to the forward-mode
differentiation framework** of Fan et al. (2019): same sensitivity equation
for a linear discrete Maxwell operator, one sensitivity solve per design
parameter, factorization/operator reuse, and exact discrete derivatives.
It is **not** a numerical reproduction of their FDTD examples, and the
spatial discretizations are not the same.
