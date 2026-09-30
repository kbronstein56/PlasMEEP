# Independent B=0 weak-form and PML audit

This derivation was written from Maxwell's equations and then compared with `assemble_anisotropic`, `pml_sx_sy`, `ElementSampler.fields`, and `guide_normal_flux`. It does not treat a passing numerical test as evidence that a sign is right.

## Phasor convention

Fields are the real part of `F(x) exp(-i ω t)`. Outgoing dependence on `+x` is `exp(+i k_x x)` with `Im(k_x) ≥ 0`. Outgoing cylindrical waves are `H_n^{(1)}`. In the code units `ε0 = μ0 = c = 1`, so `ω = k0 = 2 π f_a` and `f_a` is the Meep ordinary frequency.

## Hz equation

Take `H = (0, 0, Hz)` and an isotropic, possibly complex, relative permittivity `ε(x, y)`. Faraday and Ampere in the `e^{-iωt}` convention, with `μ = 1`, are

```
∇ × E = i ω H
∇ × H = -i ω ε E
```

The second equation gives

```
E = (i / ω) (1/ε) (∂y Hz, -∂x Hz, 0)
```

because `∇ × H = (∂y Hz, -∂x Hz, 0)`. Substituting into Faraday's `z` component, `∂x Ey - ∂y Ex = i ω Hz`, produces

```
∂x [ (1/ε) ∂x Hz ] + ∂y [ (1/ε) ∂y Hz ] + ω² Hz = 0
```

With `ρ = 1/ε` and `k0 = ω`,

```
div(ρ ∇ Hz) + k0² Hz = 0
```

or, with a source,

```
-div(ρ ∇ Hz) - k0² Hz = f
```

For a symmetric tensor `ρ` the same steps give `div(ρ ∇ Hz) + k0² Hz = 0`, where `(ρ ∇ Hz)_x = ρ_xx ∂x Hz + ρ_xy ∂y Hz` and `(ρ ∇ Hz)_y = ρ_yx ∂x Hz + ρ_yy ∂y Hz`.

## Weak form

Multiply by a test function `v`, integrate, and integrate the divergence by parts. On a boundary where the boundary term is dropped,

```
∫ (ρ ∇ Hz) · ∇ v dA - k0² ∫ Hz v dA = ∫ f v dA
```

The boundary integrand is `v (ρ ∇ Hz) · n`. Dropping it enforces the natural condition `(ρ ∇ Hz) · n = 0`.

The element matrix is therefore stiffness minus `k0²` times the consistent mass. On a P1 triangle the consistent-mass weights are `area/6` on the diagonal and `area/12` off the diagonal. Those are the exact integrals of the products of barycentric functions, not a lumped mass.

## Comparison with `assemble_anisotropic`

The loop over local basis indices uses trial index `j` and test index `i`. The gradient products are

| term in `(ρ ∇ φ_j) · ∇ φ_i` | code |
| --- | --- |
| `ρ_xx ∂x φ_j ∂x φ_i` | `rho_xx * (sy/sx) * bx[:,j] * bx[:,i]` |
| `ρ_xy ∂y φ_j ∂x φ_i` | `rho_xy * by[:,j] * bx[:,i]` |
| `ρ_yx ∂x φ_j ∂y φ_i` | `rho_yx * bx[:,j] * by[:,i]` |
| `ρ_yy ∂y φ_j ∂y φ_i` | `rho_yy * (sx/sy) * by[:,j] * by[:,i]` |

`bx[:,j]` is `∂φ_j/∂x` and `by[:,j]` is `∂φ_j/∂y`, including the signed double-area denominator. The area used in the mass is `½|twice|`. A clockwise triangle flips the sign of `twice` and of both gradient components together, so the products are unchanged. The matrix entry is `Ke - Me`. That is the weak form above, not a strong-form second derivative. The comment in the assembler writes `∂x∂x` as shorthand for the integrated gradient product.

Off-diagonal `ρ` is multiplied in the order required by `(ρ ∇ φ_j) · ∇ φ_i`. The full-complex manufactured solution is the numerical check of that index order; this paragraph is the algebraic check.

No sign discrepancy was found in the volume terms.

## Ex, Ey

Both reconstructions implement the same map. `fields_from_hz` sets `(vx, vy) = (∂y Hz, -∂x Hz)` and

```
Ex = (i/ω) (ρ_xx vx + ρ_xy vy)
Ey = (i/ω) (ρ_yx vx + ρ_yy vy)
```

`ElementSampler.fields` uses the element gradient of the P1 field and the element `ρ`, with `gy` and `-gx` in those same slots. For isotropic `ρ` this is `(Ex, Ey) = (i/ω) ρ (∂y Hz, -∂x Hz)`. The factor is `+i/ω`, not `-i/ω`.

## Poynting vector and guide power

`S = ½ Re(E × H*)`. With `H = (0, 0, Hz)`,

```
Sx = ½ Re(Ey conj(Hz))
Sy = -½ Re(Ex conj(Hz))
```

`guide_normal_flux` uses those expressions. The sign of `Sy` is required by the cross product; it is not a free choice.

For the PEC-guide mode `Hz = cos(k_y y) exp(i β x)` in vacuum, `Ey = (β/k0) Hz`, so

```
∫ Sx dy = ½ (β/k0) ∫ cos²(k_y y) dy = ½ (β/k0) (width/2)
```

when `k_y = π/width`. `guide_eh_close` uses `p_exact = 0.5 * Re(β) * (width/2) / k0`. The factor `1/2` from the cosine integral and the factor `1/2` from the Poynting definition are both present. They are not the same factor counted twice.

## PEC

On a horizontal wall, tangential `E` is `Ex ∝ ρ ∂y Hz`. The natural condition `∂n Hz = 0` is therefore `Ex = 0`. That is PEC for this polarization. Setting `Hz = 0` would constrain tangential `H`, which is PMC here. The assembler states that, and it does not apply a Dirichlet `Hz = 0` condition on metal. Metal triangles are removed by setting `ρ = 0` and dropping both stiffness and mass, so the air-side condition remains the natural Neumann condition.

## Dielectric interfaces

`Hz` is single-valued in the continuous P1 space, so tangential `H` is continuous. The natural flux `(ρ ∇ Hz) · n` is continuous across an interior edge because the boundary terms from the two sides cancel in the weak form. That flux is tangential `E`. Element-wise constant `ρ` is the correct place for a jump in `ε`. No extra interface penalty is present, and none should be.

## Loss sign

For `e^{-iωt}`, a passive conductivity enters as `ε = ε' + i σ/(ε0 ω)`, so `Im(ε) > 0`. The cold-plasma Drude value in ordinary-frequency units is

```
ε = 1 - f_p² / (f² + i f γ)
```

with `f`, `f_p`, and `γ` all in the same Meep frequency unit. The imaginary part is positive for `γ > 0`. The absorbed power density consistent with this convention is `(ω/2) Im(ε) |E|²`. The volume diagnostic uses `0.5 * k0 * Im(ε) |E|²` with `k0 = ω`. A minus sign there would describe an active medium.

`ρ = 1/ε` is used in the stiffness. The mass term stays `k0²`, not `k0² ε`. Putting `ε` in the mass and `1` in the stiffness would be the wrong polarization.

## Source

The scatterer load is the barycentric interpolant of a unit delta: the three nodes of the source triangle receive the barycentric weights. That is the consistent load for `∫ f v` with `f = δ`. It is not the coefficient `i/4` in front of `H_0^{(1)}`.

The fundamental solution of `(∇² + k0²) G = -δ` is `(i/4) H_0^{(1)}(k0 r)` in this time convention. The FEM delta and the Hankel function therefore differ by a complex scale. Every scatterer comparison in this campaign is a ratio `Hz(object) / Hz(vacuum)` against `analytic / H_0^{(1)}`. The same load is used for the object and the vacuum solve, so the scale cancels. Absolute flux of a raw delta must not be compared with an unscaled Hankel field. The guide acceptance test does not use this delta; it imposes the exact mode as Dirichlet data.

## PML, derived again

Stretch `x̃ = ∫ s_x dx'` and `ỹ = ∫ s_y dy'`, with

```
s_x = 1 + i σ_x / ω,    s_y = 1 + i σ_y / ω
```

and `σ ≥ 0`. A wave `exp(i k_x x̃)` with `k_x > 0` then decays as `exp(-k_x ∫ σ_x/ω dx)`. The imaginary part is positive. A negative imaginary stretch would amplify an outgoing wave.

The chain rule is `∂/∂x̃ = s_x^{-1} ∂/∂x`. The area element transforms as `d x̃ d ỹ = s_x s_y dx dy`. Substituting into `∫ (ρ ∇̃ u) · ∇̃ v d x̃ d ỹ` gives the weights

| tensor entry | weight in the physical-coordinate integral |
| --- | --- |
| `ρ_xx` | `s_y / s_x` |
| `ρ_xy` | `1` |
| `ρ_yx` | `1` |
| `ρ_yy` | `s_x / s_y` |
| mass `k0²` | `s_x s_y` |

The off-diagonal weights are exactly 1. The factors `s_x` and `s_y` from the Jacobian cancel the two stretch factors in `∂/∂x̃` and `∂/∂ỹ`. This is the non-obvious term. The assembler multiplies only the diagonal gradient terms by `sy/sx` and `sx/sy`, leaves `ρ_xy` and `ρ_yx` unstretched, and multiplies the mass by `s_x s_y`. That matches this table. In vacuum the isotropic reduction is `∫ (s_y/s_x) ∂x u ∂x v + (s_x/s_y) ∂y u ∂y v - k0² s_x s_y u v`.

`pml_sx_sy` evaluates `σ` at the element centroid. `σ` is zero more than `dpml` inside the box `[0, nx] × [0, ny]`, and

```
σ = σ_max (distance into the layer / dpml)²
```

on each side. The profile is quadratic and starts at zero on the inner face. The normalization is the polynomial PML value

```
σ_max = -ln(R) * (p+1) / (2 * thickness)
```

with `p = 2` and `R = 1e-15`, implemented as `-log(1e-15) * 3 / (2 * dpml)` times `PML_SIGMA_SCALE` (default 1). `ω` in `s = 1 + i σ/ω` is `2 π f_a`, the same `ω` as `k0`.

The stretch does not depend on the incidence angle. Normal incidence uses `s_y = 1` in a region with `σ_y = 0`. Oblique incidence uses both `s_x` and `s_y` when the ray has entered a corner, or only `s_x` when `σ_y = 0`. The Bloch reflection test places PML only on `+x` by shifting the mesh so the other three layers are not entered. That is a domain restriction, not a different formula.

## Result of the sign audit

Every checked minus sign and factor of `i` or `ω` in the volume form, the field reconstruction, the Poynting vector, the PEC condition, the loss sign, and the PML weights agrees with this derivation. The delta-source scale is a declared convention and cancels in the ratio tests. No discrepancy required a stop.
