# Independent magnetized-plasma tensor

Phasor convention: fields are the real part of `F exp(-i ω t)`, the same convention as the frozen B=0 FEM. `ε0 = μ0 = c = 1` in the solver. Ordinary frequency `f = ω/2π` is the Meep frequency unit. `B0` is along `+z`.

This derivation starts from the cold-fluid electron equation. It does not start from the Meep susceptibility code.

## Signs that are fixed here

- Electron charge is `-e` with `e > 0`.
- Cyclotron frequency `ω_c = e B_z / m`. It is positive when `B_z > 0`.
- Lorentz force is `q(E + v × B)` with `q = -e`.
- Time factor `exp(-i ω t)`.
- Collision drag is `-m ν v` with `ν > 0`.

## Equation of motion

```
m dv/dt = -e E - e (v × B) - m ν v
```

With `B = B_z z-hat` and `v × B = B_z (v_y, -v_x, 0)`,

```
(-iω + ν) v_x + ω_c v_y = -(e/m) E_x
(-iω + ν) v_y - ω_c v_x = -(e/m) E_y
(-iω + ν) v_z           = -(e/m) E_z
```

Polarization of the electrons is `P = n e v / (iω)`, because `J = -n e v` and `J = -iω P`. Then `ε0 χ E = P`.

## Tensor

Write `U = ω + iν` and `D = U² - ω_c²`. Define

```
η = ω_p² ω_c / (ω D)
ε_⊥ = 1 - ω_p² U / (ω D)
ε_zz = 1 - ω_p² / (ω U)
```

with `ω_p² = n e² / (ε0 m)`. The result of eliminating `v` is

```
ε_xx = ε_yy = ε_⊥
ε_xy = +i η
ε_yx = -i η
ε_zz = ε_zz
```

In ordinary frequencies `f, f_p, γ, f_c`, with `ν = 2πγ` and `ω_c = 2π f_c`, the same expressions hold with `(ω, ω_p, ν, ω_c)` replaced by `(f, f_p, γ, f_c)`. At `f_c = 0` this is the validated scalar Drude value

```
ε = 1 - f_p² / (f² + i γ f)
```

The electron cyclotron resonance sits in the circular mode `(E_x, E_y) ∝ (1, i)`. That mode's eigenvalue is `ε_⊥ - η`, whose denominator vanishes at `f = f_c`. The opposite mode `(1, -i)` resonates at `f = -f_c`. A free electron in `B_z > 0` follows `(cos ω_c t, sin ω_c t)`, which is the same `(1, i)` sense under `exp(-iωt)`.

The existing helper `gyrotropic_drude_eps_eta` matches `ε_⊥` and the magnitude of `η`, and then stores the opposite off-diagonal, `ε_xy = -iη`. That helper is not the source of this derivation. The FEM tests use `ε_xy = +iη`.

At the stored Meep Faraday point (`a = 0.028 m`, 5 GHz, `f_p = 2 GHz`, `|B| = 0.05 T`), this Lorentz tensor gives `κ(+B) = -0.078460 / a`. The stored Meep measurement is `-0.078479 / a` (0.024% in magnitude). The helper's theory value `+0.078460 / a` is the same magnitude with the opposite off-diagonal. The sign that matches both the Lorentz force and the Meep run is `ε_xy = +iη`. The helper is left unchanged; B=0 does not use the off-diagonal.

## Inverse of the xy block

```
det = ε_⊥² - η²
ρ_xx = ρ_yy = ε_⊥ / det
ρ_xy = -i η / det
ρ_yx = +i η / det
ρ_zz = 1/ε_zz
```

So `ρ_xy = -ρ_yx`, and `ρ(B) = ρ(-B)^T` because `η(-B) = -η(B)`.

## What the 2D Hz model actually solves

`H = (0, 0, H_z(x,y))` and `∂/∂z = 0`, so the wave vector lies in the xy plane, perpendicular to `B`. This is the Voigt configuration. It is not Faraday propagation along `B`.

Ampere's law gives

```
(E_x, E_y) = (i/ω) ρ (∂y H_z, -∂x H_z)
```

with the 2×2 `ρ` above. For constant `ρ_xx = ρ_yy` and `ρ_xy = -ρ_yx`, the mixed derivatives cancel and

```
div(ρ ∇ H_z) + ω² H_z = 0
```

with effective index

```
n² = 1/ρ_xx = (ε_⊥² - η²) / ε_⊥
```

`n²` depends on `η²`, so a plane-wave propagation constant is even in `B`. The polarization `E_x / H_z` is odd in `B`. Faraday rotation along `z`, `κ = (k_+ - k_-)/2` with `k_± = ω sqrt(ε_⊥ ± η)`, is a property of the tensor for `k ∥ B`. The z-invariant Hz solver cannot propagate that wave. The campaign tests `κ` analytically, and tests the Voigt operator in the FEM.

## Interface condition

`H_z` is continuous. The normal flux is

```
(ρ ∇ H_z) · n = ρ_xx ∂_n H_z + ρ_xy (1/r) ∂_φ H_z
```

on a circular interface, with the exterior values `ρ_xx = 1`, `ρ_xy = 0`. Each cylindrical harmonic stays decoupled. The scattered coefficient is a scalar Mie formula with that derivative jump. At `B = 0` the jump reduces to the validated isotropic Mie condition.
