# Gyrotropic FEM validation report

Status: **GYROTROPIC_FEM_VERIFIED_AND_VALIDATED**

The 2D Hz operator with the cold-plasma tensor solves the Maxwell problem that this formulation can represent. That problem is the Voigt configuration, `k` perpendicular to `B`. Adjoints, optimization, and the 91-bulb circulator were not started.

## Geometry

The production bulb is built in `Add_Bulb` with plasma radius `4.6 * r_bulb_inner / 6.5 = 0.230 a`. `sixport_common.r_plasma = 0.250 a` is stored on trainable-rod records and is not the material cylinder. No production geometry was changed. Details are in `GYROTROPIC_GEOMETRY_PRECHECK.md`.

## Tensor

The derivation is in `GYROTROPIC_TENSOR_DERIVATION.md`. Time dependence is `e^{-iωt}`. The electron charge is `-e`, `ω_c = e B_z / m` is positive for `B_z > 0`, and

```
ε_xx = ε_yy = ε_⊥,   ε_xy = +i η,   ε_yx = -i η
```

with collisions kept. At `B = 0` this is the validated scalar Drude permittivity. The xy inverse is analytic and matches a numerical inverse to `2e-15`.

`gyrotropic_drude_eps_eta` has the same `ε_⊥` and the same `|η|`, then stores `ε_xy = -iη`. That helper was not used as the derivation and was not edited. At the stored Faraday point the Lorentz `κ(+B)` is `-0.078460/a`. The Meep res-64 measurement is `-0.078479/a`. The helper's theory value `+0.078460/a` is the opposite off-diagonal. The sign that agrees with the Lorentz force and with that Meep run is `ε_xy = +iη`.

## What was compared

Smooth tests use a manufactured load and a homogeneous Voigt wave. Hz converges as `h^2` (order 1.999). The gradient converges as `h` (order 0.999). At `h = 0.01` the Voigt wave has Hz relative error `1.8e-5`. `β` is even in `B` because the index depends on `η^2`. `E_x / H_z` is odd in `B` and matches the reconstructed field.

A z-invariant code cannot propagate a Faraday wave along `z`. Faraday `κ = (k_+ - k_-)/2` was therefore checked as a property of the tensor, including `B → 0`, sign reversal, and the one stored Meep run. The FEM tests of this operator are the Voigt wave, the interface, the slab, and the cylinders.

Normal incidence on a slab reduces to a scalar interface with `ε_eff = (ε_⊥^2 - η^2)/ε_⊥`. Transmission is even in `B`. That is the correct Maxwell result, not a missing nonreciprocal effect. Oblique incidence keeps a `ρ_xy k_y` term. There `R(k_y, +B) = R(-k_y, -B)` exactly, while `R(k_y, +B)` differs from `R(-k_y, +B)`. The FEM field on those problems is within `0.043%`.

The production-frequency slab, where the plasma is evanescent, matches its analytic field to `1.1e-5` and matches the existing B=0 stack routine to `2e-16` in the reflection coefficient.

One, two, and three cylinders use an independent cylindrical-harmonic `T_n`. The interior boundary term is `ρ_xx ∂_r H_z + ρ_xy (1/r) ∂_φ H_z`. At `B = 0` it is the validated Mie formula. Exterior coupling is the same Graf addition used for isotropic clusters, with this `T_n` on the diagonal. FEM probe errors are below `3e-4 dB` and `0.003 deg`, and they drop on the finer mesh. Side amplitudes swap when `B` reverses. A y-flipped pair was solved separately and matches the same reference.

## Onsager and a nonreciprocal port

`ε(B) = ε(-B)^T` holds exactly for this tensor, and the same identity holds for `ρ`. On the assembled interior block the FEM matrix residual is `0`. The Green function satisfies `G(x, y; +B) = G(y, x; -B)` to `2e-15`.

Three nodal ports around one offset disk give the same statement as an S-matrix of unit nodal loads:

- `S_ij(+B)` versus `S_ji(-B)` is `6e-15`
- `S_ij(+B)` versus `S_ji(+B)` is `0.63%`, and that difference changes only from `0.644%` to `0.631%` when `h` is halved
- at `B = 0`, `S_ij` versus `S_ji` is `9e-16`

`S_ij(+B) = S_ji(+B)` is not expected and is not true for these ports.

A straight parallel-plate channel with one rod, observed through a full-width integral of `H_z`, changes by only `2e-5` when `B` reverses. That integral keeps the even mode and cancels the odd scattered field, so it is a weak nonreciprocal observable. It is not evidence that the tensor is reciprocal.

## Passivity

For `e^{-iωt}`, the dissipation density uses `(ω/2) E† ((ε - ε†)/(2i)) E`. Those eigenvalues are nonnegative. A homogeneous Voigt wave satisfies `dS_x/dx + p_abs = 0` to `1.6e-9`. In the FEM disk, absorbed power is positive, the same for both signs of `B`, and proportional to the collision rate. At collision rate `0.01` and `h = 0.01`, the circular-contour flux closes the volume integral to `0.54%`. At a much smaller collision rate the absorbed power falls below the reactive-flux sampling error, so that relative contour residual is not used as a decibel gate.

## Meep

Meep is not the ground truth for a dispersive cylinder. The comparison that is valid for Meep in this campaign is the homogeneous Faraday cell already stored at resolution 64. Magnitude and `B` reversal agree with the Lorentz tensor; the older theory sign does not. No FEM coefficient was adjusted to Meep.

## Not started

The full 91-bulb magnetized forward model, adjoints, and optimization were not started. They wait on a review of this status.
