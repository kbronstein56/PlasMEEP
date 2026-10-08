# Publication benchmark pass criteria

Frozen before the ellipse comparison in this file is computed.
These thresholds are not to be loosened after the run.

## Angle and duality convention

Hamid and Cooray, Advanced Electromagnetics, 30 Dec 2016, Table I and Fig. 3.

- Elliptic cylinder, axial ratio 2, semi-major axis along x, `k0 * a = 1`, so `k0 * b = 0.5`.
- Their TM ferrite problem is the dual of this Hz solver with `μ = 1` and the in-plane permittivity equal to their in-plane permeability.
- Tensor, as printed in the Fig. 3 caption: `ε_xx = ε_yy = 7/9`, `ε_xy = -5/18`, `ε_yx = +5/18`.
- Their Maxwell equations use the `e^{+jωt}` signs (`∇×E = -j k0 Z μ·H`). This solver uses `e^{-iωt}`. The real antisymmetric tensor is inserted as printed. The gyrotropic sign is not flipped after the comparison.
- Incidence angle `0°` is a wave traveling toward `-x`. Observation angle is the standard polar angle from `+x`: `0°` is backscatter for that incidence, `180°` is forward. This is the convention in their Section 4, where the forward width sits at `180°`.
- The published column used as the reference is their separation-of-variables column (five digits), not the older integral-equation column.

Published `σ/λ`:

| incidence | observation | published |
|----------:|------------:|----------:|
| 0° | 45° | 0.13781 |
| 0° | 90° | 0.73693 |
| 45° | 0° | 0.14538 |
| 45° | 90° | 0.62262 |
| 90° | 0° | 0.77898 |
| 90° | 45° | 0.68662 |

## Acceptance

On the finest mesh, every row above:

- relative error `|σ_FEM − σ_published| / σ_published ≤ 0.01`

Mesh sequence of at least two body-fitted resolutions. The last step, on every row:

- relative change `≤ 0.003`

The same extractor must first match an independent circular-cylinder Mie width to a relative error `≤ 0.005` at the same `k0 a = 1` and `ε = 4`. If that calibration misses, the ellipse comparison is not interpretable and is not marked PASS.

Reciprocity is reported, not used to edit the tensor: `σ(α, φ; ε)` against `σ(φ, α; ε^T)` on the pairs in their Table I.
