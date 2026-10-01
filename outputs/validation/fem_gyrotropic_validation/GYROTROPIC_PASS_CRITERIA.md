# Gyrotropic pass criteria

Written before the gyrotropic electromagnetic comparisons in this campaign.
These thresholds are not to be loosened in order to produce a pass.
Near a null, judge absolute complex amplitude and absolute power, not decibels.

## Tensor identities

- `rho @ epsilon - I` and `epsilon(B) - epsilon(-B).T` at or near floating-point roundoff, target below `1e-12` relative
- `B = 0` reproduces the validated scalar Drude value below `1e-12`

## Smooth homogeneous cases

- normalized Hz L2 ≤ 0.1%
- phase ≤ 0.1 degree where the mode is propagating and well conditioned
- manufactured solutions: L2 order approaching 2 and H1 order approaching 1 for P1 elements

## Planar and slab

- ≤ 0.05 dB on significant transmission or reflection
- ≤ 0.5 degree phase

## Scatterer and multiport

- ≤ 0.10 dB on significant channels
- ≤ 1 degree phase
- last mesh step much smaller than those allowances

## Onsager and B reversal

- constitutive and matrix residuals much smaller than the discretization error of the same mesh
- `Sij(+B) = Sji(-B)` is the relation under test
- `Sij(+B) = Sji(+B)` is not required
