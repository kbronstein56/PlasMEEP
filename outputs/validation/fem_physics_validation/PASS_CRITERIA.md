# Predeclared pass criteria

Written before the electromagnetic comparison results in this campaign.
These thresholds are not to be loosened in order to produce a pass.
If a case is resonant or a quantity is near a null, the report must say so
and judge it with absolute complex amplitude and absolute power, not with dB.

## Simple smooth analytic cases

Manufactured solutions, homogeneous media, PEC guide modes, Fresnel
interfaces, and planar slabs.

- magnitude error ≤ 0.01 dB on significant transmitted or reflected power
- phase error ≤ 0.1 degree
- normalized field error ≤ 0.1% where an L2 field norm is meaningful
- manufactured-solution orders: L2 approaching 2 and H1 approaching 1 for P1 elements

## Single cylinder and coated cylinder

- ≤ 0.05 dB on significant scattering or transmission quantities
- ≤ 0.5 degree phase
- ≤ 1% complex-field L2 error on the comparison cut

## Multiple cylinders

- ≤ 0.10 dB on significant channels
- ≤ 1 degree phase
- FEM discretization uncertainty, measured by the last two meshes, much smaller than those thresholds

## Near-null quantities

Use absolute power and absolute complex amplitude. A decibel error on a
near-zero reflection is not a failure by itself.
