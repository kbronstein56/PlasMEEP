# Comparison to Fan-group differentiable Maxwell / inverse design

Primary method reference:

Hughes, Williamson, Minkov, Fan,
“Forward-Mode Differentiation of Maxwell’s Equations,”
ACS Photonics 6, 3010–3016 (2019).

Plus the standard adjoint / inverse-design framework used throughout the
Fan group (frequency-domain and time-domain).

## What is the same (mathematical framework)

| Item | Fan-group practice | This work |
|------|--------------------|-----------|
| Discrete Maxwell residual | Linear update / `A x = b` | Frequency-domain P1 FEM `A(s) x = b` |
| Forward sensitivity | `A dx = −(dA)x` | Identical |
| Operator reuse | Same linearized operator / update | One SuperLU factor of `A` for all sources and adjoints |
| Adjoint for few objectives | One adjoint solve per scalar FoM | Six adjoint RHS (one per source) for the six-port FoM |
| Scaling | Forward ~ N_param; adjoint ~ N_obj | 16 tied (or 91) params, one scalar FoM → adjoint preferred |

## What differs (not claimed equivalent as a device)

| Item | Typical Fan demos | This circulator |
|------|-------------------|-----------------|
| Discretization | Yee FDTD or other | Anisotropic P1 FEM + PML stretch |
| Material | Often reciprocal dielectrics / permittivity pixels | Cold-plasma Lorentz gyrotropic `ε(ω,B)`, density scale `s` |
| Device | Couplers, scatterers, etc. | Six-port magnetized plasma circulator |
| Objective | Intensity / efficiency | Weighted guide-normal port powers, all six sources |
| Published digit match | N/A | **Not claimed** — different physics and mesh |

## Claim

The differentiation and inverse-design **methodology** is the same
forward/adjoint discrete-Maxwell framework. The **physical device,
constitutive model, and FEM discretization are different**. No performance
superiority versus a published Fan circulator benchmark is claimed.
