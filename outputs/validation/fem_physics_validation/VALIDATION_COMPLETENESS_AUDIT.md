# Validation completeness audit

Status of the B = 0 package: **B0_FEM_PARTIALLY_VALIDATED**.

`B0_FEM_VERIFIED_AND_VALIDATED` is withdrawn. A related test is not counted as the requested test. Thresholds are the ones frozen in `PASS_CRITERIA.md` before the electromagnetic comparisons. They were not loosened after any result.

Frozen thresholds used below:

- smooth analytic field: ≤ 0.01 dB on significant power, ≤ 0.1° phase, ≤ 0.1% normalized field L2
- one cylinder, bare or coated: ≤ 0.05 dB, ≤ 0.5° phase, ≤ 1% complex-field L2 on the cut
- several cylinders: ≤ 0.10 dB, ≤ 1° phase, and the last mesh step much smaller than those allowances
- near a null: absolute complex amplitude, not decibels

Evidence files live in this directory unless a path is given. `plasma_ring_campaign/` means `outputs/validation/fem_meep_validation/plasma_ring_campaign/`.

## 1. Method of manufactured solutions

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 1a | scalar real ρ | yes | imposed Hz | n = 8, 16, 32, 64 | no | mms.json | n=64 L2 = 4.25e-4, L2 order 1.998, H1 order 0.999 | L2 order → 2, H1 order → 1 | PASS | — |
| 1b | scalar complex ρ | yes | imposed Hz | same | no | mms.json | L2 order 1.998, H1 order 0.999 | same | PASS | — |
| 1c | diagonal anisotropic complex ρ | yes | imposed Hz | same | no | mms.json | L2 order 1.997, H1 order 0.999 | same | PASS | — |
| 1d | full complex tensor, off-diagonal | yes | imposed Hz | same | no | mms.json | L2 = 3.90e-4, L2 order 1.998, H1 order 0.999 | same | PASS | — |

## 2. Homogeneous-medium propagation

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 2a | ε = 1 | yes | plane wave | h = 0.05, 0.025, 0.0125 | no | planar_fem.json, canonical_suite.json | rel L2 1.89e-4 at h=0.025; 4.74e-5 at h=0.0125 | ≤ 0.1% field | PASS | — |
| 2b | ε = 3.8 | yes | plane wave | same | no | canonical_suite.json | Sx relative error 2.4e-4 at h=0.0125; field rel L2 8.30e-4 at h=0.025 | ≤ 0.1% field | PASS | — |
| 2c | positive complex ε | yes | ε = 2+0.3i | same | no | canonical_suite.json | Sx relative 1.3e-4; field rel L2 4.02e-4 at h=0.025 | ≤ 0.1% field | PASS | — |
| 2d | production B=0 plasma ε | yes | ε(3.85 GHz) | same | no | canonical_suite.json | field rel L2 1.39e-4 at h=0.0125; Sx relative 2.4e-4 | ≤ 0.1% field | PASS | — |
| 2e | analytical propagation constant | yes, through the field | exact kx | same | no | canonical_suite.json | field error above. A phase-slope `beta_rel` in that file is not a valid metric for complex or evanescent kx and is not used | field ≤ 0.1% | PASS | do not quote `beta_rel` for complex kx |
| 2f | analytical E/H relation | yes | Ey/Hz from the plane wave | h=0.0125 | no | conservation.json | point sample 0.34% vacuum, 0.60% plasma at h=0.0125 | ≤ 0.1% if the field bar is applied to this sample | FAIL | one more refinement of the gradient sample |
| 2g | Poynting power | yes | analytic Sx | h=0.0125 | no | canonical_suite.json | Sx relative 6.7e-5 (ε=1), 2.4e-4 (ε=3.8), 1.3e-4 (complex), 2.4e-4 (plasma) | ≤ 0.01 dB on significant power | PASS | — |

## 3. PEC analytical waveguide modes

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 3a | mode shape | yes | cos mode, natural Neumann plates | widths 0.8, 1.0, 1.4; k0 3.2, 5.0, 7.5; h down to 0.00125 | no | guide_refined.json | hardest case width=1, k0=5: rel L2 2.97e-4 at h=0.00125, 1.026e6 DOFs | ≤ 0.1% | PASS | — |
| 3b | propagation constant | yes, through the mode | exact β | same | no | guide_refined.json | same field match | ≤ 0.1% field | PASS | — |
| 3c | cutoff | yes | evanescent mode | h=0.01 | no | canonical_suite.json | rel L2 1.7e-4 class at h=0.01 on the cutoff member | ≤ 0.1% | PASS | — |
| 3d | E/H relation | yes at h=0.01 | analytic Ey/Hz | not repeated at h=0.00125 | no | canonical_suite.json | propagating EH relative from 0.30% to 3.3% at h=0.01 | ≤ 0.1% | FAIL | E/H on the h=0.00125 mode |
| 3e | analytic guide-normal power | yes | element-centered Poynting | h to 0.00125 | no | guide_refined.json | power relative 4.25e-4 at h=0.00125 (about 0.0018 dB) | ≤ 0.01 dB | PASS | — |
| 3f | multiple widths and frequencies | yes | 3 widths × 3 k0 | h=0.04, 0.02, 0.01 where the mesh divides the guide | no | canonical_suite.json | 27 guide rows | each propagating mode ≤ 0.1% at sufficient h | PASS | — |
| 3g | FEM mesh convergence | yes | same mode | h = 0.04 … 0.00125 | no | guide_refined.json | rel L2 0.212, 0.069, 0.0185, 0.00472, 0.00119, 0.000297 | observed order near 2 | PASS | — |
| 3h | Meep comparison | no | — | — | not run | — | — | ≤ 0.01 dB | UNRESOLVED | Meep parallel-plate guide |

## 4. Fresnel interfaces

Dirichlet boxes impose the analytic field on the boundary. The interior error is the FEM evidence. Extracted scattering coefficients are the driven-slab rows in section 5.

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 4a | air/quartz | yes | Fresnel Hz | h=0.025 | no | planar_fem.json | interior rel L2 within the 8.31e-4 maximum of the set | ≤ 0.1% | PASS | — |
| 4b | quartz/air | yes | Fresnel Hz, including 40° TIR with \|r\|=1 | h=0.025 | no | planar_fem.json | largest entry in the set, rel L2 8.31e-4 | ≤ 0.1% | PASS | — |
| 4c | positive complex material | yes | Fresnel Hz | h=0.025 | no | planar_fem.json | inside the same maximum | ≤ 0.1% | PASS | — |
| 4d | negative-ε production plasma | yes | Fresnel Hz | h=0.025 | no | planar_fem.json | air→plasma rel L2 3.9e-4 to 4.1e-4 | ≤ 0.1% | PASS | — |
| 4e | normal incidence | yes | Fresnel Hz | h=0.025 | no | planar_fem.json | inside the same maximum | ≤ 0.1% | PASS | — |
| 4f | at least two oblique angles | yes | 15° and 40° | h=0.025 | no | planar_fem.json | both angles in the 24-row set | ≤ 0.1% | PASS | — |
| 4g | complex r | yes | closed form | h=0.025 field | no | canonical_suite.json, planar_fem.json | analytic r stored; interior field rel L2 ≤ 8.31e-4 | field ≤ 0.1% | PASS | — |
| 4h | complex t | yes | closed form, t=1+r | same | no | same | same field bound | field ≤ 0.1% | PASS | — |
| 4i | power coefficients | yes, analytic | R+T | no FEM surface integral on the Fresnel box | no | canonical_suite.json | lossless power_sum = 1 | R+T = 1 | PASS | FEM surface integral on the oblique box |

## 5. Planar multilayer transfer matrices

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 5a | quartz slab | yes | transfer matrix; Fabry–Perot error 2.8e-17 in r | h=0.04 and 0.02 | yes, res 24 and 40 | planar_fem.json, meep_crosscheck.json | FEM T +0.0013 dB, −0.0007° at h=0.02, thickness 0.20. Meep +0.0020 dB / +0.031° at res 24; +0.0024 dB / +0.008° at res 40 | ≤ 0.01 dB, ≤ 0.1° | PASS | — |
| 5b | plasma slab | yes | transfer matrix | h=0.02 | no | planar_fem.json | thickness 0.20: T −0.0017 dB, +0.015°. Thickness 0.40: T −0.0015 dB, +0.025° | ≤ 0.01 dB, ≤ 0.1° | PASS | Meep plasma slab in this package |
| 5c | quartz/plasma/quartz | yes | transfer matrix | h=0.02 is the valid grid | no | planar_fem.json | T −0.0015 dB, +0.015° at h=0.02. h=0.04 is an invalid grid (0.10 layers are not an integer number of 0.04 cells) and is not a solver result | ≤ 0.01 dB, ≤ 0.1° | PASS | — |
| 5d | multiple thicknesses | yes | 0.20 and 0.40, and the 0.10/0.40/0.10 stack | h=0.02 | quartz only | planar_fem.json | all valid-grid rows within 0.0035 dB | ≤ 0.01 dB | PASS | — |
| 5e | multiple frequencies | yes | 3.542, 3.85, 4.158 GHz | h=0.02 | no | planar_fem.json | all within 0.0035 dB and 0.033° | ≤ 0.01 dB, ≤ 0.1° | PASS | — |
| 5f | complex R and T | yes | transfer matrix | h=0.02 | quartz Meep | planar_fem.json | quartz 0.40: R +0.0034 dB; T errors above | ≤ 0.01 dB, ≤ 0.1° | PASS | — |
| 5g | mesh-conforming interface convergence | yes | same stacks | h=0.04 fails only when faces miss the interface; h=0.02 passes | no | planar_fem.json | aligned stacks converge; off-grid q/p/q at h=0.04 is excluded | faces on element boundaries | PASS | a third slab mesh |

## 6. Bare cylindrical Mie benchmark

Probe ratios are Hz(object)/Hz(vacuum) against the analytic total/H0 ratio. Source (−4.5, 0) in the new scatterer suite. The older disk line cut used source (−6, 0).

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 6a | production radius 0.230 | yes | Mie / one-cylinder T-matrix | h=0.04, 0.02, 0.012 | prior plasma Meep, not this package | scatter_partial.json | h=0.012, 1.68e6 DOFs, worst probe outside_shell −0.0015 dB, −0.031° | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6b | smaller radius 0.12 | yes | Mie | h=0.04, 0.02 | no | scatter_partial.json | h=0.02 worst outside_shell −0.0052 dB, −0.130° | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6c | larger radius 0.40 | yes | Mie | h=0.04, 0.02 | no | scatter_partial.json | h=0.02 worst outside_shell −0.0022 dB, −0.062° | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6d | dielectric cylinder ε=3.8 | yes | Mie | h=0.04, 0.02 | yes, res 24 and 40 | scatter_partial.json, meep_crosscheck.json | FEM h=0.02 worst −0.0003 dB, −0.007°. Meep res 40 forward +8.6e-5 dB, phase about 0 | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6e | several negative-ε cases | yes | ε=−1.5 and ε=−8 | −1.5 at h=0.04 and 0.02; −8 at h=0.04 only | no | scatter_partial.json | ε=−1.5 at h=0.02 worst −0.0195 dB, −0.420°. At h=0.04 the same probe was −0.119 dB, −1.60° and fails; the finer mesh passes. ε=−8 at h=0.04 worst −0.0083 dB, −0.194° | ≤ 0.05 dB, ≤ 0.5° | PASS | second mesh for ε=−8 |
| 6f | production complex plasma | yes | ε(fs) | three meshes | prior Meep | scatter_partial.json | same as 6a | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6g | frequency sweep | analytic yes, FEM no | Mie at 3.2–6.5 GHz, 102 samples | FEM only at 3.85 GHz | no | analytic_sweeps.json | broad peak of the production cylinder at 5.30 GHz, max \|T_n\|=0.999, ε_re=−1.28. No narrow pole at 3.85 GHz | FEM ≤ 0.05 dB at several frequencies | UNRESOLVED | FEM at frequencies other than 3.85 GHz |
| 6h | forward scattering | yes | Mie | yes | dielectric Meep yes | scatter_partial.json | production h=0.012 forward inside the 6a worst-probe bound | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6i | backward scattering | yes | Mie | yes | dielectric Meep yes | scatter_partial.json | inside the same bound | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6j | side scattering | yes | Mie, +60°, −60°, 90° | yes | dielectric Meep side_90 | scatter_partial.json | inside the same bound. Side +60 and −60 agree when the geometry is y-symmetric | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 6k | near field | yes | outside-shell probe | yes | no | scatter_partial.json | production outside_shell −0.0015 dB, −0.031° at h=0.012 | ≤ 0.05 dB, ≤ 0.5° | PASS | full-field L2 map of the new solve |
| 6l | FEM mesh convergence | yes | Mie | three levels on the production disk | no | scatter_partial.json | errors fall from h=0.04 to h=0.012; finest is 0.0015 dB | last step much smaller than 0.05 dB | PASS | — |
| 6m | Meep resolution convergence | dielectric yes; plasma not in this package | Mie | — | dielectric res 24 and 40 pass. Prior plasma disk does not return to Mie | meep_crosscheck.json; prior ring campaign | dielectric forward 8.6e-5 dB at res 40. Prior unshifted plasma disk stays off Mie after a long DFT | ≤ 0.05 dB, registration not collapsed | FAIL | plasma Meep does not meet the Mie bar |
| 6n | Meep sub-cell registration | prior campaign only | Mie | — | several offsets, not averaged | prior ring-campaign records | registration spread is large; one number was not formed | report the spread | FAIL | not registration-independent for plasma |

## 7. Coated-cylinder analytic benchmark

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 7a | production plasma radius | yes | multilayer T-matrix | h=0.04, 0.02 | no | scatter_partial.json | h=0.02, 604838 DOFs, worst outside_shell −0.0037 dB, −0.081° | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 7b | vacuum gap | yes | radii 0.230 / 0.325 / 0.375 kept | same | no | fem_scatterers.py, scatter_partial.json | gap is in the mesh and in the T-matrix | geometry present | PASS | — |
| 7c | quartz shell | yes | ε=3.8 annulus | same | no | scatter_partial.json | quartz area error about 8e-7 at h_edge=0.008 | ≤ 0.05 dB | PASS | — |
| 7d | limiting cases | yes | all-air, air shell = bare core, uniform fill = bare outer cylinder | analytic | no | analytic_maxwell.py self-checks | limits at roundoff | algebraic identity | PASS | — |
| 7e | target frequency 3.85 GHz | yes | coated T-matrix | two meshes | no | scatter_partial.json | see 7a | ≤ 0.05 dB, ≤ 0.5° | PASS | — |
| 7f | nearby frequencies | analytic yes, FEM no | 3.50, 3.85, 4.20, 5.30 GHz | FEM only at 3.85 GHz | no | analytic_sweeps.json | four analytic coated responses stored | FEM ≤ 0.05 dB at nearby f | UNRESOLVED | FEM coated cylinder off 3.85 GHz |
| 7g | FEM mesh sequence | yes | coated T-matrix | h=0.04 and 0.02 | no | scatter_partial.json | finest probe −0.0037 dB, −0.081° | ≤ 0.05 dB, ≤ 0.5° | PASS | third mesh |
| 7h | Meep resolution sequence | no in this package | — | — | not run | — | — | ≤ 0.05 dB | UNRESOLVED | coated-cylinder Meep resolution sequence |
| 7i | Meep registration sequence | no | — | — | not run | — | — | spread reported separately | UNRESOLVED | coated-cylinder registration sequence |
| 7j | complex near field | yes | outside-shell probe | two meshes | no | scatter_partial.json | −0.0037 dB, −0.081° | ≤ 0.05 dB, ≤ 0.5° | PASS | full-field L2 on the new mesh |

## 8. Multiple-cylinder T-matrix

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 8a | independent implementation | yes | Bessel/Hankel code does not call the FEM assembler | not an FEM test | no | analytic_maxwell.py | one-cylinder reduction below | identity | PASS | — |
| 8b | one-cylinder reduction to Mie | yes | Mie coefficients | — | no | analytic_self_checks.json | difference 6.8e-21 | roundoff | PASS | — |
| 8c | m_max truncation | yes | m_max = 2, 3, 4, 5, 6, 8, 10, 12, 14 | — | no | analytic_sweeps.json | coated 7-bulb forward change by m_max=8 is about 2.7e-5; by m_max=12–14 below 1e-8 | observable stable | PASS | — |
| 8d | conditioning | estimated, not repaired | 2-norm of the unscaled dense matrix | — | no | analytic_sweeps.json | bare 7 at m_max=6: cond about 9e5, forward change 4e-7. Coated 7 at m_max=6: cond about 2e6. At m_max=8: cond about 2e10. At m_max=12–14: cond 1e19–1e25 while the forward probe stays put | a moderate condition number | UNRESOLVED | scaled T-matrix. The probe is stable; the raw matrix is not well conditioned past m_max=8 |
| 8e | no FEM/Meep fitting | yes | coefficients are not adjusted to FEM | — | no | analytic_maxwell.py | no scale fit. Vacuum-normalized ratios cancel source amplitude | no fitted scale | PASS | — |

## 9. Two-cylinder production-spacing benchmark

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 9a | T-matrix | yes | two-cylinder T-matrix | — | no | analytic_sweeps.json, scatter_partial.json | compared below | — | PASS | — |
| 9b | FEM convergence | yes for the axis-aligned bare pair | T-matrix | h=0.04 and 0.02 | no | scatter_partial.json | h=0.02 worst outside_shell −0.0128 dB, −0.083° | ≤ 0.10 dB, ≤ 1° | PASS | coated pair has only h=0.04: worst −0.0477 dB, −0.318°, inside the band, one mesh |
| 9c | Meep | no | — | — | not run | — | — | ≤ 0.10 dB | UNRESOLVED | two-cylinder Meep |
| 9d | orientation relative to the Yee grid | FEM of a 30° pair at one mesh; no Meep | T-matrix of the rotated pair | h=0.04 only | no | scatter_partial.json | worst outside_shell −0.0316 dB, −0.232° | ≤ 0.10 dB, ≤ 1°, and a mesh sequence | UNRESOLVED | second mesh, and any Meep orientation run |

## 10. Three-cylinder benchmark

The triangle is not mirror-symmetric across x. Side +60° and side −60° differ for that reason.

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 10a | T-matrix | yes | three-cylinder T-matrix | — | no | cluster_reference.json, scatter_partial.json | old-domain coated FEM-F forward −0.00017 dB, −0.110°. New domain h=0.02 worst −0.0088 dB, −0.178° | ≤ 0.10 dB, ≤ 1° | PASS | — |
| 10b | FEM convergence | yes | same | h=0.04 and 0.02, plus older C/F | no | scatter_partial.json | h=0.02 inside the band | last step smaller than 1° | PASS | — |
| 10c | Meep | no | — | — | not run | — | — | ≤ 0.10 dB | UNRESOLVED | three-cylinder Meep |
| 10d | multipole convergence | yes | m_max sweep includes the 3-cylinder systems | — | no | analytic_sweeps.json | forward probe stable by m_max about 8 | observable change ≪ 0.10 dB | PASS | — |

## 11. Seven-cylinder bare-plasma benchmark

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 11a | T-matrix | yes | bare seven-cylinder T-matrix | — | no | scatter_partial.json | h=0.02 worst outside_shell −0.0160 dB, −0.224° | ≤ 0.10 dB, ≤ 1° | PASS | — |
| 11b | FEM convergence | yes | same | h=0.04 and 0.02 | no | scatter_partial.json | forward phase −0.587° then −0.180°. Shell phase −0.762° then −0.224° | ≤ 0.10 dB, ≤ 1° | PASS | a third mesh so the last step is clearly below 1° |
| 11c | Meep | no | — | — | not run | — | — | ≤ 0.10 dB | UNRESOLVED | bare seven-cylinder Meep |

## 12. Seven-cylinder production coated-bulb benchmark

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 12a | multilayer T-matrix | yes | coated addition theorem | — | no | analytic_sweeps.json | truncation in 8c | observable stable | PASS | — |
| 12b | FEM-C | yes | same, source (−6, 0), domain 18×14 | h=0.030 | prior Meep | cluster_reference.json | forward −0.0139 dB, −0.61°, exterior rel L2 1.31% | ≤ 0.10 dB, ≤ 1° | PASS | — |
| 12c | FEM-F | yes | same | h=0.018 | prior Meep | cluster_reference.json | forward −0.0052 dB, −0.23°, exterior rel L2 0.50% | ≤ 0.10 dB, ≤ 1° | PASS | — |
| 12d | FEM-X | yes | same | h=0.012, 3,028,995 DOFs | prior Meep | fem_cluster7.json | forward T ratio error −0.0023 dB, −0.100°, absolute 0.00313. Backward +0.008 dB, +0.012°. Exterior rel L2 0.21% | ≤ 0.10 dB, ≤ 1° | PASS | off-axis probes on this domain |
| 12e | final-two-mesh difference | yes | FEM-F versus FEM-X | h=0.018 to 0.012 | no | fem_cluster7.json | forward step −0.0030 dB and −0.132° | much smaller than 0.10 dB and 1° | PASS | — |
| 12f | T-matrix truncation uncertainty | yes | m_max sweep | — | no | analytic_sweeps.json | forward change from m_max=8 onward below about 3e-5 | ≪ FEM residual | PASS | — |
| 12g | Meep comparison | prior campaign, not re-run | the FEM-X number is the analytic one to 0.002 dB | — | original grid and shifted grids, not averaged | prior seven-bulb long-DFT records | original-grid late \|Δ\| about 0.43–0.46, about −0.6 dB and −15° from FEM, hence from the T-matrix | ≤ 0.10 dB, ≤ 1° | FAIL | Meep does not meet the analytic cluster |
| 12h | field comparison, not one forward scalar | yes on the new domain | probes: forward, backward, +60°, −60°, gap air, outside shell | h=0.04 and 0.02 only | no | scatter_partial.json | h=0.02 forward −0.0106 dB, −0.301°; outside shell −0.0248 dB, −0.412°. Last shell-phase step is 0.77° | values inside 0.10 dB and 1°, and the last step much smaller than 1° | UNRESOLVED | third mesh. The near-probe step is not yet much smaller than 1° |

## 13. PML and domain independence

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 13a | normal propagation | yes | quartz-slab transfer matrix | one mesh, several PML settings | no | canonical_suite.json | dpml 0.6 and 1.4, sigma scale 0.5 and 2, air-pad shift −0.4: within 0.0013 dB and 0.03° of the reference | ≤ 0.01 dB | PASS | — |
| 13b | oblique propagation | no | — | — | no | — | — | ≤ 0.01 dB | UNRESOLVED | oblique PML sweep |
| 13c | cylinder radiation | only one domain | Mie bound on that domain | the Mie meshes | no | scatter_partial.json | production disk error 0.0015 dB on one domain | a padding sweep | UNRESOLVED | cylinder air-pad sweep |
| 13d | multi-cylinder radiation | only one domain | T-matrix bound | the cluster meshes | no | scatter_partial.json | seven-coated forward 0.011 dB on one new domain; 0.002 dB on the old domain | a padding sweep | UNRESOLVED | cluster padding sweep |
| 13e | air padding variation | yes for the slab | transfer matrix | one mesh | no | canonical_suite.json | pad shift −0.4 within 0.0013 dB | ≤ 0.01 dB | PASS | cylinder padding |
| 13f | PML thickness variation | yes for the slab | transfer matrix | one mesh | no | canonical_suite.json | dpml 0.6 and 1.4 within 0.0013 dB | ≤ 0.01 dB | PASS | — |
| 13g | PML strength variation | yes | transfer matrix | sigma scale 0.5 and 2 | no | canonical_suite.json | within 0.0013 dB and 0.03° | ≤ 0.01 dB | PASS | — |
| 13h | outer-domain variation | slab shift only | transfer matrix | — | no | canonical_suite.json | slab only | cylinder outer-domain sweep | UNRESOLVED | cylinder outer box |
| 13i | quantified reflection | yes as a slab residual | transfer matrix | — | no | canonical_suite.json | transmission change ≤ 0.0013 dB. Not a separately extracted reflection coefficient | ≤ 0.01 dB | PASS | dedicated reflection coefficient |

## 14. Power conservation and passivity

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 14a | lossless homogeneous | yes | analytic Sx | h=0.0125 | no | canonical_suite.json | Sx relative 6.7e-5 for ε=1 | power error ≤ 0.01 dB | PASS | — |
| 14b | lossless dielectric interface and slab | yes | R+T=1 and FEM T | h=0.02 | quartz Meep | planar_fem.json | quartz R=0.31836, T=0.68164, sum 1. FEM T within 0.003 dB | R+T=1 | PASS | — |
| 14c | dielectric cylinder | analytic contour yes | contour at r=1.2, source outside | FEM contour not computed | no | analytic_sweeps.json | net power about 1e-12 | roundoff | PASS | FEM contour |
| 14d | lossless cluster | no | — | — | no | — | — | R+T or contour at roundoff | UNRESOLVED | a lossless multi-cylinder contour |
| 14e | lossy plasma slab | yes | transfer matrix | h=0.02 field; power is analytic | no | conservation.json | thickness 0.40: \|r\|²+\|t\|²=0.99968, absorption proxy +3.2e-4 | absorption ≥ 0 | PASS | — |
| 14f | lossy plasma cylinder | analytic contour yes | same contour | FEM contour not computed | no | analytic_sweeps.json | absorption proxy +1.34e-5 | absorption ≥ 0 | PASS | FEM contour |
| 14g | lossy coated cluster | analytic contour for one coated cylinder | same | FEM contour not computed | no | analytic_sweeps.json | absorption proxy +1.81e-5 | absorption ≥ 0 | PASS | seven-bulb contour |
| 14h | absorption ≥ 0 | yes on 14e–14g | those proxies | — | no | analytic_sweeps.json, conservation.json | all three proxies positive | ≥ 0 | PASS | — |
| 14i | explicit volume absorption | no | Im(ε)\|E\|² not integrated | — | no | — | — | matches the contour proxy | UNRESOLVED | volume integral |

## 15. B = 0 reciprocity

These are Green-function pairs, G(a,b) versus G(b,a), on interior point sources. They are not a port S-matrix.

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 15a | simple guide | no port pair | — | — | no | — | — | \|Sij−Sji\| at roundoff | UNRESOLVED | guide-port Sij/Sji |
| 15b | horns only | no | — | — | no | — | — | same | UNRESOLVED | horns-only pair |
| 15c | one bulb | yes | reciprocity identity | h=0.05, 97313 nodes, coated bulb | no | reciprocity_scatter.json | 6 pairs, max relative \|Gab−Gba\| = 5.5e-14 | roundoff | PASS | a second mesh |
| 15d | small cluster | yes | same | bare pair, 97311 nodes | no | reciprocity_scatter.json | 6 pairs, max relative 3.1e-15 | roundoff | PASS | coated pair and a finer mesh |
| 15e | seven-bulb cluster | no | — | — | no | — | — | roundoff | UNRESOLVED | seven-bulb pairs |
| 15f | multiple source/receiver pairs | yes | same identity | quartz inclusion, one bulb, one pair | no | canonical_suite.json, reciprocity_scatter.json | quartz disk: 15 pairs, max relative 1.7e-15 | roundoff | PASS | — |
| 15g | direct Sij/Sji comparison | yes for the Green function on those three geometries | identity | — | no | same | numbers above | roundoff | PASS | port S-parameters |

## 16. Symmetry and invariance

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 16a | mirror | yes for y-symmetric geometries | side +60 versus side −60, and the T-matrix | the scatter meshes | no | scatter_partial.json | bare disk, dielectric, axis-aligned pair, and seven-hex side probes agree at about 1e-6 to 1e-4. The 3-bulb triangle and the 30° pair are not y-symmetric; their side difference is physical | symmetric geometries agree | PASS | — |
| 16b | rotation | yes | T-matrix of the 30° pair | h=0.04 only | no | scatter_partial.json | worst −0.0316 dB, −0.232° | ≤ 0.10 dB, ≤ 1° | PASS | second mesh |
| 16c | transformed observables | yes | each rotated geometry compared with its own analytic field | h=0.04 | no | scatter_partial.json | same number | same | PASS | — |
| 16d | no averaging used to manufacture symmetry | yes | — | — | registrations were not averaged | this audit | no left/right average was used as a result | no manufactured symmetry | PASS | — |

## 17. Source and port validation

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 17a | exact analytic guide-mode power | yes | closed-form guide power | h to 0.00125 | no | guide_refined.json | power relative 4.25e-4 | ≤ 0.01 dB | PASS | — |
| 17b | reconstructed FEM Poynting | yes | element-centered flux | same | no | guide_refined.py | one-sided nodal differences are O(h) and were not used for this number | element flux | PASS | — |
| 17c | normalization | yes | vacuum ratio cancels the point-source amplitude; guide power is absolute | — | prior straight guide 6.9e-5 dB, not re-run | scatter_partial.json, guide_refined.json | ratios and absolute guide power | source amplitude cancelled or compared absolutely | PASS | repeat of the straight-guide Meep port |
| 17d | forward/backward sign | yes | the analytic traveling wave | the guide meshes | no | guide_refined.json | field rel L2 2.97e-4 includes the sign of β | field match | PASS | — |
| 17e | monitor orientation | yes | guide-normal cut | same | no | guide_refined.py | flux is the guide-normal component | orientation stated | PASS | — |
| 17f | mesh-density independence | yes | analytic power | h sequence | no | guide_refined.json | power relative falls to 4.25e-4 | ≤ 0.01 dB at the fine mesh | PASS | — |
| 17g | monitor-sampling independence | yes | three stations on one mode | h=0.01 | no | canonical_suite.json | width=1.4, k0=5: near 0.317412495, mid 0.317412385, far 0.317412497 | stations agree | PASS | the absolute one-sided power at this h is still 1.5% high; the element flux is the converged number |
| 17h | monitor-location independence | yes | same three stations | h=0.01 | no | canonical_suite.json | relative station scatter about 3e-7 | stations agree | PASS | — |

## 18. Linear-system numerical accuracy

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 18a | residual \|\|Ax−b\|\|/\|\|b\|\| | yes | direct SuperLU | the scatter and guide solves | no | scatter_partial.json | residuals about 1e-14 to 1e-15. Example: ε=−1.5 at h=0.02, residual 6.9e-14 | near roundoff | PASS | — |
| 18b | solver tolerance | not applicable | SuperLU has no iterative tolerance | — | no | fem_validated_solver.py | direct factorization | residual above | PASS | — |
| 18c | factorization consistency | no second factorization of the same matrix | — | — | no | — | — | two factorizations agree | UNRESOLVED | repeat one factorization |
| 18d | FEM condition estimate | no | T-matrix condition is in 8d, which is a different matrix | — | no | — | — | an estimate | UNRESOLVED | FEM matrix condition |

## 19. Geometry convergence

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 19a | body-fitted circles | yes | exact disk area πr² | h=0.04, 0.02, 0.012 | no | scatter_partial.json | production plasma area relative error 7.8e-5 at h=0.012 (0.166177 versus 0.166190) | area error falling | PASS | — |
| 19b | thin quartz wall | yes | annulus area | h_edge=0.008 gives about 6.25 elements across the 0.05 wall | no | scatter_partial.json | annulus area error about 8e-7 | wall resolved | PASS | — |
| 19c | cluster interfaces | yes | same area check on every cylinder | the cluster meshes | no | fem_scatterers.py | same centroid assignment | area error | PASS | — |
| 19d | geometric error separate from field error | yes | area versus πr², reported beside the field error | three h on the production disk | no | scatter_partial.json | area 7.8e-5 at the mesh where the field error is 0.0015 dB | both numbers stated | PASS | — |

## 20. Quantitative uncertainty budget

| ID | Required benchmark | Actually run? | Independent analytic reference? | FEM refinement? | Meep comparison? | Evidence file | Numerical result | Pass criterion | Status | Missing work |
|---|---|---|---|---|---|---|---|---|---|---|
| 20a | FEM: mesh | yes | T-matrix and Mie | C/F/X and the new h sequence | — | fem_cluster7.json | seven-coated forward last step 0.003 dB and 0.13° | much smaller than 0.10 dB and 1° on that probe | PASS | near-probe third mesh, see 12h |
| 20b | FEM: geometry | yes | area error | h=0.012 | — | scatter_partial.json | plasma area relative 7.8e-5 | stated | PASS | — |
| 20c | FEM: PML | partial | slab only | slab settings | — | canonical_suite.json | ≤ 0.0013 dB on the slab | cylinder sweep missing | UNRESOLVED | 13b, 13c, 13d, 13h |
| 20d | FEM: port integration | yes | guide power | h=0.00125 | — | guide_refined.json | 4.25e-4 relative | ≤ 0.01 dB | PASS | — |
| 20e | FEM: material evaluation | yes | ε(f) at the centroid | the slab and cylinder solves | — | analytic_sweeps.json | off-grid slab at h=0.04 was excluded | faces on the grid | PASS | — |
| 20f | FEM: linear solve | yes | residual | the new solves | — | scatter_partial.json | about 1e-14 | roundoff | PASS | condition estimate |
| 20g | Meep: resolution | partial | dielectric yes; plasma no | — | dielectric res 24, 40 | meep_crosscheck.json | dielectric matches; plasma does not | ≤ 0.05 dB | FAIL | plasma resolution does not reach Mie |
| 20h | Meep: grid registration | prior only | Mie / T-matrix | — | offsets kept separate | prior campaign | spread not collapsed; seven-bulb original grid about 0.6 dB and 15° | spread reported | FAIL | — |
| 20i | Meep: runtime | prior | late DFT | — | long runs | prior campaign | long disk DFT does not return to Mie | late value near analytic | FAIL | — |
| 20j | Meep: dispersive-interface error | prior | analytic | — | plasma versus quartz | prior campaign | quartz-only ports were within about 0.08 dB of FEM; plasma interfaces are not | ≤ 0.10 dB | FAIL | — |
| 20k | Meep: DFT convergence | prior | late-time value | — | long DFT | prior campaign | plasma disk and seven-bulb do not settle on the analytic value | settled and near analytic | FAIL | — |
| 20l | Meep: port extraction | prior straight guide only | FEM port | — | one guide | prior checkpoint | about 6.9e-5 dB on the straight guide | ≤ 0.01 dB | PASS | not repeated here |
| 20m | analytic: Mie truncation | yes | coefficient tail | — | — | analytic_sweeps.json | one-cylinder reduction 6.8e-21 | roundoff | PASS | — |
| 20n | analytic: T-matrix truncation | yes | m_max sweep | — | — | analytic_sweeps.json | forward change below about 3e-5 by m_max=8 | ≪ 0.10 dB | PASS | — |
| 20o | analytic: conditioning | estimated | 2-norm | — | — | analytic_sweeps.json | cond about 2e10 at m_max=8 for coated seven, while the probe is stable | moderate cond | UNRESOLVED | scaled formulation |
| 20p | analytic: floating point | yes | Mie reduction and Graf test | — | — | analytic_self_checks.json | Mie reduction 6.8e-21; Graf worst about 1e-11 at nmax=12 | roundoff on the identity | PASS | — |

## Count of audit rows

Counted from the status cell of every data row above.

| Status | Rows |
|---|---|
| PASS | 110 |
| FAIL | 10 |
| UNRESOLVED | 24 |
| Total rows | 144 |

FAIL: 2f, 3d, 6m, 6n, 12g, 20g, 20h, 20i, 20j, 20k.

Two of those are FEM gradient samples that miss a 0.1% bar (E/H on the homogeneous wave and on the guide at h=0.01). The other eight are Meep plasma comparisons against the analytic field. No completed FEM-versus-analytic field or power row in sections 1, 4, 5, 6a–6l, 7a, 7e, 9b, 10, 11, or 12b–12e is a FAIL.

UNRESOLVED: 3h, 6g, 7f, 7h, 7i, 8d, 9c, 9d, 10c, 11c, 12h, 13b, 13c, 13d, 13h, 14d, 14i, 15a, 15b, 15e, 18c, 18d, 20c, 20o.
