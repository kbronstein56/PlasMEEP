# Overnight B=0 validation

**Status: QUARTZ_VALIDATED**

Date: 2026-09-27. No commit and no push.

The horn-only continuum comparison is unchanged: Meep moves toward the mesh-converged FEM horn solution, and a 0.1–0.2 dB finite-resolution gap is not treated as an FEM failure. Quartz-only follows that same pattern and passes. The full 91-bulb B=0 plasma problem does not. Higher-resolution Meep moves away from FEM on the through port and on P2/P6. B≠0 was not started.

Plots: `overnight_b0/quartz_vs_dx.png`, `overnight_b0/plasma_vs_dx.png`.

## 1. Horn-continuum conclusion

The converged FEM-VH horn-only solution is the reference. FEM-H+ and FEM-VH differ by at most about 0.004 dB. Meep minus FEM-VH, in dB:

| Port | 25 ppc | 35 ppc | 50 ppc |
|---|---:|---:|---:|
| P2/P6 | +0.239 | +0.192 | +0.136 |
| P3/P5 | +0.140 | +0.131 | +0.083 |
| P4 | −0.025 | −0.023 | −0.038 |

Meep moves toward FEM on P2 and P3. P4 stays inside 0.04 dB. A linear-in-dx extrapolation of the P2 error intercepts near +0.04 dB; a linear-in-dx² extrapolation intercepts near +0.11 dB. Those two models were left as an unresolved 0.04–0.11 dB band. A 60 ppc horn run would not have separated them, and it was not run. Free-space numerical dispersion is not the explanation: the FEM and Meep free-space wavenumbers agree to a few parts in 10⁴. FEM was not tuned to finite-resolution Meep.

## 2. Quartz-only FEM convergence

The operator is the canonical one, and each solve asserted it:

- exact shared Meep prism polygons, metal area 17.880 a²
- PEC elements dropped from stiffness and mass (`pec_model = excised_neumann_no_mass`)
- Meep quadratic PML, σ_max = −ln(1e−15)·3/(2d) = 25.904 at d = 2 a
- mirrored constrained triangulation (P2=P6 and P3=P5 to the printed digits)
- guide-normal Poynting flux, numerical P1 mode, no scale factor

Quartz shells use the production centers, OD 15 mm, ID 13 mm, ε = 3.8. Interiors are air. There is no plasma susceptibility.

Wall thickness: 0.050 a = 1.00 mm.

| Mesh | DOFs | Median quartz edge | Min. triangles on a radial ray through the wall |
|---|---:|---:|---:|
| FEM-H | 444,775 | 0.0159 a (0.32 mm) | 7 |
| FEM-VH | 890,959 | 0.0100 a (0.20 mm) | 9 |
| FEM-QF | 1,107,918 | 0.0060 a (0.12 mm) | 17 |
| FEM-VB | 1,417,429 | 0.0100 a (0.20 mm) | 9 |

The ray count includes triangles the sample line clips, so it is higher than a pure radial stack. The median edge is the resolution of the wall: FEM-QF puts a 0.12 mm edge in a 1.00 mm wall. Refining only the wall from FEM-VH to FEM-QF moves the ports by at most 0.006 dB. Refining the air/horn mesh from FEM-VH to FEM-VB moves them by at most 0.012 dB. FEM-VB is the converged quartz solution.

Normalized FEM-VB powers, source P1: P2 = P6 = 0.04276, P3 = P5 = 0.02483, P4 = 0.70457. The source-monitor ratio P11 = −0.893 is the unsubtracted inward flux, not a reflection coefficient.

## 3. Quartz Meep resolution series

Matched runs, source P1, receivers P1–P6, `geometry_only` (wp = 0), ε averaging on, np = 32, df = 0.10 fs, run_time = 20. One MPI job at a time.

| ppc | dx | Wall time | P11 | P12 | P13 | P14 | P15 | P16 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | 0.400 mm | 155 s | −0.8744 | 0.04785 | 0.02658 | 0.6610 | 0.02658 | 0.04785 |
| 35 | 0.286 mm | 782 s | −0.8822 | 0.04565 | 0.02558 | 0.6787 | 0.02558 | 0.04565 |
| 50 | 0.200 mm | 2631 s | −0.8878 | 0.04313 | 0.02441 | 0.6945 | 0.02441 | 0.04313 |

## 4. Quartz port errors and trends

Meep minus FEM-VB, in dB. Positive means Meep is higher.

| Port | 25 ppc | 35 ppc | 50 ppc | Trend |
|---|---:|---:|---:|---|
| P2/P6 | +0.489 | +0.284 | +0.038 | toward FEM |
| P3/P5 | +0.296 | +0.130 | −0.074 | toward FEM, then 0.07 dB past it |
| P4 | −0.277 | −0.162 | −0.063 | toward FEM |

Meep50 and the converged FEM agree to ≤ 0.08 dB on P2, P3, and P4. That passes the quartz gate (threshold about 0.15 dB). The quartz shells do not open a new discrepancy.

Linear-in-dx and linear-in-dx² fits of these three curved points disagree by several tenths of a dB (P2 intercepts −0.39 dB versus −0.07 dB). The curvature is strong, so those intercepts are not used as the continuum claim. The measured Meep50−FEM residual, 0.04–0.07 dB, is the quartz continuum estimate.

## 5. Full B=0 FEM convergence

Plasma is the uniform fp = 8 GHz, γ = 1 MHz Drude response at fs = 3.85 GHz, inside r = 4.6/6.5 × ID. Printed before the solve, and checked against `mp.Medium.epsilon(fs)`:

- ε_FEM(fs) = ε_Meep(fs) = −3.31775987 + 0.00112150 i
- complex difference ≈ 0 (imaginary part 9×10⁻¹⁹)

The plasma circle is a constrained mesh edge. On that mesh the plasma area is 15.119 a² against the analytic 15.123 a², and the quartz area matches the analytic annulus. A centroid staircase without that edge was asymmetric by about 1 dB and disagreed with the constrained mesh by several dB; those staircased numbers are in `fem_plasma_staircased.json` and are not the comparison below.

Constrained-circle FEM, normalized, source P1:

| Mesh | DOFs | P2=P6 | P3=P5 | P4 |
|---|---:|---:|---:|---:|
| FEM-VH | 916,818 | 0.13022 | 0.06934 | 0.26590 |
| FEM-VB | 1,447,959 | 0.12583 | 0.06295 | 0.30406 |
| FEM-VC | 2,114,448 | 0.12405 | 0.05938 | 0.32993 |

Symmetry is exact after the circle constraint. The levels are not converged. VB minus VH is −0.15 dB (P2), −0.42 dB (P3), +0.58 dB (P4). VC minus VB is −0.06 dB (P2), −0.25 dB (P3), +0.36 dB (P4). P4 is still rising. FEM-VC is the finest value used below, and it is not a converged continuum number.

A flat plasma slab of thickness 0.46 a, same ε, transmits within 0.09 dB of the analytic |T|² in this weak form (`plasma_slab.json`). The isotropic negative-ε operator is right on a flat interface. The device discrepancy below is not a missing factor of 1/ε.

Adding the plasma disks changes the two methods in opposite directions. From quartz-only to plasma, Meep50 P2 falls from 0.043 to 0.002 and P4 falls from 0.694 to 0.031. FEM P2 rises from 0.043 to 0.124 while P4 falls from 0.705 to 0.330.

## 6. Full B=0 Meep 25 / 35 / 50

Same formulation as the trusted column: a = 0.020, `num_mode_guide_normal`, df = 0.10 fs, prism horns, offsets 0, rotation 0, run_time = 20, uniform fp = 8 GHz, source P1, np = 32. An aborted first launch had passed ρ = 1 and therefore wp = 1 instead of fp_a. It was killed before it wrote a result. The runs below assert Drude frequency = fp_a = 0.533703 and 91 susceptibilities. P1→P2 reproduces the existing trusted files to the printed digits (25 ppc 0.01372968, 35 ppc 0.00552307, 50 ppc 0.00213008).

| | FEM-VC | Meep 25 (371 s) | Meep 35 (1298 s) | Meep 50 (4339 s) |
|---|---:|---:|---:|---:|
| P2=P6 | 0.12405 | 0.013730 | 0.005523 | 0.002130 |
| P3=P5 | 0.05938 | 0.002207 | 0.008432 | 0.002142 |
| P4 | 0.32993 | 0.318804 | 0.184630 | 0.031266 |

Meep minus FEM-VC, dB:

| Port | 25 ppc | 35 ppc | 50 ppc | Trend |
|---|---:|---:|---:|---|
| P2/P6 | −9.56 | −13.51 | −17.65 | away from FEM |
| P3/P5 | −14.30 | −8.48 | −14.43 | nonmonotonic |
| P4 | −0.15 | −2.52 | −10.23 | away from FEM |

The 25 ppc P4 agreement (0.15 dB) does not survive refinement. An existing 60 ppc file, P1 and P2 only, has P1→P2 = 0.00343, still in the small-Meep cluster. A new 60 ppc six-port run would not turn this into a shared limit. It was not run.

## 7. Continuum extrapolations

Quartz: see section 4. The measured 50 ppc residual is the estimate. The two polynomial intercepts disagree and are not adopted.

Plasma: the three Meep points are not a smooth march onto FEM. Fitting the error anyway gives large negative intercepts, and the two models disagree with each other by several dB (P4: −18.7 dB linear in dx, −11.5 dB linear in dx², rms 1.6–2.0 dB). P3’s fit rms is 2.7 dB because the sequence is nonmonotonic. Those intercepts are reported so the disagreement is visible. They are not a continuum prediction.

## 8. Second source

Not run. The P1 plasma comparison does not support an orientation check as a validation step. A P2 numerical mode exists at res 50, 70, and 100 for a later run.

## 9. Validated FEM six-port runtime

The mesh that survived validation is quartz FEM-VB, 1,417,429 DOFs. Dedicated timing after the Meep jobs finished (`quartz_vb_timing.json`):

| Stage | Time |
|---|---:|
| Mesh | 4.0 s |
| Material | 5.2 s |
| Assembly | 7.7 s |
| Factorization | 36.8 s |
| P1 solve | 0.48 s |
| Five reused LU solves | 2.58 s |
| Six-port flux integration | 0.006 s |
| Complete | 56.7 s |

Peak RSS 5.5 GB. The five extra solves measure LU backsolve cost. P3–P6 were not separate physical excitations; their backsolves reused the P2 load. A factorization of this same mesh while Meep was running earlier took 106 s. The 37 s figure is the unloaded measurement.

## 10. Final DOF count

Validated quartz mesh: **1,417,429** DOFs (FEM-VB). Finest plasma mesh, not validated: 2,114,448 DOFs (FEM-VC).

## 11. Largest finite-resolution Meep−FEM error

Quartz, against converged FEM-VB: 0.49 dB at 25 ppc on P2, shrinking to 0.074 dB at 50 ppc on P3.

Plasma, against finest FEM-VC: **17.7 dB** at 50 ppc on P2/P6. P4 at 50 ppc is 10.2 dB. Both are larger than the 25 ppc gaps.

## 12. Best estimate of continuum mismatch

Horns: 0.04–0.11 dB, unresolved between the two extrapolations, with Meep still approaching FEM at 50 ppc.

Quartz: 0.04–0.07 dB, measured at Meep 50 ppc, with FEM internally converged to about 0.01 dB.

Plasma: no shared-continuum estimate. At 50 ppc the methods differ by 10–18 dB, and Meep’s P2 and P4 are still moving away from FEM.

## 13. Measured versus extrapolated

Measured: horn Meep 25/35/50, quartz FEM four meshes, quartz Meep 25/35/50, plasma FEM three meshes, plasma Meep 25/35/50, plasma area, Drude ε(fs), and the slab |T|².

Extrapolated and not adopted: horn dx versus dx² intercepts (0.04 versus 0.11 dB), and the plasma error-intercept fits (rms > 1 dB). Quartz polynomial intercepts were computed and rejected because the three points are curved and Meep50 is already on top of FEM.

## 14. Can B=0 be considered validated?

The horn-only and quartz-only continuum problems can. Quartz Meep50 matches converged FEM to ≤ 0.08 dB, and Meep moved toward FEM as the grid was refined.

The full 91-bulb B=0 device cannot. FEM is not internally converged, and Meep’s resolution series moves away from the finest FEM result. The 25 ppc P4 agreement was a finite-resolution coincidence.

## 15. Next step for B≠0

Do not start a magnetized campaign, and do not treat circulation as validation.

The next measurement is one plasma disk, then one quartz shell with that disk, at two Meep resolutions, against this FEM operator. The flat-slab transmission already matches analytic |T|² to 0.09 dB, and the empty quartz shells already match Meep. The failure appears when the overdense disks are placed in the six-horn cavity: Meep’s side and through ports collapse with resolution, and FEM’s side ports rise. That single-disk scattering test is the isolation step. Gyrotropic tensor algebra can wait until that B=0 disk agrees.
