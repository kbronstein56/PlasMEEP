# B=0 validation campaign

**FAILED_HORNS**

The corrected horn-only FEM reproduces the trusted Meep column to 0.03 dB on P4 and 0.14 dB on P3/P5, and it is short by 0.24 dB on P2/P6. FEM-H+ (758k DOFs) and FEM-VH (1.31M DOFs) agree to 0.004 dB, so another uniform refinement will not close the gap. Symmetry splitting is 10⁻¹² dB. Quartz and plasma were not added.

Trusted Meep horns-only reference, 25 ppc, `num_mode_guide_normal`, rt=20, prism PEC, no bulbs:

| Port | Normalized power | dB |
|---|---:|---:|
| P2 = P6 | 0.03763160 | −14.244 |
| P3 = P5 | 0.08553915 | −10.678 |
| P4 | 0.69668430 | −1.570 |

## 1. PEC implementation

Validation assembly in `fem_validated_solver.py` drops every triangle whose ρ is zero from both the stiffness and the mass. Empty matrix rows, which are the metal nodes that no longer touch an air element, are pinned to 1. Those nodes are not unknowns of a Helmholtz equation. The air-side condition is the natural Neumann condition ∂Hz/∂n = 0. Dirichlet Hz = 0 is not used.

The horn boundary is the Meep prism polygon. `triangle` builds a constrained triangulation of the upper half-cell with those edges forced in, and the mesh is mirrored. The meshed metal area is 17.8800000005 a² against the polygon shoelace 17.8800000000 a².

A run aborts if the operator log does not show at least 100 excised PEC elements and σ_max ≥ 20. The FEM-VH log is:

- PEC elements: 118190
- PEC interior DOFs pinned: 48396
- PML thickness: 2.0 a
- σ_max: 25.904082296183013
- domain: 30 × 28
- model: `excised_neumann_no_mass`

## 2. PML implementation

Quadratic Meep profile with R = 10⁻¹⁵:

σ_max = −ln(10⁻¹⁵) · 3 / (2 d) = 25.904 at d = 2 a

The stretch is s = 1 + iσ/ω, with σ = σ_max (depth/d)², applied anisotropically in the weak form. The old σ_max = 2 path is gone.

## 3. Mesh symmetry

The P1–P4 axis is a mirror. The upper half is a constrained Delaunay mesh (quality angle 20°) and the lower half is the reflected triangulation, including reflected horn edges. This is one mesh, not an average of seeds.

On FEM-VH the port splits are P2−P6 = 5.8×10⁻¹³ dB and P3−P5 = 4.3×10⁻¹³ dB.

An earlier mirrored but unconstrained Delaunay, with the PEC region chosen by triangle centroids, was also symmetric to 10⁻¹² dB and still sat near −0.23 dB on P2. The residual is not left/right mesh noise.

## 4. Horn-only convergence

Polygon PEC, Meep PML, element-wise Poynting integral on the true monitor line, same numerical-mode source as Meep. Errors are FEM minus Meep in dB.

| Level | DOFs | Triangles | min / median / max edge (a) | Factor (s) | Wall for that level (s) | P2−P6 (dB) | P3−P5 (dB) | P2 | P3 | P4 | P5 | P6 |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| FEM-M | 184968 | 367210 | 5.8e−5 / 0.070 / 0.141 | 1.63 | ~5 | ~0 | ~0 | −0.269 | −0.156 | +0.052 | −0.156 | −0.269 |
| FEM-H | 386679 | 770028 | 1.7e−5 / 0.045 / 0.099 | 6.01 | ~15 | ~0 | ~0 | −0.196 | −0.061 | +0.124 | −0.061 | −0.196 |
| FEM-H+ | 757952 | 1510400 | 2.1e−6 / 0.030 / 0.071 | 13.42 | ~36 | ~0 | ~0 | −0.243 | −0.140 | +0.029 | −0.140 | −0.243 |
| FEM-VH | 1308537 | 2610206 | 5.1e−7 / 0.022 / 0.057 | 33.65 | ~85 | 5.8e−13 | 4.3e−13 | −0.239 | −0.140 | +0.025 | −0.140 | −0.239 |

FEM-VH normalized powers: P2 = P6 = 0.035613, P3 = P5 = 0.082835, P4 = 0.700743.

FEM-H+ and FEM-VH differ by at most 0.004 dB on every port. The minimum edge is a constrained-vertex sliver, not the design spacing. The median edge at FEM-VH is 0.022 a (0.44 mm).

Plots: `b0_campaign/horn_convergence.png`.

Success gate (all of P2–P6 within 0.10 dB, and symmetry within 0.10 dB): not passed. Worst port is 0.239 dB.

## 5. Quartz-only comparison

Not run. The horn gate did not pass.

## 6. Quartz mesh convergence

Not run.

## 7. Full plasma B=0 comparison

Not run.

## 8. Full-device mesh convergence

Not run.

## 9. Second source

Not run.

## 10. Corrected six-port timing

FEM-VH, horns only, one factorization, six right-hand sides, Poynting on all six ports. Measured on this machine:

| Step | Time (s) |
|---|---:|
| Mesh | 3.17 |
| PEC mask | 1.13 |
| Assembly | 7.32 |
| LU factorization | 32.67 |
| First RHS build | 0.004 |
| First solve | 0.489 |
| Each extra solve | 0.472 |
| All six port integrals | 0.006 |
| Total six-port evaluation | 47.2 |

That total is in the **30–60 s** band. It does not include a second factorization of the straight-guide incident reference. Building that reference on the same mesh costs another ~40 s. A candidate that reuses a cached incident power stays in the 30–60 s band.

The old ~11 s figure was the unvalidated operator on a coarser unstructured mesh.

## 11. Final DOF count

1,308,537 (FEM-VH). FEM-H+ at 757,952 DOFs is the same answer to 0.004 dB.

## 12. Largest Meep-vs-FEM error

0.239 dB, on P2 and on P6, FEM-VH.

## 13. RMS dB error

0.176 dB over P2–P6 at FEM-VH.

## 14. Remaining symmetry splitting

Below 10⁻¹² dB for both P2−P6 and P3−P5.

## 15. Is the FEM result converged?

The FEM ports are converged. FEM-H+ and FEM-VH agree to 0.004 dB. They are not converged onto the 25 ppc Meep values. The missing 0.24 dB is a difference between the two models, not an under-resolved mesh.

Checks that did not remove it:

- Meep rt=80 instead of rt=20 moves P2 by −0.008 dB and P4 by +0.001 dB. The reference is steady.
- Meep with subpixel averaging turned off reproduces the rt=20 powers to the printed digits. Averaging is not what the prism is doing.
- An FEM metal boundary set to the res-50 pixel staircase, which is the Yee material test, gives P2 = −0.319 dB, P3 = −0.107 dB, P4 = +0.007 dB. The staircase is farther from Meep on P2 than the exact polygon.
- After one complex scale fixed in the P1 feed, the median |Hz| error of FEM-H+ against the Meep grid is 8% in the center and 14% near the PML. The difference plot is an interference pattern in the open cavity (`b0_campaign/field_diff_FEM-H+.png`).

## 16. Exact next step

Stay on the six-horn, no-bulb problem. Do not add quartz or plasma, and do not turn on B≠0.

The next measurement is the complex Hz along the P1 feed axis and along the two diagonal cuts through the open cavity, FEM-VH against the existing Meep DFT, with the same single complex scale used for the field map. That locates the 0.24 dB as a phase error along a path. It is not another global h-refinement, and it is not a return to ρ→0-with-mass or σ_max = 2.
