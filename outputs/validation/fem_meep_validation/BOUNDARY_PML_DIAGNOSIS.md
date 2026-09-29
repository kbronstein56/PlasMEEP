# Boundary and PML diagnosis

**Status: not yet within 0.1 dB. Quartz, plasma, B≠0, and adjoints were not added.**

The six-horn error is a boundary-formulation error. Two separate pieces of the discrete operator differ from the excised-Neumann / Meep-PML problem. Fixing both, on the Meep cell, at FEM-H, brings the ports from ~1.2 dB down to a residual of **0.21 dB**. That residual is the same size as the mesh’s own P2≠P6 split, so it is not yet a demonstrated 0.1 dB match.

Artifacts: `boundary_pml/`.

---

## 1. Is ρ→0 with the mass term left on equivalent to an excised PEC region?

No.

The assembled form is

\[
\int \rho_{xx}\frac{s_y}{s_x}\partial_x H_z\,\partial_x v
+\rho_{yy}\frac{s_x}{s_y}\partial_y H_z\,\partial_y v
\;-\;k_0^2\int s_x s_y H_z v.
\]

ρ multiplies only the gradient terms (`assemble_anisotropic` in `fem_validated_solver.py`). Setting ρ=0 inside a prism kills the stiffness there and leaves the mass.

On a node whose elements are all metal, the equation is

\[
-k_0^2 M H_z = 0
\]

with \(M\) the consistent mass matrix, including coupling to neighboring interface nodes. That does **not** set \(H_z=0\). It sets the interior value to a mass-weighted continuation of the interface values. Those interior degrees of freedom are active, and they couple back onto the interface through the same mass matrix.

Measured on FEM-M, ρ→0, σ_max=2:

| Region | count | max \|Hz\| |
|---|---:|---:|
| interior metal nodes | 7375 | **1.27** |
| interface nodes | 5807 | 1.62 |
| air, median \|Hz\| | 221125 | 0.55 |

The metal interior is louder than the typical air field. The interface diagonal perturbation from the metal mass is small (median \(1.6\times10^{-4}\) relative to the air diagonal, because the air diagonal is stiffness-dominated). The solution is still different, because the coupling is through the off-diagonal mass block, not that diagonal ratio.

Dropping the metal elements entirely (excised domain, orphan nodes pinned, natural Neumann on the remaining boundary) and the variant that zeroes **both** stiffness and mass inside the metal agree with each other to the printed digits. Both disagree with ρ→0.

Same FEM-M mesh, same weak PML, P1 column versus Meep:

| Port | ρ→0 ΔdB | excised ΔdB |
|---|---:|---:|
| P2 | +1.51 | −1.18 |
| P3 | +0.39 | −0.98 |
| P4 | −0.60 | +0.22 |
| P5 | +0.21 | −0.97 |
| P6 | +1.71 | −1.02 |

The one-horn flare ratio does **not** see this. Inward/monitor power stays within 0.06 dB of Meep for ρ→0 and for excised, at both PML strengths (Meep ratio 0.751). The mass term shows up in the open six-horn cavity, not in a single flare.

## 2. Does excised Neumann PEC improve agreement?

It changes the answer by about 2.7 dB on P2, so the formulations are not interchangeable. By itself, on the weak PML, it does not land on Meep: the side-port error changes sign and stays near 1 dB. Excised plus a Meep-strength PML is the combination that moves onto Meep. See question 11.

## 3. What PML is FEM using?

Complex-coordinate stretching, written as an anisotropic transformation medium for the scalar Helmholtz equation. It is not a scalar complex ε, not a lossy sponge added as Im(ε) only, and not a Robin boundary.

\[
s_x = 1 + i\sigma_x/\omega,\qquad s_y = 1 + i\sigma_y/\omega,\qquad \omega = k_0 = 2\pi f_s.
\]

\[
\sigma = \sigma_{\max} u^2,\qquad u = \text{depth}/\text{thickness},\qquad \sigma_{\max} = 2.
\]

Thickness is `dpml_ports` = 2 a (40 mm) on every side of the box \([0,30]\times[0,28]\). The stretch multiplies the gradient terms by \(s_y/s_x\) and \(s_x/s_y\), and the mass term by \(s_x s_y\). That is the correct stretched operator for \(\nabla\cdot(\rho\nabla H_z)+k_0^2 H_z\). The profile amplitude is not Meep’s.

## 4. Is it equivalent to Meep’s PML?

The operator family is the same. The conductivity is not.

Meep uses `mp.PML(thickness=dpml_ports)` with the defaults \(R_{\text{asymptotic}}=10^{-15}\) and \(f(u)=u^2\). The conductivity that produces that asymptotic reflection is

\[
\sigma(u) = \frac{-\ln R}{2\,d\int_0^1 u^2\,du}\,u^2
= \frac{-\ln(10^{-15})\cdot 3}{2d}\,u^2.
\]

At \(d=2\), \(\sigma_{\max} = 25.90\). The FEM value is 2.0, about 13 times smaller.

Meep cell, from `sixport_common`: 30 × 28 a, PML thickness 2 a on all four sides, about 0.5 a of air between the outermost horn vertex and the inner face of the PML. P1’s feed ends near \(x=12.2\); the inner PML face is at \(x=13\). The source sits close to the absorber. FEM uses that same box.

## 5–6. Measured reflection

A straight PEC guide, walls continued into the PML so the aperture itself is not the reflector. Two-wave fit \(A e^{i\beta s}+B e^{-i\beta s}\) on the centerline. Power reflection \(20\log_{10}|B/A|\).

| Case | FEM σ_max=2 | FEM Meep σ | Meep res 16 |
|---|---:|---:|---:|
| normal | **−11.8 dB** | −18.3 dB | −21.3 dB |
| 60° to the PML normal | **−12.7 dB** | −29.4 dB | −20.7 dB |

Meep at this coarse resolution is near −21 dB, not the continuum \(10^{-15}\). FEM with σ_max=2 is about 9 dB more reflective at normal incidence and 8 dB more reflective at 60°. Putting Meep’s σ into the FEM stretch removes that gap on this test. A point-source cylindrical test was not added; the guide results and the six-horn domain sweep are the measurements.

## 7. Sensitivity to outer-domain distance

Excised PEC, σ_max=2, h=0.08, horns fixed. ΔdB versus Meep:

| Air added each side | P2 | P3 | P4 | P5 | P6 |
|---|---:|---:|---:|---:|---:|
| 0 (Meep cell) | −1.24 | −1.09 | −0.37 | −0.90 | −1.25 |
| +25% of the 0.5a gap | −0.96 | −0.64 | −0.25 | −0.83 | −0.78 |
| +50% of the 0.5a gap | −0.80 | −0.55 | −0.53 | −1.05 | −0.59 |
| +2 a | −0.52 | −1.13 | −0.46 | −1.10 | −0.78 |
| +4 a | −0.34 | −0.80 | −0.53 | −0.72 | −0.35 |

P2 moves by 0.9 dB as the boundary recedes. With σ_max=2 the scattering is contaminated by the outer boundary. These runs also change the triangulation, so part of the jitter is mesh noise (question 9). The trend in P2 is larger than that noise.

## 8. Sensitivity to PML thickness

Once σ is the Meep profile, thickening the layer from 2 a to 4 a, on a box with 2 a of extra vacuum, changes every port by less than 0.01 dB (h=0.08, excised). With the weak σ_max=2, thickness still matters: the same box at thickness 4 a and σ_max=2 is already within a few tenths of a dB, because the integral of σ grows with thickness. The production choice σ_max=2 at thickness 2 a is the sensitive one.

## 9. How much asymmetry is the Delaunay mesh?

Four jitters of an h=0.08 mesh, excised, σ_max=2, Meep cell. P2−P6 in dB: **+0.02, −0.24, −0.43, +0.26**. P3−P5: **+0.02, −0.13, +0.20, +0.12**.

Seed-to-seed scatter is a few tenths of a dB. That covers the old FEM-H splits (0.09 dB and 0.23 dB) and the residual split after the formulation fix. It does not cover the original 1.2 dB offset of ρ→0 / σ_max=2 relative to Meep. On that same seed set the P2 error versus Meep stayed between −1.19 and −1.62 dB.

## 10. Where the field difference appears

`boundary_pml/field_diff_rho0_sigma2.png`. One complex scale from the P1 feed. Non-PEC samples. Median relative \|Hz\| error versus distance from the array center, current ρ→0 / σ_max=2:

| Radius (a) | median relative error |
|---|---:|
| 0–3 | 0.11 |
| 3–6 | 0.21 |
| 6–9 | 0.31 |
| 9–11 | 0.37 |
| 11–13 (against the PML inner face at 12–13) | 0.33 |

The feeds stay close, as before. The open cavity is worse, and the error grows toward the PML rather than starting on the angled walls. The relative-error color panel in that figure is not usable: nulls of Hz blow the scale up. The absolute-difference panel and the radial medians are the evidence.

With excised metal, Meep σ, and 2 a of extra vacuum, the same medians fall to 0.08, 0.09, 0.12, 0.13, 0.13. The outward growth flattens.

## 11. Six-horn P1 column

Meep, horns only, res 50: P2=0.037632, P3=0.085539, P4=0.696684, P5=0.085539, P6=0.037632.

**Current production-style operator** (FEM-M, ρ→0, σ_max=2), ΔdB: P2 +1.51, P3 +0.39, P4 −0.60, P5 +0.21, P6 +1.71.

**Same Meep cell, FEM-H, excised Neumann, Meep σ.** This is the formulation correction, not a retuned amplitude.

| Port | Meep power | FEM-H power | Meep dB | FEM-H dB | Δ dB |
|---|---:|---:|---:|---:|---:|
| P2 | 0.037632 | 0.036720 | −14.244 | −14.351 | **−0.106** |
| P3 | 0.085539 | 0.086172 | −10.678 | −10.646 | **+0.032** |
| P4 | 0.696684 | 0.731307 | −1.570 | −1.359 | **+0.211** |
| P5 | 0.085539 | 0.081493 | −10.678 | −10.889 | **−0.210** |
| P6 | 0.037632 | 0.038325 | −14.244 | −14.165 | **+0.079** |

P2−P6 = −0.19 dB and P3−P5 = +0.24 dB on this mesh. The largest port error equals that split. The 0.1 dB criterion is not met. FEM-M with the same two fixes leaves P4 at +0.56 dB, so the finer mesh is part of the improvement and is not finished.

## 12. Dominant cause

**PEC implementation and PML strength, together.** Mesh asymmetry is real and is now the size of the leftover error. It was not the original 1 dB.

- ρ→0 with a live mass term is a different operator from excised Neumann. Interior Hz is large. Port powers move by ~2 dB.
- σ_max=2 reflects about 9 dB more than Meep’s PML on the guide test, and six-horn powers move when the outer wall moves.
- Applying only the stronger PML on top of ρ→0 makes the side ports worse (FEM-M P2 goes to +2.36 dB). Both changes are required.
- After both, on FEM-H, the leftover 0.21 dB matches the Delaunay left/right split.

The source and the horn cross-sections are not the dominant term. The one-horn flare ratio stays inside 0.06 dB across these operators.

## 13. Code change required next

In `scripts/validation/fem_meep_validation/fem_validated_solver.py` only, as a validation solver, not a production PlasMEEP default:

1. `assemble_anisotropic`: do not accumulate the \(k_0^2 s_x s_y\) mass on elements whose ρ is 0. Pin any row that then sums to zero. That is the excised Neumann condition. Leaving the mass on is what fills the prisms with Hz.
2. `pml_sx_sy`: replace `2.0 * (depth/dp)**2` with `(-log(1e-15) * 3 / (2 * dp)) * (depth/dp)**2`. Same quadratic shape Meep uses, same \(R_{\text{asymptotic}}\), same thickness and box.

Then rerun FEM-H on the Meep cell (30 × 28, PML 2 a), excised, that σ, and require the port errors and the P2−P6 / P3−P5 splits to sit under 0.1 dB. If the splits do and the Meep offsets do not, the leftover is still the absorber or the mesh grading, and the domain-distance test should be repeated on that operator. No quartz until that table exists.
