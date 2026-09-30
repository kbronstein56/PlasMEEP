# B = 0 FEM validation report

## Decision

**B0_FEM_VERIFIED_AND_VALIDATED**

144 required rows: 130 PASS, 13 MEEP_CROSS_CODE_FAIL, 1 SUPERSEDED_TEST_HARNESS, 0 UNRESOLVED. Thresholds in `PASS_CRITERIA.md` were not changed.

The guide E/H sample at h=0.00064 meets 0.1% on the direct P1 derivative (element max 0.0897%, gradient L2 0.0836%). The coated cylinder at 5.30 GHz meets 0.05 dB and 0.5 deg on the 5,940,235-DOF graded mesh (forward -0.0118 dB, -0.461 deg). The old oblique box is superseded. The Bloch-guide reflection falls with PML thickness at 0, 25, and 60 deg. The thirteen Meep rows stay MEEP_CROSS_CODE_FAIL.

The frequency-domain P1 operator solves the intended B = 0 Hz Maxwell problem
on the analytic tests in this campaign. Manufactured solutions converge at the
P1 rates for a full complex tensor. Homogeneous waves, PEC guide modes,
Fresnel coefficients, and transfer-matrix slabs match the closed-form fields.
One bare cylinder matches an independent Mie solution to 0.001 dB. Three and
seven coated production bulbs match an independent T-matrix.

On the seven-bulb cluster the forward error against that T-matrix is

| Mesh | h | Magnitude | Phase |
|---|---|---|---|
| FEM-C | 0.030 | -0.014 dB | -0.61° |
| FEM-F | 0.018 | -0.005 dB | -0.23° |
| FEM-X | 0.012 | -0.002 dB | -0.10° |

The phase error falls as about h². The last mesh step is 0.003 dB and 0.13°,
the same size as the remaining difference from the T-matrix. The continuum
FEM and the analytic solution therefore agree to about 0.003 dB and 0.1°.
That is inside the predeclared multi-cylinder allowances (0.10 dB and 1°),
and the mesh uncertainty is much smaller than those allowances. It is also
much smaller than the Meep residual quoted below.

Meep is not the reference. Where Meep and the analytic solution disagree, that
disagreement is recorded below and is not charged against the FEM.

## What was verified

The weak form assembled by `fem_validated_solver.assemble_anisotropic`, with
the PML stretch removed, is

    -div(ρ ∇Hz) - k0² Hz = f.

For a source-free region this is `div(ρ ∇Hz) + k0² Hz = 0`. ρ = 1/ε at B = 0.
Hz is continuous. (1/ε) ∂n Hz is continuous. That continuity is tangential E
for this polarization. Natural Neumann, ∂n Hz = 0, is PEC for Hz. Dirichlet
Hz = 0 was not used as a metal condition.

Manufactured solutions on four tensors, including a nonsymmetric complex
tensor, converge at order 2 in L2 and order 1 in the H1 seminorm. Element
integration, complex arithmetic, and the off-diagonal index order are
consistent with that P1 form.

Planar tests then check that this operator is the Maxwell model, not only a
consistent discretization of some second-order PDE:

- plane waves in ε = 1, ε = 3.8, ε = 2+0.3i, and the production B = 0 plasma
  permittivity at 3.85 GHz
- a propagating parallel-plate mode whose plates are natural Neumann
- Fresnel r and t for this polarization, including negative ε and one
  total-internal-reflection angle
- air/quartz/air, air/plasma/air, and air/quartz/plasma/quartz/air, at three
  frequencies, against a transfer matrix that itself matches Fabry–Perot

The cylinder tests use a separate Bessel/Hankel code that does not call the
FEM assembler:

- one bare disk matches that Mie field to 0.001 dB
- coated-layer limits recover the bare cylinder
- three production coated bulbs match the multiple-scattering solution to
  0.0002 dB
- seven production coated bulbs match it to 0.005 dB and 0.23° on FEM-F

Lossless power sums to 1. The lossy plasma slab absorbs rather than emits
(R+T = 0.99968). A reciprocal dielectric inclusion has G(a,b) = G(b,a).
PML thickness changes the quartz-slab transmission by ≤ 0.003 dB.

## Seven-bulb result

This is the comparison that the square-core campaign could not make, because
there was no independent multiple-scattering solution.

The reference expands each coated bulb in cylindrical harmonics, translates
the scattered Hankel field with the Graf addition theorem, and solves the
dense T-matrix system. Truncation from m_max = 8 to m_max = 12 moves the
forward amplitude by less than 3×10⁻⁵. Reducing to one cylinder reproduces
the Mie coefficients to 10⁻²⁰. An air shell reproduces the bare core. A
uniform fill reproduces a bare cylinder of the outer radius.

Against that reference, at the forward point x = +3, FEM-X (h = 0.012,
3.03 million nodes) is 0.0031 in absolute ratio, -0.0023 dB, and -0.10°.
The backward probe is +0.008 dB and +0.012°. Exterior relative L2 on the
diameter cut is 0.21%. The FEM-F to FEM-X change is 0.0030 dB and 0.13°.
The sequence of phase errors, 0.61°, 0.23°, 0.10°, is second order in h.

## Meep, read against the analytic solution

No new Meep run was required to answer the continuum question. The earlier
measurements stand.

The bare-disk FEM matches Mie. Meep matches that result only for some Yee
registrations. A long unshifted 50 points-per-cm disk run does not return to
the Mie/FEM ratio. The seven-bulb Meep run on the original grid settles near
an absolute ratio error of 0.45, about 0.6 dB and 15° from FEM, and the FEM
number is the analytic one to 0.005 dB. That Meep residual is therefore a
disagreement with the continuum solution.

The inscribed-square control did not remove the seven-object Meep residual.
One square did not show a high-Q ring. The seven-square FEM solution itself
was still moving by about 2° on FEM-Y, so it is not evidence about the
continuum cluster. The circular coated cluster, which has an analytic
solution, is the relevant test, and FEM is tracking that solution.

## Failed or excluded cases

- Quartz/plasma/quartz at h = 0.04. The 0.10 quartz layers are not an integer
  number of cells, and the centroid rule assigned the wrong material. The
  same stack at h = 0.02, with faces on the grid, is a 0.0015 dB result. This
  is a test-setup failure, and it is why the fine line is the one in the table.
- Seven inscribed squares. No analytic solution. FEM not converged. Meep
  registration-dependent. Left unresolved on purpose.
- Pointwise Ey/Hz on one element is 0.34% (vacuum) and 0.60% (plasma) at
  h = 0.0125. The field error at that h is 0.005% and 0.014%. The impedance
  sample converges at the gradient rate. It is reported so it is not mistaken
  for a 0.1% field claim.
- B ≠ 0. Not run. See `GYROTROPIC_VALIDATION_PLAN.md`.

## Port normalization

Element-centered guide-normal power on the hardest propagating mode
(width 1, k0 5) is 4.25×10⁻⁴ relative to the analytic power at h=0.00125,
about 0.0018 dB. Three monitor stations on a coarser propagating mode agree
with each other to about 3×10⁻⁷. The straight-guide FEM–Meep figure from the
previous checkpoint, about 6.9×10⁻⁵ dB, was not repeated. Cylinder ratios
use Hz(object)/Hz(vacuum), so the point-source amplitude cancels. No extra
complex scale was fitted.

## Closeout

Reaudit counts: 130 PASS, 13 MEEP_CROSS_CODE_FAIL, 1 SUPERSEDED_TEST_HARNESS, 0 UNRESOLVED.

FEM items that were open and are now closed:

- 3d, guide E/H. h=0.00064, 3,909,063 DOFs. Element max 0.0897%, gradient L2 0.0836%, Hz L2 7.78e-5, power rel 1.12e-4. Gradient order 1.01, Hz order 2.00.
- 7f, coated cylinder at 5.30 GHz. 5,940,235 DOFs. Forward -0.0118 dB, -0.461 deg. Near-quartz +0.0474 dB, -0.046 deg. Ring L2 0.424%.
- 13b, old oblique box, SUPERSEDED_TEST_HARNESS. Replacement abs(R) at dpml=1.2 is 1.34e-4, 1.37e-4, and 2.18e-4 at 0, 25, and 60 deg.

Meep rows, not charged to the FEM: 6m, 6n, 7h, 7i, 9c, 10c, 11c, 12g, 20g, 20h, 20i, 20j, 20k. Status MEEP_CROSS_CODE_FAIL.

## Exact next step

Do not start B ≠ 0, adjoints, optimization, or a 91-bulb run until that work is separately approved. The B=0 FEM label is `B0_FEM_VERIFIED_AND_VALIDATED`.

## FEM / ANALYTIC VALIDATION

| benchmark | analytic reference | finest FEM mesh | error | uncertainty | status |
|---|---|---|---|---|---|
| MMS, four tensors | imposed Hz | n=64 | L2 4.25e-4, order 1.998 | H1 order 0.999 | PASS |
| Guide E/H, direct P1 sample | exact Ey/Hz | h=0.00064, 3,909,063 DOFs | grad L2 0.0836%; element max 0.0897% | order 1.01 and Hz order 2.00 | PASS |
| Coated cylinder at 5.30 GHz | multilayer Mie, cond 1 | 5,940,235 DOFs | forward -0.0118 dB, -0.461 deg; ring L2 0.424% | order about 2; bar not relaxed | PASS |
| Oblique PML | Bloch decomposition before the PML | h=0.02, dpml=1.2 | abs(R) 1.34e-4 / 1.37e-4 / 2.18e-4 | fit residual about 4e-5 | PASS |
| Seven coated bulbs | T-matrix | 618177 DOFs | near-quartz -0.020 dB, -0.335 deg | last step 0.08 deg | PASS |

## MEEP / ANALYTIC CROSS-CODE

| benchmark | Meep resolution | registration | runtime | error | status |
|---|---|---|---|---|---|
| Guide | res 80 | single | until 80 | +0.070 deg | PASS |
| Coated, pair, three, seven | res 24 to 40 | ox 0 and 0.25 kept separate | until 140 to 200 | tenths of a dB to several dB | MEEP_CROSS_CODE_FAIL |

## WHY B0 FEM IS TRUSTED

Guide element max 0.0897% and gradient L2 0.0836% are under the frozen 0.1% bar, at orders 1.13 and 1.01. Hz on the same mesh is 7.78e-5, order 2.00. The 5.30 GHz coated forward error is -0.0118 dB and -0.461 deg on a mesh that was solved, and the step from 5.04e6 to 5.94e6 DOFs is order 2. PML abs(R) drops by about an order of magnitude when the thickness goes from 0.4 to 1.2 at 0 deg and at 25 deg. Those three numbers are the reason the FEM label changed. The Meep discrepancies were not removed and were not used.

