# B = 0 FEM validation report

## Decision

**B0_FEM_PARTIALLY_VALIDATED**

The earlier `B0_FEM_VERIFIED_AND_VALIDATED` label stays withdrawn. The
item-by-item record is `VALIDATION_COMPLETENESS_AUDIT.md`, recomputed by
`reaudit_closeout.py`: 144 required rows, 128 PASS, 16 FAIL, 0 UNRESOLVED.

The homogeneous E/H sample now meets 0.1% (vacuum 0.084% at h=0.003125,
plasma 0.077% at h=0.0015625). The guide gradient does not: at h=0.0008
the element median is 0.0987%, the element max is 0.115%, and the gradient
L2 is 0.105%. The coated cylinder at 5.30 GHz is −0.102 dB and −4.73° at
h=0.012. The oblique PML box does not isolate reflection. Those three are
the FEM failures. The other thirteen FAIL rows are Meep against the analytic
continuum. Details are in `B0_VALIDATION_CLOSEOUT.md`.

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

Reaudit counts: 128 PASS, 16 FAIL, 0 UNRESOLVED.

FEM FAIL rows:

- 3d, guide E/H. h=0.0008, element max 0.115%, gradient L2 0.105%. Class C.
- 7f, coated cylinder at 5.30 GHz. Forward −0.102 dB, −4.73° at h=0.012. Class C.
- 13b, oblique PML box. rel L2 about 2.4%, insensitive to the PML knobs. Class F.

Meep FAIL rows, not charged to the FEM: 6m, 6n, 7h, 7i, 9c, 10c, 11c, 12g, 20g, 20h, 20i, 20j, 20k.

Bare-cylinder frequency sweep, coated frequencies other than 5.30 GHz, pair orientations, the locally refined seven-bulb probes, cylinder and pair PML sweeps, lossless contour power, the absorption identity, Green reciprocity, and the repeated factorizations pass.

Solve counts. The pre-closeout package stored 275 analytic evaluations, 187 FEM solves, and 8 Meep simulations. This closeout added 77 recorded FEM factorizations, 45 T-matrix condition solves plus the cluster solves paired with the new scatter meshes (about 92 analytic evaluations), and 36 Meep simulations (2 guide resolutions and 17 scatter cases, each with a vacuum run and an object run). Totals for the B=0 package: about 367 analytic evaluations, 264 FEM solves, and 44 Meep simulations. The earlier plasma-ring campaign is additional evidence for the Meep FAIL rows and is not in the 44.

## Exact next step

Do not start B ≠ 0, adjoints, optimization, or a 91-bulb run. The label stays
`B0_FEM_PARTIALLY_VALIDATED` until 3d and 7f meet the frozen thresholds on a
finer mesh, and until the oblique PML test is replaced by a measurement that
actually tracks the PML. Meep plasma agreement is a separate cross-code item.
