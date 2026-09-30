# B = 0 FEM validation report

## Decision

**B0_FEM_PARTIALLY_VALIDATED**

The earlier `B0_FEM_VERIFIED_AND_VALIDATED` label is withdrawn. The item-by-item
record is `VALIDATION_COMPLETENESS_AUDIT.md`: 144 required rows, 110 PASS,
10 FAIL, 24 UNRESOLVED. A related test was not used to pass a missing one.

The ten FAIL rows are the E/H point samples (0.34% vacuum and 0.60% plasma
at h=0.0125; guide E/H 0.30% to 3.3% at h=0.01) and the Meep plasma
comparisons. Completed FEM field and power comparisons against Mie,
transfer matrices, and the T-matrix are inside the frozen thresholds.
The twenty-four UNRESOLVED rows are listed at the end of this report.

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

## What is still open

FAIL:

- E/H point sample, homogeneous wave, 0.34% (vacuum) and 0.60% (plasma) at h=0.0125
- E/H on the guide at h=0.01, 0.30% to 3.3% depending on the mode
- Meep plasma disk versus Mie: resolution and registration, from the earlier campaign
- Meep seven coated bulbs versus the T-matrix: about 0.6 dB and 15° on the original grid

UNRESOLVED:

- Meep parallel-plate guide, coated cylinder, two cylinders, three cylinders, bare seven cylinders
- FEM frequency sweep of the bare and coated cylinders (analytic sweep exists)
- rotated-pair and coated-pair mesh sequences (one mesh each, already inside 0.10 dB and 1°)
- seven-coated near-probe mesh uncertainty (h=0.02 is inside 1°, last step 0.77°)
- oblique, cylinder, and multi-cylinder PML sweeps
- lossless-cluster contour and a volume absorption integral
- guide-port Sij, horns-only reciprocity, and seven-bulb reciprocity
- a second factorization of one FEM matrix, and an FEM condition estimate
- a scaled T-matrix (the forward probe is stable; the raw 2-norm condition is not)

## Exact next step

Do not start B ≠ 0, adjoints, optimization, or a 91-bulb run. The gyrotropic
file is a plan only. The B = 0 package stays partial until the unresolved
rows are run and the failed rows either meet the frozen thresholds or are
explicitly accepted as known limits.
