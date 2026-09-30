# Gyrotropic validation plan

This is a plan for B ≠ 0. It is not a result, and it is not an authorization
to run the production magnetized device or any adjoint.

The B = 0 package is `B0_FEM_PARTIALLY_VALIDATED`. The earlier verified label was withdrawn. Do not start this ladder until that package is accepted. Each rung below uses an analytic or algebraic reference first. Meep is a later comparison, after that reference exists.

## Ladder

1. Exact B → 0 recovery. At zero bias the tensor ε must reduce to the scalar
   already validated here, including the production fp, γ, and fs. The
   reduction is a sample-point identity, not a solve.

2. Sign relation. ε(+B) and ε(-B) must satisfy the expected transpose
   relation for the Hall terms. State the relation from the cold-plasma
   dielectric tensor and check it at sample points before any solve.

3. ρ = ε⁻¹. Invert the 2×2 tensor analytically at the same sample points and
   compare with the code that builds ρ. This is the quantity the weak form
   actually multiplies into the stiffness.

4. Homogeneous gyrotropic eigenmode. In a uniform magnetized plasma, Maxwell’s
   equations have a known plane-wave dispersion for propagation perpendicular
   and parallel to B. Impose that mode as a manufactured or Dirichlet field
   and require the same P1 rates already seen at B = 0. This is the first
   solve, and it does not need a port or a bulb.

5. Propagation and phase. Where the dispersion relation gives a real or
   complex β, compare phase per meter and the E/H ratio, as in the B = 0
   homogeneous test.

6. Faraday rotation. A slab or a long uniform region, with the analytic
   rotation angle as the reference. Sweep frequency, |B|, the sign of B, and
   the collision rate. +B and -B must rotate in opposite directions by the
   same amount when the medium is otherwise unchanged.

7. Planar magnetoplasma interface and slab, if the Fresnel or transfer-matrix
   algebra closes for this tensor. If it does not close cleanly, say so and
   do not replace it with a Meep number.

8. Energy and passivity. For the collisionless Hermitian tensor, check the
   power identity appropriate to that tensor. With collisions, absorption
   must not come out negative.

9. Onsager–Casimir. For the same geometry,
   S_ij(+B) = S_ji(-B).
   This is the magnetized replacement for the B = 0 reciprocity test already
   passed. Do it on a small scatterer before any horn or array.

10. One small nonreciprocal scatterer. A single cylinder or slab in a bias
    field, with the analytic mode or the Onsager check as the reference.
    Refine the mesh until the +B/-B residual is smaller than the effect
    being claimed.

11. Meep, only after those analytic checks. A disagreement is a finding, not
    a reason to edit the FEM tensor.

12. Full-device B reversal last. It is not a validation of the tensor. It is
    a device-level confirmation after the tensor, the interface, and a single
    scatterer have already passed.

## Out of scope until that ladder passes

Adjoints, shape optimization, Bayesian searches, 91-bulb magnetized runs,
and any change to the production bias path.
