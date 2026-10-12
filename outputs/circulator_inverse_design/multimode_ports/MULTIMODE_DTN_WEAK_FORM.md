# Modal DtN / Robin port weak form

Time convention: **e^{-iωt}**. Starting from the validated Hz weak form after integration by parts:

    ∫_Ω ∇H_z·∇v − k₀² H_z v  +  ∫_{∂Ω} (∂_n H_z) v  = 0

## First-order absorbing (outgoing) condition

For a single propagating mode m with outward normal **n̂** and real β_m > 0, an outgoing wave satisfies:

    ∂_n H_z − i β_m H_z = 0

Moving to the left-hand side adds the boundary Galerkin term:

    + ∫_Γ i β_m H_z v ds

**Frozen implementation sign:** the code uses **−i β_m** on the stiffness (Robin) block:

    + ∫_Γ (−i β_m) H_z v ds

This matches the straight-guide validation (TE0/TE1 transmission |S| ≈ 1, |R| ≲ 0.006) with the `fem_validated_solver` / local unity-stretch bulk operator.

Multi-mode port: sum over retained modes (propagating + optional evanescent):

    + Σ_m ∫_Γ (−i β_m) (H_z · φ_m)(v · φ_m) ds   (Galerkin rank-one blocks on port DOFs)

## Prescribed incident mode

Prescribe a single mode m with **field** amplitude a_{m,f} and profile φ_m (L²-normalized on Γ).

At the **left** end (x = 0, **n̂** = −x̂), incident wave travels **+x** into the domain: H_inc = a φ_m e^{+iβx}.

The inhomogeneous boundary condition on the total field is equivalent to a weak-form RHS:

    b_p += ∫_Γ 2 i β_m a_{m,f} φ_m N_p ds

At the **right** end (x = L, **n̂** = +x̂), incident travels **−x** into the domain with the same RHS coefficient pattern when the outward normal is used consistently.

**Unit power drive:** choose a_{m,f} = 1 / √(β_m/(2ω)) so |â_m|² = 1.

## Relation to extraction

Modal amplitudes are recovered with the same (H_z, ∂_n H_z) projection as `modal_s_core.extract_ab`, using φ_m instead of the flat fundamental profile.

## Full91 placement

Physical device unchanged. The DtN boundary is placed at the **outer end of the straight feed** (farther from the flare than the legacy monitor), optionally after a modest **numerical lead extension** documented separately from throat width.

PML at domain edges may remain for non-port radiation; port modes are **not** terminated by thin feed PML alone.
