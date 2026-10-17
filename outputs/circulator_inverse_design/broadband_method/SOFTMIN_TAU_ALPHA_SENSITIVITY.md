# Soft-min τ / α sensitivity (fixed spectra)

No FEM solves. Synthetic J(f) shapes over 11 samples (proxy for 4–6 GHz).

## Meanings

- **τ**: how closely soft-min approximates worst-frequency J. Small τ → near hard min; large τ → near mean.
- **α**: weight on mean vs soft-min. α=0 pure soft-min; α=1 pure mean.

FD differentiation PASS (rel ≈ 5.7e-11) is **not** a reason to freeze values.

## Interpretation

| (τ, α) | softmin−worst (sharp dip) | w on worst | notes |
|--------|--------------------------:|-----------:|-------|
| (0.02, 0.10) | 0.048 | 0.91 | tightest worst-tracking in scored box |
| (0.02, 0.25) | 0.048 | 0.77 | same soft-min, more mean weight |
| (0.05, 0.25) prototype | larger soft gap | lower w_worst | smoother; FD-checked; kept for now |

**Keep prototype τ=0.05, α=0.25** until a real dense J(f) spectrum from a trustworthy
SF design is available. FD PASS alone does not justify the values; this study shows
τ≈0.02 tracks isolated dips more aggressively and may be preferable once feature
widths are measured.

**Meanings (for Jesse):**
- **τ** — how closely soft-min ≈ worst-frequency J (small → hard min).
- **α** — mean vs worst weighting (0 → pure soft-min; 1 → pure mean).

## Candidate table (scored box)

| τ | α | dip soft−worst | dip w_worst | flat eff_n |
|--:|--:|---------------:|------------:|-----------:|
| 0.02 | 0.1 | 0.0480 | 0.909 | 11.00 |
| 0.02 | 0.25 | 0.0480 | 0.773 | 11.00 |

Artifact: `softmin_tau_alpha_sensitivity.json`.
