# Circulator parameter bounds

Frozen before production optimization.

- `s ∈ [0.05, 1.8]`
- `f_p^2 = s f_p,ref^2` with `f_p,ref = 8 GHz`
- Analytic `∂ε/∂s` remains finite at `s=0` (depends on `f_p,ref^2`, not `1/√s`)
- In-range samples: no permittivity-determinant singularity; dissipation matrix PSD

s=0 is mathematically ok (fp^2=0) but excluded from the closed interval to avoid optimizer pile-up on a vacuum corner; s_lo=0.05 keeps a weak plasma. s_hi=1.80 stays below the exploratory s=2 map and away from denser regimes not exercised in the validated Jacobian campaign. Production s=1 is interior.

**BOUNDS_AUDIT: PASS**
