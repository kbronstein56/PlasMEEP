# Absorption fix validation

Root cause: `material_absorption` used elementwise `Im(ε)`; corrected to
`D=(ε−ε†)/(2i)`. Gate remains **≤5%**. Frozen DtN S package unchanged.

**ABSORPTION_FIX: PASS**

| case | GATE / PASS | max residual | P_mat mean |
|------|-------------|-------------:|-----------:|
| homogeneous_lossy_Voigt_operator | True | — | — |
| uniform_FULL91_C | PASS | 0.006705259086471292 | 0.010449828929032907 |
| optimized_FULL91_C | PASS | 0.004632797125813132 | 0.09117095760522854 |
| one_cylinder | None | — | — |
