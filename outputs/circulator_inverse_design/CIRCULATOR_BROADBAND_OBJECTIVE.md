# Broadband objective (frozen before the broadband run)

Warm-start from the best single-frequency sixfold design.

## Frequency grid

    f ∈ {3.75, 3.80, 3.85, 3.90, 3.95} GHz

Same `B = +0.05 T`, same six sources, same C6 tying, same bounds
`s ∈ [0.05, 1.80]`, same power weights.

## Normalization

Reuse the **same frozen** `P_ref` from the single-frequency uniform design
at 3.85 GHz (not recomputed per frequency). This keeps the scale fixed for
comparison across the band.

## Objective

At each frequency `f_n` evaluate the single-frequency objective `J(f_n)`.

    J_bb = mean_n J(f_n) + 0.5 * min_n J(f_n)

The min term is a robust penalty so one excellent point cannot mask a
disastrous band edge. The subgradient of `min` uses the active (worst)
frequency.

Gyrotropic tensors are recomputed at every frequency. One factorization
per frequency; six forward + six adjoint RHS per frequency.
