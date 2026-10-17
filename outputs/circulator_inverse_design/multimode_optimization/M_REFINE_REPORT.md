# M-mesh C6 refinement (warm-start from opt_uniform_C)

Frozen objective / weights / bounds / ports / projected measure / S definition — unchanged.
Efficient adjoint: 1 LU + 6 TE0 forward + 6 combined adjoint.

- Start: transferred C optimum on M → J=0.4875, T_des=0.6946, TE1=0.0224
- Best after maxiter=8: J=0.7229, T_des=0.8187, TE1=0.0265
- IL=0.87 dB, iso=20.30 dB, RL=15.97 dB
- through=0.0010, other=0.0462, R=0.0253
- elapsed=1.47 h, n_fev=17, nit=8

| metric | C opt (on C) | C→M transfer | M refined |
|--------|-------------:|-------------:|----------:|
| J | 0.708 | 0.487 | 0.723 |
| T_des | 0.804 | 0.695 | 0.819 |
| IL dB | 0.95 | 1.58 | 0.87 |
| iso dB | 19.0 | 15.35 | 20.30 |
| RL dB | 16.6 | 11.48 | 15.97 |
| TE1 | 0.019 | 0.0224 | 0.0265 |
| through | 0.002 | 0.0084 | 0.0010 |
| other | 0.049 | 0.1137 | 0.0462 |
| R | 0.022 | 0.0712 | 0.0253 |

**M recovered** from transferred J=0.487 → 0.723 (above C-mesh J=0.708).
Artifact: `opt_uniform_C_warm_M.json`.
