# Robustness and 91-variable refinement

## Robustness (sixfold peak design on M)

| perturbation | J | t_des | iso (dB) |
|--------------|--:|------:|---------:|
| nominal | 5.380 | 4.358 | 15.83 |
| s×0.95 | −2.15 | 0.51 | −3.14 |
| s×1.05 | 2.98 | 2.49 | 11.49 |
| f=3.80 GHz | 1.76 | 1.60 | 17.76 |
| f=3.90 GHz | −0.44 | 0.23 | −0.44 |
| B=0.045 T | 0.15 | 0.79 | 13.75 |
| B=0.055 T | −0.03 | 1.54 | 11.28 |
| noise 2% (3 trials) | 4.8–5.1 | 4.05–4.18 | 14.1–15.4 |

Characterization only — objective not retuned.

## Free-91 refinement

Warm start = sixfold peak. After 8 L-BFGS iterations on M:

| | sixfold | free-91 |
|--|--------:|--------:|
| J | 5.380 | 5.476 |
| t_des avg | 4.358 | 4.481 |
| t_des min | 4.356 | 4.478 |
| desired std (6 sources) | ~0 | 0.003 |

Added freedom gives ~2% objective gain without breaking six-port
equivalence. Prefer sixfold for simplicity unless the gain matters.
