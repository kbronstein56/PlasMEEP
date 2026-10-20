# Crossed-refinement plan (fixed design `opt_uniform_C_warm_M`)

No reoptimization. Objective / ports / S definition unchanged.
Frozen M→F bars: `MESH_TRANSFER_CRITERION.md`.

## Numerical error sources (separated)

| knob | meaning |
|------|---------|
| `h_iface` | plasma / quartz circle resolution |
| `h_air` | cavity / bulk air |
| `h_horn` | horn wall / feed region |
| `h_pml` | PML (held with air family when refining F-side) |

## Mesh matrix

### M-centric one-knob (new)

| grade | h_iface | h_horn | h_air | h_pml | role |
|-------|--------:|-------:|------:|------:|------|
| **M** | 0.008 | 0.025 | 0.040 | 0.08 | baseline |
| **Mi** | **0.005** | 0.025 | 0.040 | 0.08 | finer plasma/quartz only |
| **Ma** | 0.008 | 0.025 | **0.030** | 0.08 | finer air only |
| **Mh** | 0.008 | **0.016** | 0.040 | 0.08 | finer horn only |
| **Mia** | **0.005** | 0.025 | **0.030** | 0.08 | iface + air |
| **F** | 0.005 | 0.016 | 0.030 | 0.06 | combined reference |

### F-centric existing (already on disk; used as finer references)

| grade | h_iface | h_horn | h_air | h_pml | role | factorable? |
|-------|--------:|-------:|------:|------:|------|-------------|
| **R** | 0.004 | 0.016 | 0.055 | 0.11 | iface finer than F, air coarser | yes (~4.9M) |
| **A** | 0.005 | 0.016 | 0.022 | 0.06 | F iface/horn, finer air | yes (~5.1M) |
| **Q** | 0.005 | 0.016 | 0.018 | 0.06 | finest air, F iface | yes (~6.1M) |
| **X2** | 0.0042 | 0.016 | 0.022 | 0.06 | crossed iface+A air | yes (~6.2M) |
| X / S | — | — | — | — | SuperLU abort | **skip** |

## Memory discipline

Prior SuperLU abort band ≈ 6.6–7.2M nodes (X, S). Do not factor those.
Estimate LU ≈ 8e−6 GiB/DOF; keep free RAM ≫ estimated LU.

## Metrics per mesh (fixed \(s\))

J, T_des, IL, isolation, RL, TE1, through, other TE0, R, P_out, P_material,
power residual. Plus node count, factor time, RSS / est. LU memory.

## Interpretive questions

1. Does Mi (iface alone) recover most of M→F ΔJ / Δiso / ΔRL?
2. Does Ma (air alone) matter?
3. Is Mh (horn) negligible?
4. Is Mia ≈ F (horn refinement secondary)?
5. Do A/Q/X2 move beyond F (is F itself converged)?
