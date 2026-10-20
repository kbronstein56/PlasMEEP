# MESH_CONVERGENCE_MREFINED_DESIGN

Fixed design throughout: `opt_uniform_C_warm_M.json`  
Frequency 3.85 GHz, B = +0.05 T, projected TE0+TE1 direct-DtN.  
**No reoptimization.** Frozen bars: `MESH_TRANSFER_CRITERION.md`.

## Executive verdict

| Gate | Status |
|------|--------|
| MULTIMODE_FINE_MESH | **EVALUATION_COMPLETE** |
| MESH_CONVERGENCE (M→F frozen bars) | **FAIL** |
| Power balance (fixed abs.) | **PASS** on all factorable grades (rem ≤ 0.004) |
| F as continuum reference | **NOT fully converged** vs finest air (Q) |
| Q ↔ X2 mutual agreement | **PASS** all frozen bars (independent continuum bracket) |

**Do not quote C-mesh or M-mesh FoMs as final physics.**

---

## 1. Design history (preserved)

| Stage | Mesh | J | T_des | iso dB | RL dB | TE1 |
|-------|------|--:|------:|-------:|------:|----:|
| C optimum | C | 0.708 | 0.804 | 19.0 | 16.6 | 0.019 |
| C→M transfer | M | 0.487 | 0.695 | 15.4 | 11.5 | 0.022 |
| **M refined** | M | **0.723** | **0.819** | **20.30** | **15.97** | **0.0265** |
| M→F transfer | F | 0.670 | 0.795 | 18.41 | 13.24 | 0.0203 |

M→F exceeds frozen bars on **J**, **isolation**, **RL** (T_des / TE1 / through / other within bars).

---

## 2. Crossed-refinement matrix

### Spacing targets

| grade | h_iface | h_horn | h_air | h_pml | nodes | role |
|-------|--------:|-------:|------:|------:|------:|------|
| M | 0.008 | 0.025 | 0.040 | 0.08 | 1.85M | baseline |
| Mh | 0.008 | **0.016** | 0.040 | 0.08 | 1.92M | horn only |
| Ma | 0.008 | 0.025 | **0.030** | 0.08 | 2.31M | air only |
| Mi | **0.005** | 0.025 | 0.040 | 0.08 | 3.55M | plasma/quartz only |
| Mia | **0.005** | 0.025 | **0.030** | 0.08 | 4.02M | iface+air |
| F | 0.005 | 0.016 | 0.030 | 0.06 | 4.14M | combined reference |
| R | 0.004 | 0.016 | 0.055 | 0.11 | 4.90M | iface↑ air↓ (adversarial) |
| A | 0.005 | 0.016 | 0.022 | 0.06 | 5.06M | F + finer air |
| Q | 0.005 | 0.016 | 0.018 | 0.06 | 6.05M | finest air |
| X2 | 0.0042 | 0.016 | 0.022 | 0.06 | 6.23M | finest iface @ A air |

Skipped (SuperLU abort history): **S**, **X**.

### Factorization resources

| grade | est LU GiB | factor_s | rem (power) |
|-------|----------:|---------:|------------:|
| M | 14.8 | 30 | 0.0028 |
| Mh | 15.4 | 30 | 0.0037 |
| Ma | 18.5 | 44 | 0.0027 |
| Mi | 28.4 | 74 | 0.0012 |
| Mia | 32.2 | 90 | 0.0013 |
| F | 33.1 | 93 | 0.0027 |
| R | 39.2 | 113 | 0.0019 |
| A | 40.5 | 132 | 0.0018 |
| Q | 48.4 | 180 | 0.0032 |
| X2 | 49.9 | 174 | 0.0022 |

All factorable grades stayed well under free RAM (~245 GiB).

---

## 3. Fixed-design metrics

| grade | J | T_des | IL dB | iso dB | RL dB | TE1 | through | other | R | P_out | P_mat |
|-------|--:|------:|------:|-------:|------:|----:|--------:|------:|--:|------:|------:|
| M | 0.7229 | 0.8187 | 0.869 | 20.30 | 15.97 | 0.0265 | 0.0010 | 0.0462 | 0.0253 | 0.925 | 0.072 |
| Mh | 0.7177 | 0.8206 | 0.859 | 20.28 | 15.97 | 0.0333 | 0.0011 | 0.0462 | 0.0254 | 0.934 | 0.069 |
| Ma | 0.6976 | 0.8085 | 0.923 | 19.43 | 14.04 | 0.0218 | 0.0012 | 0.0590 | 0.0395 | 0.939 | 0.058 |
| Mi | 0.7039 | 0.8111 | 0.909 | 20.06 | 14.54 | 0.0273 | 0.0009 | 0.0538 | 0.0352 | 0.936 | 0.062 |
| Mia | 0.6709 | 0.7968 | 0.986 | 18.40 | 13.29 | 0.0211 | 0.0018 | 0.0673 | 0.0469 | 0.946 | 0.053 |
| F | 0.6698 | 0.7953 | 0.995 | 18.41 | 13.24 | 0.0203 | 0.0018 | 0.0677 | 0.0475 | 0.944 | 0.053 |
| R | 0.6739 | 0.7971 | 0.985 | 17.58 | **18.67** | **0.0438** | 0.0037 | 0.0426 | 0.0136 | 0.915 | 0.083 |
| A | 0.6372 | 0.7817 | 1.070 | 17.23 | 12.71 | 0.0209 | 0.0026 | 0.0767 | 0.0536 | 0.950 | 0.048 |
| Q | 0.6199 | 0.7726 | 1.121 | 16.81 | 12.43 | 0.0202 | 0.0031 | 0.0813 | 0.0572 | 0.951 | 0.046 |
| X2 | 0.6326 | 0.7798 | 1.080 | 17.06 | 12.67 | 0.0208 | 0.0028 | 0.0779 | 0.0541 | 0.951 | 0.047 |

Machine-readable: `cross_refinement_metrics.json`.

---

## 4. Error-source identification

### M→F ΔJ attribution (from M)

| refinement | ΔJ | fraction of M→F |
|------------|---:|----------------:|
| Mh (horn only) | −0.005 | **10%** |
| Ma (air only) | −0.025 | **48%** |
| Mi (iface only) | −0.019 | **36%** |
| Mia (iface+air) | −0.052 | **98%** |
| M→F total | −0.053 | 100% |

**Mia ≈ F** on every FoM (all frozen bars PASS). Horn refinement beyond M is negligible for J / iso / RL once iface+air match F.

### Answers to campaign questions

1. **Does finer plasma/quartz account for most of M→F?**  
   Partially (~36% of ΔJ). **Air is comparable or larger (~48%).** Combined iface+air explains essentially all of it.

2. **Does air refinement matter?**  
   **Yes — strongly**, especially for isolation, RL, reflection, and other leakage.

3. **Are isolation / RL more mesh-sensitive than T_des?**  
   **Yes.** M→F: ΔT_des = −0.023 (within bar 0.04); Δiso = −1.90 dB and ΔRL = −2.73 dB (fail bars). Same pattern on further F→Q air refinement.

4. **Does TE1 converge differently from TE0?**  
   TE1 is **stable from F onward** (~0.020–0.021). TE0-derived J / iso / RL continue to move with air. (R’s TE1=0.044 with coarsened air is an outlier — do not use R as a reference.)

5. **Localized high-field rods?**  
   Node-count jump on Mi (1.85M → 3.55M) shows plasma/quartz resolution dominates DOFs; FoM shifts track both iface and air, consistent with interface fields **and** cavity standing-wave / port path resolution — not horn walls alone.

### Beyond F (is F continuum?)

| pair | frozen-bar result | note |
|------|-------------------|------|
| F → A | **PASS** | first air step from F |
| F → X2 | **PASS** | finer iface @ A air |
| F → Q | **FAIL** (iso only, Δ=−1.60 dB) | finest air moves iso past bar |
| A → Q | **PASS** | air family converging |
| A → X2 | **PASS** | iface step at fixed A air tiny |
| **Q → X2** | **PASS** | finest-air ↔ finest-iface agree |

**F is not a fully converged continuum reference** for this design: further cavity-air refinement (Q) still moves isolation outside the frozen bar vs F. However the two finest factorable adversarial meshes **Q and X2 agree within all frozen bars**, providing independent continuum evidence.

---

## 5. Convergence classification (frozen bars)

Reference for “converged” here = mutual agreement of the finest factorable pair **Q ↔ X2**.

| metric | class | evidence |
|--------|-------|----------|
| TE1 | **CONVERGED** | F/A/Q/X2 all ~0.020–0.021; Q↔X2 Δ=6e−4 |
| through | **CONVERGED** | small absolute values; within bars across F/A/Q/X2 |
| T_des | **TRENDING → near CONVERGED** | Q↔X2 Δ=0.007 (≪ 0.04); still soft drift F→Q |
| J | **TRENDING** | Q↔X2 within bar (Δ=0.013) but F→Q Δ=−0.050 at bar edge; air still mattered |
| IL | **TRENDING** | tracks T_des |
| RL | **TRENDING** | Q↔X2 Δ=0.24 dB (≪ 1.5); F→Q Δ=−0.81 |
| isolation | **NOT_CONVERGED vs F; TRENDING in continuum bracket** | F→Q fails bar; Q↔X2 PASS (Δ=0.25 dB) |
| other TE0 / R | **TRENDING** | rise under air refinement; Q↔X2 close |

**MESH_CONVERGENCE overall: FAIL** relative to the original M→F claim.  
**Continuum bracket (Q, X2): usable for quoting physics with uncertainty.**

---

## 6. Recommended single-frequency performance

Quote the **Q–X2 continuum bracket** (not M, not C). Midpoint ± half-range:

| quantity | Q | X2 | recommended | numerical uncertainty |
|----------|--:|---:|-------------:|----------------------:|
| J | 0.620 | 0.633 | **0.626** | ±0.007 |
| T_des | 0.773 | 0.780 | **0.776** | ±0.004 |
| IL (dB) | 1.121 | 1.080 | **1.10** | ±0.02 |
| isolation (dB) | 16.81 | 17.06 | **16.9** | ±0.13 |
| RL (dB) | 12.43 | 12.67 | **12.55** | ±0.12 |
| TE1 | 0.0202 | 0.0208 | **0.0205** | ±0.0003 |
| through | 0.0031 | 0.0028 | **0.0030** | ±0.0002 |
| other TE0 | 0.081 | 0.078 | **0.080** | ±0.002 |
| R | 0.057 | 0.054 | **0.056** | ±0.002 |

Power residual on Q/X2: ~0.002–0.003 (PASS).

**Caveat:** Uncertainty is the Q↔X2 half-range under the frozen bars, **not** a formal continuum proof that no further air/iface refinement would move FoMs. S/X remain unfactorable. Isolation remains the most fragile metric.

For comparison only (not physics quotes):

| | M-opt on M | M-opt on F | continuum (Q–X2) |
|--|----------:|----------:|-----------------:|
| J | 0.723 | 0.670 | **0.626 ± 0.007** |
| iso dB | 20.3 | 18.4 | **16.9 ± 0.13** |

---

## 7. Implications

- Remaining M→F drift is **physical unresolved mesh error**, dominated by **cavity air + plasma/quartz interface**, not optimizer / horn resolution.
- More M-mesh L-BFGS would not fix continuum bias.
- Next mesh work (if needed): factorable air refinement beyond Q **or** iface beyond X2 without crossing SuperLU ~6.6M abort — not more design iterations on M.
- Broadband optimization should use continuum-aware expectations (lower J / iso than M quotes) and keep adaptive sampling; channel count for 4–6 GHz is already OK.

## Artifacts

- Plan: `CROSS_REFINEMENT_PLAN.md`
- Metrics: `cross_refinement_metrics.json`
- Mesh stats: `cross_refinement_mesh_stats.json`
- S matrices: `S12_crossref_<grade>_Bp.npz`
- Grades added: `Mi`, `Ma`, `Mh`, `Mia` in `full91_geometry.py`
