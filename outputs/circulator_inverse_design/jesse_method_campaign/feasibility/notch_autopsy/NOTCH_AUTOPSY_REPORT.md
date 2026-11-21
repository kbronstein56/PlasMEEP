# Feasibility notch autopsy (forward-only)

**Stage-6: NOT STARTED / AWAITING AUTHORIZATION**
**No new optimization.**

Modal powers are averages over the six TE0 source directions, extracted from the
production 12×12 eigenmode-port S matrix (`|S_{mode,port; TE0 source}|^2`).
`P_out` is the sum of all outgoing TE0+TE1 modal powers (power-budget closure vs 1).

Classification heuristic uses absolute power fractions:
- **A_TE1_conversion** — large total TE1
- **B_reflection** — large R0+R1
- **C_wrong_TE0_port** — large rev0+thr0+other0
- **D_mixed** — two channels both large
- **E_power_budget_inconsistent** — |1−P_out| > 0.15

## Numerical modal-power table

| design | B | f [GHz] | role | A2 | des0 | rev0 | wrong0 | R0 | R1 | TE1 | P_out | miss | class |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| A_uniform | 0.050 | 4.609375 | crossband_probe | 0.744 | 0.088 | 0.003 | 0.067 | 0.652 | 0.022 | 0.134 | 0.942 | 0.058 | B_reflection |
| A_uniform | 0.050 | 4.909375 | crossband_probe | 2.344 | 0.103 | 0.001 | 0.022 | 0.578 | 0.043 | 0.270 | 0.972 | 0.028 | B_reflection |
| A_uniform | 0.050 | 4.971875 | control:best_in_band | 3.791 | 0.103 | 0.000 | 0.031 | 0.598 | 0.037 | 0.158 | 0.890 | 0.110 | B_reflection |
| A_uniform | 0.050 | 5.000000 | control:band_center | 1.156 | 0.160 | 0.006 | 0.100 | 0.580 | 0.009 | 0.093 | 0.932 | 0.068 | B_reflection |
| A_uniform | 0.050 | 5.056250 | control:near_soft_side | -0.153 | 0.125 | 0.014 | 0.142 | 0.534 | 0.003 | 0.018 | 0.818 | 0.182 | E_power_budget_inconsistent |
| A_uniform | 0.050 | 5.088625 | dither:soft_min | -0.687 | 0.174 | 0.031 | 0.255 | 0.261 | 0.053 | 0.116 | 0.806 | 0.194 | E_power_budget_inconsistent |
| A_uniform | 0.050 | 5.089625 | dither:soft_min | -0.683 | 0.160 | 0.028 | 0.242 | 0.309 | 0.045 | 0.101 | 0.811 | 0.189 | E_power_budget_inconsistent |
| A_uniform | 0.050 | 5.090625 | notch:soft_min | -0.642 | 0.134 | 0.021 | 0.210 | 0.340 | 0.037 | 0.087 | 0.770 | 0.230 | E_power_budget_inconsistent |
| A_uniform | 0.050 | 5.091625 | dither:soft_min | -0.516 | 0.127 | 0.017 | 0.218 | 0.370 | 0.033 | 0.085 | 0.800 | 0.200 | E_power_budget_inconsistent |
| A_uniform | 0.050 | 5.092625 | dither:soft_min | -0.623 | 0.102 | 0.014 | 0.198 | 0.423 | 0.021 | 0.074 | 0.796 | 0.204 | E_power_budget_inconsistent |
| A_uniform | 0.050 | 5.164062 | crossband_probe | -6.279 | 0.012 | 0.026 | 0.505 | 0.133 | 0.046 | 0.249 | 0.899 | 0.101 | C_wrong_TE0_port |
| B_stage5 | 0.050 | 4.609375 | crossband_probe | 1.747 | 0.315 | 0.011 | 0.088 | 0.167 | 0.291 | 0.334 | 0.904 | 0.096 | D_mixed |
| B_stage5 | 0.050 | 4.784375 | control:left_of_notch | 2.464 | 0.330 | 0.006 | 0.051 | 0.064 | 0.040 | 0.443 | 0.888 | 0.112 | A_TE1_conversion |
| B_stage5 | 0.050 | 4.907375 | dither:A2_local_min | -3.952 | 0.047 | 0.052 | 0.329 | 0.002 | 0.294 | 0.507 | 0.885 | 0.115 | D_mixed |
| B_stage5 | 0.050 | 4.908375 | dither:A2_local_min | -4.655 | 0.034 | 0.057 | 0.345 | 0.007 | 0.253 | 0.490 | 0.877 | 0.123 | D_mixed |
| B_stage5 | 0.050 | 4.909375 | notch:A2_local_min | -4.964 | 0.030 | 0.061 | 0.349 | 0.021 | 0.207 | 0.471 | 0.872 | 0.128 | D_mixed |
| B_stage5 | 0.050 | 4.910375 | dither:A2_local_min | -4.649 | 0.036 | 0.063 | 0.341 | 0.046 | 0.160 | 0.447 | 0.870 | 0.130 | D_mixed |
| B_stage5 | 0.050 | 4.911375 | dither:A2_local_min | -3.960 | 0.050 | 0.064 | 0.328 | 0.079 | 0.123 | 0.417 | 0.874 | 0.126 | D_mixed |
| B_stage5 | 0.050 | 4.967187 | control:best_in_band | 6.586 | 0.371 | 0.000 | 0.065 | 0.270 | 0.100 | 0.188 | 0.893 | 0.107 | B_reflection |
| B_stage5 | 0.050 | 4.967187 | control:right_of_notch | 6.586 | 0.371 | 0.000 | 0.065 | 0.270 | 0.100 | 0.188 | 0.893 | 0.107 | B_reflection |
| B_stage5 | 0.050 | 4.978125 | notch:high_R_local_min | -3.625 | 0.011 | 0.003 | 0.005 | 0.656 | 0.064 | 0.147 | 0.819 | 0.181 | E_power_budget_inconsistent |
| B_stage5 | 0.050 | 5.162062 | dither:A2_local_min | -2.910 | 0.065 | 0.041 | 0.222 | 0.128 | 0.103 | 0.349 | 0.763 | 0.237 | E_power_budget_inconsistent |
| B_stage5 | 0.050 | 5.163062 | dither:A2_local_min | -3.990 | 0.042 | 0.044 | 0.237 | 0.146 | 0.100 | 0.337 | 0.762 | 0.238 | E_power_budget_inconsistent |
| B_stage5 | 0.050 | 5.164062 | notch:A2_local_min | -4.915 | 0.030 | 0.052 | 0.245 | 0.153 | 0.105 | 0.332 | 0.760 | 0.240 | E_power_budget_inconsistent |
| B_stage5 | 0.050 | 5.165062 | dither:A2_local_min | -4.826 | 0.030 | 0.062 | 0.247 | 0.148 | 0.115 | 0.332 | 0.757 | 0.243 | E_power_budget_inconsistent |
| B_stage5 | 0.050 | 5.166062 | dither:A2_local_min | -4.353 | 0.039 | 0.072 | 0.245 | 0.136 | 0.127 | 0.336 | 0.754 | 0.246 | E_power_budget_inconsistent |
| C_B010_Bs5 | 0.050 | 4.607375 | Bswap_0p05:dither | -1.649 | 0.044 | 0.005 | 0.201 | 0.401 | 0.005 | 0.336 | 0.982 | 0.018 | D_mixed |
| C_B010_Bs5 | 0.050 | 4.608375 | Bswap_0p05:dither | -1.817 | 0.042 | 0.006 | 0.199 | 0.404 | 0.005 | 0.337 | 0.982 | 0.018 | D_mixed |
| C_B010_Bs5 | 0.050 | 4.609375 | Bswap_0p05:notch | -1.981 | 0.040 | 0.006 | 0.198 | 0.407 | 0.006 | 0.337 | 0.982 | 0.018 | D_mixed |
| C_B010_Bs5 | 0.050 | 4.610375 | Bswap_0p05:dither | -2.142 | 0.038 | 0.006 | 0.196 | 0.410 | 0.006 | 0.338 | 0.982 | 0.018 | D_mixed |
| C_B010_Bs5 | 0.050 | 4.611375 | Bswap_0p05:dither | -2.300 | 0.036 | 0.007 | 0.195 | 0.413 | 0.007 | 0.339 | 0.982 | 0.018 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.565625 | control:left_of_notch | -3.219 | 0.036 | 0.018 | 0.130 | 0.422 | 0.137 | 0.364 | 0.953 | 0.047 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.607375 | dither:A2_le_m5_notch | -5.297 | 0.042 | 0.194 | 0.362 | 0.274 | 0.014 | 0.293 | 0.970 | 0.030 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.608375 | dither:A2_le_m5_notch | -5.322 | 0.040 | 0.186 | 0.346 | 0.308 | 0.012 | 0.277 | 0.971 | 0.029 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.609375 | notch:A2_le_m5_notch | -5.333 | 0.038 | 0.177 | 0.330 | 0.342 | 0.011 | 0.263 | 0.973 | 0.027 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.610375 | dither:A2_le_m5_notch | -5.333 | 0.037 | 0.168 | 0.313 | 0.374 | 0.009 | 0.249 | 0.974 | 0.026 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.611375 | dither:A2_le_m5_notch | -5.324 | 0.036 | 0.160 | 0.298 | 0.405 | 0.008 | 0.236 | 0.975 | 0.025 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.684375 | control:right_of_notch | -0.470 | 0.124 | 0.016 | 0.125 | 0.486 | 0.013 | 0.235 | 0.970 | 0.030 | B_reflection |
| C_B010_Bs5 | 0.100 | 4.706250 | notch:A2_secondary_min | -4.793 | 0.019 | 0.029 | 0.089 | 0.794 | 0.018 | 0.082 | 0.985 | 0.015 | B_reflection |
| C_B010_Bs5 | 0.100 | 4.787500 | control:best_in_band | 3.691 | 0.232 | 0.001 | 0.406 | 0.043 | 0.027 | 0.254 | 0.935 | 0.065 | D_mixed |
| C_B010_Bs5 | 0.100 | 4.996875 | notch:high_R_local_min | -4.141 | 0.003 | 0.000 | 0.001 | 0.817 | 0.156 | 0.171 | 0.992 | 0.008 | B_reflection |

## Dither smoothness (±1, ±2 MHz)

- **A_uniform** parent 5.090625 GHz: A2=['-0.687', '-0.683', '-0.642', '-0.516', '-0.623'], TE1=['0.116', '0.101', '0.087', '0.085', '0.074'], wrong0=['0.255', '0.242', '0.210', '0.218', '0.198'], Rtot=['0.315', '0.354', '0.376', '0.403', '0.444']; smooth/continuous redistribution
- **B_stage5** parent 4.909375 GHz: A2=['-3.952', '-4.655', '-4.964', '-4.649', '-3.960'], TE1=['0.507', '0.490', '0.471', '0.447', '0.417'], wrong0=['0.329', '0.345', '0.349', '0.341', '0.328'], Rtot=['0.297', '0.260', '0.228', '0.206', '0.202']; smooth/continuous redistribution
- **B_stage5** parent 5.164062 GHz: A2=['-2.910', '-3.990', '-4.915', '-4.826', '-4.353'], TE1=['0.349', '0.337', '0.332', '0.332', '0.336'], wrong0=['0.222', '0.237', '0.245', '0.247', '0.245'], Rtot=['0.231', '0.247', '0.258', '0.263', '0.263']; smooth/continuous redistribution
- **C_B010_Bs5** parent 4.609375 GHz: A2=['-5.297', '-5.322', '-5.333', '-5.333', '-5.324'], TE1=['0.293', '0.277', '0.263', '0.249', '0.236'], wrong0=['0.362', '0.346', '0.330', '0.313', '0.298'], Rtot=['0.287', '0.320', '0.352', '0.384', '0.413']; smooth/continuous redistribution
- **C_B010_Bs5 B=0.05 swap** parent 4.609375: A2=['-1.649', '-1.817', '-1.981', '-2.142', '-2.300'], TE1=['0.336', '0.337', '0.337', '0.338', '0.339'], wrong0=['0.201', '0.199', '0.198', '0.196', '0.195'], Rtot=['0.406', '0.409', '0.413', '0.416', '0.420']

## Physical interpretation

### 1. Are the narrow failures TE1, reflection, wrong-port, or mixed?

**Primarily mixed (D), with a clear channel ranking by case — not a single pure mechanism.**

| Site | A2 | des0 | wrong0 | Rtot (R0/R1) | TE1 | P_out | Dominant read |
|------|---:|-----:|-------:|-------------:|----:|------:|---|
| C @4.609375 (B=0.10) | −5.33 | 0.038 | **0.330** | **0.352** (0.342/0.011) | 0.263 | 0.973 | **D: reflection ≈ wrong TE0** (rev0=0.177); TE1 secondary |
| B @4.909375 | −4.96 | 0.030 | **0.349** | 0.228 (0.021/**0.207**) | **0.471** | 0.872 | **D: TE1 + wrong TE0**; TE1-heavy, R mostly R1 |
| B @5.164062 | −4.92 | 0.030 | 0.245 | 0.258 | 0.332 | 0.760 | mixed + larger unaccounted power |
| C @4.706 / @4.997 | −4.79 / −4.14 | ≈0 | ≈0 | **0.81–0.97** | low–mod | ≈0.99 | **B: reflection resonances** (separate family) |
| A soft min @5.090625 | −0.64 | 0.134 | 0.210 | 0.376 | 0.087 | 0.770 | mild; not an A2≤−5 corpse |

Power-budget note: class **E** (`|1−P_out|>0.15`) can reflect **physical plasma absorption** as well as numerics. Prefer reading E alongside the modal splits rather than as “solver failure.” Band-C A2≤−5 notch has **good** closure (P_out≈0.97) → physical modal redistribution, not an accounting glitch.

Dithers at ±1/±2 MHz are **smooth/continuous** at all primary notches → narrow but physical resonances, not single-bin FEM spikes.

### 2. Does Band A avoid the resonance family, or place it outside the band?

**Both, depending on which resonance.**

- At Band-C notch freq 4.609375 and Band-B 4.909375, the **same A-uniform q** stays mild (A2=+0.74 and +2.34) — those specific deep failures are **design-dependent**, not hard geometry zeros present for every q.
- At **5.164062 GHz** (just above the 5.1 GHz Band-A edge), the same A-uniform q collapses: **A2=−6.28**, class **C_wrong_TE0_port** (wrong0=0.505, des0=0.012).

So Band A’s in-band success is partly **lucky band placement relative to a wrong-port failure just above 5.1 GHz**, not proof that this density pattern is globally free of the failure family.

### 3. Does Band B put the first such resonance at/near the band edge?

**No — the soft floor is mid-band.**

Deepest Band-B minimum is at **4.909 GHz**, ~159 MHz inside 4.75–5.25, with TE1≈0.47 and wrong0≈0.35. A second comparable minimum sits near **5.164 GHz** (also interior). The usable-band limitation is an **in-band modal tradeoff**, not edge clipping.

### 4. Does increasing B move/split the same resonances, or change modal character?

**Qualitative change on this fixed q, not a simple rigid shift.**

Same Band-C optimum q at 4.609375 GHz:
- **B=0.10 (native):** A2=−5.33; wrong0=0.330, rev0=0.177, R0=0.342, TE1=0.263 — isolation/routing collapse.
- **B=0.05 (swap):** A2=−1.98; wrong0=0.198, rev0=0.006, R0=0.407, TE1=0.337 — more reflective, far less reverse TE0.

B=0.10 dither shows a flat A2≈−5.33 bottom; B=0.05 dither is smoother and shallower without the same reverse-port catastrophe. Higher B here **amplifies wrong-port content at this notch** rather than merely translating an unchanged resonance. (This does not contradict the earlier feasibility result that a *reoptimized* B=0.10 design can beat a B=0.05 Stage-5 basin on Band C overall — it shows bias changes modal character, so B is not a free “more is better” knob on a fixed density pattern.)

### 5. Smallest justified next physical change

**Recommended smallest next test (still no Stage-6 / no new long opt):**

1. **Primary:** a **feed / TE1-cutoff or mode-filter geometry probe** (or equivalent horn-width change already supportable in the repo), because the Band-B soft floor is **TE1-dominated (≈0.47)** and Band-C’s deep notch still carries TE1≈0.26 alongside wrong-port/reflection.
2. **Co-primary:** treat **wrong-port routing** as first-class — any geometry/B experiment must track rev0/thr0/other0, not only TE1.
3. **Not sufficient alone:** denser Stage-6 sampling, or “just raise B” on a frozen q.
4. **Initialization/basin** remains real for narrow bands (A Stage-5 vs uniform) but is secondary once these modal channels are visible.

Field maps: `C_B010_Bs5_field_notch_vs_good.png` — notch (4.609 GHz) vs best-in-band (4.7875 GHz) at identical |Hz| scale. Notch shows much higher reflection (R≈0.34 vs 0.04) with comparable TE1≈0.25–0.26; consistent with a mixed reflective/wrong-routing resonant state rather than pure TE1 conversion.


### Field-map note (Band C)

At identical |Hz| scale, the 4.609 GHz notch shows a more speckled/high-intensity central lattice response and strong standing-wave structure in the driven arm (consistent with R≈0.34), whereas 4.7875 GHz is cleaner through-routing with R≈0.04. TE1 remains ~0.25–0.26 in both frames — the notch is **not** explained by TE1 alone; it looks like a mixed reflective / wrong-routing resonant state of the fixed array+horn geometry.

Stage-6 remains: **NOT STARTED / AWAITING AUTHORIZATION.**
