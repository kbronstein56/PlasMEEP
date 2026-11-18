# Circulator feasibility study

**Status:** COMPLETE (queued ladder + B/DOF probes finished; no further expensive opts)
**Stage-6:** **NOT STARTED / AWAITING AUTHORIZATION**
**Physics package (frozen):** FULL91-C, magnetized plasma, C6/16 unless noted, density bounds [0.05, 1.8],
A2_hybrid τ=0.05 α=0.25, validated fast kernel, n_workers=12, OMP/OpenBLAS/MKL=1.

Companion artifacts:
- `CIRCULATOR_FEASIBILITY_BASELINE.md`
- `CIRCULATOR_PHYSICAL_LIMITS_CONTEXT.md`
- `feasibility/FEASIBILITY_COMPARISON_TABLE.md` (+ `.json`)
- `feasibility/FEASIBILITY_NOTCH_REPORT.md`
- `feasibility/feasibility_worstA2_vs_bandwidth.png`
- `feasibility/feasibility_A2_overlays_ABCDE.png`
- `feasibility/feasibility_A2_overlays_Bfield_C.png`
- `feasibility/feasibility_A2_overlays_DOF_B.png`

---

## Executive answer

**What limits broadband performance?** Evidence from the completed ladder points primarily to a
**required-bandwidth / resonant-geometry limit** of the present 91-element + C6 + B=0.05 T package,
with a secondary, real but incomplete assist from magnetic bias. C6→C3 DOF expansion and density-bound
saturation are **not** supported as the dominant cause. Denser objective sampling alone is already falsified
by Stages 2–5.

- **~0.2 GHz (Band A):** convincingly achievable — uniform start audit A2≈-0.6416, IL≈9.1 dB, iso≈8.1 dB, zero A2≤−5.
- **~0.5 GHz (Band B):** transitional soft floor — best C6 audit A2≈-4.9639 with n(A2≤−5)=0; IL/iso still poor at the worst site.
- **≥1.0 GHz (Bands C/D/E):** not demonstrated achievable — deep narrow notches return (C audit A2≈-9.6490; E true-worst A2≈-8.7411; D audit A2≈-13.6664).

**Stage-6 on the current 4–6 GHz formulation is not justified by this study.** Next effort should change
the scientific ask (narrower band and/or higher justified B / geometry), not add another dense full-band L-BFGS stage.

---

## Master comparison table

See `feasibility/FEASIBILITY_COMPARISON_TABLE.md` for the full per-case table (band, B, symmetry/n_var,
init, nit/n_eval, J_BB, audit spacing, worst A2/f/IL/iso, n≤−5/≤−7, bound counts).

### Snapshot (dense-audit authoritative)

| case | BW [GHz] | B | DOF | init | J_BB | audit A2 | f_worst | IL | iso | n≤−5 | n≤−7 | lo/hi |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| A_4p9_5p1_C6_B0.0500_stage5 | 0.200 | 0.050 | C6/16 | Stage-5 q_best | -2.7183 | -7.5401 | 5.0344 | 28.1 | -2.2 | 2 | 1 | 1/0 |
| A_4p9_5p1_C6_B0.0500_uniform | 0.200 | 0.050 | C6/16 | uniform q=1 | 0.0270 | -0.6416 | 5.0906 | 9.1 | 8.1 | 0 | 0 | 0/0 |
| B_4p75_5p25_C3_B0.0500_Bs5 | 0.500 | 0.050 | C3/31 | file:B_4p75_5p25_C6_B0.0500_stage5_opt.json | -2.5568 | -5.0208 | 4.9094 | 15.4 | -3.1 | 1 | 0 | 0/0 |
| B_4p75_5p25_C3_B0.0500_stage5 | 0.500 | 0.050 | C3/31 | Stage-5 q_best | -4.1865 | -5.0008 | 5.2266 | 14.1 | -5.5 | 2 | 0 | 16/5 |
| B_4p75_5p25_C6_B0.0500_Auni | 0.500 | 0.050 | C6/16 | file:A_4p9_5p1_C6_B0.0500_uniform_opt.json | -5.6099 | -11.4619 | 5.0328 | 30.1 | -16.1 | 55 | 20 | 0/0 |
| B_4p75_5p25_C6_B0.0500_stage5 | 0.500 | 0.050 | C6/16 | Stage-5 q_best | -2.5840 | -4.9639 | 4.9094 | 15.2 | -3.0 | 0 | 0 | 0/0 |
| C_4p5_5p5_C6_B0.0500_stage5 | 1.000 | 0.050 | C6/16 | Stage-5 q_best | -5.5296 | -9.6490 | 5.0844 | 29.8 | -10.2 | 72 | 12 | 1/0 |
| C_4p5_5p5_C6_B0.1000_Bs5 | 1.000 | 0.100 | C6/16 | file:B_4p75_5p25_C6_B0.0500_stage5_opt.json | -3.9924 | -5.3331 | 4.6094 | 14.2 | -6.6 | 21 | 0 | 1/14 |
| C_4p5_5p5_C6_B0.1000_stage5 | 1.000 | 0.100 | C6/16 | Stage-5 q_best | -3.5315 | -8.7517 | 4.5266 | 25.5 | -9.4 | 5 | 1 | 2/6 |
| D_4p25_5p75_C6_B0.0500_Bs5 | 1.500 | 0.050 | C6/16 | file:B_4p75_5p25_C6_B0.0500_stage5_opt.json | -6.4713 | -13.6664 | 5.1641 | 33.3 | -23.8 | 55 | 26 | 0/0 |
| E_4p0_6p0_C6_B0.0500_stage5 | 2.000 | 0.050 | C6/16 | Stage-5 campaign (production) | -6.0560 | -8.7411 | 4.7344 | 25.0 | -11.0 | 100 | 27 | 0/0 |

---

## 1. Bandwidth scaling (A→B→C→D→E)

| Band | classification | notes |
|---|---|---|
| A 0.2 GHz uniform | robustly feasible (tested) | best positive control |
| A 0.2 GHz Stage-5 start | notch-limited (dense audit) | bad basin on a solvable band |
| B 0.5 GHz Stage-5 start | borderline / soft floor near A2≈−5 | soft floor; 1.5625 MHz refine agreed |
| B 0.5 GHz A-uniform start | notch-limited (dense audit) | same BW, bad basin → collapse |
| C 1.0 GHz @0.05 | notch-limited (dense audit) | many narrow deep notches |
| D 1.5 GHz @0.05 | notch-limited (dense audit) | Band-B warm start |
| E 2.0 GHz Stage-5 | notch-limited (dense audit) | production full-band |

**Threshold:** performance changes from *robustly feasible* at **~0.2 GHz** to a *soft ~A2≈−5 floor*
at **~0.5 GHz**, then to *deep multi-notch failure* by **~1.0 GHz** and remains failed at 1.5–2.0 GHz.
The A-uniform→B collapse at fixed physics shows bandwidth pressure is not only an optimizer quirk.

Plot: `feasibility/feasibility_worstA2_vs_bandwidth.png`, overlays `feasibility_A2_overlays_ABCDE.png`.

---

## 2. B-field sensitivity (Band C, dense audits)

Doubling B from 0.05→0.10 T on Band C **improves** the dense-audit floor for the better start (0.05: A2=-9.6490, n≤−5=72; 0.10 Band-B warm start: A2=-5.3331, n≤−5=21), but **does not eliminate** A2≤−5 structure — notches remain, with changed depth/count/location. The Stage-5-start 0.10 T run remains deeply notched (A2=-8.7517). Conclusion: increased magnetization helps the Pareto tradeoff but is **not sufficient alone** for clean 1 GHz.

Plot: `feasibility/feasibility_A2_overlays_Bfield_C.png`.

---

## 3. Design DOF / symmetry (Band B, C3 vs C6)

On Band B, C6/16 Stage-5 start reaches audit A2=-4.9639 (n≤−5=0). C3/31 from Stage-5: A2=-5.0008 (n≤−5=2, bounds lo/hi=16/5). C3 from C6 Band-B: A2=-5.0208. Extra DOFs do **not** materially beat the C6 soft floor; one C3 start saturates many bounds without better audit performance. Evidence argues **against** C6/16 as the primary 0.5–2 GHz limiter.

Plot: `feasibility/feasibility_A2_overlays_DOF_B.png`.

---

## 4. Notch character

Full region list: `feasibility/FEASIBILITY_NOTCH_REPORT.md`.

Pattern across bandwidth:
- Narrow-band success (A uniform) shows **no** A2≤−5 regions.
- Borderline Band B shows at most shallow/touching structure near A2≈−5 without a swarm of deep notches.
- Wider bands (C/E and D if failed) add **more** deep, still **narrow** notches distributed across the band
  rather than one broad stopband — consistent with migrating/newly admitted resonant zeros.

---

## 5. Optimization-grid artifact vs design limitation

- `A_4p9_5p1_C6_B0.0500_stage5`: HIDDEN FAIL; opt_worst=-3.3549, audit_worst=-7.5401, gap=-4.1852
- `A_4p9_5p1_C6_B0.0500_uniform`: grid≈audit; opt_worst=-0.5877, audit_worst=-0.6416, gap=-0.0539
- `B_4p75_5p25_C3_B0.0500_Bs5`: HIDDEN FAIL; opt_worst=-3.5531, audit_worst=-5.0208, gap=-1.4677
- `B_4p75_5p25_C3_B0.0500_stage5`: grid≈audit; opt_worst=-4.9680, audit_worst=-5.0008, gap=-0.0327
- `B_4p75_5p25_C6_B0.0500_Auni`: HIDDEN FAIL; opt_worst=-7.2591, audit_worst=-11.4619, gap=-4.2027
- `B_4p75_5p25_C6_B0.0500_stage5`: HIDDEN FAIL; opt_worst=-3.5714, audit_worst=-4.9639, gap=-1.3925
- `C_4p5_5p5_C6_B0.0500_stage5`: HIDDEN FAIL; opt_worst=-6.9993, audit_worst=-9.6490, gap=-2.6497
- `C_4p5_5p5_C6_B0.1000_Bs5`: grid≈audit; opt_worst=-5.3069, audit_worst=-5.3331, gap=-0.0262
- `C_4p5_5p5_C6_B0.1000_stage5`: HIDDEN FAIL; opt_worst=-4.7420, audit_worst=-8.7517, gap=-4.0097
- `D_4p25_5p75_C6_B0.0500_Bs5`: HIDDEN FAIL; opt_worst=-8.5171, audit_worst=-13.6664, gap=-5.1493
- `E_4p0_6p0_C6_B0.0500_stage5`: HIDDEN FAIL; opt_worst=-7.5854, audit_worst=-8.7411, gap=-1.1557

Interpretation: Stage-5-style full-band designs remain vulnerable to inter-node notches (H1 partial),
but the ladder also contains cases where dense audits **agree** with a poor floor (C, E) and a case
where a different init solves a narrow band the Stage-5 basin cannot (A). Together this means sampling
artifacts are real **and** insufficient to explain the bandwidth collapse.

---

## 6. Hypotheses (H1–H6)

| H | Verdict | Confidence |
|---|---|---|
| H1 sampling/opt only | **Weak as sole explanation** — denser grids failed; A works with same optimizer; Auni→B collapses | medium-high |
| H2 limited DOF/C6 | **Not supported as primary** — C3 ≉ better on Band B | medium |
| H3 limited nonreciprocity | **Partial** — B=0.10 helps but leaves notches on 1 GHz | medium |
| H4 density bounds | **Unlikely primary** — best A/B solutions largely interior; saturated C3 did not win | medium-high |
| H5 geometry/modal | **Leading** for ≥1 GHz ask under current passive platform | medium |
| H6 combination | **Likely** — bandwidth×geometry primary, B secondary, basin tertiary | medium-high |

---

## 7. Smallest next diagnostics (no Stage-6)

1. **Narrow-band physics audit of a Band-A success vs Band-C failure:** at the Band-A uniform q_best and Band-C 0.05 q_best, extract port S-parameters + field energy / modal indicators on a fine local grid around each A2≤−5 notch (forward-only). Goal: decide whether failures are transmission zeros / trapped resonances of the fixed 91-element+horn geometry.
2. **Single-frequency continuous minimax probe (no new multi-day opt):** on Band C at B=0.05, take the densest audit minima and run a few dozen forward evaluations with ±1–2 MHz dither around each notch under fixed q — confirm notch width/depth are physical resonances, not FEM/noise artifacts.
3. **Only if (1–2) show movable zeros:** a *small* justified-B ladder on Band B (already nearly clean at 0.05) at one additional documented B between 0.05 and 0.10, forward+short reopt, to map whether the soft A2≈−5 floor lifts with gyrotropy before any Stage-6.
4. **Do not** launch Stage-6 on 4–6 GHz with the current formulation: denser sampling already failed through 6.25 MHz / 321 freqs; feasibility ladder shows collapse by ~1 GHz required BW.

---

## 8. Final recommendations

1. **Convincing BW of current package:** ~**200 MHz** centered near 5 GHz; ~**500 MHz** only at a weak soft floor.
2. **Degradation with BW:** yes — abrupt between 0.2 and 1.0 GHz.
3. **4–6 GHz:** not demonstrated achievable in the tested design space.
4. **Most important constraint:** required **bandwidth** interacting with **fixed geometry / modal structure**; B helps secondarily; DOF/bounds/optimizer secondary.
5. **One change most likely to help:** **narrow the target band** (or retarget usable sub-bands) **or** raise B within documented hardware limits **and** reconsider geometry — not denser Stage-6 sampling.
6. **Next major effort:** **E (narrow target)** and/or **C (justified higher B)** and/or **D (geometry)** — **not A (Stage-6 on current 4–6 GHz)**.

Stage-6 remains: **NOT STARTED / AWAITING AUTHORIZATION.**

