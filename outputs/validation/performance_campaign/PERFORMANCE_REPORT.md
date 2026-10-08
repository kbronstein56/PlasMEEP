# PlasMEEP High-Resolution Performance Campaign

**Date:** 2026-09-13  
**Branch context:** `agent/eigenmode-ports` (benchmark-only; no production default changes)  
**Question:** Can a scientifically trustworthy high-resolution full 91-bulb forward solve run in ≤10 minutes (ideally 3–5)?

---

## CURRENT baselines (one-source, full PMM, B=0, np=32)

| ppc | res | Cells (30×28×res²) | Steps (T≈409) | Device wall | s/step | Cell-updates/s |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 50 | 2.10 M | 40934 | **6.11 min** | 8.96 ms | 2.35×10⁸ |
| 50 | 100 | 8.40 M | 81869 | **~72 min** | 52.8 ms | 1.59×10⁸ |
| 60 | 120 | 12.1 M | 98242 | **~130 min** | 79.4 ms | 1.52×10⁸ |

Setup (`set_epsilon`) is seconds; **≥99.6% of device wall is FDTD stepping**.

### TARGET speedups required

| Goal | 50 ppc factor vs 72 min | 60 ppc factor vs 130 min |
|---|---:|---:|
| ≤10 min | ≥7.2× | ≥13× |
| 5 min | ≥14× | ≥26× |
| 3 min | ≥24× | ≥43× |

### Timestep budget at 50 ppc (measured s/step)

| Goal | Max steps | Max Meep T | vs current T=409 |
|---|---:|---:|---|
| 3 min | 3406 | **17** | need ×24 shorter |
| 5 min | 5677 | **28** | ×14 shorter |
| 10 min | 11354 | **57** | ×7.2 shorter |

**Conclusion:** With present s/step, a 50 ppc solve finishing in ≤10 min can afford only **T ≲ 57**. The production Gaussian (`df=0.10 fs`, `T_src≈389`) is mathematically incompatible with that budget unless the pulse is abandoned or radically shortened **and** still accurate.

---

## Phase 0 — Where the time goes

**Hypothesis confirmed:** High-res slowdown is **mostly expected FDTD workload** (cells × steps ~ res³), not an unexplained CPU-idle problem.

- 25→50: work ×8.0, wall ×11.8 (CUP/s drops to 0.68× — cache/MPI overhead, not idle cores)
- 50→60: work ×1.73, wall ×1.80 (nearly ideal)
- FDTD fraction ≈ 99.6–99.9% of device time
- np=32 already optimal; 48/64 slower

WSL exposes **1 NUMA node**, ~252 GiB RAM, AVX2 present. Runtime/NUMA retuning cannot deliver 7–20×.

Domain: horns/feed dominate (~45% area), PML ~26%, PMM ~11%, empty air ~18%. Conservative trims (air buffer / thinner PML) give only **~7–15%** cell reduction — not a path to 7×.

---

## Phase 1 — Short pulse + longer ringdown (`df=0.20 fs`)

| post-source | T_end | Device (25ppc) | ΔP11 vs 0.10 | ΔP1→P2 vs 0.10 |
|---:|---:|---:|---:|---:|
| 20 (prior) | 215 | 3.25 min | −1.27 dB | −8.18 dB |
| 40 | 235 | 3.49 min | −2.42 dB | −4.46 dB |
| 80 | 275 | 4.06 min | −0.39 dB | +0.16 dB |
| 120 | 315 | 4.63 min | +0.54 dB | +0.04 dB |
| df=0.15 rt40 | 300 | 4.46 min | −0.85 dB | −0.99 dB |

**Finding:** Longer ringdown **partially recovers** transmission (P1→P2 can reach ~0.04 dB), but P11 remains ~0.4–0.5 dB off and the curve is non-monotonic.  

**50 ppc projection** for the more accurate cases: **~60–80 min** — **fails ≤10 min**.  

**Verdict:** Reject as high-resolution speed solution (accuracy incomplete + projected runtime too high).

Plot: `phase1/ringdown_error_vs_T.png`

---

## Phase 2 — CW / frequency-domain

### `mp.Simulation.solve_cw` (Meep 1.30)

| Material | Result |
|---|---|
| Vacuum (complex fields) | Works (CG converges) |
| **Drude** | **`RuntimeError: dispersive materials are not yet implemented for solve_cw`** |
| **Gyrotropic Drude** | **Same error** |

**Verdict:** `solve_cw` **cannot** solve the PlasMEEP plasma model (B=0 or B≠0). Not a production path.

### ContinuousSource FDTD (same spatial mode, CW envelope)

| Cycles | until | Device 25ppc | ΔP11 | ΔP1→P2 | Est. 50ppc |
|---:|---:|---:|---:|---:|---:|
| 40 | 156 | 2.37 min | −2.73 dB | (sign failure) | ~11 min* |
| 80 | 311 | 4.66 min | −1.20 dB | −5.28 dB | ~42 min |

\*Optimistic projection if inaccurate short-CW were accepted — it is **not** accurate.

**Verdict:** Reject. No promotion to 50 ppc.

---

## Phase 3 — Domain trim

Max credible cell cut ~13% (`nx,ny` 30×28 → 28×26). **Not meaningful** for a 7× target. No 50 ppc trim run.

---

## Phase 4 — Runtime / build

- pymeep 1.30.1 + MPICH 4.3.2 + FFTW MPI (conda-forge, generic march)
- WSL: **1 NUMA node** → no multi-node interleave win
- AVX2 available; conda build is not `-march=native`
- Plausible native-build gain: **≲1.2–1.5×**, not 7×
- **Do not rebuild Meep** in this campaign

---

## Phase 5 — High-res experiments

**No acceleration candidate survived low-res accuracy + ≤10 min projection.**  
Ran one **instrumented trusted 50 ppc** one-source / multi-receiver profile (production `df=0.10`, `rt=20`):

| Metric | Measured |
|---|---:|
| Device wall | **4346 s = 72.44 min** |
| Setup (Python+geometry) | 4.07 s |
| s/step | 53.09 ms |
| Cell-updates/s | 1.582×10⁸ |
| Steps | 81869 |
| P1→P2 (normalized) | 0.002130081229278762 (**bit-match** to prior P1P2 study S21) |

Confirms prior ~72 min one-source baseline; transmission extraction is stable. (On-port P11 with 6-monitor flux subtraction can be signed/nonphysical; use 2-port study S11=0.279 as the trusted reflection reference.)

Artifacts: `phase5/trusted_ppc50_profile.json`, `.log`

---

## Phase 6 — B≠0

**Skipped** (no ≥3× accurate B=0 win).

---

## Phase 7 — Inverse-design architecture (using measured numbers)

Six-port candidate = **6 source solves** (current API).

| Resolution | 1-source (np32) | 6× sequential | 2× np32 concurrent | 4× np16 concurrent |
|---|---:|---:|---:|---:|
| 25 ppc | 6.11 min | 36.7 min | 18.3 min | ~13.6 min |
| 50 ppc | ~72 min | ~7.2 h | ~3.6 h | ~2.7 h |
| 60 ppc | ~130 min | ~13 h | ~6.5 h | ~4.8 h |

**Frequency-domain multi-RHS reuse:** Not applicable — `solve_cw` unsupported for Drude. No factorization to amortize across ports in Meep’s installed API.

**25 ppc + concurrency** can approach **~14 min/candidate** with the **trusted** Gaussian — near but not under a hard 10 min six-port budget without accuracy sacrifice.

---

## Ranked options by realistic impact

| Rank | Option | Speed | Accuracy | B≠0 | Hit ≤10 min @50? |
|---:|---|---|---|---|---|
| — | Production Gaussian 0.10 | baseline | trusted | yes | **No** |
| 1 | Optimize at 25–35 ppc + validate ranking | large (workload) | needs correlation study | yes | N/A (different res) |
| 2 | Concurrent 4×np16 at opt res | ~2.7× vs sequential | same physics | yes | only at low ppc |
| 3 | Short Gaussian + long ringdown | ~1.3× at fixed res if accurate | partial; not enough | yes | **No** (~60–80 min) |
| 4 | ContinuousSource FDTD | similar T-scaling | **poor** in tests | yes | No if accurate |
| 5 | Domain trim | ≤1.15× | untested | yes | No |
| 6 | NUMA/build tuning | ≤1.5× optimistic | n/a | yes | No |
| ✗ | `solve_cw` | would be large if it worked | — | — | **Blocked by Drude** |

---

## Explicit answers

1. **Why 1–2 h at 50–60 ppc?** ~8–12 M cells × ~8–10×10⁴ steps at ~50–80 ms/step; FDTD-dominated; Gaussian needs T≈409.
2. **Software bug vs workload?** Predominantly **expected FDTD workload**. Mild CUP/s drop with res; not idle-CPU.
3. **Can CURRENT TD algorithm hit ≤10 min @50?** **Not with trusted 0.10 fs pulse.** Need T≲57; production pulse needs T≈409.
4. **3–5 min @50?** **No** with current algorithm + accuracy.
5. **Would CW change that?** **Installed `solve_cw` cannot run Drude plasma.** ContinuousSource FDTD did not give accurate S-params in the affordable-T regime.
6. **Fastest accurate full-91 run this campaign?** Trusted **25 ppc / 6.11 min** (scaling). Among experiments, no new method was both faster **and** accurate to ≲0.1 dB on P11 & P1→P2.
7. **Fastest accurate high-res?** **50 ppc = 72.44 min / source** (this campaign’s instrumented profile; matches prior).
8. **Factor still needed @50 for 10 / 5 min?** ~**7× / 14×** beyond trusted TD, **or** drop resolution for optimization.
9. **Implement next in production?**  
    - Keep high-res for **validation / formulation**.  
    - Build **resolution-ranking correlation** study (25/35/50).  
    - Optimize at lowest ppc that preserves ranking.  
    - Use **4×np16 concurrency** for six-port candidates.  
    - Do **not** ship wider Gaussian or ContinuousSource without accuracy recovery.
10. **What res can meet target and still deserve validation?** **~25–35 ppc** is the realistic opt band for ≤10–15 min six-port candidates; must prove ranking vs 50+.
11. **Ranking experiment:** Fixed design set (20–40 ρ); evaluate objective at 25/30/35/40/50; Spearman + scatter vs 50; gate on rank agreement of top-K.
12. **3–5 min per six-port candidate?** **Not realistic** at 50+ ppc. At 25 ppc with concurrency, **~14 min** is realistic; 3–5 min would need further T-cut **with** proven accuracy (not demonstrated).
13. **≤10 min six-port?** **Not at 50 ppc with trusted TD.** **Borderline/plausible at 25 ppc** only with concurrency (~14 min today — still slightly over 10).
14. **If no:** Need either (a) **opt at lower validated ppc**, (b) a **true frequency-domain solver that supports gyrotropic Drude** (not Meep `solve_cw`), or (c) **much more hardware** (many nodes) — single-socket MPI is saturated.

---

## Artifacts

- `phase0/accounting.json`
- `phase1/` ringdown JSONs + `ringdown_error_vs_T.png`
- `phase2/cw_compat_v2.json`, `cw_fdtd_cyc*_ppc25.json`, `CONCLUSION.json`
- `phase3/domain_audit.json`
- `phase4/runtime_audit.txt`
- `phase5/` trusted 50 ppc profile (when complete)

**No production defaults modified. No commits. No 75/100 ppc runs.**
