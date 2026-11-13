# Stage-5 frequency-sampling strategy — **NOT AUTHORIZED / DO NOT LAUNCH**

Prepared after Stage-4 completion + dense audit + notch localization + 48-frequency
high-worker parallel scaling. **Stage-5 remains NOT STARTED / AWAITING AUTHORIZATION.**

This is **not** another “add 5.08125 and 4.725” enrichment proposal.

---

## 1. Localized Stage-4 holes (fixed `q_best`, no optimization)

Source: `BB_STAGE4_NOTCH_LOCALIZATION.json` (6.25 MHz + 3.125 MHz where required).

### 5.08125 GHz (true dense-audit worst)

| Quantity | 6.25 MHz | 3.125 MHz refine |
|----------|---------:|-----------------:|
| Center | **5.08125** | **5.08125** (unchanged) |
| A2 | **−11.04** | **−11.04** |
| IL / isolation / T_des | 32.61 dB / −13.59 dB / 1.20e−2 | same at center |
| Deep core (A2 < −7) | 5.08125–5.0875 (**6.25 MHz**) | 5.08125–5.090625 (**9.4 MHz**) |
| Bad core (A2 < −5) | ~5.075–5.1375 (broader pedestal + secondary lobe) | ~5.072–5.094 (**22 MHz**) |
| Half-depth width | 43.8 MHz (inflated by secondary structure) | **9.4 MHz** |
| Structure | 1 dominant sharp minimum; secondary lobe near 5.12 | 3.125 confirms true floor at 5.08125 |

**3.125 MHz was required** to resolve the core; the deepest sample remains 5.08125.

### ~4.725 GHz feature → true center **4.73125 GHz**

| Quantity | 6.25 MHz |
|----------|---------:|
| Center | **4.73125** |
| A2 | **−10.15** |
| IL / isolation / T_des | 27.66 dB / −14.65 dB / 1.75e−3 |
| Deep core (A2 < −7) | 4.725–4.7375 (**12.5 MHz**) |
| Bad core (A2 < −5) | 4.725–4.74375 (**18.75 MHz**) |
| Half-depth width | **18.75 MHz** |
| Structure | **single** local minimum; 4.775 sits on the repaired shoulder of this pedestal |

**3.125 MHz not required** — 6.25 MHz already resolves a unique minimum with depressed neighbors.

### Width comparison across stages

| Feature | Deep-core width | Bad-core width | Notes |
|---------|----------------:|---------------:|-------|
| Stage-3 **5.8125** | **≈6.25 MHz** | ≈12.5 MHz | prior localization |
| Stage-4 **5.08125** | **≈6–9 MHz** | ≈22 MHz | 3.125 confirms |
| Stage-4 **4.73125** | **≈12.5 MHz** | ≈18.8 MHz | 6.25 adequate |

**Narrowest reliably observed bad-feature (deep-core) width: ≈ 6.25 MHz.**

Implication: any objective grid whose **maximum unsampled gap > 6.25 MHz** can again hide a corpse of the observed class between nodes (exactly the Stage-2→3→4 failure mode).

---

## 2. High-worker parallel scaling (48-frequency probe)

Source: `BROADBAND_PARALLEL_SCALING_48FREQ.json`  
Design: Stage-4 `q_best`; A2_hybrid τ=0.05 α=0.25; FULL91-C; OMP/OPENBLAS/MKL=1.

Probe: **48** frequencies spanning 4–6 GHz (keeps all tested worker counts busy).

| n_workers | batches | wall (s) | speedup vs 12w | ≈equiv cores | peak RAM (GiB) | J/grad match |
|----------:|--------:|---------:|---------------:|-------------:|---------------:|:------------:|
| **12** | 4 | **1094.7** | 1.00× | 11.6 | **44.7** | exact |
| **24** | 2 | **1038.8** | 1.05× | 23.6 | **86.3** | exact |
| **32** | 2 | 1063.1 | 1.03× | 26.6 | 111.7 | exact |
| **48** | 1 | 1083.9 | 1.01× | 47.3 | 163.5 | exact |

### Interpretation

- Numerical equivalence holds at all tested worker counts (`rel_err_J` = `rel_err_grad` = 0).
- **Wall time is essentially flat from 12 → 48 workers.**
- Per-frequency solve time inflates under concurrency (~264 s mean at 12w → ~1062 s at 48w): the job is **memory-/bandwidth-bound**, not worker-starved.
- Peak RAM scales roughly linearly with worker count (~3.7 GiB / worker).
- **Useful parallelism saturates near 12 workers** for this FEM adjoint stack. 24 workers buys ~5% wall at ~2× RAM; 32–48 buy nothing in wall and burn RAM.

**Recommended production worker count for denser Stage-5 grids: `n_workers = 12`**  
(`omp_threads_per_worker = 1`). Optional: 24 only if RAM budget is free and a few percent wall matters.

**Increased frequency parallelism does *not* make a scientifically dense grid “free.”**  
It keeps 12 workers busy on more batches; cost still scales ≈ `ceil(n_freq/12) × ~274 s` per BB eval.

Uncontended batch model used below: **`t_batch ≈ 274 s`** at 12 workers (from 48-freq / 4-batch measurement and Stage-4 production ≈16 min/eval at 47 freqs).

---

## 3. Strategy comparison

Reference: Stage-4 completed **74** BB objective/gradient evaluations in ~22.6 h wall at 47 freqs / 12 workers.  
Estimates below use **74 evals × ceil(n_freq/12) × 274 s** at `n_workers=12`.

| Strategy | n_freq | max gap | batches @12 | ≈wall / BB eval | ≈total @74 evals | Can a ~6–12 MHz notch hide? |
|----------|-------:|--------:|------------:|----------------:|-----------------:|------------------------------|
| **A. Adaptive isolated enrichment** (Stage-4 set + local 5.08/4.73 neighborhoods) | ~55 | still **25–100 MHz** in unenriched regions | 5 | ~23 min | ~28 h | **Yes — high** (new holes appear elsewhere) |
| **B. Full-band 12.5 MHz** | **161** | **12.5 MHz** | 14 | ~64 min | ~79 h | **Yes — material** (gap = 2× narrowest deep core) |
| **C. 25 MHz backbone + local 6.25/12.5 densification** of known failure neighborhoods | ~120–150 | **>6.25 MHz** outside densified bands | 10–13 | ~46–59 min | ~56–73 h | **Yes — residual** wherever densification is incomplete |
| **D. Deterministic max-gap rule ≤ 6.25 MHz over 4–6 GHz** (= full 6.25 MHz grid) | **321** | **6.25 MHz** | 27 | ~123 min | ~152 h (~6.3 d) | **Lowest** among deterministic fixed grids matching observed width |

### Strategy notes

- **A** repeats Stages 2–4. Soft-min will again overweight sampled nodes and migrate failures off-grid.
- **B** is much safer than A but still **scientifically insufficient** against a 6.25 MHz deep core (a notch can sit between 12.5 MHz nodes with both neighbors only mildly depressed — Stage-3 5.8125 between 5.800/5.825).
- **C** is a useful compromise only if densified bands cover *all* vulnerable physics; Stage-4 showed new holes at previously quiet frequencies (5.08). Hard to certify completeness a priori.
- **D** is the unique deterministic rule implied by “max unsampled gap ≤ narrowest observed deep-core width.”

Changing τ/α is **not** a substitute for D: soft-min weights only exist on evaluated frequencies.

---

## 4. Scientific diagnosis

1. **Root cause:** primarily **unsampled frequency intervals**, not a bug in A2_hybrid’s algebra. Soft-min *amplifies* the sampling defect by concentrating weight on the worst *sampled* nodes, allowing the continuum to form corpses between them.

2. **Would changing τ/α fix an off-grid failure?** **No.** τ/α only reshape the soft-min over the discrete objective set. A frequency that is never evaluated never enters soft-min, never gets an adjoint contribution, and cannot be repaired by the optimizer.

3. **Are dense-audit failures real?** **Yes.** Forward-only localization reproduces deep, smooth, repeatable notches (unique minima with depressed neighbors; 3.125 MHz confirms 5.08125). IL ~27–33 dB and strongly negative isolation are physically consistent resonant failures, not solver noise.

4. **Required sampling to claim 4–6 GHz broadband performance:**  
   **maximum objective gap ≤ ~6.25 MHz** (narrowest observed deep core). Preferentially the full uniform 6.25 MHz grid.

5. **Is full 12.5 MHz enough?** **No** for the observed notch class. Selected regions (and, safely, the whole band) need **6.25 MHz**. 3.125 MHz is only needed for *characterization*, not necessarily as the full objective spacing once 6.25 resolves the floor.

6. **Can more workers make the safer grid practical?** **Only partially.** 12→48 workers do **not** reduce wall time for this FEM; they raise RAM. The safer grid is practical on this machine **as a multi-day, checkpointed run** (~6 days at 74 evals), not as a same-day Stage-4-length job.

---

## 5. ONE recommended Stage-5 strategy — **AWAITING AUTHORIZATION**

### Recommendation: **Strategy D — full-band 6.25 MHz fixed grid**

**Why this and not “add two frequencies”:** Stages 2–4 proved isolated enrichment is structurally insufficient. The measured deep-core width (~6.25 MHz) implies a **band-wide max-gap constraint**, which is exactly the uniform 6.25 MHz grid.

### Exact frequency-generation rule

```python
def stage5_freqs_GHz(f0=4.0, f1=6.0, df=0.00625):
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]
# n_freq = 321, max gap = 6.25 MHz, deterministic, L-BFGS-safe
```

Reproducible helper: `scripts/.../stage5_frequency_set.py` (prints the list; does not optimize).

| Item | Value |
|------|------:|
| **n_freq** | **321** |
| max gap | **6.25 MHz** |
| objective / τ / α / bounds | **unchanged** A2_hybrid / 0.05 / 0.25 / [0.05, 1.8] |
| C6 / FULL91-C / FEM / adjoint | **unchanged** |
| optimizer | L-BFGS-B, maxiter=30 (or as authorized) |
| warm start | Stage-4 `q_best` |
| **n_workers** | **12** |
| omp / OpenBLAS / MKL | **1 / 1 / 1** |
| batches / BB eval | **27** (`ceil(321/12)`) |
| **≈wall / BB eval** | **~123 min (~2.05 h)** |
| **≈peak RAM** | **~45–50 GiB** (scale from 44.7 GiB @48-freq/12w) |
| **≈total @74 evals** | **~152 h (~6.3 days)** |
| checkpoint/resume | required (same pattern as Stages 3–4) |

Optional accelerator (not required): `n_workers=24` → similar total wall (~150 h) at ~90 GiB RAM — **not recommended** as default.

### Scientific justification vs Stages 2–4

| Stage | n_freq | max gap | Outcome |
|------:|-------:|--------:|---------|
| 2 | sparse + few enrichments | ≫12.5 MHz | holes at 5.75/5.825 |
| 3 | 30 | 25–50 MHz away from enrichments | repaired 5.75/5.825; **new** 5.8125 & 4.775 |
| 4 | 47 | still ≫6.25 MHz off local ladders | improved targets; **new** 5.08125 & 4.73125 |
| **5 (proposed)** | **321** | **6.25 MHz everywhere** | every 6.25 MHz interval is an objective node |

This is materially less vulnerable because soft-min can no longer park a 6–12 MHz corpse in a silent gap: every such interval contributes to J and ∂J.

### Frozen scientific package (unchanged)

- A2_hybrid, τ=0.05, α=0.25  
- C6, bounds [0.05, 1.8], FULL91-C  
- validated FEM + adjoint, port/B conventions, mesh/solver tolerances  
- deterministic frequency set (no random/jitter during L-BFGS)

---

## Status

```
STAGE-5: NOT STARTED / AWAITING AUTHORIZATION
```

Do **not** launch Stage-5 until explicitly authorized.  
Do **not** treat isolated enrichment of 5.08125 / 4.73125 as an adequate Stage-5 plan.
