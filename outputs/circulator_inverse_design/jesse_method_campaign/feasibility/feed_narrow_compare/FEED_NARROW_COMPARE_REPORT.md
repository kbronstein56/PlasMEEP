# Forward-only narrowed-feed comparison

**Stage-6: NOT STARTED / AWAITING AUTHORIZATION**
**No optimization.** Fixed plasma-density vectors q.

## Feed geometry and analytic cutoffs

- Original clear width: **40.000 mm** (2.0000 lu) → TE1 cutoff **3.7474 GHz**
- Narrowed clear width: **23.061 mm** (1.1530 lu) → TE1 cutoff **6.5000 GHz**
- Formula: `f_c,m = m * c / (2 W_m), W_m = W_lu * a; TE0 cutoff-free`
- Target: Place TE1 cutoff at 6.5 GHz so TE1 is evanescent over Band B (≤5.25 GHz) and Band C (≤5.5 GHz) with margin.

Runtime: 56.2 min for 382 forward evals.

## Watch-frequency table

| design | f_watch | feed | A2 | des0 | wrong0 | Rtot | TE1 | P_out | TE1_prop |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|
| B_stage5 | 4.909375 | original | -4.964 | 0.030 | 0.349 | 0.228 | 0.471 | 0.872 | True |
| B_stage5 | 4.909375 | narrowed | -0.602 | 0.160 | 0.450 | 0.049 | 0.004 | 0.661 | False |
| B_stage5 | 5.164062 | original | -4.915 | 0.030 | 0.245 | 0.258 | 0.332 | 0.760 | True |
| B_stage5 | 5.164062 | narrowed | -5.586 | 0.039 | 0.474 | 0.139 | 0.004 | 0.653 | False |
| C_B010_Bs5 | 4.609375 | original | -5.333 | 0.038 | 0.330 | 0.352 | 0.263 | 0.973 | True |
| C_B010_Bs5 | 4.609375 | narrowed | -6.520 | 0.020 | 0.358 | 0.600 | 0.004 | 0.979 | False |
| C_B010_Bs5 | 4.706250 | original | -4.793 | 0.019 | 0.089 | 0.813 | 0.082 | 0.985 | True |
| C_B010_Bs5 | 4.706250 | narrowed | -3.866 | 0.041 | 0.167 | 0.779 | 0.007 | 0.987 | False |
| C_B010_Bs5 | 4.996875 | original | -4.141 | 0.003 | 0.001 | 0.973 | 0.171 | 0.992 | True |
| C_B010_Bs5 | 4.996875 | narrowed | -8.124 | 0.002 | 0.014 | 0.896 | 0.001 | 0.912 | False |

## Resonance classifications under feed narrowing

### B_stage5
- **4.909375 GHz → 5_weakened_with_TE1_drop**: A2 -4.96→-0.60, TE1 0.471→0.004 (ΔA2=4.36, ΔTE1=-0.468, Δdes=0.130, Δwrong=0.102, ΔR=-0.179, shift=3.1 MHz)
- **5.164062 GHz → 2_TE1_to_wrong_TE0**: ΔTE1=-0.328, Δwrong0=0.229 (ΔA2=-0.67, ΔTE1=-0.328, Δdes=0.009, Δwrong=0.229, ΔR=-0.119, shift=0.0 MHz)

### C_B010_Bs5
- **4.609375 GHz → 3_TE1_to_reflection**: ΔTE1=-0.259, ΔR=0.247 (ΔA2=-1.19, ΔTE1=-0.259, Δdes=-0.018, Δwrong=0.028, ΔR=0.247, shift=-34.4 MHz)
- **4.706250 GHz → mixed_or_persistent**: ΔA2=0.93, ΔTE1=-0.076, Δdes=0.022, Δwrong=0.078, ΔR=-0.034 (ΔA2=0.93, ΔTE1=-0.076, Δdes=0.022, Δwrong=0.078, ΔR=-0.034, shift=-3.1 MHz)
- **4.996875 GHz → mixed_or_persistent**: ΔA2=-3.98, ΔTE1=-0.170, Δdes=-0.001, Δwrong=0.013, ΔR=-0.077 (ΔA2=-3.98, ΔTE1=-0.170, Δdes=-0.001, Δwrong=0.013, ΔR=-0.077, shift=-3.1 MHz)

## Scientific interpretation

### Method note
Narrowed-feed S uses TE0-only propagating channels (`n_prop=1`) with TE1 retained only as an
evanescent Robin load; TE1 row amplitudes are filled by the same DtN projection formula used in
production. This avoids an unphysical `|S|^2≫1` artifact that appears if imaginary-β TE1 is treated
as a propagating S channel on the narrowed aperture. Original-feed evaluations remain the production
12×12 path (`n_prop=2`).

### 1–6 classification summary

| Resonance | Outcome class | What happened |
|---|---|---|
| Band-B ~4.909 | **5 weakened (+ TE1→des0)** | A2 −4.96→−0.60; TE1 0.47→0.00; des0↑; R↓; wrong0 slightly↑; shift ~3 MHz |
| Band-B ~5.164 | **2 TE1→wrong TE0** | TE1 eliminated but wrong0 0.25→0.47; A2 slightly worse (−4.92→−5.59); no shift |
| Band-C ~4.609 | **3 TE1→reflection** | TE1 gone; R 0.35→0.60; A2 worse (−5.33→−6.52); local min shifts **−34 MHz** and deepens |
| Band-C ~4.706 | **mixed/persistent** | Already reflection-dominated (R~0.81); mild A2 improvement, still bad |
| Band-C ~5.0 | **mixed/persistent** | Reflection family persists / A2 worsens; TE1 was not the driver |

### Distinguishing the three hypotheses

1. **TE1 as a true bandwidth bottleneck** — **supported for the Band-B soft floor (~4.909 GHz)** only: suppressing TE1 converts power into desired TE0 and largely removes the soft A2≈−5 feature.
2. **TE1 as an outlet for a deeper cavity resonance** — **supported for Band-B 5.164 and Band-C 4.609**: TE1 collapses but failures remain or worsen via wrong-port or reflection redistribution; C 4.609 also **shifts and deepens**, consistent with a geometry-sensitive cavity/reflection resonance that used TE1 as one leakage channel.
3. **Feed geometry as a contributor to the resonant spectrum** — **yes**: changing only the feed width (fixed q) moves/reshapes several features. Feed is not inert. But narrowing is not a universal cure.

### Is a narrowed-feed optimization justified?

**Conditionally yes — as a small Band-B experiment only, not Stage-6.**

- The Band-B soft floor is the one place where TE1 suppression clearly helps A2 at fixed q. A short re-optimization of q under the narrowed feed on Band B (~0.5 GHz) is the smallest justified follow-up.
- Band-C / reflection-family notches are **not** fixed by TE1 cutoff; a narrowed-feed opt on 1–2 GHz is unlikely to be sufficient and is **not** justified as the next major step.
- Do **not** start Stage-6. Do **not** treat feed narrowing as proof that 4–6 GHz becomes feasible.

**Recommended next (when authorized):** one forward-validated Band-B L-BFGS under narrowed feed (same q bounds/B=0.05), then dense audit — stop and compare before any wider-band work.

No optimization was launched after this comparison.

Stage-6 remains: **NOT STARTED / AWAITING AUTHORIZATION.**

## Artifacts

- `feed_geometry_cutoffs.json`
- `feed_narrow_results.json`
- `feed_narrow_classifications.json`
- `B_stage5_spectra_orig_vs_narrow.png`
- `C_B010_Bs5_spectra_orig_vs_narrow.png`
- `C_B010_field_4p609_orig_vs_narrow.png`
