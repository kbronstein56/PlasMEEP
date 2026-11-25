# Targeted Band-B Middle-Ring Density Reoptimization

**Stage-6 status: NOT STARTED. No 1 GHz / 2 GHz campaign.**

Start: Stage-5 `q` with `q[4] -= 0.15`. Free orbits: [4, 6, 5, 8, 3]. Surrogate objective: softmin hybrid **J_BB** (τ=0.05, α=0.25) on a protective Band-B frequency set. Candidate selection by hard worst A2 on that set; acceptance by dense 1.5625 MHz audit vs baseline.

## Dense-audit comparison (1.5625 MHz)

| Design | Worst A2 | Worst f (GHz) | Mean A2 | #≤−5 | #≤−7 | min Tdes | mean Tdes | mean R | mean TE1 | mean Pout |
|--------|---------:|--------------:|--------:|-----:|-----:|---------:|---------:|-------:|---------:|----------:|
| Original Stage-5 | -4.964 | 4.909375 | -0.288 | 0 | 0 | 0.009 | 0.142 | 0.284 | 0.384 | nan |
| q[4]−0.15 (no opt) | -8.616 | 4.909375 | -0.022 | 3 | 2 | 0.004 | 0.147 | 0.291 | 0.378 | nan |
| Targeted local opt | -8.889 | 5.198437 | 0.341 | 14 | 4 | 0.002 | 0.139 | 0.319 | 0.362 | nan |

## Critical finding: protective-grid aliasing

The optimizer selected a design with protective-grid worst A2 = **−4.879 at 5.200 GHz**.
The dense 1.5625 MHz audit of the **same** `q` finds A2 = **−8.889 at 5.198437 GHz**
(1.6 MHz away). At exactly 5.200 GHz the dense audit matches −4.879. So the surrogate
frequency set undersampled a newly created razor notch.

Additionally, the mandated start `q[4]−0.15` alone worsens the dense 4.909 GHz floor
from −4.964 to **−8.616** — the localization win at ~5.16 did not survive full-band
worst-case scrutiny.

## Optimization summary

- Selection rule: `best_by_protective_worst_A2`
- L-BFGS accepted iters: 12, evals: 31, wall: 74.3 min
- Optimizer message: STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT
- Δq (free) vs q4−0.15 start: q[4]=-0.2334, q[6]=-0.2456, q[5]=-0.2453, q[8]=+0.1949, q[3]=-0.2373
- Protective-grid worst at best: -4.879048076702787 @ 5.2 GHz

## Scientific questions

### Q1 — Did targeted reopt improve worst-case Band-B A2?
**no**
- orig worst=-4.964 @4.9094
- q4−0.15 worst=-8.616 @4.9094
- targeted worst=-8.889 @5.1984
- Δ vs orig = -3.925

### Q2 — Suppress original resonance or merely shift it?
**shifted_and_or_reshaped**
- new worst frequency 5.198437 GHz (orig 4.909375)

### Q3 — Did another resonance become the new limiting feature?
**yes**
- worst moved by 289.1 MHz

### Q4 — Can 0.5 GHz band become substantially more robust than −4.96?
**no — protective improvement was aliasing; dense floor worsened**
- achieved dense floor -8.889

### Q5 — Density-controllable vs geometry needed?
**density_can_reshape_spectrum_but_worst_case_not_reliably_improved_by_this_protocol**

### Q6 — Smallest justified next experiment?
**Do not expand to 16-var or Stage-6. Smallest next experiment: one forward-only single-orbit probe (e.g. orbit08_plus or orbit05_minus from localization) with a dense local audit around both 4.91 and any new soft floors — OR stop density tweaks and treat the −4.96 Band-B floor as the current validated limit for C6/40 mm/B=0.05.**
- Do not start Stage-6 or 1–2 GHz optimization from this result.

Do **not** expand to all 16 variables or launch Stage-6 from this script.

