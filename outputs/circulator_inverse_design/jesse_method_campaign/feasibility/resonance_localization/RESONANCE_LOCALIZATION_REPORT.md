# Band-B ~5.16 GHz Resonance Localization & Sensitivity

**Stage-6 status: NOT STARTED. No optimization launched.**

Baseline: Band-B C6 Stage-5 `q`, B = 0.05 T, original 40 mm feed. Forward-only. Evals ≈ 738, wall ≈ 125.2 min.

## Executive verdict

- The ~5.16 GHz feature is an **array_core_cavity_plasma_hybrid**: field energy brightens in the array core at the notch; plasma |Hz|² fraction rises ~4.0× vs the left control, with **middle-ring** bulbs dominating plasma participation.
- It is **not** feed-dominated (consistent with the feed-width sweep). Center bulb alone is negligible; outer ring can move the notch but with severe Band-B collateral damage.
- Δq = ±0.15 does **not** push the resonance outside 4.75–5.25 GHz (needs ≳86 MHz up-shift).
- **Smallest justified next change:** `orbit04_minus` (Δf=-31.6 MHz, ΔA2_min=+1.46, mean sentinel ΔA2=-1.66). Supports a **targeted middle-ring density reopt**, not a geometry redesign yet.

## Resonance characterization

- Dense-audit minimum (1.5625 MHz): **5.164062 GHz**
- Controls: left 5.059375 GHz, right 5.226562 GHz, best-in-band 4.967187 GHz

| tag | f (GHz) | A2 | Tdes | R | TE1 | plasma |Hz|² frac | dominant ring |
|-----|--------:|---:|-----:|--:|----:|----------------------:|---------------|
| notch_5p16 | 5.164062 | -4.915 | 0.030 | 0.153 | 0.332 | 0.126 | middle |
| ctrl_left | 5.059375 | 3.998 | 0.279 | 0.158 | 0.355 | 0.032 | outer |
| ctrl_right | 5.226562 | 2.671 | 0.223 | 0.055 | 0.490 | 0.088 | middle |
| ctrl_best_BandB | 4.967187 | 6.586 | 0.371 | 0.270 | 0.188 | 0.087 | outer |

**Spatial character:** `array_core_cavity_plasma_hybrid`
- Plasma energy fraction at notch: 0.126 (ctrl_left 0.032, contrast ×3.96); dominant ring: **middle**
- Ring fractions: center=0.002, inner=0.240, middle=0.569, outer=0.189
- Top orbits: q[8](r=3.61)=0.141, q[6](r=3.00)=0.138, q[3](r=2.00)=0.108, q[9](r=3.61)=0.075, q[2](r=1.73)=0.075

## Orbit / ring catalog

| ring | orbits | radii (a) |
|------|--------|-----------|
| center | [0] | 0.000 |
| inner | [1, 2, 3] | 1.000, 1.732, 2.000 |
| middle | [4, 5, 6, 7, 8, 9, 10] | 2.646, 2.646, 3.000, 3.464, 3.606, 3.606, 4.000 |
| outer | [11, 12, 13, 14, 15] | 4.359, 4.359, 4.583, 4.583, 5.000 |

Perturbation amplitude: Δq = ±0.15 (clipped to [0.05,1.8]).

## Sensitivity ranking

### Top by |Δf_res|

| rank | name | Δf (MHz) | ΔA2_min | sens f (MHz/Δq) | sentinel min ΔA2 |
|-----:|------|---------:|--------:|----------------:|-----------------:|
| 1 | `ring_outer_plus` | +43.44 | -11.260 | +289.6 | -6.91 |
| 2 | `orbit09_minus` | +43.44 | +0.843 | -289.6 | -2.87 |
| 3 | `orbit12_minus` | +43.44 | +2.119 | -289.6 | -7.05 |
| 4 | `orbit14_minus` | +43.44 | +3.375 | -289.6 | -11.49 |
| 5 | `ring_outer_minus` | -31.56 | +2.371 | +210.4 | -5.61 |
| 6 | `orbit02_plus` | -31.56 | +0.482 | -210.4 | -6.54 |
| 7 | `orbit04_minus` | -31.56 | +1.458 | +210.4 | -1.92 |
| 8 | `orbit06_plus` | -31.56 | +2.186 | -210.4 | -3.94 |
| 9 | `orbit06_minus` | -31.56 | +1.373 | +210.4 | -2.70 |
| 10 | `orbit10_plus` | -31.56 | -4.294 | -210.4 | -5.53 |

### Top by shallowness (ΔA2_min)

| rank | name | ΔA2_min | Δf (MHz) | sentinel min ΔA2 |
|-----:|------|--------:|---------:|-----------------:|
| 1 | `orbit14_minus` | +3.375 | +43.44 | -11.49 |
| 2 | `orbit03_plus` | +2.839 | -19.06 | -3.34 |
| 3 | `ring_outer_minus` | +2.371 | -31.56 | -5.61 |
| 4 | `orbit06_plus` | +2.186 | -31.56 | -3.94 |
| 5 | `orbit12_minus` | +2.119 | +43.44 | -7.05 |
| 6 | `orbit08_plus` | +1.709 | +0.00 | -2.99 |
| 7 | `ring_middle_plus` | +1.589 | -12.50 | -8.74 |
| 8 | `orbit04_minus` | +1.458 | -31.56 | -1.92 |
| 9 | `orbit05_minus` | +1.379 | +0.00 | -0.59 |
| 10 | `orbit06_minus` | +1.373 | -31.56 | -2.70 |

## Scientific questions

### Q1 — Localized or global?
**array_core_cavity_plasma_hybrid**

### Q2 — Strongest frequency control?
- Best ring: `ring_outer_plus` Δf=+43.44 MHz
- Best orbit: `orbit09_minus` Δf=+43.44 MHz (r=3.6055512754639896)

### Q3 — Strongest depth control?
- Best ring: `ring_outer_minus` ΔA2=+2.371
- Best orbit: `orbit08_plus` ΔA2=+1.709

### Q4 — Can a small perturbation move it outside Band B?
**not_with_delta_q_0p15**
- Distance to band edges: −414.1 MHz / +85.9 MHz

### Q5 — Geometry redesign vs targeted density reopt?
**targeted_density_reopt_plausible**
- `orbit04_minus`: Δf=-31.6 MHz, ΔA2=+1.46, sentinel min ΔA2=-1.9172833331107908

### Q6 — Single smallest justified next design change?
**`orbit04_minus`** — Δf=-31.56 MHz, ΔA2_min=+1.458, mean sentinel ΔA2=-1.6568522048736043
- Highest combined score among forward-only ±0.15 C6/ring density perturbations for moving/shallowing ~5.16 GHz with limited Band-B damage.
- Physically: lower the C6 density scale on middle-ring orbit r≈2.646 a (orbit index 4). Related helpful directions: decrease q[6] or increase q[3]/q[8] with care — always check Band-B sentinels.

Do **not** auto-launch optimization. Stage-6 remains **NOT STARTED**.

