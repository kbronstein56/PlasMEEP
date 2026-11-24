# Band-B Fixed-q Feed-Width Sweep

**Stage-6 status: NOT STARTED**

Fixed Band-B C6 Stage-5 optimized `q` at B = 0.05 T. No q optimization.
Dense forward audit over 4.65–5.35 GHz (Band B = 4.75–5.25 GHz + guard).
Evals: 525. Wall time: 73.8 min.

## Executive verdict

**Best fixed-q feed for worst-case Band-B A2 is still the original 40.00 mm**
(worst A2 = −4.96 at 4.909 GHz). Every intermediate and narrowed width
**degrades** the Band-B worst point (down to −6.73 at 35 mm).

Narrowing does what the prior TE1 test predicted at the old 4.91 GHz notch —
that family softens and redshifts — but the ~5.16 GHz family stays inside
Band B, deepens, and reclassifies TE1 → wrong-TE0 (or reflection at mid-widths).
TE1 is suppressible; worst-case broadband A2 is not improved by feed width alone.

**Recommendation: do not launch a narrowed/intermediate-feed Band-B q reoptimization
on the strength of this fixed-q sweep.** Stage-6 remains not started.

## Analytical cutoffs (solver modal convention)

| W (mm) | W (lu) | TE1 fc (GHz) | TE2 fc (GHz) | Note |
|-------:|-------:|-------------:|-------------:|------|
| 40.00 | 2.0000 | 3.747 | 7.495 | original production; TE1 open throughout Band B |
| 35.00 | 1.7500 | 4.283 | 8.565 | TE1 open throughout Band B |
| 30.00 | 1.5000 | 4.997 | 9.993 | TE1 cutoff inside Band B |
| 28.00 | 1.4000 | 5.353 | 10.707 | TE1 modestly above Band-B upper edge |
| 26.50 | 1.3250 | 5.656 | 11.313 | extra TE1 cutoff margin over Band B |
| 25.00 | 1.2500 | 5.996 | 11.992 | extra TE1 cutoff margin over Band B |
| 23.06 | 1.1530 | 6.500 | 13.001 | prior aggressive 6.5 GHz cutoff; extra TE1 cutoff margin over Band B |

Formula: `f_c,m = m · c / (2 W_m)`, `W_m = W_lu · a`, `a = 20 mm`.

## Worst-case Band-B dense A2 vs width

| W (mm) | TE1 fc | Worst A2 | Worst f (GHz) | Mean A2 | Mean TE1 | Mean wrong TE0 | Mean R | Worst class |
|-------:|-------:|---------:|--------------:|--------:|---------:|---------------:|-------:|-------------|
| 40.00 | 3.747 | -4.964 | 4.9094 | -0.696 | 0.378 | 0.106 | 0.472 | `dominated_by_TE1` |
| 35.00 | 4.283 | -6.730 | 5.1500 | -0.529 | 0.457 | 0.116 | 0.538 | `dominated_by_reflection` |
| 30.00 | 4.997 | -5.969 | 4.9750 | 0.170 | 0.115 | 0.301 | 0.512 | `dominated_by_reflection` |
| 28.00 | 5.353 | -6.178 | 5.1641 | 0.091 | 0.015 | 0.202 | 0.361 | `dominated_by_wrong_TE0` |
| 26.50 | 5.656 | -5.952 | 5.1641 | -0.004 | 0.012 | 0.199 | 0.381 | `dominated_by_wrong_TE0` |
| 25.00 | 5.996 | -5.779 | 5.1641 | -0.044 | 0.005 | 0.196 | 0.384 | `dominated_by_wrong_TE0` |
| 23.06 | 6.500 | -5.582 | 5.1641 | -0.111 | 0.004 | 0.195 | 0.399 | `dominated_by_wrong_TE0` |

## Resonance family tracking

### near_4p91

| W (mm) | f (GHz) | Δf from seed (MHz) | A2 | TE1 | wrong TE0 | R | class |
|-------:|--------:|-------------------:|---:|----:|----------:|--:|-------|
| 40.00 | 4.9094 | -0.6 | -4.964 | 0.471 | 0.349 | 0.228 | `dominated_by_TE1` |
| 35.00 | 4.8250 | -85.0 | -5.133 | 0.677 | 0.003 | 1.216 | `dominated_by_reflection` |
| 30.00 | 4.7625 | -147.5 | -2.006 | 0.007 | 0.200 | 0.652 | `dominated_by_reflection` |
| 28.00 | 4.7625 | -147.5 | -1.703 | 0.002 | 0.194 | 0.654 | `dominated_by_reflection` |
| 26.50 | 4.7625 | -147.5 | -1.584 | 0.002 | 0.192 | 0.657 | `dominated_by_reflection` |
| 25.00 | 4.7625 | -147.5 | -1.572 | 0.001 | 0.188 | 0.662 | `dominated_by_reflection` |
| 23.06 | 4.7625 | -147.5 | -1.542 | 0.001 | 0.183 | 0.670 | `dominated_by_reflection` |

### near_5p16

| W (mm) | f (GHz) | Δf from seed (MHz) | A2 | TE1 | wrong TE0 | R | class |
|-------:|--------:|-------------------:|---:|----:|----------:|--:|-------|
| 40.00 | 5.1641 | +4.1 | -4.915 | 0.332 | 0.245 | 0.258 | `dominated_by_TE1` |
| 35.00 | 5.1500 | -10.0 | -6.730 | 0.448 | 0.080 | 0.517 | `dominated_by_reflection` |
| 30.00 | 5.1625 | +2.5 | -5.200 | 0.226 | 1.036 | 0.193 | `dominated_by_wrong_TE0` |
| 28.00 | 5.1641 | +4.1 | -6.178 | 0.022 | 0.513 | 0.100 | `dominated_by_wrong_TE0` |
| 26.50 | 5.1641 | +4.1 | -5.952 | 0.012 | 0.500 | 0.118 | `dominated_by_wrong_TE0` |
| 25.00 | 5.1641 | +4.1 | -5.779 | 0.005 | 0.487 | 0.126 | `dominated_by_wrong_TE0` |
| 23.06 | 5.1641 | +4.1 | -5.582 | 0.004 | 0.474 | 0.139 | `dominated_by_wrong_TE0` |

## Scientific questions

### Q1 — Intermediate width TE1 suppress without 23 mm penalties?
**Answer: no — TE1 can be suppressed, but not without a worse Band-B worst-case (wrong-TE0 / reflection penalties at ~5.16 GHz)**

Criterion for promising: mean TE1 ≪ 40 mm **and** worst Band-B A2 competitive with 40 mm. Mean-channel cosmetics alone do not count.
- **No** intermediate width meets the promising criterion. Widths with TE1 cutoff above Band B do suppress TE1, but they deepen the ~5.16 GHz family into a wrong-TE0 (or reflection) notch that is *worse* than the original 40 mm worst-case A2.
  - W=28.00 mm: TE1_suppressed=True, worst A2=-6.178 (40 mm baseline -4.964), mean wrong=0.202, mean R=0.361
  - W=26.50 mm: TE1_suppressed=True, worst A2=-5.952 (40 mm baseline -4.964), mean wrong=0.199, mean R=0.381
  - W=25.00 mm: TE1_suppressed=True, worst A2=-5.779 (40 mm baseline -4.964), mean wrong=0.196, mean R=0.384

Note: some intermediate widths improve *mean* Band-B A2 while degrading the *worst* point — the feasibility metric is worst-case.

### Q2 — Do resonance frequencies move smoothly with W?
**Answer: partial — ~5.16 GHz family is nearly pinned in frequency while deepening/reclassifying; ~4.91 GHz family softens and redshifts (tracker may latch onto a shallower neighbor once the TE1 notch collapses)**
- ~4.91: 4.9094 → 4.7625 GHz (Δ = -146.9 MHz); A2 -4.96 → -1.54
- ~5.16: 5.1641 → 5.1641 GHz (Δ = +0.0 MHz); A2 -4.92 → -5.58

### Q3 — Can feed width move the dominant bad resonance outside Band B?
- Family softens/exits: 4.91=True, 5.16=False
- near_4p91 at W=30.00: f=4.7625, A2=-2.01, softened
- near_4p91 at W=28.00: f=4.7625, A2=-1.70, softened
- near_4p91 at W=26.50: f=4.7625, A2=-1.58, softened
- near_4p91 at W=25.00: f=4.7625, A2=-1.57, softened
- near_4p91 at W=23.06: f=4.7625, A2=-1.54, softened

### Q4 — Clear candidate for **worst-case dense A2 over entire Band B**?
**Best fixed-q width: W = 40.00 mm** (worst A2 = -4.964 at 4.9094 GHz).
vs 40 mm: worst A2 = -4.964 at 4.9094 GHz (Δ = +0.000).
vs 23.06 mm: worst A2 = -5.582.
Improves entire-Band-B worst (not merely the old 4.909 min): **False**

### Q5 — Justify a narrowed/intermediate-feed Band-B reoptimization?
**Recommendation: NO**
- No narrowed/intermediate width improves worst Band-B A2 over the original 40 mm feed.
- Fixed-q evidence: narrowing trades the TE1-dominated ~4.91 notch for a deeper wrong-TE0/reflection-dominated ~5.16 notch. A joint W+q reopt would need to *overcome* that conversion, not merely retune around the old 4.91 min. That is a larger ask than a 'small' Band-B calculation is likely to settle.

Do **not** auto-launch that optimization. Stage-6 remains **NOT STARTED**.

## Plots

- `A2_heatmap_f_vs_W.png`
- `A2_overlay_all_widths.png`
- `worst_A2_and_metrics_vs_W.png`
- `resonance_family_tracks.png`
- `TE0_des_overlay.png`
- `TE0_wrong_overlay.png`
- `R_tot_overlay.png`
- `TE1_overlay.png`
- `P_out_overlay.png`
- `TE0_des_heatmap.png`
- `TE0_wrong_heatmap.png`
- `R_tot_heatmap.png`
- `TE1_heatmap.png`

