# Horn continuum comparison

**The 0.24 dB gap at 25 ppc is largely Meep’s finite grid. A leftover offset below about 0.11 dB is unresolved. Quartz stays out.**

Six-horn, no bulbs, source P1, `num_mode_guide_normal`, df = 0.10 fs, run_time = 20, prism PEC. FEM is the converged polygon model (excised Neumann metal, Meep PML, mirrored constrained mesh). Error in the tables below is Meep minus FEM, in dB.

## 1. Complex-field paths

One complex scale was fixed on the P1 feed, within 0.35 a of the source, and applied to every cut. FEM-VH against the saved 25 ppc Meep Hz grid.

Plot: `horn_continuum/path_phase.png`. Curves: `horn_continuum/path_*.npz`.

| Cut | Complex rel. L2 | Magnitude rel. L2 | Mean phase (deg) | Where it leaves the feed |
|---|---:|---:|---:|---|
| P1 axis, source → throat → center → P4 | 0.088 | 0.031 | −3.9 | 6.1 a from the source, past the throat |
| +60° through the cavity, toward P2 | 0.087 | 0.043 | −3.6 | near the far end of the cut |
| −60° mirror, toward P6 | 0.087 | 0.043 | −3.6 | same, to 10⁻¹² |
| Transverse line x = 0, \|y\| ≤ 6 | 0.111 | 0.053 | −5.3 | the whole central cut sits at this level |

Along the P1 axis the error changes by region:

| Segment | Complex L2 | Magnitude L2 | Mean phase (deg) |
|---|---:|---:|---:|
| Feed, source to throat | 0.035 | 0.023 | −0.9 |
| Just past the throat | 0.093 | 0.041 | −2.7 |
| Central cavity | 0.089 | 0.033 | −4.0 |
| P4 feed | 0.116 | 0.019 | −6.4 |

The feed matches. Magnitude error stays near 2–4% after the throat. The phase offset grows by about 5° from the feed to the P4 feed and then sits there. A single phase slope fitted through the whole axis is not a wavenumber: standing-wave nulls make that fit jump. On the central-cavity segment the two phase slopes agree to 0.003 a⁻¹. The +60° and −60° cuts are identical, so this is not a left/right bug.

Reading of the port error: scattering amplitude at the few-percent level, plus a few degrees of phase picked up crossing the open cavity. The source region and the PML ends are not where the solutions separate. Free-space k, below, is far too close to explain a 5° path error by numerical dispersion.

## 2. Free-space phase constant

Homogeneous air, point Hz source, phase fit on the outgoing cylindrical wave at fs. k0 = 1.613801 a⁻¹. FEM uses a structured mesh with h = 0.04 a. Meep uses resolution 50 (25 ppc).

| Solver | Direction | k_num | k_num/k0 − 1 |
|---|---|---:|---:|
| continuum k0 | — | 1.613801 | 0 |
| FEM | axial | 1.614250 | +2.8×10⁻⁴ |
| FEM | 45° | 1.614929 | +7.0×10⁻⁴ |
| Meep 25 ppc | axial | 1.614590 | +4.9×10⁻⁴ |
| Meep 25 ppc | 45° | 1.615335 | +9.5×10⁻⁴ |

Meep minus FEM is about 2×10⁻⁴ of k0 on axis and 2.5×10⁻⁴ on the diagonal. Over a 15 a path that is ~0.3° of phase. Both solvers sit slightly above k0, which is the expected Hankel-fit bias plus a small grid dispersion. This does not account for the horn-path phase or the 0.24 dB port gap.

## 3. Meep resolution table

Same horns-only problem. 25 ppc is the existing reference. 35 and 50 ppc were run here. Incident power is the cached straight-guide value at the same resolution (P_inc = 2909.72, 2909.79, 2910.05). No compatible 35/50 ppc horns-only six-port result was already on disk.

35 and 50 ppc were run with Meep subpixel averaging off. At 25 ppc, averaging on and off produced the same powers to the printed digits, so that switch is not a new geometry.

| | dx (mm) | res | P2 | P3 | P4 | P5 | P6 | Wall (s) | Ranks |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Meep 25 ppc | 0.400 | 50 | 0.037632 | 0.085539 | 0.696684 | 0.085539 | 0.037632 | 182 | 16 |
| Meep 35 ppc | 0.286 | 70 | 0.037223 | 0.085367 | 0.697043 | 0.085367 | 0.037223 | 915 | 32 |
| Meep 50 ppc | 0.200 | 100 | 0.036743 | 0.084435 | 0.694668 | 0.084435 | 0.036743 | 2493 | 48 |
| FEM-VH | — | — | 0.035613 | 0.082835 | 0.700743 | 0.082835 | 0.035613 | 47 (six-port) | 1 |

P2 = P6 and P3 = P5 at every Meep resolution, to the printed digits.

## 4. FEM-H+ / FEM-VH

| Level | DOFs | P2 | P3 | P4 | P2−P6 (dB) |
|---|---:|---:|---:|---:|---:|
| FEM-H+ | 757952 | 0.035586 | 0.082818 | 0.701322 | ~0 |
| FEM-VH | 1308537 | 0.035613 | 0.082835 | 0.700743 | 6×10⁻¹³ |

H+ minus VH, in dB: P2 −0.0033, P3 −0.0009, P4 +0.0036. No finer FEM mesh was run. The side-port change is already far below the Meep resolution trend.

## 5. Error versus resolution

Error = 10 log10(P_Meep / P_FEM-VH).

| Port | 25 ppc | 35 ppc | 50 ppc |
|---|---:|---:|---:|
| P2, P6 | +0.239 | +0.192 | +0.136 |
| P3, P5 | +0.140 | +0.131 | +0.083 |
| P4 | −0.025 | −0.023 | −0.038 |

Plot: `horn_continuum/ports_vs_resolution.png`.

## 6. Direction

| Port | As Meep is refined |
|---|---|
| P2, P6 | Moves toward FEM-VH. 0.239 → 0.192 → 0.136 dB. Still above 0.10 dB at 50 ppc. |
| P3, P5 | Moves toward FEM-VH. 0.140 → 0.131 → 0.083 dB. Under 0.10 dB at 50 ppc. |
| P4 | Stays within 0.04 dB of FEM-VH. The 50 ppc point is 0.015 dB farther than 25 ppc. |

## 7. Extrapolation

Three points. A straight line in dx and a straight line in dx² both fit, and they do not agree on the intercept.

| Port | Intercept if error ∝ ends as a + b·dx | Intercept if a + b·dx² |
|---|---:|---:|
| P2, P6 | +0.037 dB | +0.111 dB |
| P3, P5 | +0.038 dB | +0.078 dB |
| P4 | −0.046 dB | −0.037 dB |

The linear-dx fit has the smaller residual, but three points and two parameters cannot choose the model. No single continuum value is reported.

60 ppc was not run. At resolution 120 the two models predict P2 errors of about 0.12 dB and 0.13 dB. That gap is too small to decide which intercept is real. 75 and 100 ppc were not run.

## 8. Uncertainties

FEM discretization uncertainty, taken as the H+ to VH change: ≤ 0.004 dB on P2–P6.

Meep discretization uncertainty is the unresolved part of the trend. At 50 ppc the P2 error is still 0.14 dB and still falling. The two extrapolations place the infinite-resolution P2 offset anywhere from 0.04 dB to 0.11 dB. That spread is the Meep continuum uncertainty on the side ports. P4 is already inside 0.04 dB.

## 9. What the 0.24 dB was

Most of the 25 ppc side-port gap is Meep’s finite grid: refining Meep moves P2 and P3 toward the converged FEM value, and free-space k agrees to a few parts in 10⁴. The field cuts put the remaining difference in the open-cavity scattering, not in the source, the PML, or a wrong phase constant.

A formulation offset of order 0.05–0.11 dB is still allowed by the extrapolation, so the continuum limit is not closed. It is not a locked 0.24 dB operator error.

## 10. Quartz

Do not add quartz yet. P2 at 50 ppc is still 0.14 dB from FEM-VH, and the continuum intercept is model-dependent. The next Meep point that would actually separate a 0.04 dB intercept from a 0.11 dB intercept is finer than 60 ppc, which this task does not run.
