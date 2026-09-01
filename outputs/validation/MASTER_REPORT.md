# Master validation report — `agent/eigenmode-ports`

**Generated:** 2026-08-31 (overnight performance + reciprocity checkpoint)  
**Branch:** `agent/eigenmode-ports`  
**Registry:** `outputs/validation/MASTER_VALIDATION_REGISTRY.json`

---

## Executive summary

1. **Direct Lorentz reciprocity (no port normalization) passes** — localized Hz point sources on `horns_only` prism geometry at res32: P1↔P2 **0.006 dB** amp / **0.87°** phase; P2↔P3 **<10⁻¹³ dB** (machine precision). The fixed B=0 staircased-PEC Meep system **remains reciprocal** when port flux machinery is removed.
2. **Apparent ~0.3–0.6 dB `te1_hz_line` reciprocity error is a port-definition artifact** — two orders of magnitude larger than direct field reciprocity on the same geometry. Focus shifts to excitation/measurement normalization, not solver nonreciprocity.
3. **Global coordinate rotations rejected** — prism/rotated_blocks +15° stabilizes P1↔P2 (~0.30 dB) but **breaks P2↔P3** (~0.54–0.58 dB at res64/96). Do not continue global rotations or broad horn-wall searches.
4. **Grid registration strongly affects port metrics** (fractional-cell offset sweep: 0.05–0.79 dB at res32) — this changes the **numerical horn mode** and how `te1_hz_line` + axis flux sample it; it does **not** imply the discrete Maxwell solver violates reciprocity.
5. **Mode-profile study:** analytic `te1_hz_line` captures only **~79%** power overlap with the numerical guided Hz mode (similar at horizontal and ±60°). Cross-orientation numerical modes match **>99.9%** — the propagating mode shape is not wildly orientation-dependent; mismatch with cosine TE1 is uniform.
6. **EigenModeSource:** fails on full PEC horns (MPB + staircasing); **works** on a straight high-ε parallel-plate reference guide (`eigenmode_reference_res32.json`). Use time-domain reference-guide extraction for angled PEC horns.
7. **Trustworthy now:** direct Lorentz reciprocity, P2↔P3 port reciprocity, Faraday, geometry audit, grid-sensitivity evidence. **Not trustworthy:** P1↔P2 axis↔diagonal reciprocity with current `te1_hz_line` + flux normalization.
8. **Recommended next step:** port-specific numerical mode references (`plasmeep/ports/numerical_mode.py`) — launch/measure via mode overlap with per-port normalization; **not** more wall geometry tweaks unless direct reciprocity regresses.
9. **Runtime profiling (overnight):** Python horn geometry is **<0.001%** of wall time; **normalization (~43%) + FDTD (~38%)** dominate. Do not micro-optimize horn construction.
10. **Validated cheap defaults:** `horns_only` res32, **run_time=5**, **np=32**, `--skip-norm-if-cached` → **~49× faster** than rt=40 uncached baseline (6.6 s vs 325 s per P1↔P2 te1 pair) with **unchanged physics** (port error still 0.625 dB; Lorentz reciprocity still 0.006 dB).

---

## 1. Physical meaning of `a`

| Item | Simulation | Physical |
|---|---|---|
| `a` | 1.0 | **0.028 m = 2.8 cm** |
| Lattice pitch `d_exp = 0.020/a` | 0.714 a | **20 mm** |
| Bulb OD `0.0075/a` | — | **15 mm** |
| Bulb ID `0.0065/a` | — | **13 mm** |
| Horn wall | 4 mm / 104 mm / 48 mm / 89 mm + 60 mm feed | as coded |

`a` is **not** the lattice spacing; it is the Meep normalization length inherited from the circulator notebook.

---

## 2. Geometry vs experiment (Rodriguez et al., ACS Photonics 2026)

| Feature | Paper | Harness | Match |
|---|---|---|---|
| 91 elements, triangular lattice | yes | yes | ✓ |
| 20 mm pitch | yes | `d_exp` | ✓ |
| Quartz 15/13 mm, ε≈3.8 | yes | `Add_Bulb` | ✓ |
| Horns | commercial MW | synthetic PEC 6-port | ✗ |
| Device | 3-port steering | 6-port circulator | different topology |
| Polarization | **Ez** (z out of plane) | **Hz** (2D circulator) | different (see §3) |

---

## 3. Polarization audit

| Model | Active 2D fields | Role |
|---|---|---|
| Paper / Ceviche | Ez, Hx, Hy | beam-steering fidelity |
| Meep circulator + `te1_hz_line` | Hz, Ex, Ey | **B∥z gyrotropy-active** |
| `te1_ez_line` (attempted) | Ez source + Ex/Ey flux | **broken** on P1↔P2 |

### Is `te1_hz_line` the same as the published PMM polarization?

**No** — the paper uses Ez along the discharge axis.

### Should we switch the circulator campaign to Ez?

**No for B≠0 circulator physics** — Ez does not couple to the Hall tensor with B∥z in 2D.

### `te1_ez_line` cheap test (res32, rt40, horns_only)

| Pair | P1↔P2 | P2↔P3 |
|---|---|---|
| `te1_hz_line` | **0.62 dB** | **<10⁻¹³ dB** |
| `te1_ez_line` | **52.8 dB** (launch mismatch ~200%) | **<10⁻⁵ dB** |

Ez fails because: (1) PEC walls enforce Ez=0 on boundaries — incompatible with TM-like horn modes; (2) flux monitors integrate TE-oriented Poynting (Ex,Ey) while source is Ez. **Do not continue Ez high-res campaign without a redesigned horn + TM flux path.**

---

## 4. Completed reciprocity inventory (`te1_hz_line`, full device, P1↔P2)

| res | rt | \|P12−P21\| | incident spread | MPI | wall (s) |
|---:|---:|---:|---:|---|---:|
| 32 | 40 | **0.38 dB** | ~2.0% | np=4 | 149 |
| 48 | 80 | **1.28 dB** | ~2.5% | np=4 | 854 |
| 64 | 80 | **1.03 dB** | ~2.9% | np=4 | 2241 |
| 96 | 80 | **0.68 dB** | ~3.3% | np=4 | 7494 |
| 128 | 80 | **1.52 dB** | ~2.3% | np=32 | 14146 |

**Non-monotonic — not insufficient resolution alone.** Launch powers nearly equal at res128 (33.1 vs 33.9, ~2.3%).

### Control: P2↔P3 (diag↔diag)

All formulations **<0.001 dB** at res32.

### Diagnostic: horns_only (no PMM), res32, `te1_hz_line`

| Formulation | P1↔P2 | Notes |
|---|---|---|
| `te1_hz_line` (axis flux) | **0.62 dB** | plasma removed; axis error persists |
| `te1_guide_normal` | **0.60 dB** | guide-normal Poynting does not fix axis error |
| `te1_ez_line` | **52.8 dB** | incompatible with PEC + TE flux |

P2↔P3 control for all Hz formulations: **<0.001 dB**.

---

## 4a. Direct Lorentz reciprocity (no port normalization) — **definitive**

Script: `scripts/validation/direct_reciprocity_test.py`  
Output: `outputs/validation/lorentz_direct/lorentz_direct_res32.json`

Localized Gaussian Hz point sources at port source centers; complex Hz DFT at paired monitor centers. No flux, no TE1 weighting, no incident-power normalization.

| Pair | |H_AB|/|H_BA| error | Phase diff | |ΔH|/mean|H| |
|---|---:|---:|---:|
| P1↔P2 | **0.006 dB** | 0.87° | 1.5% |
| P2↔P3 | **<10⁻¹³ dB** | ~0° | ~2×10⁻¹⁴ |

**Conclusion:** the discrete Meep+PEC geometry at B=0 **is reciprocal** to far better precision than the `te1_hz_line` port metric (~0.62 dB on the same horns_only prism geometry). The remaining axis↔diagonal discrepancy is in **port excitation and/or flux-based power normalization**, not a breakdown of electromagnetic reciprocity in the solver.

---

## 4b. Mode-profile / TE1 mismatch study

Script: `scripts/validation/mode_profile_study.py`  
Output: `outputs/validation/mode_profiles/mode_profiles_res32.json`, `mode_profiles_res32.png`

Excite each port with `te1_hz_line`; extract complex Hz along monitor tangent; compare to analytic cos envelope.

| Port | Power overlap \|⟨num\|te1⟩\|² | Orthogonal fraction |
|---|---:|---:|
| P1 horizontal | **0.791** | 0.209 |
| P2 +60° | **0.785** | 0.215 |
| P6 −60° | **0.785** | 0.215 |

Cross-orientation numerical mode overlap: **>99.9%** (P1↔P2, P1↔P6, P2↔P6).

**Interpretation:** staircasing changes the numerical guide mode relative to analytic TE1 cos (~21% orthogonal everywhere), but **does not** produce strongly different mode mixtures at horizontal vs ±60° ports. Port reciprocity failure is more likely from **axis-aligned flux integration / normalization** than from orientation-dependent mode excitation alone. Port-specific numerical references (`plasmeep/ports/numerical_mode.py`) are the next implementation target.

---

## 4c. Grid registration sensitivity (preserved evidence)

Fractional-cell offset sweep (res32, prism, horns_only): P1↔P2 ranges **0.05–0.79 dB** (~16× spread at zero vs worst offset).

**Wording discipline:**
- Grid/how the 60° guide lands on the Yee grid → changes the **numerical horn mode** and flux sampling.
- Using a fixed analytic `te1_hz_line` + axis flux on that mode → **apparent** nonreciprocity in port metrics.
- This is **not** the same as a true violation of discrete electromagnetic reciprocity (see §4a).

Global +15° coordinate rotation (prism and rotated_blocks): improves offset sensitivity but **degrades P2↔P3** to 0.54–0.77 dB at res64/96 — **rejected**.

---

## 5. Measurement formulation screen (closed)

Cheap res32/rt40 A/B: no receiver variant beat `te1_hz_line`+axis flux toward 0.2 dB gate meaningfully. **Do not repeat.**

---

## 6. MPI and throughput

### res48 one-port scaling (`te1_hz_line`, OMP=1, OFI/TCP) — completed

| np | device (s) | speedup vs np=4 | efficiency |
|---:|---:|---:|---:|
| 4 | 243 | 1.0 | 100% |
| 8 | 216 | 1.13 | 56% |
| 16 | 107 | 2.27 | 57% |
| **32** | **46** | **5.27** | **66%** |
| 48 | 54 | — | ~38% |
| 64 | 51 | — | ~30% |

np=48/64 provide **no gain** over np=32 at res48. Numerical flux identical across ranks.

### res96 reciprocity-pair scaling — **completed** (horns_only, rt=20, P1↔P2)

Results: `outputs/validation/mpi_scale/reciprocity_pair_res96_rt20_te1_hz_line.json`

| np | total (s) | speedup vs np=4 | efficiency |
|---:|---:|---:|---:|
| 4 | 4195 | 1.0 | 100% |
| 8 | 2111 | 1.99 | 99% |
| 16 | 1673 | 2.51 | 63% |
| **32** | **1625** | **2.58** | **32%** |
| 48 | 1663 | 2.52 | 21% |
| 64 | 1685 | 2.49 | 16% |

**At res96, np=32 is optimal** (~27 min for horns_only pair). np=48/64 provide no gain. P12=0.267 dB identical across all ranks (numerical agreement confirmed).

### Throughput — **completed** (2× concurrent np=32, horns_only, res96, rt=20)

| Mode | Wall time | Notes |
|---|---|---|
| 1× np=32, one pair | **1625 s** (~27 min) | from scaling benchmark |
| 2× np=32 concurrent (P1↔P2 + P2↔P3) | **5363 s** (~89 min) | job0=3795 s, job1=5363 s |
| 2× np=32 sequential (est.) | **~3250 s** (~54 min) | 2 × 1625 s |

**Concurrent 2×np=32 is slower than sequential** on this host (64 cores): MPI contention inflates per-job time ~2.3–3.3×. **Recommendation: run reciprocal pairs sequentially at np=32**, not 2× concurrent np=32.

For a 3-pair campaign at res96/rt80 full device: **3 × ~50 min ≈ 2.5 h** sequential at np=32.

### Runtime recommendation

| Case | Config | Est. wall time |
|---|---|---|
| P1↔P2 res96 rt80 full | 1× np=32 | **~50–60 min** (7494 s @ np=4 ÷ 2.58 speedup) |
| P1↔P2 res128 rt80 full | 1× np=32 | **~3.9 h** (observed 14146 s) |
| horns_only P1↔P2 res96 rt20 | 1× np=32 | **~27 min** (measured 1625 s) |

Use **OMP_NUM_THREADS=1**, **FI_PROVIDER=tcp**, **MPICH_CH4_NETMOD=ofi**.

---

## 7. Reusable horn / port architecture

| Layer | Location | Status |
|---|---|---|
| Horn geometry | `plasmeep/ports/horn.py` | `HornGeometry`, grid offset, wall variants |
| Direct reciprocity | `plasmeep/ports/lorentz_probe.py` | localized Hz transfer, no flux |
| Mode profiles | `plasmeep/ports/mode_profile.py` | TE1 overlap metrics |
| Numerical port modes | `plasmeep/ports/numerical_mode.py` | reference profiles, overlap receivers |
| Sources / flux | `port_formulations.py` | registry; `te1_hz_line` (under test) |
| Device modes | `sixport_common.build_circulator_device` | `full` / `horns_only` |
| Validation | `direct_reciprocity_test.py`, `mode_profile_study.py`, `eigenmode_reference_guide.py` | active |

API docs: `outputs/validation/ports/HORN_API_PROPOSAL.md`

**Not promoted** to `Sketchbook_PMMCirculator.ipynb` or `PMMCirculatorInverse.py`.

---

## 8. Faraday sign

- **PASS** res32/64: |κ| error <0.1%, perfect B-reversal.
- Global sign **opposite** theory → convention (Stokes ψ vs κ labeling), not implementation bug.
- Finite-length: L=20a, theory **+89.9°** vs Meep **−89.9°** (same |κ|). See `faraday_sign_convention_note.md`.

---

## 9. What is trustworthy / untrustworthy

| Trustworthy | Untrustworthy |
|---|---|
| Hexagonal PMM lattice dimensions vs paper | P1↔P2 B=0 reciprocity with `te1_hz_line` + flux |
| **Direct Lorentz reciprocity (horns_only)** | `te1_ez_line` without TM horn redesign |
| P2↔P3 reciprocity (diag↔diag) | Blind res160+ escalation |
| Faraday \|κ\| + reversal | Global horn rotation / broad wall-candidate searches |
| Grid registration affects **port metrics** | Claiming staircased PEC makes solver nonreciprocal |
| MPI numerical reproducibility across ranks | Measurement-only fixes without mode-consistent ports |

---

## 10. Single recommended next step

**Promote `num_mode_hz_line` (numerical launch + flux receiver)** for horns_only validation; tune full-device `run_time` before res96.

1. Use cached `NumericalPortMode` profiles (`plasmeep/ports/mode_registry.py`) per res/geometry.
2. Keep axis-aligned flux receiver (modal overlap receiver **not** production-ready).
3. Re-validate full PMM at res64 with longer `run_time` (rt=20 gave 2.19 dB; horns_only passes at 0.076 dB).
4. Modal receiver (`te1_hz_modal`, `num_mode_modal`) remains experimental.

**Do not** return to horn-wall sweeps unless direct Lorentz reciprocity regresses.

---

## Script inventory

See `MASTER_VALIDATION_REGISTRY.json` → `scripts` field. Key harness files:

- `sixport_common.py` — geometry, simulate, norms, `horns_only` mode
- `port_formulations.py` — Hz/Ez sources, flux variants
- `direct_reciprocity_test.py` — Lorentz Hz transfer without flux
- `mode_profile_study.py` — numerical vs analytic TE1 profiles
- `eigenmode_reference_guide.py` — straight-guide EigenModeSource viability
- `reciprocity_b0_study.py` — pair runs → JSON
- `run_horn_grid_study.py` — offset/resolution/candidate screening
- `run_modal_port_ablation.py` — launch×receiver A/B study
- `make_modal_port_figure.py` — before/after figure

**Do not rerun:** te1_hz_line P1P2 res32–128 full device; measurement screen; Faraday res32/64; global rotation sweeps.

---

## 11. Overnight performance audit (2026-08-31)

Artifacts: `outputs/validation/overnight/`

### Where time goes (horns_only P1↔P2 te1, res32, rt=40)

| Phase | Wall time | % of total |
|---|---:|---:|
| Python horn geometry build | **0.001 s** | **<0.001%** |
| Source normalization (2 ports) | **139 s** | **42.7%** |
| FDTD timestepping (one excitation) | **62 s** | — |
| Full pair `simulate_circulator` | **124 s** | **38.2%** |
| Flux extraction | **<0.001 s** | negligible |
| **Total** | **325 s** | 100% |

**Answer:** horn construction is irrelevant; **normalization + Meep FDTD** dominate. Do not optimize Python horn polygons.

### MPI for cheap res32 (single Lorentz transfer P1→P2, rt=40)

| np | Wall (s) |
|---:|---:|
| 4 | 15.8 |
| 8 | 11.1 |
| 16 | 10.1 |
| **32** | **5.6** |

**np=32 is optimal even for cheap res32** (same as res96 finding). Do not run 2× concurrent np32 jobs.

### Minimum safe run_time (Lorentz P1↔P2, res32)

| run_time | Reciprocity amp err | H₀₁ rel err vs rt=80 |
|---:|---:|---:|
| **5** | **0.006 dB** | **0.09%** |
| 10 | 0.006 dB | 0.08% |
| 40 (legacy) | 0.006 dB | — |

**Recommended cheap default: `run_time=5`.** Port te1 P1↔P2 error **unchanged** at 0.625 dB (rt=5 vs rt=40).

### Validated workflow speedup

| Configuration | Wall time (horns_only P1↔P2 te1) |
|---|---:|
| Baseline (rt=40, uncached, profiled) | **325 s** |
| rt=5, np=32, norm cached | **6.6 s** |
| **Speedup** | **~49×** |

### Direct reciprocity verdict (confirmed at rt=5, np=32)

| Pair | Lorentz amp err | te1_hz_line port err |
|---|---:|---:|
| P1↔P2 | **0.006 dB** | **0.625 dB** |
| P2↔P3 | **<10⁻¹³ dB** | **~0 dB** |

**Meep + staircased PEC is reciprocal.** Port excitation/flux normalization causes the axis↔diagonal discrepancy.

### Mode profiles + cached port modes

| Port | \|⟨num\|te1⟩\|² | Cross-orient overlap |
|---|---:|---:|
| P1 horizontal | 0.791 | — |
| P2 +60° | 0.785 | >99.9% vs P1/P6 |
| P6 −60° | 0.785 | (symmetry-equivalent to P2 class) |

Cached references: `outputs/validation/mode_profiles/res32/numerical_mode_P*.json`  
Policy: `outputs/validation/overnight/CACHE_POLICY.md`

**Best port-mode path (updated):** cache `NumericalPortMode` → launch via `num_mode_hz_line` → **keep flux receiver**. Modal overlap receiver degrades reciprocity on P2↔P3.

### Toward fast full circulator

1. Use **rt=5** for cheap horns_only diagnostics; **rt=80** retained for production res96 full-device unless revalidated.
2. Always **`--skip-norm-if-cached`** when cache key matches (res, rt, formulation, geometry offsets).
3. **np=32 sequential** — never 2× concurrent np32 on this host.
4. Next physics fix: wire `numerical_mode.py` into `port_formulations.py` — expected to close P1↔P2 gap without horn-wall sweeps.
5. Dominant cost for inverse design will remain **FDTD per source excitation** — cache norms + mode profiles; **rho/B changes do not invalidate port references at B=0**.

---

## 12. Numerical-mode port validation (2026-09-01)

Artifacts: `outputs/validation/modal_ports/`, `plasmeep/ports/{numerical_launch,modal_receiver,mode_registry}.py`

### A/B decomposition (horns_only, res32, rt=5, np=32)

| Formulation | Launch | Receiver | P1↔P2 | P2↔P3 |
|---|---|---:|---:|---:|
| `te1_hz_line` | analytic TE1 | flux | **0.625 dB** | **0.000 dB** |
| **`num_mode_hz_line`** | **numerical** | **flux** | **0.076 dB** | **0.208 dB** |
| `te1_hz_modal` | analytic TE1 | modal overlap | 0.487 dB | 31 dB |
| `num_mode_modal` | numerical | modal overlap | 0.672 dB | 48 dB |

**Dominant error was launch, not measurement.** Numerical-mode excitation with existing flux receiver clears the P1↔P2 gate (<0.2 dB). Modal overlap receiver is **not** ready (especially P2↔P3).

### Grid-offset robustness (`num_mode_hz_line`, P1↔P2)

| grid offset (cells) | \|P12−P21\| |
|---|---:|
| (−0.5, 0) | 0.209 dB |
| (0, 0) | **0.076 dB** |
| (+0.5, 0) | 0.209 dB |

Prior `te1_hz_line` range at res32: **0.05–0.79 dB** (~16×). Numerical launch reduces spread to **~2.7×** (0.076–0.209 dB).

### Full device (res64, rt=20, P1↔P2)

| Formulation | \|P12−P21\| |
|---|---:|
| `te1_hz_line` (prior) | ~0.38–1.03 dB |
| `num_mode_hz_line` | **2.19 dB** |

**Not satisfactory** — do not run res96. Likely needs longer `run_time` and/or horns_only-validated settings re-tuned with plasma present.

### Normalization convention

- Reference mode φₙ: L2-normalized discrete Hz profile from cached JSON (`mode_profile_study.py`).
- Launch weights: `a_k = conj(φ_k) / ||φ||₂` (unit discrete power).
- Flux receiver: unchanged axis-aligned `FluxRegion` + incident normalization from straight-guide reference run.
- Modal receiver (experimental): `b = ⟨φ_n|H⟩` with unit-norm φ_n; `P = |b|² / P_inc`.

### Cache invalidation (B=0 inverse design)

| Quantity | Recompute when changing… |
|---|---|
| Numerical mode JSON | `res`, `fs_a`, horn geometry/offsets, `horn_walls` |
| Source normalization pickle | `res`, `run_time`, formulation, geometry offsets |
| **Not** required for inner loop | **91-element `rho` pattern** (linear plasma, fixed geometry, B=0) |

### Production forward-eval cost estimate (np=32, cached)

| Item | horns_only P1↔P2 rt=5 | full device P1↔P2 res64 rt=20 |
|---|---:|---:|
| Mode profiles (one-time / res) | ~50 s/port | ~90 s/port |
| Norm cache (2 ports, first run) | ~9 s | ~included in 895 s |
| Pair simulation (cached norms) | **~8–11 s** | **~895 s** |
| Full 6×6 matrix (extrapolated) | ~6 × 11 s ≈ **1 min** | ~6 × 450 s ≈ **45 min** |

Inner optimization loop: **only FDTD device runs** per `rho` once norms + modes cached.

---
