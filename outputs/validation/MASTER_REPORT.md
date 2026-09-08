# Master validation report — `agent/eigenmode-ports`

**Generated:** 2026-09-03 (physical units cleanup + 50 points/cm port validation)  
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
| `a` | 1.0 | **0.020 m = 2.0 cm (lattice pitch; was 0.028 m until §19)** |
| Lattice pitch `d_exp = 0.020/a` | **1.0 a** | **20 mm** |
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

**Complete full-device res64 runtime convergence** for corrected `num_mode_hz_line`. Horns_only gates now pass (P1↔P2 0.076 dB, P2↔P3 ~0 dB after per-port profile fix). **Do not promote as default** until full-device reciprocity is acceptable at demonstrated-converged `run_time`.

1. Cached profiles: one JSON per port (`audit_mode_cache.py` verifies res + tangent alignment).
2. Keep flux receiver (modal overlap remains experimental).
3. Full-device rt convergence: `run_full_device_rt_convergence.py` (prior rt=20 → 2.19 dB likely under-converged).

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
| `te1_hz_line` | analytic TE1 | flux | **0.625 dB** | **~0 dB** |
| **`num_mode_hz_line`** (after per-port fix) | **numerical** | **flux** | **0.076 dB** | **~0 dB** |
| `te1_hz_modal` | analytic TE1 | modal overlap | 0.487 dB | 31 dB |
| `num_mode_modal` | numerical | modal overlap | 0.672 dB | 48 dB |

**Dominant error was launch, not measurement.** Modal overlap receiver remains experimental.

### P2↔P3 regression root cause (fixed 2026-09-02)

`canonical_profile_port` incorrectly mapped **P3 → P2 profile** (and P5 → P6, P4 → P1), assuming ±60° symmetry classes share one JSON file. On the hexagon:

| Port | angle | old profile | tangent_dot | result |
|---|---:|---|---:|---|
| P2 | +60° | P2 | 1.0 | OK |
| P3 | +120° | P2 (wrong) | **0.5** | **0.208 dB** reciprocity error |
| P6 | −60° | P6 | 1.0 | OK |

**Fix:** one cached profile per port (`numerical_mode_P{1..6}.json`); `audit_mode_cache.py` logs file path + `res_match` + `tangent_dot` at launch. P3/P4/P5 profiles extracted at res32 and res64.

### Cache audit (res64 full-device)

All six ports load `outputs/validation/mode_profiles/res64/numerical_mode_P*.json` with `stored_res=64`, `tangent_dot=1.0`. **No res32 interpolation.** Prior full-device P1↔P2 result (2.19 dB at rt=20) used correct P1/P2 references; failure is **not** a cache-resolution mismatch.

### Grid-offset robustness (`num_mode_hz_line`, P1↔P2)

| grid offset (cells) | \|P12−P21\| |
|---|---:|
| (−0.5, 0) | 0.209 dB |
| (0, 0) | **0.076 dB** |
| (+0.5, 0) | 0.209 dB |

Prior `te1_hz_line` range at res32: **0.05–0.79 dB** (~16×). Numerical launch reduces spread to **~2.7×**.

### Full device (res64, P1↔P2) — runtime convergence **STOPPED**

| run_time | `num_mode_hz_line` \|P12−P21\| | wall (s) | MPI |
|---:|---:|---:|---|
| 10 | 3.46 dB | 2992 | serial (np=1) |
| 20 | **2.19 dB** | 1914 | serial |
| 40 | 5.64 dB | 3180 | serial |
| 60/80 | *not run* | — | campaign stopped |

**Runtime does not converge** toward horns_only (0.076 dB). Non-monotonic behavior rules out "just needs longer rt". Early serial runs lacked `mpirun -np 32`; orchestrator fixed in `mpi_runner.py` — **do not re-run rt60/80** to chase convergence.

**Next diagnostic:** direct Lorentz reciprocity on full B=0 device (uniform rho, res64, rt=20, np=32). **Result: 1.48 dB amp error** (`lorentz_direct_full_res64_rt20_P1P2_mpi.json`) vs horns_only 0.006 dB. `num_mode_hz_line` port metric: 2.19 dB (serial). Full plasma geometry degrades field reciprocity; port formulation adds ~0.7 dB more. **Not purely a port-normalization artifact** like horns_only.

### Normalization convention

- Reference mode φₙ: L2-normalized discrete Hz profile from cached JSON (`mode_profile_study.py`).
- Launch weights: `a_k = conj(φ_k) / ||φ||₂` (unit discrete power).
- Flux receiver: unchanged axis-aligned `FluxRegion` + incident normalization from straight-guide reference run.
- Modal receiver (experimental): `b = ⟨φ_n|H⟩` with unit-norm φ_n; `P = |b|² / P_inc`.

### Cache invalidation (B=0 inverse design)

| Quantity | Recompute when changing… | **Not** invalidated by `rho` alone |
|---|---|---|
| Numerical mode JSON | `res`, `fs_a`, horn geometry/offsets, `horn_walls`, **per-port orientation** | ✓ at B=0 |
| Source normalization pickle | `res`, `run_time`, formulation, geometry offsets | ✓ at B=0 |
| Device FDTD response | every `rho` / B evaluation | — |

Inner optimization loop: **only FDTD device runs** per `rho` once norms + modes cached.

### Production forward-eval cost estimate (np=32, cached)

| Item | horns_only P1↔P2 rt=5 | full device P1↔P2 res64 |
|---|---:|---:|
| Mode profiles (one-time / res, 6 ports) | ~50 s/port | ~90 s/port |
| Norm cache (2 ports, first run) | ~9 s | ~270 s (rt=10) |
| Pair simulation (cached norms) | **~8–11 s** | **~450–900 s** (rt-dependent) |
---

## 13. Overnight direct-Lorentz diagnostic (2026-09-02)

Artifacts: `outputs/validation/lorentz_overnight/` (`overnight_summary.json`, `comparison_table.json`)

Harness: `run_lorentz_overnight_campaign.py`, `direct_reciprocity_test.py`, `mpi_runner.py` (np=32 verified)

### Phase 1 — Test audit

| Item | Value |
|---|---|
| Source | `mp.Source` + `GaussianSource`, **Hz**, amplitude **1.0** (real) |
| Receiver | complex **Hz** DFT at monitor point, frequency `fs_a` |
| Reciprocity relation | **H_AB = H_BA** (passive linear reciprocal TE, B=0) |
| Timing | `until_after_sources=run_time` (fixed, no decay stop) |
| Fields large enough? | **Yes** at res32 (\|H\|≈0.97–1.5); res64 \|H\|≈0.09–0.11 still >1e-3 |

### Answers to morning questions

1. **Is 1.48 dB real or test artifact?** **Real** — observable is documented; fields are not near null at res32. At res64 fields are smaller but still above threshold; error is larger (1.48 dB), not a ratio-of-noise artifact.

2. **Does direct reciprocity converge with runtime?** **No** (res32, full active): 0.445 → 0.452 → 0.553 → 0.606 dB for rt=5/10/20/40. Stabilizes ~0.45–0.6 dB, not 0.006 dB.

3. **Where does failure first appear?** **91-bulb quartz geometry alone → 0.13 dB**. Active Drude 1–7 bulbs → ~0.07–0.08 dB. **91-bulb active Drude → 0.55 dB** (res32). Partial-array tests are position-dependent (non-monotonic with index-ordered bulb count).

4. **Does quartz/simple dielectric pass?** **Mostly** — geometry_only 0.13 dB; dielectric_fill 0.085 dB at 91 bulbs.

5. **Does scalar Drude pass?** **Small arrays yes** (~0.08 dB for 1 bulb); **full 91-bulb array no** (0.55 dB res32).

6. **Does GyrotropicDrude(B=0) pass?** **Identical to scalar Drude** — forced gyrotropic path gives bit-identical results at 1 and 91 bulbs. **Not a GyrotropicDrude B=0 implementation bug.**

7. **Scalar vs gyrotropic difference?** **None observed** — at B=0 library already uses `DrudeSusceptibility`; forced gyrotropic matches exactly.

8. **Fields numerically meaningful?** **Yes** at res32. Res64 \|H\|≈0.1 — monitorable but lower dynamic range.

9. **Smallest reproducible failing case?** **91-bulb active plasma array** at res32 (0.55 dB). Geometry-only 91 bulbs is milder (0.13 dB).

10. **Supported hypothesis?** Cumulative **dispersive 91-element PMM array** breaks approximate reciprocity of the localized Hz Green's function test at finite resolution — **not** port normalization, **not** GyrotropicDrude B=0, **not** insufficient runtime. res64 amplifies error (1.48 dB). Port metric adds ~0.7 dB more on top.

11. **Fastest validated config?** **res32, np=32, rt=20, P1↔P2 direct Lorentz** (~10–15 s/case). Full ladder completes in ~15 min.

### Material-complexity ladder (res32, rt=20, np=32)

| Case | amp err (dB) | \|H_AB\| | \|H_BA\| |
|---|---:|---:|---:|
| horns_only | **0.006** | 0.966 | 0.965 |
| geometry_91 (wp=0) | 0.129 | 1.227 | 1.209 |
| dielectric_91 | 0.085 | 1.440 | 1.426 |
| Drude 1 bulb | 0.079 | 1.095 | 1.085 |
| Drude 7 bulbs | 0.074 | 1.523 | 1.511 |
| Drude 91 bulbs | **0.553** | 1.419 | 1.332 |
| Gyrotropic B=0 91 | **0.553** (same) | — | — |
| full active **res64** | **1.476** | 0.109 | 0.092 |

### Decomposition (unchanged)

| Layer | P1↔P2 error |
|---|---:|
| Direct Lorentz (full res64) | ~1.5 dB |
| `num_mode_hz_line` port metric | ~2.2 dB |
| Port-only excess | ~0.7 dB |

**Do not promote `num_mode_hz_line`.** §14 supersedes the §13 root-cause conclusion: point-probe failure is a sampling artifact, not FDTD nonreciprocity. Horns_only numerical launch (0.076 dB) remains valid for empty-horn validation only.

---

## 14. Matched discrete-overlap reciprocity (2026-09-02)

Artifacts: `outputs/validation/discrete_reciprocity/` (`campaign_summary.json`, `comparison_table.json`)

Harness: `run_discrete_reciprocity_campaign.py`, extended `direct_reciprocity_test.py`, `lorentz_probe.py` (`evaluate_reciprocity_pair`, `build_hz_grid_patch`)

### Method

Replaces point-source → point-probe with **matched 3×3 Hz Yee-grid patches** (half_width=1, uniform weights):

- Sources: `mp.Source` per patch DOF at exact cell centers `(i+0.5)/res`
- Receiver: block DFT weighted overlap on the same DOFs (no continuous interpolation)

### Point vs matched discrete (P1↔P2, rt=20, np=32)

| Geometry / material | Point test | Matched discrete test |
|---|---:|---:|
| horns_only (res32) | **0.006 dB** | **~0 dB** |
| quartz 91 bulbs (res32) | 0.129 dB | **~0 dB** |
| scalar Drude 7 bulbs (res32) | 0.074 dB | **~0 dB** |
| scalar Drude 91 bulbs (res32) | **0.553 dB** | **~0 dB** |
| scalar Drude 91 bulbs (**res64**) | **1.476 dB** | **~0 dB** |

### Bulb-count scaling (point vs discrete, res32)

| n bulbs | Point (dB) | Discrete (dB) |
|---:|---:|---:|
| 0 (horns) | 0.006 | ~0 |
| 1 | 0.079 | ~0 |
| 7 | 0.074 | ~0 |
| 19 | 0.503 | ~0 |
| 37 | −0.163 | ~0 |
| 91 | **0.553** | **~0** |

Point error grows erratically with scatterer count; discrete error stays at machine precision throughout.

### Grid-location audit

Legacy point probes are **off the Hz Yee grid** (fractional offsets up to **±0.5 cells**). Example P1 (91-bulb case):

| Site | fractional offset (cells) |
|---|---|
| P1 source | (−0.24, −0.50) |
| P1 monitor | (+0.33, −0.50) |
| P2 source | (+0.13, +0.29) |
| P2 monitor | (+0.42, −0.46) |

Meep interpolates continuous source/monitor coordinates; in a 91-element scattering environment this breaks the symmetric discrete Green-function test even though the underlying FDTD system remains reciprocal.

### Conclusion

**Does the 91-element B=0 Drude system violate discrete reciprocity?** **No.**

**Was the ~0.55–1.48 dB point-test failure real?** **No** — it was caused by **how the reciprocal response was sampled** (off-grid point interpolation amplified by multiple scattering), not by GyrotropicDrude, scalar Drude physics, or dispersive FDTD nonreciprocity at B=0.

### Next steps

1. Use `evaluate_reciprocity_pair(..., observable='discrete')` as the canonical direct reciprocity check.
2. Revisit full-device `num_mode_hz_line` port development — the ~1.5 dB “fundamental” direct-Lorentz term is removed.
3. Port receivers should likewise be anchored to Yee-grid DOFs where possible.

---

## 15. Six-port formulation gates (2026-09-02)

Artifacts: `outputs/validation/port_gates/` (`campaign_summary.json`, `grid_offset_P0P1_res32.json`)

Harness: `run_port_gates_campaign.py`, `run_grid_offset_port_gates.py`, `reciprocity_b0_study.py --discrete-control`

Direct-Lorentz / material investigation **closed** (§14). Legacy point-probe metric retained only as documented sampling-artifact lesson.

### P2↔P3 regression — resolved

Per-port numerical mode caches (no symmetry sharing) give **P2↔P3 ≈ 0 dB** for all `num_mode_*` formulations at horns_only res32/rt5/np32. The prior 0.208 dB failure was incorrect P3←P2 profile sharing.

### Horns_only A/B decomposition (res32, rt=5, np=32, zero grid offset)

| Launch + receiver | P1↔P2 | P2↔P3 | discrete FDTD control |
|---|---:|---:|---:|
| te1_hz_line + axis flux | 0.625 dB | ~0 dB | ~0 dB |
| num_mode_hz_line + axis flux | 0.076 dB | ~0 dB | ~0 dB |
| **num_mode_guide_normal** + guide-normal flux | **0.061 dB** | ~0 dB | ~0 dB |
| num_mode_yee_sdotn + Yee S·n integration | 0.070 dB | ~0 dB | ~0 dB |

**Gates at zero offset:** P1↔P2 ≤0.2 dB ✓, P2↔P3 ≤0.05 dB ✓.

**Best candidate:** `num_mode_guide_normal` — orientation-consistent guide-normal flux fixes axis↔diagonal asymmetry vs axis-aligned `num_mode_hz_line`; marginally better P1↔P2 than hz_line.

### Grid-offset robustness (P1↔P2, ±0.5 cells) — **not yet passing**

| offset (cells) | num_mode_hz_line | num_mode_guide_normal | discrete control |
|---|---:|---:|---:|
| (0, 0) | 0.076 | **0.061** | ~0 |
| (±0.5, 0) | 0.209 | 0.197 | ~0 |
| (0, ±0.5) | 0.85–0.87 | 0.93–0.96 | ~0 |
| (±0.5, ±0.5) | 1.02–1.05 | 1.07–1.10 | ~0 |

Port metric remains grid-sensitive; **discrete FDTD reciprocity stays at machine precision at all offsets** — confirming errors are port sampling/integration, not physics.

### Full device (res32, rt=20, np=32, uniform rho, B=0)

| Pair | `num_mode_guide_normal` port | discrete FDTD control |
|---|---:|---:|
| P1↔P2 | **0.166 dB** | ~0 dB |
| P2↔P3 | **0.003 dB** | ~0 dB |

Major improvement vs pre-fix full-device `num_mode_hz_line` (~2.19 dB at res64/rt20 with contaminated axis flux). Port residual on P1↔P2 remains; FDTD layer is clean.

### Cache architecture (inverse-design prep)

| Artifact | Valid when rho changes? | Valid when B/res/grid offset changes? |
|---|---|---|
| Horn geometry | yes | per res/offset/rotation |
| Numerical mode JSON (per port) | yes | per res/frequency/offset |
| Incident normalization pickle | yes | per formulation/offset/rt |
| Yee receiver weights (future) | yes | per res/offset |
| PMM FDTD solve | **no** — must rerun | per rho/B/res |

### Promotion status

**Do not promote to production yet.** Horns_only zero-offset gates pass; grid-offset sensitivity and full-device P1↔P2 (~0.17 dB) need further Yee-anchored launch/receiver work before res64 confirmation and runtime sweep.

### Next steps

1. Yee-snap numerical launch + guide-normal flux jointly (reduce grid-offset spread).
2. Full-device rt convergence at res32 with `num_mode_guide_normal`.
3. res64 confirmation once res32 port metric is satisfactory.

---

## 16. Six-hour port block (2026-09-02)

Artifacts: `outputs/validation/port_block/` (`block_summary.json`, `followup_summary.json`)

Harness: `run_port_block_campaign.py`, `num_mode_yee_guide_normal`, `audit_mode_symmetry.py`

### P2↔P3 symmetry — confirmed fixed

`mode_symmetry_audit.json`: P2 and P3 use **separate** cached profiles (`tangent_dot=1.0` each). horns_only P2↔P3 ≈ **0 dB** for all `num_mode_*` formulations.

### Best formulations

| Role | Formulation | Notes |
|---|---|---|
| **Source** | per-port numerical mode (`make_numerical_mode_sources`) | 64-sample complex Hz profiles |
| **Receiver (primary)** | **guide-normal flux** (`num_mode_guide_normal`) | orientation-consistent ∫S·n̂ |
| **Receiver (grid-x robust)** | Yee-snapped guide-normal (`num_mode_yee_guide_normal`) | launch + flux DOFs on Hz Yee grid |

### Horns_only gates (res32, rt=5, np=32)

| Formulation | P1↔P2 | P2↔P3 | discrete FDTD |
|---|---:|---:|---:|
| num_mode_guide_normal | **0.061 dB** | ~0 | ~0 |
| num_mode_yee_guide_normal | 0.105 dB | ~0 | ~0 |

### Grid-offset P1↔P2 (±0.5 cells)

| offset | guide_normal | yee_guide_normal |
|---|---:|---:|
| (±0.5, 0) | 0.197 dB | **0.115 dB** |
| (0, ±0.5) | 0.93–0.96 dB | 0.80–0.83 dB |
| discrete control (all) | ~0 | ~0 |

Yee snapping **halves** x-offset port error but does not fix y-offset sensitivity.

### Full device (B=0, uniform rho, np=32)

**`num_mode_guide_normal` rt sweep (res32):**

| rt | P1↔P2 port |
|---:|---:|
| 5 | 0.301 dB |
| 10 | 0.314 dB |
| **20** | **0.166 dB** |
| 40 | 0.224 dB |

**Minimum safe runtime: rt=20** (non-monotonic above rt=20).

| Test | port reciprocity | discrete FDTD |
|---|---:|---:|
| res32 rt=20 P1↔P2 | **0.166 dB** ✓ | ~0 |
| res32 rt=20 P2↔P3 | **0.003 dB** ✓ | ~0 |
| **res64 rt=20 P1↔P2** | **2.90 dB** ✗ | ~0 |

Discrete reciprocity passes at all resolutions; **res64 port metric regression is a measurement/discretization issue**, not plasma physics.

### MPI wall time (np=32, cached norm, res32, rt=20)

| Item | Time |
|---|---:|
| 2-port reciprocity study (norm+device+discrete) | **~25 s** |
| Device solve only (2×2 ports) | ~12 s |
| **Per source excitation (est.)** | **~6 s** |
| 6-source forward matrix (est.) | **~36 s** |

**Cached (rho-invariant):** numerical modes, horn geometry, incident normalization, flux quadrature weights.

**Per rho evaluation:** PMM FDTD solve + flux extraction only (~6 s/source at res32).

**Est. optimization throughput:** ~100–140 forward evals/hour (3–4 active ports) to ~60/hour (full 6×6), excluding res64.

### Production promotion verdict

**NOT ready for production promotion.**

Blockers:
1. res64 full-device P1↔P2 port reciprocity **2.9 dB** (discrete still ~0).
2. Grid **y-offset** ±0.5 cell still ~0.8–1.0 dB port error.
3. ~4.5% incident-power mismatch between P1/P2 on full device.

**Acceptable for:** res32 horns validation, res32 full-device B=0 preliminary optimization prototyping with `num_mode_guide_normal` + discrete FDTD control.

---

## 17. res64 axis-port diagnosis (2026-09-02 continuation)

Artifacts: `outputs/validation/port_block/res64_horns_summary.json`, `full_num_mode_yee_sdotn_P0P1_g0_0_res64_rt*.json`, `run_res64_followup.py`

### horns_only res64 gates (rt=5, np=32)

| Formulation | P1↔P2 | P2↔P3 | discrete FDTD |
|---|---:|---:|---:|
| num_mode_guide_normal | **0.248 dB** | ~0 | ~0 |
| num_mode_yee_sdotn | **0.251 dB** | ~0 | ~0 |

Horns_only reciprocity is acceptable at res64; axis pair is marginally above the 0.2 dB gate. P1/P2 incident-power mismatch ~**7.8%** (axis vs +60° ports).

### Full device res64 — axis P1↔P2 fails for all receivers

| Formulation | rt | port P1↔P2 | discrete FDTD | flux subtraction |
|---|---:|---:|---:|---|
| num_mode_guide_normal | 20 | **2.90 dB** | ~0 | **broken** (positive diagonal) |
| num_mode_yee_sdotn | 20 | **2.89 dB** | ~0 | N/A (DFT S·n) |
| num_mode_yee_sdotn | 40 | **3.49 dB** | ~0 | — |
| num_mode_yee_sdotn | 60 | **4.02 dB** | ~0 | — |

**Key findings:**

1. **`load_minus_flux_data` fails on full device at res64** for guide-normal flux (diagonal reflection entries go positive). horns_only and res32 full device subtract correctly.
2. **Switching to DFT S·n (`num_mode_yee_sdotn`) does not fix reciprocity** — same ~2.9 dB axis error. The failure is not receiver-specific.
3. **Longer runtime makes res64 worse** (non-convergence in the port metric sense; discrete control stays ~0).
4. **Discrete Hz-patch reciprocity remains machine precision** at all res64 full-device points tested.
5. **Measured port transmission at res64 full (~0.3%) is ~35× smaller than discrete |g| (~11.5%)** — modal launch/measurement severely under-couples on axis with PMM present at res64.
6. **Diagonal pair P2↔P3 at res64 full device passes** (`full_num_mode_guide_normal_P1P2_g0_0_res64_rt20`: **0.0006 dB**). Failure is **axis-specific** (P1↔P2), not general res64.

### Full device res32 — yee_sdotn confirms guide_normal

| Formulation | rt | P1↔P2 port | discrete FDTD |
|---|---:|---:|---:|
| num_mode_guide_normal | 20 | **0.166 dB** | ~0 |
| num_mode_yee_sdotn | 20 | **0.158 dB** | ~0 |

Both receivers valid at res32 full device.

### Grid-offset P1↔P2 (horns_only, res32, rt=5) — yee_sdotn sweep

| offset | yee_sdotn | yee_guide_normal (§16) |
|---|---:|---:|
| (0, 0) | **0.070 dB** | 0.105 dB |
| (±0.5, 0) | 0.180 dB | **0.115 dB** |
| (0, ±0.5) | 0.87–0.90 dB | 0.80–0.83 dB |

Yee snapping on launch+DFT helps zero-offset; **y-offset sensitivity remains ~0.9 dB** for all tested receivers.

### Receiver recommendation by resolution

| Resolution | Source | Receiver | Status |
|---|---|---|---|
| **res32** | per-port numerical mode | **guide-normal flux** (`num_mode_guide_normal`) | **validated** (full rt=20: 0.166 dB P1↔P2) |
| res64 horns_only | per-port numerical mode | guide-normal or yee_sdotn | ~0.25 dB P1↔P2 |
| res64 full device | per-port numerical mode | any tested flux/DFT | **blocked** (~3 dB P1↔P2) |

### Production promotion (updated)

**Still NOT ready** for production promotion at res64. **Ready for res32 optimization prototyping** with `num_mode_guide_normal`, rt=20, np=32, discrete reciprocity control.

Remaining blockers:
1. res64 full-device axis P1↔P2 port metric (~3 dB; P2↔P3 passes at 0.0006 dB).
2. Grid y-offset ±0.5 cell on y-axis (~0.9 dB at res32).
3. P1/P2 incident-power mismatch (~4.5% res32, ~7.8% res64).

**Acceptable for:** res32 horns + full-device B=0 optimization prototyping with `num_mode_guide_normal`, rt=20, np=32.

---

## 18. Optimization forward-eval cost (res32, np=32, cached norm)

Measured on this workstation with `mpi_runner.py`, `OMP_NUM_THREADS=1`, `FI_PROVIDER=tcp`, `MPICH_CH4_NETMOD=ofi`.

| Item | Wall time |
|---|---:|
| Norm cache build (one-time, 2 ports) | ~7–90 s (formulation-dependent) |
| Full-device 2-port study (cached norm + discrete control) | **~25 s** |
| Device FDTD only (2×2 ports, rt=20) | **~12 s** |
| **Per source excitation** | **~6 s** |
| 6-source forward 6×6 matrix (est.) | **~36 s** |

**Cached across rho (invariant):** numerical mode profiles, horn geometry, incident normalization, flux/DFT quadrature weights.

**Per rho evaluation (inner loop):** PMM FDTD solve + port extraction only.

**Estimated throughput (res32, rt=20, cached norm):**
- 3–4 active ports (typical circulator objective): **~100–140 evals/hour**
- Full 6×6 characterization: **~60 evals/hour**

**Optimizer needs:** typically 3 forward solves per rho (one per active input port), not full 6×6 — budget **~20 s/rho** at res32 with warm cache.

---

## 19. Physical units cleanup + high-resolution port validation (2026-09-02)

### Why `a = 28 mm` was used before

The production notebook (`Sketchbook_PMMCirculator.ipynb`, cell 4) sets `a = 0.028 m` with the comment *"1 Meep length unit = 2.8 cm"*. This was an **arbitrary Meep normalization length**, not the paper lattice constant. Physical device dimensions were always entered as SI lengths divided by `a`:

- `d_exp = 0.020 / a` → **20 mm** lattice pitch regardless of `a`
- bulb / horn sizes likewise via `X_m / a`

So **`a = 28 mm` did not represent the paper geometry** as the lattice constant; it was a legacy normalization choice. The physical lattice pitch was still exactly 20 mm.

### Cleanup: `a` = lattice constant

`scripts/validation/sixport_common.py` now sets:

| Item | Before | After cleanup |
|---|---|---|
| Normalization length `a` | 0.028 m (arbitrary) | **0.020 m = 20 mm lattice pitch** |
| `d_exp` in a-units | 0.714 | **1.0** |
| `fs_a` (3.85 GHz) | 0.3596 c/a | **0.2568 c/a** |
| `fp_a` (8.00 GHz) | 0.7472 c/a | **0.5337 c/a** |
| Domain `nx×ny` (a-units) | 23×21 | **30×28** (same physical extent) |

All **SI physical dimensions are unchanged**; only normalized Meep coordinates and nondimensional frequencies rescale.

### Physical geometry audit

| Quantity | Paper (Rodriguez et al. 2026) | Code before `a` cleanup | After cleanup |
|---|---:|---:|---:|
| Lattice center-to-center pitch | 20 mm | 20 mm (`0.020/a`, `a=28 mm`) | **20 mm** (`d_exp=1.0 a`) |
| Quartz tube OD | 15 mm | 15 mm | **15 mm** |
| Quartz tube ID | 13 mm | 13 mm | **13 mm** |
| Quartz wall thickness | 1 mm | 1 mm (via OD−ID) | **1 mm** |
| Quartz ε_r | ≈3.8 | 3.8 | **3.8** |
| Horn aperture | — (commercial) | 104 mm | **104 mm** |
| Horn throat | — | 48 mm | **48 mm** |
| Horn flare depth | — | 89 mm | **89 mm** |
| Straight feed | — | 60 mm | **60 mm** |
| Plasma radius (model) | — | 5 mm | **5 mm** |

Six-port topology is unchanged; only the normalization convention is now physically transparent.

### New user-facing resolution interface

**Primary:** `points_per_cm` with `dx_mm = 10 / points_per_cm`.

Every validation JSON now reports:

- `points_per_cm`
- `dx_mm`
- internal Meep `res` (= `round(a_m / dx_m)`)
- normalization length `a_m`

Implementation: `scripts/validation/physical_units.py`; `reciprocity_b0_study.py` accepts `--points-per-cm` (preferred) or legacy `--res`.

| points/cm | dx (mm) | Meep `res` (`a=20 mm`) | Old Meep res64 equiv. ppc (`a=28 mm`) |
|---:|---:|---:|---:|
| 50 | 0.200 | **100** | ~22.9 |
| 75 | 0.133 | **150** | — |
| 100 | 0.100 | **200** | — |

Low-resolution Meep res32/res64 campaigns are **debugging history only**; substantive validation starts at **≥50 points/cm**.

### High-resolution B=0 port campaign (50 points/cm)

Harness: `run_hires_port_validation.py` → `reciprocity_b0_study.py`, formulation **`num_mode_guide_normal`**, full 91-bulb device, `rt=20`, **np=32**, matched discrete Yee-grid reciprocity control.

Artifacts: `outputs/validation/hires_ports/`

| Item | P1↔P2 | P2↔P3 |
|---|---|---|
| points/cm | 50 | 50 |
| dx_mm | 0.20 | 0.20 |
| Meep `res` | 100 | 100 |
| `a` | 0.020 m | 0.020 m |
| MPI ranks | 32 | 32 |
| Wall time | 3.74 h | 3.73 h |
| Discrete FDTD reciprocity | **7.5×10⁻¹³ dB** | **1.6×10⁻¹³ dB** |
| Port reciprocity \|Pij−Pji\| | **2.259 dB** | **4.8×10⁻¹³ dB** |
| Raw T (dB) | P12 −24.46 / P21 −26.72 | P23 = P32 = −25.36 |
| Incident power | P1 2910 / P2 2859 | P2 = P3 = 2859 |
| Incident mismatch | 1.8% | ~0 |
| Flux-subtraction diagonal | P1 +0.279 / P2 +0.015 (**unhealthy**) | P2/P3 +0.015 (**unhealthy**) |
| Numerical-mode TE1 overlap | P1 0.792 / P2 0.790 | same P2/P3 caches |

**Does the port formulation behave correctly when the spatial discretization is genuinely fine?**

- **FDTD control: yes.** Matched discrete Hz-patch reciprocity remains machine precision at 50 points/cm. Do not reopen plasma/Drude/Faraday debugging.
- **Diagonal ports (P2↔P3): yes.** Port metric is machine precision, same as at low res.
- **Axis port (P1↔P2): no.** **2.26 dB** remains. Old Meep res64 (~22.9 points/cm, `a=28 mm`) was **~2.9 dB**; finer physical grid did **not** remove the axis-port failure.
- **Flux subtraction is unhealthy** (positive reflection diagonals) on both pairs at this formulation, so the 2.26 dB number is a port-extraction problem, not a discrete Maxwell failure.

The apparent res64 axis-port problem **persists at genuinely fine physical resolution**. It is not merely an intermediate-grid pathology of old res64.

### 75 / 100 points/cm — not started

Estimate from 50 points/cm (cell count ∝ ppc², ~3.74 h wall per 2-port case at np=32):

| Target | dx | Meep `res` | Cell scale vs 50 ppc | Est. wall / pair | Est. both pairs |
|---:|---:|---:|---:|---:|---:|
| 75 ppc | 0.133 mm | 150 | 2.25× | **~8.4 h** | **~17 h** |
| 100 ppc | 0.10 mm | 200 | 4× | **~15 h** | **~30 h** |

75/100 were not launched: 50 points/cm already answers the high-res question (axis-port metric still fails; discrete control still passes). A 100 ppc confirmation is not justified until the P1 flux-subtraction / axis-port formulation is fixed.

### Preserved findings (unchanged)

- Matched discrete B=0 FDTD reciprocity ≈ machine precision (now also at 50 points/cm)
- Off-grid point-probe reciprocity failure was a sampling artifact
- Numerical horn launch substantially improved the port source
- P2↔P3 has generally behaved much better than P1↔P2 — **and remains so at 50 points/cm**

---

## 20. S-parameter reciprocity mismatch — overnight decision tree (2026-09-04/05)

### Phase 1 audit — how S21/S12 are computed

Artifacts: `outputs/validation/sparam_audit/phase1_sparam_path_audit.{md,json}`

| Quantity | Formula | Uses `load_minus_flux_data`? |
|---|---|---|
| Incident Pk | `\|flux\|` on straight-feed reference oriented as port k | no |
| S21-like | `flux(P2 \| drive P1) / incident_P1` | **no** |
| S12-like | `flux(P1 \| drive P2) / incident_P2` | **no** |
| S11-like | source-port flux **after** subtract / incident | **yes** |

**Critical finding:** `load_minus_flux_data` is called **only on the source-port monitor**. It cannot create S21≠S12. Transmission off-diagonals are always raw receive flux / incident. Unhealthy positive reflection diagonals are a **separate** reflection-subtraction bug.

Reference roles:
- P1 full-device run uses **P1** reference only; P2 uses **P2** only — no cross-port reference sharing.
- Norm cache keys: `res`, `run_time`, `formulation`, grid/monitor offsets, horn walls, rotation, port subset.
- Mode caches: per-port JSON under `mode_profiles/res{N}/` (P1/P2/P3 at res100 present; cross-port numerical-mode power overlap >99.99%).
- Reference vs device monitors use the **same** `monitor_center_for_port` / flux-region builders; only the geometry content differs.

### Overnight experiments (completed except abandoned monitor-shift follow-up)

Harness: `dual_receiver_diagnostic.py` — one FDTD run measures **both** guide-normal flux and numerical-mode overlap (no subtraction on modal path).

Campaign: `run_sparam_reciprocity_overnight.py` → `outputs/validation/sparam_audit/`

Planned sequence:
1. horns-only P1↔P2 / P2↔P3 at 50 ppc (flux vs modal; sub vs nosub)
2. full PMM P1↔P2 / P2↔P3 at 50 ppc (dual receive)
3. full P1↔P2 nosub confirmation (expect identical transmission to sub)
4. horns-only monitor moved +2 / +5 cells into feed

The monitor-shift follow-up was intentionally not resumed after the first unfinished case.

### Overnight results (auto-filled)

| Case | Device | Settings | Flux Δ dB | Modal Δ dB | Discrete dB | Flux-sub healthy | Wall |
|---|---|---|---:|---:|---:|---|---:|
| `horns_dual_P1P2_ppc50_rt20_nosub` | horns_only | skip_sub=True mon=0.0 | 0.0102 | 0.5856 | -1.446e-14 | True | 9731s |
| `horns_dual_P1P2_ppc50_rt20_sub` | horns_only | skip_sub=False mon=0.0 | 0.0102 | 0.5856 | -1.446e-14 | True | 19214s |
| `horns_dual_P2P3_ppc50_rt20_nosub` | horns_only | skip_sub=True mon=0.0 | 1.776e-14 | 6.75e-14 | 1.929e-15 | True | 9735s |
| `horns_dual_P2P3_ppc50_rt20_sub` | horns_only | skip_sub=False mon=0.0 | 1.776e-14 | 6.75e-14 | 1.929e-15 | True | 19159s |
| `full_dual_P1P2_ppc50_rt20_nosub` | full | skip_sub=True mon=0.0 | 2.259 | 1.095 | 7.541e-13 | True | 17541s |
| `full_dual_P1P2_ppc50_rt20_sub` | full | skip_sub=False mon=0.0 | 2.259 | 1.095 | 7.541e-13 | True | 22073s |
| `full_dual_P2P3_ppc50_rt20_sub` | full | skip_sub=False mon=0.0 | 4.832e-13 | 1.918e-13 | 1.639e-13 | True | 22186s |

### Decision-tree status

- Horns-only P1↔P2: flux Δ=0.0102 dB, modal Δ=0.5856 dB
- Full PMM P1↔P2: flux Δ=2.259 dB, modal Δ=1.095 dB
- **Interpretation:** horns nearly OK; full PMM amplifies a small port error.
- **Interpretation:** both flux and modal fail → launch/grid/mode definition, not flux subtraction.
- Full P1↔P2 sub vs nosub flux Δ: 2.259 vs 2.259 (expect nearly identical transmission)

### Direct Poynting-vector integral (`S·n`) at 50 points/cm

Receiver used: `num_mode_yee_sdotn` in `reciprocity_b0_study.py`, which samples complex `Ex`, `Ey`, `Hz` on a DFT slab and integrates

\[
P_{\text{port}} \approx \int \tfrac12 \mathrm{Re}(E \times H^*) \cdot \hat{n}\; dl
\]

along the horn cross-section (`extract_dft_sdotn_power()` in `port_formulations.py`). This uses **no** `load_minus_flux_data`, no reference-field subtraction on transmission, and no modal-overlap receiver.

Existing 50 ppc dual-receiver JSONs did **not** save `Ex/Ey`, so this test required fresh Meep runs.

| Method | Geometry | Pair | S21/S12-like reciprocity error |
|---|---|---|---:|
| matched discrete Yee-grid control | horns_only | P1↔P2 | **1.45×10⁻¹⁴ dB** |
| existing guide-normal flux | horns_only | P1↔P2 | **0.0102 dB** |
| direct Poynting `S·n` integral | horns_only | P1↔P2 | **0.0861 dB** |
| matched discrete Yee-grid control | full PMM | P1↔P2 | **7.54×10⁻¹³ dB** |
| existing guide-normal flux | full PMM | P1↔P2 | **2.259 dB** |
| direct Poynting `S·n` integral | full PMM | P1↔P2 | **2.314 dB** |
| matched discrete Yee-grid control | full PMM | P2↔P3 | **1.64×10⁻¹³ dB** |
| existing guide-normal flux | full PMM | P2↔P3 | **4.83×10⁻¹³ dB** |
| direct Poynting `S·n` integral | full PMM | P2↔P3 | **1.03×10⁻¹² dB** |

Result: **No.** Simply integrating `S·n` across the horn does **not** recover `S21 = S12` for the axis/diagonal `P1↔P2` pair at `B=0`. It preserves reciprocity for `P2↔P3`, but `P1↔P2` remains **~2.31 dB** on the full PMM and **~0.086 dB** on horns-only.

This strengthens the current diagnosis:
- the failure is **not** caused by `load_minus_flux_data`;
- it is **not** specific to modal-overlap normalization;
- it survives a direct field-based `S·n` receiver;
- the full PMM amplifies a smaller axis-port measurement asymmetry already present in horns-only.

### Phase 7 (modes) — preliminary

At 50 points/cm / res100: P1/P2/P3 TE1 power overlaps ≈0.79; P1↔P2 and P2↔P3 numerical-mode shape overlaps **>0.99999**. Stale incompatible modes are unlikely.

---

## 21. Port-physics audit and reciprocity-paired diagnostic (2026-09-07)

### Phase 1 — read-only audit

#### What exact source distribution launches P1 / P2 / P3?

Current `num_mode_*` launches are built from per-port cached `NumericalPortMode` profiles:

- cache load: `plasmeep/ports/mode_registry.py` → `get_numerical_mode()`
- source weights: `plasmeep/ports/numerical_mode.py` → `source_amplitudes()`
- Meep source placement: `plasmeep/ports/numerical_launch.py` → `make_numerical_hz_sources()`

Mathematically, for cached complex transverse samples `φ_k` on offsets `s_k`, the source uses:

\[
J_k \propto \overline{\phi_k}
\]

as distributed `mp.Hz` point sources along the local port tangent. This is **not** an eigenmode source from MPB; it is a discrete impressed `Hz` current pattern derived from a previously extracted numerical line profile.

For the snapped variants (`num_mode_yee_*` / `num_mode_grid_*`), continuous points are first snapped to nearest `Hz` Yee sites and duplicate snaps are merged. Thus:

- P1 and P2 do **not** land on identical DOF layouts after rotation.
- They are **separate numerical profiles**, not one mathematically rotated object reused at two angles.

#### Are P1 and P2 the same properly normalized electromagnetic mode, just rotated?

**No.** They are separate cached per-port profiles with matching qualitative shape but different discrete support:

- `canonical_profile_port()` is identity in `plasmeep/ports/mode_registry.py`
- each port loads its own JSON (`numerical_mode_P1.json`, `numerical_mode_P2.json`, ...)
- snapped P1 and P2 use different Yee-site patterns because the Cartesian lattice treats the horizontal and 60° cross-sections differently

So they are “similar numerical source profiles” rather than one exact rotated discrete mode.

#### How are amplitude, phase, and power normalized?

There are three layers:

1. **Mode-vector normalization**  
   `NumericalPortMode.normalized_field()` L2-normalizes the cached complex `Hz` samples.

2. **Source normalization**  
   `source_amplitudes(target_power=1)` returns `conj(φ)` scaled so `Σ|J_k|² = 1` in discrete source-weight norm.

3. **Reference-run normalization**  
   `scripts/validation/sixport_common.py` → `normalize_port()` runs an isolated straight-guide reference for each port and defines incident power from the chosen receiver:
   - flux formulations: `|Φ_inc|`
   - direct Poynting (`num_mode_yee_sdotn`): `|∫ S·n dl|`
   - modal formulations: `|⟨φ,H_inc⟩|²`

Device “S-parameters” are then power ratios:

\[
P_{ij} = \frac{\text{receiver}_i(\text{field from source }j)}{\text{incident power of source }j}
\]

with dB reported as `10 log10(Pij)`. These are **not** complex scattering amplitudes.

#### Do P1 and P2 use the same Yee-grid source DOFs?

Not after snapping. The axis-aligned P1 aperture and the 60° P2 aperture collapse onto different snapped Hz-site patterns. That means even when the same conceptual chord width and cached offsets are used, the actual discrete launch operator differs by orientation.

#### What receiver is paired with each source formulation?

- `num_mode_guide_normal`: total guide-normal flux across a weighted chord
- `num_mode_yee_sdotn`: direct field-based `∫ (1/2) Re(E×H*)·n dl`
- `num_mode_modal`: `Hz`-only overlap `⟨φ,H⟩`

None of the first two is a mode projector. They measure **total net power crossing the plane**, not the amplitude of one designated port mode.

The third (`num_mode_modal`) is closer to an adjoint modal receiver, but it still uses only scalar `Hz` samples, not a full electromagnetic forward/backward mode pair with Lorentz-orthogonality normalization.

#### Is the current receiver plane a clean single-mode port plane?

Geometrically, yes; electromagnetically, not necessarily.

From `plasmeep/ports/horn.py`:

- flare depth = **89 mm**
- straight feed length = **60 mm**
- legacy receiver plane = **18 mm** from throat (`monitor_frac = 0.30`)
- source plane = **42 mm** from throat (`source_frac = 0.70`)
- clear PEC guide width = **40 mm**

At `fs = 3.85 GHz`, that 40 mm guide is only just above first cutoff (`fc1 ≈ c / (2W) ≈ 3.75 GHz`). So there is approximately **one propagating mode**, but the monitor plane is only **0.45 guide widths** and about **0.053 guided wavelengths** from the throat. Evanescent / higher-order reactive content from the flare and PMM transition can therefore still contaminate the field there.

#### Does total `∫ S·n dl` measure the quantity needed for an S-parameter?

Not in general. It measures **total net real power** crossing that section. A true modal S-parameter requires the complex amplitude of a designated port mode. Those coincide only when the cross-section is effectively single-mode and contamination-free.

#### Exact reciprocity issue

For a reciprocal multiport electromagnetic system, `S21 = S12` requires:

1. the two ports be defined as reciprocal forward/backward modes of the same linear operator;
2. excitation and readout use the corresponding adjoint modal pairing;
3. normalization be performed with the same modal convention at both ports.

Our current `num_mode_hz_line + flux` and direct `S·n` receivers fail this because the source is a discrete `Hz` current pattern while the receiver is **total power through a plane**, not the adjoint modal functional.

Current `num_mode_modal` improves the situation but still loses reciprocity mathematically because it projects only `Hz` on a continuous sampled line, not the exact snapped DOFs of the launched source and not the full electromagnetic modal bilinear form.

This is why the matched discrete reciprocity control can be machine precision while the extracted port powers differ by 2.31 dB:

- the matched discrete test uses the **same Yee DOFs and identical weights** as both source and receiver;
- the S-parameter path adds a non-matching measurement operator (flux, total `S·n`, or incomplete modal projection).

### Phase 2 — chosen single best diagnostic

Choice: **C (both)**.

The implemented diagnostic is:

- source: exact snapped Yee-grid modal source on the actual port line
- receiver: adjoint `Hz` modal overlap on the **same snapped Yee DOFs**
- plane: moved deeper into the straight feed, from 18 mm to **30 mm** past the throat (`feed_frac = 0.50`)

Implementation:

- new shared helper: `plasmeep/ports/grid_modal.py`
- new formulation: `num_mode_grid_modal_clean`

Hypothesis:
- the old measurement loses reciprocity because it combines an orientation-dependent discrete source with a non-adjoint or contaminated receiver plane.

Support:
- if `horns_only` P1↔P2 at 50 ppc becomes very tight, then the main remaining error was the old source/receiver definition, not the solver.

Falsification:
- if `horns_only` remains noticeably nonreciprocal, then even this paired `Hz`-only modal definition is still missing the operation needed to match the reciprocal electromagnetic response, and a fuller E/H Lorentz-mode pairing is the next likely step.

### Phase 3 — horns-only gate at 50 points/cm

Case: `paired_modal_horns_P1P2_ppc50_rt20.json`

| Method | P1→P2 | P2→P1 | Difference dB |
|---|---:|---:|---:|
| matched discrete control | — | — | **1.45×10⁻¹⁴** |
| direct Poynting `S·n` (previous) | -14.260 | -14.346 | **0.0861** |
| new `num_mode_grid_modal_clean` | **+4.691** | **-4.691** | **9.382** |

Raw facts:

- incident normalization P1: **1.0119×10^209**
- incident normalization P2: **3.4356×10^208**
- incident mismatch: **2.945×**
- normalized off-diagonal powers became exact reciprocals: **2.945** and **0.3395**

This is a mathematically useful failure. The new paired-`Hz` modal receiver/source did **not** produce a valid modal power normalization even in `horns_only`.

Interpretation:

1. The underlying horns-only field solve is still reciprocal (matched discrete control remains machine precision).
2. The new receiver/source pairing is closer in DOF matching, but the quantity

\[
| \langle \phi, H \rangle |^2
\]

with `φ` a scalar `Hz` line profile is **not** a power-orthogonal electromagnetic port norm.
3. Therefore the reference “incident power” extracted from that functional is not invariant across P1/P2 orientation and can become arbitrarily ill-scaled.
4. Because the horns-only gate failed, a full PMM confirmation was **not** run.

### Current best diagnosis

The single most likely remaining cause is now more specific:

- not `load_minus_flux_data`
- not total Poynting integration alone
- not Meep reciprocity
- not merely “wrong monitor location”

Instead, the missing mathematical ingredient is a **true electromagnetic forward/backward port-mode bilinear form**. A scalar `Hz` overlap, even on exact snapped DOFs and a cleaner plane, is still not the Lorentz-reciprocal adjoint required for a physically normalized S-matrix.

The likely next step is a port definition based on a full E/H reference mode and the appropriate Lorentz overlap functional, rather than `Hz`-only Euclidean projection.

## 22. Meep eigenmode-coefficient audit for full E/H ports (2026-09-07)

### Phase 1 — established Meep mathematics

Meep's documented mode-decomposition path (`Simulation.get_eigenmode_coefficients`, `Mode_Decomposition`) does **not** use a scalar `Hz` norm. It uses a full electromagnetic bilinear form over a port cross-section `S`:

\[
\left\langle \psi,\psi' \right\rangle
= \int_S \left[
\mathbf{E}^*(\rho)\times\mathbf{H}'(\rho)
+ \mathbf{E}'(\rho)\times\mathbf{H}^*(\rho)
\right]\cdot\hat{\mathbf{n}}\,dA
\]

where `\hat{n}` is the port normal and `\psi=(E_\parallel,H_\parallel)` contains the **tangential electric and magnetic fields together**. For our 2D `Hz` polarization on a line port this reduces to

\[
\left\langle \psi,\psi' \right\rangle
= \int \left[E_t^* H_z' + E_t' H_z^*\right]\,dl
\]

with `E_t = \mathbf{E}\cdot\hat{\mathbf{t}}` on the guide cross-section.

Meep separates forward/backward amplitudes by projecting onto the forward and backward traveling versions of the same mode. For a reciprocal uniform guide, the tangential fields satisfy:

- `E_t^- = E_t^+`
- `H_z^- = -H_z^+`

So if we define

\[
I_1 = \int E_{t,m}^* H_z\,dl,\qquad
I_2 = \int E_t H_{z,m}^*\,dl
\]

then the Meep-style coefficients are

\[
\alpha^+ = \frac{1}{2}(I_1 + I_2),\qquad
\alpha^- = \frac{1}{2}(I_1 - I_2)
\]

after the reference mode has been normalized to **unit modal power**.

### Modal power normalization

Per Meep's documentation, the eigenmodes are normalized so that the real power flux of the reference mode is unity:

\[
P_m = \Re\int_S \mathbf{E}_m^*\times\mathbf{H}_m\cdot\hat{\mathbf{n}}\,dA = 1
\]

In the present 2D line-port reduction:

\[
P_m = \Re\int E_{t,m}^* H_{z,m}\,dl = 1
\]

With that normalization, the coefficient powers satisfy

\[
|\alpha^\pm|^2 = P^\pm
\]

which is exactly the quantity needed for a power-normalized modal S-parameter.

### Why this is the right adaptation target

This mathematics is directly reusable with a **numerical reference mode** extracted from a time-domain straight PEC guide, even though MPB/EigenModeSource is not suitable for the staircased horn. The required ingredients are:

1. full complex `Ex`, `Ey`, and `Hz` for the numerical guide mode,
2. a fixed port line with known tangent/normal,
3. a consistent quadrature/interpolation rule for the measured fields and the reference mode,
4. unit-power normalization of the reference mode before using the overlap.

The previous `Hz`-only `|<\phi,Hz>|^2` diagnostic failed because it discarded `E_t`, so it could not represent Meep's forward/backward electromagnetic bilinear form or its power normalization.

### Phase 2/3 — numerical reference-mode audit and reusable implementation

New reusable helper: `plasmeep/ports/em_mode.py`  
Gate driver: `scripts/validation/eh_modal_horns_gate.py`  
Output: `outputs/validation/eh_modal/horns_only_P1P2_ppc50_rt20.json`

Implementation choices:

- store full complex `Ex`, `Ey`, and `Hz` samples on the port line;
- reduce `Ex/Ey` to tangential `E_t = E\cdot\hat{t}`;
- normalize the reference mode by `P_m = Re ∫ E_t^* H_z dl = 1`;
- compute forward/backward amplitudes by
  - `α+ = 0.5 (I1 + I2)`
  - `α− = 0.5 (I1 − I2)`
  - `I1 = ∫ E_{t,m}^* H_z dl`
  - `I2 = ∫ E_t H_{z,m}^* dl`
- explicitly reject NaN/Inf/absurdly large field magnitudes and non-positive mode powers.

The reference runs at `50 points/cm` produced **finite, sensible** modal quantities for both ports:

| Port | `|α_inc,+|²` | `|α_inc,-|²` | chosen incident power |
|---|---:|---:|---:|
| P1 | `1.98×10⁻⁴` | **2910.18** | **2910.18** |
| P2 | `19.56` | **2812.47** | **2812.47** |

Interpretation:

- the inward (`α−`) component is dominant for both reference runs, as expected;
- the incident modal powers are of the same order but **not identical** (~3.4% mismatch);
- the coefficients are finite, unlike the invalid `Hz`-only exact-DOF test.

### Phase 4 — horns-only hard gate (`P1↔P2`, 50 ppc, `np=32`, `rt=20`)

| Method | `S21` | `S12` | `|S21/S12|` mismatch |
|---|---:|---:|---:|
| matched discrete Yee-grid control | — | — | **1.45×10⁻¹⁴ dB** |
| existing guide-normal flux | — | — | **0.0102 dB** |
| direct Poynting `S·n` integral | — | — | **0.0861 dB** |
| full E/H modal overlap (new) | `0.02765` | `0.03087` | **0.4779 dB** |

Raw transmitted modal powers:

- `P1 -> P2`: `|α_out,+|² = 80.465`
- `P2 -> P1`: `|α_out,+|² = 86.808`

Raw source-port decomposition in the device runs:

- `P1` source case: `|α_src,-|² = 2892.32`, `|α_src,+|² = 32.87`
- `P2` source case: `|α_src,-|² = 2826.55`, `|α_src,+|² = 42.40`

### Horns-only verdict

This full E/H formulation is **mathematically defensible and numerically finite**, but it does **not** restore near-reciprocal horns-only `P1↔P2` behavior at `50 points/cm`.

What improved:

- the result is no longer an invalid normalization catastrophe;
- the overlap now uses the correct electromagnetic bilinear form rather than scalar `Hz` projection;
- forward/backward amplitudes and incident modal powers are explicitly separated.

What did **not** improve enough:

- horns-only `P1↔P2` remains **0.478 dB**, only modestly better than the prior `Hz`-only modal receiver result (**0.5856 dB**), and much worse than both
  - guide-normal flux (**0.0102 dB**), and
  - direct horn-aperture `S·n` (**0.0861 dB**).

Most likely remaining cause after this gate:

The receiver math is now much closer to Meep's eigenmode-coefficient formulation, so the next blockage is more likely in the **reference/source pairing and modal-plane equivalence**, not in the bilinear form itself. In particular:

1. the time-domain straight-guide reference mode used for `P1` and `P2` is not yielding equal incident modal powers;
2. the launched horn fields in the full horns geometry still contain a sizeable opposite-going / non-mode component at the measurement plane (`|α_src,+|²` is not small);
3. the chosen clean plane is still not acting like a single-mode guide cross-section for this reciprocal decomposition.

Therefore, the full E/H overlap is a **necessary correction**, but by itself it is **not sufficient** to make the present horn source/plane pair behave like a true reciprocal eigenmode port.
