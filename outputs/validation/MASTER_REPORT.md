# Master validation report — `agent/eigenmode-ports`

**Generated:** 2026-08-31 (port-mode vs true reciprocity checkpoint)  
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

**Implement port-specific numerical mode launch/measurement** on `horns_only` cheap tests:

1. Extract reference Hz profiles per orientation (`numerical_mode.py`; reuse ±60° symmetry where justified).
2. Launch via numerical mode weights; receive via mode overlap (not total axis-aligned flux alone).
3. Normalize each port to equal incident modal power.
4. Re-test P1↔P2 reciprocity — target <0.2 dB with mode-consistent ports.

**Do not** continue global coordinate rotations or broad horn-wall searches unless direct Lorentz reciprocity regresses.

EigenModeSource: viable only on simplified straight guides (`eigenmode_reference_guide.py`); use time-domain reference extraction for angled PEC horns.

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
- `validation_registry.py` — consolidate JSON artifacts

**Do not rerun:** te1_hz_line P1P2 res32–128 full device; measurement screen; Faraday res32/64; global rotation sweeps.
