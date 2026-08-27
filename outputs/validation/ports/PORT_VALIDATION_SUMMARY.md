# Port Validation Summary — agent/eigenmode-ports

**Date:** 2026-08-27 (resolution / MPI / polarization cycle)  
**Branch:** `agent/eigenmode-ports`  
**Launcher (fixed):** `te1_hz_line` + axis-aligned flux (harness only; not production)

---

## Executive conclusion

1. **Measurement reformulation is closed** — no alternative receiver beat `te1_hz_line`+axis flux toward 0.2 dB at cheap res; not escalated further.
2. **High-resolution P1↔P2 does not converge to the 0.2 dB gate** — strongly **non-monotonic** through res=128; res96 local minimum (0.68 dB), **res128 regresses to 1.52 dB**.
3. **Grid refinement alone does not explain the residual** as a simple discretization error that vanishes with Δx→0. Next hypothesis: **horn/PEC geometry, mode purity, or model–experiment polarization mismatch** — not more FluxRegion variants or blind res sweeps.
4. **MPI at np=32** (FI_PROVIDER=tcp, OMP_NUM_THREADS=1) makes res≈100+ practical (~5× faster than np=4 at res48 smoke).

---

## What is `a`?

| Quantity | Value |
|---|---|
| Code constant | `a = 0.028` **m** |
| **Physical Meep unit** | **2.8 cm** (not the 20 mm lattice pitch) |
| Lattice pitch in code | `d_exp = 0.020/a` → **20 mm** = 0.714 a-units |

Frequencies nondimensionalize as \(f_a = f\,a/c\) via PMMI (same as notebook).

---

## Geometry vs beam-steering paper

Paper: 91 elements, triangular lattice, six bulbs/side, **20 mm** pitch, quartz **OD 15 mm / ID 13 mm / ε≈3.8**.

| Feature | Paper | Meep harness | Match |
|---|---|---|---|
| N, lattice, pitch | 91, triangular, 20 mm | yes | ✓ |
| Quartz OD/ID/ε | 15/13 mm, 3.8 | `r_outer=7.5 mm`, `r_inner=6.5 mm`, ε=3.8 | ✓ |
| Horns | 3 commercial horns | 6 synthetic PEC flares + feeds | **different** |
| Device | 3-port steering | 6-port circulator | **different** |

Full audit: `GEOMETRY_POLARIZATION_AUDIT.md`.

---

## Polarization (components, not TE/TM labels)

| Model | Active E/H in 2D |
|---|---|
| Paper / Ceviche | **Ez** (z out of page); plots \|Ez\|² |
| Meep baseline + `te1_hz_line` | **Hz** source; **Ex, Ey** in Poynting receivers |
| Gyrotropy with **B ∥ z** | Couples **Ex ↔ Ey**; **Hz** is the magnetically active 2D polarization |

**Answer:** Current circulator model is **not** the same polarization as the published beam-steering FDFD setup (Ez vs Hz). For **B≠0 circulator physics in 2D**, Hz is **required**; matching paper Ez is a separate experimental-fidelity task. **Did not block** this resolution study.

---

## Measurement cycle (completed — no further variants)

Cheap res32/rt40 P1↔P2 (dB): baseline+axis **0.18**, TE1+axis **0.38**, TE1+guide-normal **0.34**, TE1+DFT S·n **0.32**, dense guide-normal **0.35**. None materially improved; **stopped**.

---

## High-resolution reciprocity convergence (`te1_hz_line`, B=0, P1↔P2)

Fixed: TE1 launch + axis flux. Incident mismatch stays **~2–3%** on excited ports.

| res | Δx (mm) | px/cm | px / 20 mm lattice | \|P12−P21\| (dB) | incident spread | wall (h) | MPI |
|---:|---:|---:|---:|---:|---:|---:|---|
| 32 | 0.875 | 11.4 | 22.9 | **0.384** | 2.0% | 0.04 | np=4 |
| 48 | 0.583 | 17.1 | 34.3 | 1.282 | 2.5% | 0.24 | np=4 |
| 64 | 0.438 | 22.9 | 45.7 | 1.033 | 2.9% | 0.62 | np=4 |
| 96 | 0.292 | 34.3 | 68.6 | **0.682** | 3.3% | 2.08 | np=4 |
| 128 | 0.219 | 45.7 | 91.4 | **1.519** | 2.3% | 3.93 | np=32 |

**Verdict:** Error **does not** monotonically decrease toward 0.2 dB. Local minimum at res96; **res128 worse than res64**. This is **not** consistent with a pure resolution/discretization error washing out at ~0.22 mm grid spacing.

Plot: `outputs/validation/advisor_update/fig5_te1_reciprocity_convergence.png`  
JSON: `outputs/validation/ports/te1_hz_line_convergence.json`

---

## MPI scaling (representative device port, res=48, rt=30)

Env: `OMP_NUM_THREADS=1`, `FI_PROVIDER=tcp`, `MPICH_CH4_NETMOD=ofi`, `UCX_TLS=tcp,self`  
(WSL default UCX shm fails at ≥16 ranks unless OFI/TCP is set.)

| ranks | device (s) | speedup vs np=4 | efficiency | flux (identical) |
|---:|---:|---:|---:|---|
| 4 | 243.4 | 1.00 | 100% | 0.21533 |
| 8 | 215.6 | 1.13 | 56% | 0.21533 |
| 16 | 107.1 | 2.27 | 57% | 0.21533 |
| 32 | 46.2 | **5.27** | **66%** | 0.21533 |
| 48 | 54.1 | 4.50 | 38% | 0.21377 |
| 64 | 51.2 | 4.75 | 30% | 0.21377 |

**Recommended production config:** **np=32**, OMP=1, FI_PROVIDER=tcp. np=48/64 show no gain at res48.

Artifacts: `outputs/validation/mpi_scale/mpi_scale_summary_res48_te1_hz_line_ofi.json`

---

## Validation script inventory (no new standalone scripts)

| Script | Role |
|---|---|
| `sixport_common.py` | Shared geometry, simulate, norm cache, `physical_resolution_report` |
| `port_formulations.py` | Formulation registry (launcher/measurement) |
| `reciprocity_b0_study.py` | P1↔P2 / pair reciprocity + JSON |
| `incident_power_orientation.py` | Six-port incident uniformity |
| `run_port_validation.py` | A/B orchestrator |
| `run_reciprocity_campaign.py` | Overnight batch driver |
| `mpi_scale_oneport.py` | MPI benchmark (extend, don't duplicate) |
| `faraday_benchmark.py` | Isolated gyrotropy check (preserved PASS) |
| `make_advisor_figures.py` | Plots from existing JSON |

**Reusable code to promote later:** horn builder + explicit-component sources/monitors → see `HORN_API_PROPOSAL.md`.

---

## Faraday (secondary, preserved)

PASS at res32/64; |κ| error &lt;0.1%. Global κ sign **opposite** theory with perfect B-reversal — likely **convention** (Stokes ψ vs circular-basis labeling), not magnitude bug. Note: `outputs/validation/faraday_sign_convention_note.md`.

---

## Gate status

| Criterion | Status |
|---|---|
| P2↔P3 ≲ 0.05 dB | **PASS** |
| P1↔P2 ≲ 0.2 dB improving with res | **FAIL** (non-monotonic; 0.68 dB @96, 1.52 dB @128) |
| Incident within few % | **PASS** (TE1) |
| Resolution explains residual | **FAIL** (no convergence to gate) |
| Production untouched | **PASS** |

**Overall: FAIL** — launcher/incident fixed; reciprocity residual not cleared by measurement or grid refinement to res128.

---

## Recommended next step

1. **Stop** measurement A/B and blind resolution escalation.
2. Investigate **horn/PEC feed redesign** (straight reference guide, de-embedded segment) and/or **3D / Ez experimental-fidelity** path if matching paper hardware.
3. Keep **np=32 OFI/TCP** for any remaining high-res diagnostics.
4. Optional: finite-length Faraday degree report (secondary).

Campaign logs: `outputs/validation/ports/pair_P1P2_te1_hz_line_res128_rt80.json`, `measurement_cycle_verdict.json`, `GEOMETRY_POLARIZATION_AUDIT.md`.
