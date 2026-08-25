# Overnight Summary — B=0 Six-Port Reciprocity Residual (FINAL)

**Branch:** `agent/sixport-reciprocity`  
**Campaign status:** complete (`campaign_ressweep` EXIT:0; incident-power study complete)  
**Wall clock (res-sweep alone):** ~6.3 h at `np=4`  
**Artifacts:** `outputs/validation/reciprocity/`

---

## 1. Specialist panel (pre-run)

Independent reviews agreed before new Meep runs:

| Specialist | Main finding |
|---|---|
| Physics | Faraday benchmark PASS → residual unlikely magnetized-plasma physics at B=0 |
| Numerics | Canonical pair structure → Cartesian rotated-feed discretization |
| Ports/measurements | Global norm bug unlikely; orientation-dependent source/monitor sampling likely |
| Code | Use notebook-faithful harness; production S-param API unfinished |

Canonical notebook matrix (res=64, rt=80, B=0): mean |Pij−Pji|≈**0.75 dB**, max≈**1.64 dB**.

Pair anatomy (notebook):

- **~1.64 dB:** axis↔diagonal (P1↔P2, P1↔P6, P3↔P4, P4↔P5)
- **~1.19 dB:** other axis↔diagonal
- **~0 dB:** diagonal↔diagonal and P1↔P4 (opposite axis)

Ports: P1/P4 along ±x; P2/P3/P5/P6 at 60° — poorly aligned with the Yee grid.

---

## 2. Infrastructure

| Artifact | Role |
|---|---|
| `scripts/validation/sixport_common.py` | Notebook-faithful geometry + `simulate_circulator` |
| `scripts/validation/reciprocity_b0_study.py` | CLI B=0 study / JSON out |
| `scripts/validation/run_reciprocity_campaign.py` | Overnight orchestrator |
| `scripts/validation/mpi_scale_oneport.py` | MPI scaling |
| `scripts/validation/incident_power_orientation.py` | Axis vs diagonal launched power |

Commits: `b96aaf9`, `0e57ee5`.

---

## 3. MPI scaling (res=32, rt=30, one port)

| ranks | device time (s) | speedup |
|------:|----------------:|--------:|
| 1 | 122.6 | 1.00 |
| 2 | 72.5 | 1.69 |
| 4 | 58.3 | 2.10 |
| 8 | 46.3 | 2.65 |

**Choice used:** `np=4`.

---

## 4. Final discriminating results

### A. Orientation contrast (decisive)

| Case | res | rt | max \|Pij−Pji\| | incident |
|---|---:|---:|---:|---|
| P1↔P2 (axis–diag) | 32 | 40 | **0.181 dB** | P1=32.74, P2=29.66 |
| P2↔P3 (diag–diag) | 32 | 40 | **0.0003 dB** | P2≡P3 |

Residual appears only for axis↔diagonal pairs.

### B. Runtime (rules out “just run longer”)

| Case | res | rt | max \|Pij−Pji\| |
|---|---:|---:|---:|
| P1↔P2 | 32 | 40 | 0.181 dB |
| P1↔P2 | 32 | 80 | **0.486 dB** |

Longer runtime increases measured residual. Not a fix.

### C. Resolution sweep on worst pair P1↔P2 (rt=80) — COMPLETE

| res | max \|Pij−Pji\| | T12 / T21 (dB) | wall (device+norm) | notes |
|---:|---:|---|---:|---|
| 32 | 0.486 dB | −10.17 / −10.66 | ~4.5 min | absolute T not yet notebook-like |
| 48 | 0.452 dB | −29.45 / −29.00 | ~30.6 min | absolute T unstable |
| **64** | **1.636 dB** | **−23.00 / −24.64** | ~75.9 min | **matches notebook max (1.6356 dB) and T values** |
| **96** | **0.098 dB** | −19.74 / −19.65 | ~264.7 min | sharp drop |

Harness fidelity at production resolution is confirmed: res=64 reproduces the notebook P1↔P2 residual and transmission levels.

### D. Incident power orientation bias — COMPLETE

| res | axis/diag | max/min | axis mean | diag mean |
|---:|---:|---:|---:|---:|
| 32 | 1.104 | 1.104 | 32.735 | 29.656 |
| 48 | 1.136 | 1.136 | 33.182 | 29.216 |
| 64 | 1.112 | 1.112 | 32.751 | 29.444 |
| 96 | 1.104 | 1.104 | 32.746 | 29.675 |

Within each orientation class, ports agree to machine precision (P1≡P4, P2≡P3≡P5≡P6).  
Across classes, launched power stays **~10–14% unequal at all resolutions** — does **not** converge toward 1.

---

## 5. Reassessment of root cause (full sweep)

### What the full sweep changes

Earlier interim conclusion (“naive resolution will not help”) is **too strong**.

With the complete P1↔P2 curve:

1. Residual is **non-monotonic** in resolution: ~0.45–0.49 (coarse) → **1.64 peak at res=64** → **0.098 at res=96**.
2. Production/notebook resolution sits near a **bad local residual**, not a converged floor.
3. Going to res=96 reduces the worst-pair residual by **~17×** relative to res=64.
4. Incident axis/diag bias remains ~10% even when reciprocity residual falls to ~0.1 dB, so the residual is **not** just a constant launched-power scale factor leaking into normalized T.

### Confirmed root-cause class

**Orientation-dependent Cartesian port inequivalence** (Hz point-line sources / flux monitors on axis vs 60° feeds) remains the best explanation:

- notebook pair anatomy,
- P1↔P2 ≫ P2↔P3,
- persistent axis vs diagonal incident-power split,
- Faraday material path already PASS,
- B=0 true nonreciprocity ruled out.

### Refined ranking

1. **Primary — rotated-port discretization / sampling inequivalence on the Yee grid**  
   Still the cause class. Evidence stronger after full sweep because res=64 exactly recovers the notebook defect and res=96 shows the defect is numerically reducible.

2. **Secondary — under-resolved absolute transmission at res≤64**  
   Absolute T/R levels jump between 32/48/64/96; the reciprocity error peaks where the notebook was run. Convergence of the *observable* needs higher resolution or better ports, not only more runtime.

3. **Contributor — power/flux observables rather than complex modal S**  
   Can amplify orientation impurity, but cannot explain near-zero diag↔diag residuals.

4. **Ruled out — insufficient run_time as the fix**
5. **Ruled out — PlasMEEP B=0 material / Faraday sign bug**
6. **Ruled out — physical nonreciprocity at B=0**

---

## 6. What remains unresolved

- Exact split of residual among source sampling, monitor flux geometry, and mode impurity.
- Whether eigenmode sources/monitors eliminate the stubborn ~10% incident axis/diag bias.
- Whether res=96 residual (~0.1 dB) is acceptable as a temporary validation floor, or ports must be redesigned first.
- Full 6×6 at res=96 (only P1↔P2 was swept; notebook pattern predicts other axis↔diag pairs behave similarly).
- +B/−B transpose residual after B=0 floor is reduced.

---

## 7. Recommended next validation gate

**Do not** tune plasma / `rho` / B to force reciprocity.

**Preferred scientific next step (ports redesign):**

1. Replace Hz point-line feeds with orientation-invariant excitation/monitoring (eigenmode source + mode-coefficient monitors preferred).
2. Re-test **P1↔P2 vs P2↔P3** at fixed res (64 and/or 96).
3. Target:
   - diag↔diag ≲ 0.05 dB  
   - axis↔diag ≲ 0.2 dB and improving under refinement  
   - incident max/min → ~1 within a few percent
4. Then full 6×6 B=0, then +B/−B transpose.

**Interim numerical option (if ports redesign is deferred):**

- Treat **res=96** as the minimum resolution for trusting reciprocity-sensitive claims with the *current* sources: worst-pair residual falls to ~0.1 dB.
- Do **not** treat res=64 notebook residual (~1.64 dB) as an irreducible physics floor.

**Scientific choice for you:** accept ~0.1 dB floor at res=96 temporarily, or force eigenmode ports before further six-port validation.

---

## 8. Bottom line

The completed campaign confirms the B=0 residual is an **orientation-dependent numerical/port artifact**, faithfully reproduced by the harness at res=64 (matching the notebook’s **1.636 dB** P1↔P2 error), and **strongly reduced at res=96 (~0.098 dB)**. Incident launched-power bias between axis and diagonal ports remains ~10% at all resolutions. Faraday material validation stands. Next work is ports redesign and/or adopting higher resolution—not plasma-parameter tuning or more blind B=0 physics speculation.
