# FEM vs Meep decision

**QUARTZ_VALIDATED**

Date: 2026-09-27. Overnight record: `OVERNIGHT_B0_VALIDATION.md`.

Horn-only and quartz-only FEM and Meep approach the same continuum. Quartz Meep at 50 ppc matches the converged FEM-VB solution (1,417,429 DOFs) to ≤ 0.08 dB on P2, P3, and P4, after Meep moved toward FEM from 25 to 50 ppc. The full 91-bulb B=0 plasma problem does not. Meep P4 falls from 0.319 at 25 ppc to 0.031 at 50 ppc while the finest FEM value is 0.330, and Meep P2 falls from 0.014 to 0.002 while FEM sits at 0.124. B≠0 was not started.

The 2026-09-26 horn gate below used 25 ppc Meep as the reference and called a 0.24 dB gap a failure. The later resolution series showed that gap shrinking. That update is in `HORN_CONTINUUM_COMPARISON.md` and section 1 of the overnight record.

---

# FEM vs Meep decision — history through the boundary diagnosis

**Earlier status: six-horn residual was 0.21 dB on an unstructured mesh. That number is superseded by the campaign above.**

Date: 2026-09-26. Boundary diagnosis: `BOUNDARY_PML_DIAGNOSIS.md`.

ρ→0 with the mass term left on is not an excised Neumann PEC, and the FEM PML conductivity is about 13× weaker than Meep’s. Both show up in the open cavity. On FEM-H, excised metal plus Meep’s σ brings the ports to within 0.21 dB, which is the size of the mesh’s own P2≠P6 split. Quartz and B≠0 stay stopped.

---

## Stage 0 — what the prior work actually showed

| Item | Status | Evidence |
|---|---|---|
| 1. Same 91-bulb centers | DONE AND TESTED | Same `Rod_Array` centers; center delta 0 m (`stage2/geometry_compare.json`) |
| 2. Quartz OD/ID/ε=3.8 | DONE AND TESTED | Same radii as `Add_Bulb`; OD 15 mm, ID 13 mm, wall 1 mm |
| 3. B=0 plasma ε | IMPLEMENTED BUT NOT VALIDATED | Same Drude formula as Faraday code; not compared cell-by-cell to Meep’s Yee ε |
| 4. Horn geometry | IMPLEMENTED BUT NOT VALIDATED | Prism vertices shared, but FEM ρ=0 mask ≠ solid prisms (area Jaccard **0.32**) |
| 5–8. Source/monitor planes, mode, tangent | IMPLEMENTED BUT NOT VALIDATED | Coordinates and mode JSON loaded; no field-line match vs Meep DFT |
| 9–11. Ex/Ey, Poynting, ∫S·n̂ | DONE AND TESTED on a straight guide | See Stage 3. Not validated on the 91-bulb device |
| 12–14. Incident norm, subtraction, six-port powers | IMPLEMENTED BUT NOT VALIDATED | Numbers exist; they disagree with Meep |
| 15. Mesh convergence | DONE AND TESTED — **not converged** | FEM-M vs FEM-H P14 differs by **2.0 dB** |
| 16. Full ρ=ε⁻¹ | IMPLEMENTED BUT NOT VALIDATED on device | Code + algebra tests only |
| 17. B≠0 device plasma | NOT IMPLEMENTED as a Meep comparison | |
| 18. B reversal | IMPLEMENTED BUT NOT VALIDATED | Algebraic unit test only (`phase9`), not a Faraday FEM run |
| 19. B=0 Meep vs FEM | DONE AND TESTED — **failed** | max \|ΔdB\| = **9.32** |
| 20. B≠0 Meep vs FEM | NOT IMPLEMENTED | Stopped by the B=0 gate |

Prior campaign notes that called P14 “close” used a \|Hz\|² proxy or a wall mask that is not Meep’s prism solid. Those are not acceptance results.

---

## 1. Authoritative Meep reference

`stage1/authoritative_meep_P1_reference.json`

Trusted 25 ppc, df=0.10 fs, B=0, `num_mode_guide_normal`, prism, offsets 0, rt=20, fp=8 GHz.

| Quantity | Value | Origin |
|---|---:|---|
| P_inc (P1) | 2909.721 | production norm cache |
| P11 | 0.168267 | prior P1P2 study (flux subtraction intact) |
| P1→P2 | 0.01372968 | this all-six run; Δ vs prior study **4×10⁻¹¹** |
| P1→P3 | 0.00220668 | all-six run |
| P1→P4 | 0.318804 | all-six run |
| P1→P5 | 0.00220668 | all-six run |
| P1→P6 | 0.01372968 | all-six run |

P11 was **not** remeasured in the all-six run: pickled `FluxData` arrays are all zeros, so `load_minus_flux_data` did nothing and the raw diagonal was −0.86. Transmission ports do not need that subtraction. P2 confirms the run is the same physics as the trusted pair study.

---

## 2. Geometry comparison

`stage2/geometry_compare.json`, overlay `stage2/geometry_overlay.png`.

- 91 centers, OD, ID, plasma radius `4.6/6.5·ID`, quartz ε: **identical by construction**. Delta = 0 m.
- Domain 0.60 m × 0.56 m, PML 40 mm: same parameters.
- Horn **solids** are not the same mask. Distance-to-polyline ρ=0 vs point-in-prism: area Jaccard **0.32**. Expected prism area 0.00715 m²; polyline mask 0.00814 m²; polygon mask 0.00747 m².
- Capture radius vs Meep half-thickness differs by **0.20 mm**.

Replacing the mask with the true prism interiors on FEM-M made the port **distribution worse** (below). The bug is not a 0.2 mm radius tweak.

---

## 3. Simple-guide port validation — PASS

`stage3/straight_guide_compare.json`

Parallel-plate guide, clear width = circulator feed, cosine Hz line, Gaussian df=0.10 fs in Meep, monochromatic FEM.

Transmission T = P(downstream) / P(reference), each solver normalized to itself:

| Solver | T (dB) |
|---|---:|
| Meep res40 | +0.00062 |
| FEM h=0.05a | +0.00069 |
| FEM − Meep | **+0.000069 dB** |

Both guides send positive power downstream and negative power upstream. Ex/Ey reconstruction, S = ½ Re(E × H*), and ∫ S·n̂ are consistent with Meep on this problem, well inside 0.1 dB.

---

## 4. B=0 full device — FAIL

Best prior FEM-H (polyline ρ=0 walls, numerical mode, Poynting, field-subtraction P11) vs the authoritative column (`phase6/b0_comparison_table.json`):

| Port | Meep | FEM | FEM−Meep (dB) |
|---|---:|---:|---:|
| P11 | 0.1683 | 0.0984 | **−2.33** |
| P12 | 0.01373 | 0.02603 | **+2.78** |
| P13 | 0.002207 | 0.01848 | **+9.23** |
| P14 | 0.3188 | 0.3608 | **+0.54** |
| P15 | 0.002207 | 0.01885 | **+9.32** |
| P16 | 0.01373 | 0.02531 | **+2.66** |

max \|ΔdB\| = **9.32**. RMS = **5.66**.

This is not a global scale error. Ratios to P14 are already wrong: Meep P13/P14 = 0.0069; FEM ≈ 0.051 (~**+8.7 dB** on the ratio). B=0 symmetry is present on both sides (P2≈P6, P3≈P5), so the error is a symmetric scattering difference, not a swapped port.

True prism walls on FEM-M (`stage5/polygon_walls_FEM-M.json`), ratios to P4 vs Meep:

| Port | ratio error (dB) |
|---|---:|
| P12 | +10.2 |
| P13 | +14.7 |
| P14 | 0 (by construction of the ratio) |
| P15 | +14.0 |
| P16 | +10.7 |

Correcting the wall solid did not move FEM toward Meep.

---

## 5. Mesh convergence

`phase8/mesh_convergence.json` (polyline walls, same Poynting path):

| Grade | DOFs | h_crit | vs Meep max \|ΔdB\| | M−H on P14 |
|---|---:|---|---:|---:|
| FEM-M | 234307 | 0.4 mm | ~9.9 | — |
| FEM-H | 661890 | 0.2 mm | 9.32 | **−2.03 dB** |

FEM-H is not converged. The 9 dB side-port error is already present on FEM-M and barely moves, so another interface refinement is unlikely to be the whole explanation. No finer mesh was run.

---

## 6–8. Gyrotropic tensor, Faraday, B≠0 device

Not accepted.

- Tensor convention in code: ε = [[ε⊥, −iη], [+iη, ε⊥]], ρ = ε⁻¹, from `gyrotropic_drude_eps_eta`.
- B=0 → η=0 and B → −B off-diagonal reversal: algebraic tests passed (`phase9/gyrotropic_unit_tests.json`).
- No FEM Faraday rotation vs the trusted Meep Faraday benchmark.
- No B≠0 full-device comparison (gate: B=0 already >1 dB).

---

## 9–10. Timing

Not reported as a validated candidate time. The straight guide and the failed device solves are fast; that does not change the physics result. Adjoint `A^H` reuse was measured earlier on an unvalidated matrix and is not a success criterion here.

---

## 11. Where it disagrees

| Test | Result |
|---|---|
| Straight-guide ∫S·n̂ | **Agrees** (7×10⁻⁵ dB) |
| Bulb centers / quartz radii | **Agree** (0 m) |
| Horn solid vs FEM mask | **Disagree** (Jaccard 0.32) |
| Full-device port distribution | **Disagrees by up to 9.3 dB**, and true prism fill made it worse on FEM-M |
| Mesh M vs H | **Still moving** (~2 dB on P14) |

The flux formula is not the bug. The multi-port horn + lattice solution is a different scattering problem than Meep. Likely contributors, not yet separated: (1) ρ=0 triangles are not Meep PEC prisms (staircasing, gaps, thickness), (2) FEM-H is not mesh-converged, (3) soft nearest-node Hz sources vs Meep Yee sources inside a real horn. No amplitude fudge was applied.

---

## 12. Horn localization (this step)

`horn_localization/HORN_LOCALIZATION_REPORT.md`

The Meep prism vertices were exported and FEM now consumes that file. The old Jaccard 0.32 was an edge-tube mask around those same quads. With the solid prisms:

- One horn, flare ratio FEM−Meep = **0.032 dB**.
- Six horns, no bulbs, FEM-H vs Meep: P2 **+1.25 dB**, P3 **+0.60 dB**, P4 **−0.88 dB**, P5 **+0.37 dB**, P6 **+1.15 dB**.
- Feed, throat, and monitor |Hz| agree to about 1–2.5% on the inner 80% of each line. The center line between the horns is the first clear divergence (relative L2 **0.11**).
- Horn-only P14 is still moving with mesh refinement (0.656 → 0.614 → 0.569) and that motion is away from Meep (0.697).

Quartz, B=0 plasma, and the 91-bulb device were not rerun. The previous full-device FEM-H table is not a corrected result.

**Status: FAILED**
