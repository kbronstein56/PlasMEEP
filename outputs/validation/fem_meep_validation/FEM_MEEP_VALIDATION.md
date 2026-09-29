# FEM ↔ Meep Validation Report

**Verdict: NOT VALIDATED**

Date: 2026-09-22  
Work tree: `outputs/validation/fem_meep_validation/`  
Reusable code: `scripts/validation/fem_meep_validation/`

The ~11 s FEM-H solver is **not** yet trustworthy as an inverse-design inner loop.
B=0 port powers do not match trusted Meep within the required band after implementing
true Poynting receivers, numerical-mode sources, and TE-consistent metal walls.
Per the autonomous gate (**stop B≠0 if B=0 exceeds ~1 dB**), full-device magnetized
validation was **not** claimed.

---

## 1. Exact Meep source definition

Formulation: `num_mode_guide_normal` (see `phase1/MEEP_PORT_SPEC.md`).

- Cached `NumericalPortMode` at res50 / fs ≈ 3.85 GHz.
- Complex Hz samples on offsets along port tangent through `source_center`.
- Amplitudes \(a_k = \sqrt{P_0}\,\overline{\phi_k}\) with \(\sum |a_k|^2 = P_0\).
- Meep `GaussianSource(frequency=fs_a, fwidth=0.10·fs_a)` point Hz sources.
- Time convention: \(e^{-i\omega t}\).

## 2. Exact Meep receiver / power definition

- Monitor at `monitor_center_for_port`, outward normal \(n̂ =\) `effective_port_dir`.
- Discrete \(\displaystyle P_{\mathrm{raw}} = \int S\cdot n̂\,ds\) with 31-point trapezoid,
  span \(0.96\times\) clear_width.
- \(S_x = \tfrac12\mathrm{Re}(E_y H_z^*)\), \(S_y = -\tfrac12\mathrm{Re}(E_x H_z^*)\).
- Incident: straight-feed reference → \(P_{\mathrm{inc}}=|P_{\mathrm{raw,ref}}|\).
- Source port: `load_minus_flux_data` (field subtraction), then \(T = P_{\mathrm{raw}}/P_{\mathrm{inc}}\).

Trusted Meep reference (this campaign):

| Quantity | Value | Source |
|---|---|---|
| P11 | 0.168267 | prior `full_num_mode_guide_normal_P1P2_ppc25_rt20.json` |
| P1→P2 | 0.01372968 | this P1→all6 run (bit-match prior) |
| P1→P3 | 0.00220668 | this P1→all6 run |
| P1→P4 | 0.318804 | this P1→all6 run |
| P1→P5 | 0.00220668 | this P1→all6 run |
| P1→P6 | 0.01372968 | this P1→all6 run |
| P_inc | 2909.721 | production norm cache |

Artifacts: `phase2/meep_P1_all6_ppc25_df0.10.json`, `phase2/trusted_prior_P1P2_ppc25.json`.

**Note:** Pickled Meep `FluxData.E/H` arrays are all zeros in this environment, so
`load_minus_flux_data` from disk cache is a no-op (P11 printed raw ≈ −0.86).
Transmission ports do not need subtraction and are trusted. P11 taken from the prior
validated P1P2 study.

## 3. How FEM reproduces each piece

| Meep piece | FEM implementation |
|---|---|
| Numerical-mode Hz line source | Same offsets/amps scattered onto nearest mesh nodes (`fem_validated_solver.py`) |
| Guide-normal flux | Reconstruct \(E=(i/\omega)\rho\,\mathrm{curl}\,H\); \(P=\int S\cdot n̂\,ds\) |
| Incident norm | Straight feed with ρ→0 sidewalls; \(P_{\mathrm{inc}}=\|P_{\mathrm{raw,ref}}\|\) |
| Source-port subtraction | Field subtraction on monitor line: \((E,H)_{\mathrm{dev}}-(E,H)_{\mathrm{ref}}\) |
| Metal horns | **ρ→0** in wall triangles (ε→∞). **Not** Dirichlet Hz=0 (that is PMC for TE) |
| Plasma B=0 | Isotropic \(\rho=1/\varepsilon\) from Meep Drude at fs |
| Plasma B≠0 | Full \(\rho=\varepsilon^{-1}\) with Meep Faraday tensor (unit-tested; device not validated) |

Critical bug found and fixed during diagnosis: earlier FEM-H “agreement”
(`FEM-H_level3_direct.json`, P14≈0.33) used **no metal enforcement** (nnz=4.63M).
Adding Dirichlet Hz=0 (nnz=4.37M) destroyed transmission (P14≈0.016). Correct TE PEC
is Neumann / ρ→0, restoring nnz=4.63M and P14≈0.36.

## 4. B=0 port-by-port comparison (FEM-H best path)

Method: ρ→0 walls + numerical mode + Poynting + field-sub P11.  
Table: `phase6/b0_comparison_table.json`.

| Port | Meep | FEM | Meep dB | FEM dB | FEM−Meep dB | rel. err |
|---|---:|---:|---:|---:|---:|---:|
| P11 | 0.1683 | 0.0984 | −7.74 | −10.07 | **−2.33** | −42% |
| P1→P2 | 0.01373 | 0.02603 | −18.62 | −15.85 | **+2.78** | +90% |
| P1→P3 | 0.002207 | 0.01848 | −26.56 | −17.33 | **+9.23** | +738% |
| P1→P4 | 0.3188 | 0.3608 | −4.96 | −4.43 | **+0.54** | +13% |
| P1→P5 | 0.002207 | 0.01885 | −26.56 | −17.25 | **+9.32** | +754% |
| P1→P6 | 0.01373 | 0.02531 | −18.62 | −15.97 | **+2.66** | +84% |

## 5. Maximum / RMS dB error

- **max |ΔdB| = 9.32 dB** (P1→P5)
- **RMS |ΔdB| = 5.66 dB**
- Opposite port P14 alone: **0.54 dB** (acceptable in isolation)
- Band label: **inconsistent (>1 dB)**

## 6. FEM mesh-convergence (M vs H, same physics path)

`phase8/mesh_convergence.json`

| Grade | DOFs | h_crit | FEM−Meep max \|ΔdB\| | M−H ΔdB (P14) |
|---|---:|---|---:|---:|
| FEM-M | 234307 | 0.4 mm | ~9.9 dB | — |
| FEM-H | 661890 | 0.2 mm | 9.32 dB | **−2.03 dB** |

FEM-M and FEM-H still differ by ~2 dB on P14 under the corrected BC.
**0.2 mm local spacing is not proven Meep-equivalent** from convergence.

## 7. Is FEM-H physically validated at B=0?

**No.**

## 8. Gyrotropic ε tensor used (Meep / Faraday convention)

\[
\varepsilon =
\begin{bmatrix}
\varepsilon_\perp & -i\eta \\
+i\eta & \varepsilon_\perp
\end{bmatrix}
\]

with \(\varepsilon_\perp,\eta\) from `faraday_benchmark.gyrotropic_drude_eps_eta`
(ordinary-frequency Meep Drude, bias along z).

## 9. ρ = ε⁻¹ implementation

\[
\rho = \varepsilon^{-1},\quad
\rho_{xx}=\varepsilon_{yy}/\det,\;
\rho_{xy}=-\varepsilon_{xy}/\det,\;
\rho_{yx}=-\varepsilon_{yx}/\det,\;
\rho_{yy}=\varepsilon_{xx}/\det.
\]

Assembled in anisotropic stretched weak form (reduces to isotropic \(1/\varepsilon\) at B=0).

## 10. B-reversal / unit-test result

`phase9/gyrotropic_unit_tests.json`: **all_pass = true**

- B=0 → η=0 → isotropic ✓  
- B→−B reverses off-diagonals of ε and ρ ✓  
- Faraday |κ| reverses with B ✓  
- Tensor signs match Meep Faraday convention ✓  

Device-level B≠0 **not** compared (B=0 gate failed).

## 11. B≠0 full-device Meep vs FEM

**Not run** (stopped by >1 dB B=0 gate).  
No matching trusted magnetized six-port Meep JSON was located at comparable settings;
Faraday homogeneous benchmark remains the tensor sign reference.

## 12. Nonreciprocal behavior agreement

**Not demonstrated** on the full device.

## 13. Final validated first-RHS time (corrected path)

On FEM-H (661890 DOFs, ρ→0 walls):

| Step | Time |
|---|---:|
| Factor | ~5.7 s |
| First solve | ~0.19 s |
| Factor + first solve | ~5.9 s |
| Full wall (mat+asm+factor+solve+PP) | ~15–20 s (PP ~3 s for E/S) |

(Prior ~11 s campaign number used wrong Dirichlet PMC walls.)

## 14. Final validated six-port candidate time

Extra RHS solve ≈ **0.18 s** (measured).  
Six solves after one factor: **~6.7 s** solve-only.  
With material/assemble/postprocess once: estimate **~20–25 s** per (ρ,ω) candidate
depending on port PP cost — not optimized (per instructions).

## 15. Number of DOFs

**661890** (FEM-H), nnz **4631312** with ρ→0 walls.

## 16. Is ~0.2 mm local mesh adequate?

**Not based on convergence.** FEM-M vs FEM-H still differ by ~2 dB on P14; Meep
errors remain ≫1 dB on side/diagonal ports at FEM-H.

## 17. Trustworthy for inverse-design inner loop?

**No.**

## 18. What still prevents production use

1. **B=0 Meep mismatch** up to ~9 dB on P13/P15; ~2.7 dB on P12/P16; ~2.3 dB on P11.
2. **Mesh not converged** (M↔H ~2 dB on P14).
3. **CW FEM vs pulsed Meep DFT** may still contribute residual error.
4. **Wall geometry** is polyline ρ→0 tubes, not exact Meep prism solids.
5. **Meep FluxData pickle broken** here (zeros) — in-process norm required for P11.
6. **No trusted magnetized full-device Meep column** located for B≠0 FEM check.
7. Port postprocessing (nodal E) is a large fraction of wall time.

---

## Frequency / objective (Phase 13)

`circulator_objective` in `sixport_common` operates on a **single-frequency** 6×6
power matrix at `fs_a`. FEM at fixed ω evaluating dispersive ε(ω) is the correct
model for that objective. Broadband would require rebuild/refactor per frequency
(A changes with ω); six sources at one ω reuse one factorization.

---

## Autonomous gate outcome

| Gate | Result |
|---|---|
| Real port power implemented | Yes (Poynting) |
| B=0 matches Meep adequately | **No** (max 9.3 dB) |
| Full gyrotropic B≠0 | Code + unit tests only |
| Magnetized device vs Meep | **Not run** |
| Mesh convergence acceptable | **No** |

**Final label: NOT VALIDATED**

Do **not** select READY FOR INVERSE-DESIGN INTEGRATION.

---

## Key artifacts

| Path | Content |
|---|---|
| `phase1/MEEP_PORT_SPEC.md` | Exact Meep math |
| `phase2/meep_P1_all6_ppc25_df0.10.json` | Meep P1→all6 (+ trusted P11) |
| `phase2/numerical_mode_P1_res50.json` | Mode profile for FEM |
| `phase6/b0_comparison_table.json` | B=0 table |
| `phase7/` | Diagnosis (Dirichlet bug, ρ→0, field-sub) |
| `phase8/mesh_convergence.json` | FEM-M vs H |
| `phase9/gyrotropic_unit_tests.json` | Tensor unit tests |
| `phase12/timing_summary.json` | Timing |
| `scripts/validation/fem_meep_validation/` | Validation-only solvers |
