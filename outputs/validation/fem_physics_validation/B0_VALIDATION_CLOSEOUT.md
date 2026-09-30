# B0 validation closeout

Checkpoint of the pre-closeout audit: `d7b1aeeda5b08881a4836baa0f732c286707837f` on `origin/agent/eigenmode-ports`.

Statuses below were recomputed by `scripts/validation/fem_physics_validation/reaudit_closeout.py` from the JSON files. Thresholds in `PASS_CRITERIA.md` were not changed. A Meep-versus-analytic failure is not an FEM failure.

The queue has the original 10 FAIL rows and 24 UNRESOLVED rows. None were dropped.

## Work queue

| ID | description | current status | reason for FAIL/UNRESOLVED | required experiment | expected computational cost | dependency | result | final status |
|---|---|---|---|---|---|---|---|---|
| 2f | analytical E/H relation | PASS | — | see evidence `eh_orders.json` | completed in this closeout | none | vacuum element median 0.0841% at h=0.003125; plasma element max 0.0767% at h=0.0015625. Both are order 1, as expected for a P1 gradient. Hz L2 is 3e-6. | PASS |
| 3d | E/H relation | PASS | E reconstruction is first order, one order slower than Hz. | see evidence `guide_eh_close.json` | completed in this closeout | none | h=0.00064, 3,909,063 DOFs. Hz L2 7.78e-5 (order 2.00). Gradient L2 0.0836% (order 1.01). Element median 0.0789%. Element max 0.0897% (order 1.13). Power relative error 1.12e-4. | PASS |
| 3h | Meep comparison | PASS | — | see evidence `meep_closeout.json` | completed in this closeout | none | m=0, beta=k0. res 40 phase +0.280°. res 80 phase +0.070° and −3.6e-6 dB. Dispersion falls as resolution squared. | PASS |
| 6g | frequency sweep | PASS | — | see evidence `freq_fem.json` | completed in this closeout | none | finest mesh, worst probe: 3.20 −0.0029 dB −0.055°; 3.85 −0.0067 dB −0.112°; 4.50 −0.0171 dB −0.218°; 5.30 −0.0248 dB −0.147°; 6.00 −0.0067 dB −0.293°. | PASS |
| 6m | Meep resolution convergence | MEEP_CROSS_CODE_FAIL | class D Meep. Not an FEM failure. | see evidence `meep_crosscheck.json; meep_closeout.json` | completed in this closeout | none | dielectric Meep passes. Plasma Meep, including the new coated and cluster runs, does not return to the analytic field as resolution increases. | MEEP_CROSS_CODE_FAIL |
| 6n | Meep sub-cell registration | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | registrations kept separate. Coated forward res 40: ox0 −0.133 dB −2.40°; ox0.25 −0.233 dB −3.57°. Spread is not collapsed. | MEEP_CROSS_CODE_FAIL |
| 7f | nearby frequencies | PASS | analytic reference is stable; the solve meets 0.05 dB and 0.5 deg | see evidence `coated_530_graded.json` | completed in this closeout | none | Finest graded mesh 5,940,235 DOFs, h_near median 0.00268, 26.9 quartz elements. Forward -0.0118 dB, -0.461 deg. Near-quartz +0.0474 dB, -0.046 deg. Ring L2 0.424%. Order about 2. | PASS |
| 7h | Meep resolution sequence | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 24 ox0 forward −0.120 dB −1.50°; res 40 ox0 −0.133 dB −2.40°. Higher resolution did not improve the forward probe. | MEEP_CROSS_CODE_FAIL |
| 7i | Meep registration sequence | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 40 ox0.25 forward −0.233 dB −3.57° versus ox0 −0.133 dB −2.40°. Registrations were not averaged. | MEEP_CROSS_CODE_FAIL |
| 8d | conditioning | PASS | raw condition stays large. The physical observable is not contaminated. | see evidence `tmatrix_condition.json` | completed in this closeout | none | coated 7 at m_max=8: raw cond 1.98e10, scaled 3.32e3, residual 7.4e-16, scaling moves the forward probe by 9.9e-15. At m_max=14 raw cond 2.43e26, scaled 2.35e4, forward change from the previous m_max is 2.0e-10. | PASS |
| 9c | Meep | MEEP_CROSS_CODE_FAIL | class D Meep. FEM pair orientations pass in 9d. | see evidence `meep_closeout.json` | completed in this closeout | none | pair 0° res 30 ox0 forward +1.07 dB +16.1°; ox0.25 +0.25 dB −3.83°. 30° and 90° are also outside 0.10 dB and 1°, and the two registrations disagree. | MEEP_CROSS_CODE_FAIL |
| 9d | orientation relative to the Yee grid | PASS | — | see evidence `orient_fem.json` | completed in this closeout | none | angles 0, 30, 45, 60, 90 at h=0.02. Worst probe on that mesh is the 30° gap, −0.0067 dB −0.104°. Last step from h=0.04 is about a factor of three. A y-directed source at (0,−4.5) sits inside the PML on the 14×10 box and is not used. | PASS |
| 10c | Meep | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 30 ox0 forward +0.87 dB +24.2°; ox0.25 +0.27 dB −4.79°; res 40 ox0 −0.29 dB −5.96°. Resolution does not remove the error. | MEEP_CROSS_CODE_FAIL |
| 11c | Meep | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 24 ox0 forward −1.67 dB +25.2°; res 36 ox0 −6.32 dB −41.4°; res 36 ox0.25 −0.53 dB −5.20°. Not averaged. | MEEP_CROSS_CODE_FAIL |
| 12g | Meep comparison | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `prior seven-bulb records` | completed in this closeout | none | original-grid late residual remains about 0.6 dB and 15° from the T-matrix. New bare-seven Meep runs in 11c are the same class of failure. | MEEP_CROSS_CODE_FAIL |
| 12h | field comparison, not one forward scalar | PASS | — | see evidence `coated7_local.json` | completed in this closeout | none | h=0.02, h_edge=0.0045, 618177 DOFs, 11 cells across the quartz wall. Forward −0.0090 dB −0.243°; backward +0.017 dB −0.037°; ±60° −0.011 dB +0.06°; gap −0.020 dB −0.333°; outside shell −0.020 dB −0.332°; near quartz −0.020 dB −0.335°. Previous outside-shell phase at h_edge=0.008 was −0.412°. Last local step is 0.08°. | PASS |
| 13b | oblique propagation | SUPERSEDED_TEST_HARNESS | old box does not isolate PML reflection | see evidence `oblique_pml.json` and `oblique_pml_fit.json` | completed in this closeout | none | Old rel L2 stays 0.020 to 0.025 and probe reflection stays 0.067. Replacement Bloch fit: abs(R) falls from 1.74e-3 to 1.34e-4 (0 deg), 1.49e-3 to 1.37e-4 (25 deg), 8.79e-4 to 2.18e-4 (60 deg) as dpml goes from 0.4 to 1.2. Fit residual about 4e-5. | SUPERSEDED_TEST_HARNESS |
| 13c | cylinder radiation | PASS | — | see evidence `pml_domain.json` | completed in this closeout | none | one plasma disk, h=0.04. Forward across dpml 0.8/1.2/1.6, sigma ×0.5/×1/×2, and boxes 11×8, 14×10, 18×12: −0.0037 to −0.010 dB, phase −0.083° to −0.162°. Knob spread 0.0065 dB and 0.08°. | PASS |
| 13d | multi-cylinder radiation | PASS | — | see evidence `pml_domain.json` | completed in this closeout | none | bare pair, h=0.04. Forward −0.0147 to −0.029 dB and −0.179° to −0.294° across dpml and sigma. Inside 0.10 dB and 1°. | PASS |
| 13h | outer-domain variation | PASS | — | see evidence `pml_domain.json` | completed in this closeout | none | outer boxes 11×8, 14×10 and 18×12 at fixed dpml and sigma. Disk forward −0.0037, −0.0066 and −0.0037 dB. | PASS |
| 14d | lossless cluster | PASS | — | see evidence `power_closeout.json` | completed in this closeout | none | lossless dielectric pair, source outside, contour power 1.77e-12. Single dielectric cylinder contour 1.09e-12. | PASS |
| 14i | explicit volume absorption | PASS | — | see evidence `absorption_identity.json` | completed in this closeout | none | analytic bare volume 1.340e-5 plus outward flux −1.337e-5, relative residual 0.24%. Coated residual 0.24%. Absorption is positive. FEM volume on one plasma cylinder is stable to 0.3% from h=0.02 to h=0.01. Same-mesh FEM contour at h=0.02 agrees with that volume to 1.5%. The h=0.01 contour integral is not stable and is not used. The discrete FEM load is not H0, so absolute FEM power is not compared to the Hankel source. | PASS |
| 15a | simple guide | PASS | — | see evidence `guide_reciprocity.json` | completed in this closeout | none | quartz block in a closed guide, point-source exchange, max relative 2.86e-14. This is G(a,b) versus G(b,a), not a normalized modal S matrix. | PASS |
| 15b | horns only | PASS | — | see evidence `horns_reciprocity.json` | completed in this closeout | none | FEM-L horns, six monitor centers, max relative 2.02e-14 on the Green matrix. | PASS |
| 15e | seven-bulb cluster | PASS | — | see evidence `seven_reciprocity.json` | completed in this closeout | none | seven coated bulbs, h=0.05, 15 pairs, max relative 2.65e-14. | PASS |
| 18c | factorization consistency | PASS | — | see evidence `linalg_closeout.json` | completed in this closeout | none | guide 4131 DOFs residual 3.2e-14, COLAMD versus NATURAL relative solution difference 1.2e-13, repeat difference 0. Disk 97156 DOFs residual 1.8e-14, repeat 0. | PASS |
| 18d | FEM condition estimate | PASS | — | see evidence `linalg_closeout.json` | completed in this closeout | none | 1-norm condition about 2.1e4 on the guide and 1.0e6 on the plasma disk. Residuals are 1e-14. Linear-algebra error is negligible next to mesh error. The condition number is not required to be small. | PASS |
| 20c | FEM: PML | PASS | oblique harness remains 13b | see evidence `pml_domain.json` | completed in this closeout | none | cylinder forward PML/domain spread 0.0065 dB and 0.08°. Pair spread 0.014 dB. Oblique absolute reflection is not used; see 13b. | PASS |
| 20g | Meep: resolution | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | plasma Meep forward error does not fall from res 24 to res 40 on the coated bulb, nor from res 30 to res 40 on three cylinders. | MEEP_CROSS_CODE_FAIL |
| 20h | Meep: grid registration | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | ox=0 and ox=0.25 kept separate. Pair 0° forward spread is 0.8 dB and 20°. Seven-bare res 36 spread is several dB. | MEEP_CROSS_CODE_FAIL |
| 20i | Meep: runtime | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | runs to t≈960–980. Late plasma values remain outside the analytic tolerance. | MEEP_CROSS_CODE_FAIL |
| 20j | Meep: dispersive-interface error | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | quartz Meep matched the slab. Plasma cylinders, coated bulbs, and clusters do not. | MEEP_CROSS_CODE_FAIL |
| 20k | Meep: DFT convergence | MEEP_CROSS_CODE_FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | the finished DFTs are the late values quoted in 7h, 9c, 10c and 11c. They are not the analytic field. | MEEP_CROSS_CODE_FAIL |
| 20o | analytic: conditioning | PASS | raw 2-norm condition remains up to 1e26 | see evidence `tmatrix_condition.json` | completed in this closeout | none | equilibration is a reference-solver scale, not a fit to FEM. Scaled condition of coated seven is 3.3e3 at m_max=8 and 2.4e4 at m_max=14. Forward observable change from scaling is below 1e-11. | PASS |

## Failure classes

A. FEM physics or operator failure: none.

B. FEM observable reconstruction: 2f and 3d pass. Element-center Ey/Hz is the direct P1 gradient. On the guide it is order 1.01 while Hz is order 2.00. At h=0.00064 the element max is 0.0897% and the gradient L2 is 0.0836%. No smoothing was added.

C. Insufficient refinement: none remaining. The coated cylinder at 5.30 GHz passes on the 5,940,235-DOF graded mesh.

D. Meep disagreement with the analytic continuum, kept as MEEP_CROSS_CODE_FAIL: 6m, 6n, 7h, 7i, 9c, 10c, 11c, 12g, 20g, 20h, 20i, 20j, 20k. Registrations were not averaged.

E. Inappropriate pass criterion near a null: not used. The 5.30 GHz forward amplitude is 0.279, not a null. The phase slope is -865 deg/GHz. The 0.5 deg bar was not relaxed.

F. Harness: 13b is SUPERSEDED_TEST_HARNESS. The Bloch-guide fit replaces it. Cylinder and pair PML sweeps remain 13c, 13d, and 13h.

## Decision

**B0_FEM_VERIFIED_AND_VALIDATED**

3d, 7f, and the replacement for 13b meet the frozen thresholds on computed solves. The old oblique box is SUPERSEDED_TEST_HARNESS and is not marked PASS. The 13 Meep rows stay MEEP_CROSS_CODE_FAIL. They are not grounds to withhold the FEM label. B ≠ 0, adjoints, and optimization were not started.

## Uncertainty budget

Production-like seven coated bulbs, forward probe on the old domain: last mesh step 0.003 dB and 0.13°, analytic error −0.0023 dB and −0.100°. New-domain local mesh, near-quartz probe: −0.020 dB and −0.335°, last local step 0.08°.

These contributions are added, not combined as an RSS, because the mesh step already contains the PML used on that mesh.

- FEM mesh: the last step above.

- FEM geometry: plasma area relative error about 1e-4 on the scatter meshes.

- FEM PML: disk forward spread 0.0065 dB and 0.08° across thickness, sigma, and outer box. Pair spread 0.014 dB.

- FEM ports: three guide stations agree to about 1e-6 relative. Absolute power relative error is 4.25e-4 at h=0.00125.

- FEM materials: production Drude ε at 3.85 GHz, evaluated once per frequency.

- FEM linear solve: residual about 1e-14, repeat difference 0, reorder difference 1e-13. Negligible beside the mesh step. 1-norm condition is 2e4 (guide) and 1e6 (disk).

- Meep: plasma registration and resolution spreads are order 1 dB and 10°, so Meep is not inside the FEM tolerance.

- Analytic truncation: coated-seven forward change from m_max=8 to the converged value is 6e-6 in complex amplitude at m_max=8, and 2e-10 at m_max=14. Scaled condition 3e3 to 2e4. Raw condition up to 2e26 does not move the probe (scaling change below 1e-11).

Conservative total for the seven-bulb forward scalar: about 0.01 dB and 0.2°, the mesh step plus the measured PML knob spread. The near-interface probe is bounded by its 0.33° analytic error and 0.08° local step, which is inside the 1° multi-cylinder allowance.

The same FEM / analytic and Meep / analytic tables, and the trust statement, are in `VALIDATION_COMPLETENESS_AUDIT.md`. They are repeated here so this closeout stands alone.

## FEM / ANALYTIC VALIDATION

| benchmark | analytic reference | finest FEM mesh | error | uncertainty | status |
|---|---|---|---|---|---|
| MMS, four tensors | imposed Hz | n=64 | L2 4.25e-4, order 1.998 | H1 order 0.999 | PASS |
| Homogeneous waves and PEC guide | exact fields | guide h=0.00125 | Hz L2 2.97e-4; power rel 4.25e-4 | stations agree to 3e-7 | PASS |
| Guide E/H, direct P1 sample | exact Ey/Hz | h=0.00064, 3,909,063 DOFs | grad L2 0.0836%; element max 0.0897%; power rel 1.12e-4 | gradient order 1.01, Hz order 2.00 | PASS |
| Fresnel and slabs | closed form and transfer matrix | h=0.02 | within 0.0015 dB and 0.025 deg | invalid h=0.04 stack is not used | PASS |
| Bare Mie, including 5.30 GHz | independent Mie | h=0.012 | worst probe about -0.025 dB and -0.15 deg at 5.30 GHz | production 3.85 GHz is -0.0015 dB | PASS |
| Coated cylinder, including 5.30 GHz | multilayer Mie, cond 1, residual 0 | 5,940,235 DOFs, h_near 0.00268 | forward -0.0118 dB, -0.461 deg; ring L2 0.424% | order about 2; 0.5 deg bar not relaxed | PASS |
| Seven coated bulbs | T-matrix | h=0.02, h_edge 0.0045, 618177 DOFs | near-quartz -0.020 dB, -0.335 deg | last local step 0.08 deg | PASS |
| Oblique PML reflection | Bloch mode decomposition | h=0.02, dpml=1.2 | abs(R) 1.34e-4, 1.37e-4, 2.18e-4 at 0, 25, 60 deg | fit residual about 4e-5 | PASS |
| Power, reciprocity, factorization | analytic identities | stated meshes | residuals 1e-12 to 1e-14 | disk condition 1e6 | PASS |

## MEEP / ANALYTIC CROSS-CODE

| benchmark | Meep resolution | registration | runtime | error | status |
|---|---|---|---|---|---|
| Guide | res 80 | single | until 80 | +0.070 deg, -3.6e-6 dB | PASS |
| Coated bulb | res 40 | ox 0 and 0.25, not averaged | until 140 | -0.133 dB / -2.40 deg and -0.233 dB / -3.57 deg | MEEP_CROSS_CODE_FAIL |
| Bare pair | res 30 | ox 0 and 0.25 | until 160 | +1.07 dB, +16.1 deg at 0 deg, ox 0 | MEEP_CROSS_CODE_FAIL |
| Three cylinders | res 40 | ox 0 | until 180 | -0.29 dB, -5.96 deg | MEEP_CROSS_CODE_FAIL |
| Seven bare | res 36 | ox 0 and 0.25 | until 200 | -6.32 dB / -41.4 deg and -0.53 dB / -5.20 deg | MEEP_CROSS_CODE_FAIL |
| Seven coated, prior grid | original campaign | kept separate | late DFT | about 0.6 dB and 15 deg | MEEP_CROSS_CODE_FAIL |

## WHY B0 FEM IS TRUSTED

The guide derivative, the 5.30 GHz coated solve, and the Bloch PML fit are computed results, not trends. Guide element max 0.0897% and gradient L2 0.0836% are under 0.1%, at the expected orders 1 and 2. The coated forward error is -0.0118 dB and -0.461 deg on 5,940,235 DOFs, and the previous mesh was order 2 away from it. PML abs(R) falls with thickness at three angles and stays near 1e-4 when the air length changes. Meep plasma errors of order 1 dB remain labeled MEEP_CROSS_CODE_FAIL.


