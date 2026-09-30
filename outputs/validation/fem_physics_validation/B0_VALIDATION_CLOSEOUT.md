# B0 validation closeout

Checkpoint of the pre-closeout audit: `d7b1aeeda5b08881a4836baa0f732c286707837f` on `origin/agent/eigenmode-ports`.

Statuses below were recomputed by `scripts/validation/fem_physics_validation/reaudit_closeout.py` from the JSON files. Thresholds in `PASS_CRITERIA.md` were not changed. A Meep-versus-analytic failure is not an FEM failure.

The queue has the original 10 FAIL rows and 24 UNRESOLVED rows. None were dropped.

## Work queue

| ID | description | current status | reason for FAIL/UNRESOLVED | required experiment | expected computational cost | dependency | result | final status |
|---|---|---|---|---|---|---|---|---|
| 2f | analytical E/H relation | PASS | — | see evidence `eh_orders.json` | completed in this closeout | none | vacuum element median 0.0841% at h=0.003125; plasma element max 0.0767% at h=0.0015625. Both are order 1, as expected for a P1 gradient. Hz L2 is 3e-6. | PASS |
| 3d | E/H relation | FAIL | class C. About h=0.00064 would put the max under 0.1% if the order stays 1. Not an operator error. | see evidence `eh_orders.json` | completed in this closeout | none | at h=0.0008, 2.50e6 DOFs: element median 0.0987%, element max 0.115%, grad L2 0.105%, Hz L2 1.22e-4. Median meets 0.1%. Max and grad L2 do not. Order of the gradient is 1. | FAIL |
| 3h | Meep comparison | PASS | — | see evidence `meep_closeout.json` | completed in this closeout | none | m=0, beta=k0. res 40 phase +0.280°. res 80 phase +0.070° and −3.6e-6 dB. Dispersion falls as resolution squared. | PASS |
| 6g | frequency sweep | PASS | — | see evidence `freq_fem.json` | completed in this closeout | none | finest mesh, worst probe: 3.20 −0.0029 dB −0.055°; 3.85 −0.0067 dB −0.112°; 4.50 −0.0171 dB −0.218°; 5.30 −0.0248 dB −0.147°; 6.00 −0.0067 dB −0.293°. | PASS |
| 6m | Meep resolution convergence | FAIL | class D Meep. Not an FEM failure. | see evidence `meep_crosscheck.json; meep_closeout.json` | completed in this closeout | none | dielectric Meep passes. Plasma Meep, including the new coated and cluster runs, does not return to the analytic field as resolution increases. | FAIL |
| 6n | Meep sub-cell registration | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | registrations kept separate. Coated forward res 40: ox0 −0.133 dB −2.40°; ox0.25 −0.233 dB −3.57°. Spread is not collapsed. | FAIL |
| 7f | nearby frequencies | FAIL | class C. Estimated h≈0.004, about 9e6 DOFs, to reach 0.5°. Threshold not changed. Production 3.85 GHz still passes. | see evidence `freq_fem.json; coated_530_local.json` | completed in this closeout | none | 3.50, 3.85 and 4.20 GHz pass at h=0.02. At 5.30 GHz, h=0.012 and h_edge=0.004, 1.68e6 DOFs: forward −0.102 dB, −4.73°. Wall refinement from h_edge=0.008 to 0.004 moved the phase from −12.3° only to −10.9°. Global h dominates, order about 2. | FAIL |
| 7h | Meep resolution sequence | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 24 ox0 forward −0.120 dB −1.50°; res 40 ox0 −0.133 dB −2.40°. Higher resolution did not improve the forward probe. | FAIL |
| 7i | Meep registration sequence | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 40 ox0.25 forward −0.233 dB −3.57° versus ox0 −0.133 dB −2.40°. Registrations were not averaged. | FAIL |
| 8d | conditioning | PASS | raw condition stays large. The physical observable is not contaminated. | see evidence `tmatrix_condition.json` | completed in this closeout | none | coated 7 at m_max=8: raw cond 1.98e10, scaled 3.32e3, residual 7.4e-16, scaling moves the forward probe by 9.9e-15. At m_max=14 raw cond 2.43e26, scaled 2.35e4, forward change from the previous m_max is 2.0e-10. | PASS |
| 9c | Meep | FAIL | class D Meep. FEM pair orientations pass in 9d. | see evidence `meep_closeout.json` | completed in this closeout | none | pair 0° res 30 ox0 forward +1.07 dB +16.1°; ox0.25 +0.25 dB −3.83°. 30° and 90° are also outside 0.10 dB and 1°, and the two registrations disagree. | FAIL |
| 9d | orientation relative to the Yee grid | PASS | — | see evidence `orient_fem.json` | completed in this closeout | none | angles 0, 30, 45, 60, 90 at h=0.02. Worst probe on that mesh is the 30° gap, −0.0067 dB −0.104°. Last step from h=0.04 is about a factor of three. A y-directed source at (0,−4.5) sits inside the PML on the 14×10 box and is not used. | PASS |
| 10c | Meep | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 30 ox0 forward +0.87 dB +24.2°; ox0.25 +0.27 dB −4.79°; res 40 ox0 −0.29 dB −5.96°. Resolution does not remove the error. | FAIL |
| 11c | Meep | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | res 24 ox0 forward −1.67 dB +25.2°; res 36 ox0 −6.32 dB −41.4°; res 36 ox0.25 −0.53 dB −5.20°. Not averaged. | FAIL |
| 12g | Meep comparison | FAIL | class D Meep. | see evidence `prior seven-bulb records` | completed in this closeout | none | original-grid late residual remains about 0.6 dB and 15° from the T-matrix. New bare-seven Meep runs in 11c are the same class of failure. | FAIL |
| 12h | field comparison, not one forward scalar | PASS | — | see evidence `coated7_local.json` | completed in this closeout | none | h=0.02, h_edge=0.0045, 618177 DOFs, 11 cells across the quartz wall. Forward −0.0090 dB −0.243°; backward +0.017 dB −0.037°; ±60° −0.011 dB +0.06°; gap −0.020 dB −0.333°; outside shell −0.020 dB −0.332°; near quartz −0.020 dB −0.335°. Previous outside-shell phase at h_edge=0.008 was −0.412°. Last local step is 0.08°. | PASS |
| 13b | oblique propagation | FAIL | class F harness. This box does not isolate PML reflection. Cylinder and pair sweeps are 13c and 13d. | see evidence `oblique_pml.json` | completed in this closeout | none | physical-region rel L2 stays 0.020 to 0.025 while dpml, sigma scale, and length change. Probe reflection magnitude stays 0.067. The error does not track the PML. The first Dirichlet-through-PML run (rel L2 0.26) is discarded. | FAIL |
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
| 20g | Meep: resolution | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | plasma Meep forward error does not fall from res 24 to res 40 on the coated bulb, nor from res 30 to res 40 on three cylinders. | FAIL |
| 20h | Meep: grid registration | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | ox=0 and ox=0.25 kept separate. Pair 0° forward spread is 0.8 dB and 20°. Seven-bare res 36 spread is several dB. | FAIL |
| 20i | Meep: runtime | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | runs to t≈960–980. Late plasma values remain outside the analytic tolerance. | FAIL |
| 20j | Meep: dispersive-interface error | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | quartz Meep matched the slab. Plasma cylinders, coated bulbs, and clusters do not. | FAIL |
| 20k | Meep: DFT convergence | FAIL | class D Meep. | see evidence `meep_closeout.json` | completed in this closeout | none | the finished DFTs are the late values quoted in 7h, 9c, 10c and 11c. They are not the analytic field. | FAIL |
| 20o | analytic: conditioning | PASS | raw 2-norm condition remains up to 1e26 | see evidence `tmatrix_condition.json` | completed in this closeout | none | equilibration is a reference-solver scale, not a fit to FEM. Scaled condition of coated seven is 3.3e3 at m_max=8 and 2.4e4 at m_max=14. Forward observable change from scaling is below 1e-11. | PASS |

## Failure classes

A. FEM physics or operator failure: none. Coated-cylinder phase at 5.30 GHz moves toward the analytic value as h decreases.

B. FEM observable reconstruction: 2f is resolved. Element-center Ey/Hz is the P1 gradient. It is order 1, while Hz is order 2. One-sided differences are slower (vacuum 0.252% at h=0.003125). Nodal recovery is faster (vacuum 0.021% at that mesh). No smoothing was added to the production solver.

C. Insufficient refinement: 3d and 7f. Guide gradient L2 is 0.105% and the element max is 0.115% at h=0.0008, against a frozen 0.1% bar. The median is 0.0987%. Coated cylinder at 5.30 GHz is −0.102 dB and −4.73° at h=0.012. Production 3.85 GHz passes.

D. Meep disagreement with the analytic continuum: 6m, 6n, 7h, 7i, 9c, 10c, 11c, 12g, 20g, 20h, 20i, 20j, 20k. Registrations were not averaged. Higher resolution did not systematically approach the analytic field.

E. Inappropriate pass criterion near a null: not used.

F. Harness: 13b. The oblique box rel L2 stays near 2.4% and does not move with PML thickness or sigma. Cylinder and pair PML sweeps are 13c, 13d, and 13h, and those pass. Other harness items that were not scored as rows: `beta_rel` for complex kx, the first Dirichlet-through-PML oblique run, the `eh_orders.py` power strip, and a y-directed source placed inside the PML.

## Decision

**B0_FEM_PARTIALLY_VALIDATED**

Mandatory FEM items that do not pass: 3d, 7f, and 13b. 13b is a harness limit, not a cylinder-solver error. 3d and 7f are refinement shortfalls with the observed P1 rates, not sign errors. Meep FAIL rows do not change the FEM label. B ≠ 0, adjoints, and optimization were not started.

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

