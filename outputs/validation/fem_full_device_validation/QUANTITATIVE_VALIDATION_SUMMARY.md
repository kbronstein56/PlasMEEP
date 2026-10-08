# Quantitative versus structural validation

Jesse’s question: the earlier “130 PASS” count mixed two kinds of rows. This note separates them. A quantitative row has an independent reference, a measured error, and a threshold that was frozen before the number was used. A structural row checks a sign, a symmetry, a convention, or a procedure, and it does not by itself bound a port observable.

The historical B0 reaudit is unchanged: 144 rows, 130 PASS, 13 `MEEP_CROSS_CODE_FAIL`, 1 superseded harness, 0 unresolved. Of those 144 result cells, 119 contain a numeral (108 PASS, 11 Meep failures) and 25 point at another row or state a procedure (22 PASS, 2 Meep failures, 1 superseded). The 22 procedural PASS rows are not extra physics. Meep failures are cross-code disagreements. They are not charged as FEM physics failures, and Meep is not the continuum reference.

Thresholds below are the ones already frozen in `PASS_CRITERIA.md` and `GYROTROPIC_PASS_CRITERIA.md`. Nothing here loosens them.

## B0, quantitative

| Problem | Reference | Observable | Threshold | Measured | Mesh | Status |
|---------|-----------|------------|-----------|----------|------|--------|
| Manufactured solution, real and complex ρ, including a full tensor | imposed Hz | L2 and H1 order | L2 → 2, H1 → 1 | L2 order 1.997–1.998, H1 order 0.999, L2 3.9e-4–4.3e-4 at n=64 | n=64 | PASS |
| PEC parallel-plate guide | exact mode | field L2, power | ≤ 0.01 dB, field ≤ 0.1% | power relative 4.25e-4; Hz L2 2.97e-4; finest E/H element max 0.090% | up to 3.91e6 DOF | PASS |
| Fresnel interfaces | closed-form r, t | interior field | ≤ 0.1% | max relative L2 8.31e-4 | h=0.025 | PASS |
| Quartz, plasma, and coated slabs | transfer matrix | T, phase | ≤ 0.01 dB, ≤ 0.1° | worst about 0.0015 dB and 0.025° | h=0.02 | PASS |
| Bare plasma cylinder | Mie | probe dB and phase | ≤ 0.05 dB, ≤ 0.5° | worst −0.0015 dB, −0.031° | 1.68e6 DOF | PASS |
| Coated bulb at 3.85 GHz | multilayer T-matrix | probe dB and phase | ≤ 0.05 dB, ≤ 0.5° | worst −0.0037 dB, −0.081° | 6.05e5 DOF | PASS |
| Coated bulb at 5.30 GHz | multilayer Mie | forward probe | ≤ 0.05 dB, ≤ 0.5° | −0.0118 dB, −0.461° | 5.94e6 DOF | PASS |
| Two and three cylinders, seven-coated forward | T-matrix | probe dB and phase | ≤ 0.10 dB, ≤ 1° | worst retained passes are a few 0.01 dB and a few 0.1° | h=0.02 and FEM-X | PASS |
| PML slab and oblique Bloch fit | transfer matrix / fit | reflection | slab ≤ 0.01 dB; fit residual reported | slab ≤ 0.0013 dB; oblique \|R\| down to 1.3e-4 | one slab mesh; h=0.02 fit | PASS |
| Quartz reciprocity | G(a,b)=G(b,a) | reaction | roundoff | max relative 1.7e-15 | h=0.025 | PASS |
| Linear solve | — | \|Ax−b\|/\|b\| | roundoff | ~1e-14 | scatter solves | PASS |

Meep plasma cylinders, coated bulbs, and clusters remain `MEEP_CROSS_CODE_FAIL` against the analytic field. That is a finite-grid disagreement, kept separate from the FEM rows.

## Gyrotropic, quantitative

| Problem | Reference | Observable | Threshold | Measured | Mesh | Status |
|---------|-----------|------------|-----------|----------|------|--------|
| Tensor at B=0 vs scalar Drude | scalar formula | ε entries | roundoff | 1.78e-15 | — | PASS |
| ε(B)=ε(−B)^T, ρ ε = I | algebra | residual | roundoff | 0 and 7.5e-16 | 75 cases | PASS |
| Manufactured solution | imposed Hz | L2, H1 order | → 2 and → 1 | 1.999 and 0.999 | n=64 | PASS |
| Homogeneous Voigt wave | analytic β and polarization | Hz L2 | field ≤ 0.1% | 1.84e-5 at h=0.01 | h=0.01 | PASS |
| Oblique interface | analytic R | field and R(k_y,+B)=R(−k_y,−B) | field ≤ 0.1% | field ≤ 4.3e-4; identity 0 | graded | PASS |
| Normal slab, including 3.85 GHz | stack | field | ≤ 0.1% | 1.17e-4 and 1.08e-5 | graded | PASS |
| Plane-wave dissipation | dS_x/dx + p_abs = 0 | residual | reported | 1.6e-9 | analytic | PASS |
| One gyrotropic cylinder | independent T_n | probe dB, phase | ≤ 0.05 dB, ≤ 0.5° | 1.0e-4 dB, 7.6e-4°; B=0 vs Mie 6e-17 | h=0.02 | PASS |
| Two and three cylinders | Graf T-matrix | probe dB, phase | ≤ 0.10 dB, ≤ 1° | worst 2.6e-4 dB, 2.5e-3° | h=0.02 and finer | PASS |
| Three nodal ports | S_ij(+B)=S_ji(−B) | reaction | roundoff vs discretization | 6e-15 | h=0.0125 | PASS |
| Contour absorption | outward flux + volume absorption | relative balance | the recorded 0.54% | 0.54% at h=0.01 | h=0.01 | PASS |

## Full device, quantitative

Frozen port target for a significant channel, carried from the full-device campaign and not loosened: ≤ 0.10 dB and ≤ 1 degree between the last two acceptable meshes. Weak channels use absolute power.

| Problem | Reference | Observable | Threshold | Measured | Mesh | Status |
|---------|-----------|------------|-----------|----------|------|--------|
| Geometry | production builders | bulb centers | exact | max difference 3.55e-15 | — | PASS |
| B=0 operator bridge | scalar assembly on the same mesh | \|A−A_scalar\|/\|A\| | roundoff | 0 on mesh C | 8.23e5 nodes | PASS |
| B=0 reaction symmetry | S = S^T | max \|S−S^T\| | roundoff | 2.5e-14 on F | 4.14e6 nodes | PASS |
| B=0 port convergence | mesh F versus mesh M | normalized port power and monitor phase | ≤ 0.10 dB, ≤ 1° | through +0.226 dB and +1.03°; next −0.31 dB and +1.24°; adjacent −0.14 dB and +1.87° | M 1.85e6, F 4.14e6 | FAIL |
| Air-versus-interface diagnostic | mesh R versus F | same ports | not a pass gate | R refines quartz to 0.004 a and coarsens air to 0.055 a. Through power moves from 0.424 back to 0.388, toward mesh M. The open B=0 error is not cured by interface refinement alone | R 4.90e6 | diagnostic, not a pass |
| B≠0 port convergence | F versus M at ±0.05 T | same | ≤ 0.10 dB, ≤ 1° | every receiving port ≤ 0.083 dB and ≤ 0.87° | M, F | PASS |
| Reaction Onsager | S_ij(+B)−S_ji(−B) | max and RMS | discretization, not 0.1 dB | max 1.17e-13 on F; < 8e-14 on M at four \|B\| | F 4.14e6 | PASS |
| Same-B nonreciprocity | S−S^T at +B | relative size | a nonzero result is the physics | 0.108 of the largest reaction entry | F | PASS |
| Absorption sign | (ε−ε†)/(2i) | P_abs | ≥ 0 | 1.865e-3 at B=0; 6.997e-4 at both ±0.05 T | F | PASS |
| Horn power plus absorption | not a closed surface | residual / driven flux | not claimed as balance | 0.11 at B=0, 0.23 at ±0.05 T | F | UNRESOLVED |
| B continuation | samples 0 to 0.05 T | P2−P6 and absorption | smooth, contrast → 0 at B=0 | contrast −2e-6, +8.8e-3, +3.8e-3, +7.4e-4, +3.7e-4; absorption stays positive | M | PASS |
| Frequency continuation | 0.95 to 1.05 f0 at +0.05 T | port power | smooth tensor; steep device response is allowed | Re ε_xx monotone; through power 0.042 at 3.754 GHz and 0.0079 at 3.850 GHz | M | PASS as a sweep, operating point is steep |

## Counts

| Class | PASS | FAIL | UNRESOLVED |
|-------|-----:|-----:|-----------:|
| B0 quantitative result cells | 108 | 0 FEM; 11 Meep cross-code | 0 in the frozen reaudit |
| B0 structural / procedural cells | 22 | 0 FEM; 2 Meep statements | 0; 1 superseded harness |
| Gyrotropic quantitative rows in the matrix | 19 | 0 | 0 |
| Gyrotropic structural rows | 3 | 0 | 0 |
| Full-device quantitative rows above | 8 | 1 | 1 |
| Full-device structural rows | 2 (sixfold layout, PEC = Neumann) | 0 | 0 |

The full-device FAIL is the B=0 port step between meshes M and F. The full-device UNRESOLVED item is the horn sum, which is not a closed energy balance. Both block `FULL91_GYROTROPIC_FORWARD_VALIDATED`. Differentiability is not started while either remains.
