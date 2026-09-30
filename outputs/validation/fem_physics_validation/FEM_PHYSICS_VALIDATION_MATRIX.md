# FEM physics validation matrix

Status: **B0_FEM_PARTIALLY_VALIDATED**.

The row-by-row completeness record is `VALIDATION_COMPLETENESS_AUDIT.md`
(144 rows: 128 PASS, 16 FAIL, 0 UNRESOLVED). This file is the numeric
benchmark table. Thresholds are frozen in `PASS_CRITERIA.md`. The closeout
that produced those counts is `B0_VALIDATION_CLOSEOUT.md`.

Analytic results are the reference. A Meep disagreement is not charged
against an FEM row that already matches that analytic result.

Complex and phase errors below are FEM minus analytic on vacuum-normalized
Hz ratios, unless the row says otherwise. No global complex scale was fitted.

## Benchmark table

| Benchmark | Parameters tested | Analytic reference | FEM meshes | Meep cases | Complex / magnitude error | Phase error | Power error | Field L2 | Conservation residual | Reciprocity residual | Numerical uncertainty | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MMS, real ρ | unit square, k0=1.7, n=8…64 | imposed Hz | 4 | — | nodal max 3.0e-4 at n=64 | — | — | L2 4.25e-4, order 1.998 | — | — | H1 order 0.999 | PASS |
| MMS, complex scalar ρ | same | imposed Hz | 4 | — | nodal max 3.2e-4 | — | — | L2 order 1.998 | — | — | H1 order 0.999 | PASS |
| MMS, diagonal complex ρ | same | imposed Hz | 4 | — | nodal max 3.5e-4 | — | — | L2 order 1.997 | — | — | H1 order 0.999 | PASS |
| MMS, full complex ρ | off-diagonal tensor | imposed Hz | 4 | — | nodal max 2.2e-4 | — | — | L2 3.90e-4, order 1.998 | — | — | H1 order 0.999 | PASS |
| Homogeneous ε=1 | plane wave | exact kx, Sx | h=0.05, 0.025, 0.0125 | not run | — | — | Sx rel 6.7e-5 | rel L2 4.74e-5 | — | — | E/H sample 0.34% | PASS field; FAIL E/H vs 0.1% |
| Homogeneous ε=3.8 | plane wave | exact | same | not run | — | — | Sx rel 2.4e-4 | rel L2 8.30e-4 at h=0.025 | — | — | — | PASS |
| Homogeneous ε=2+0.3i | plane wave | exact | same | not run | — | — | Sx rel 1.3e-4 | rel L2 4.02e-4 at h=0.025 | — | — | `beta_rel` in the JSON is not valid for complex kx | PASS |
| Homogeneous plasma ε(fs) | 3.85 GHz | exact | same | not run | — | — | Sx rel 2.4e-4 | rel L2 1.39e-4 | — | — | E/H sample 0.60% | PASS field; FAIL E/H vs 0.1% |
| PEC guide, hardest mode | width 1, k0 5 | exact mode and power | h=0.04 down to 0.00125 | not run | — | — | power rel 4.25e-4 | rel L2 2.97e-4, 1.026e6 DOFs | — | — | order about 2; E/H at h=0.01 is 0.30–3.3% | PASS field and power; FAIL E/H; Meep UNRESOLVED |
| PEC guide matrix | 3 widths × 3 k0 | exact | 27 rows | not run | — | — | stations agree to ~3e-7 | cutoff rel L2 ~1.7e-4 at h=0.01 | — | — | — | PASS |
| Fresnel set | air/quartz, quartz/air, complex, plasma; 0°, 15°, 40° | closed-form r, t | h=0.025, 24 rows | not run | r, t imposed on the boundary | — | analytic R+T=1 | max interior rel L2 8.31e-4 | — | — | one mesh | PASS |
| Quartz slab | thickness 0.20 and 0.40; 3 frequencies | transfer matrix | h=0.02 | res 24 and 40 | T +0.0013 dB | −0.0007° | R+T=1 | — | 0 | — | Meep +0.0024 dB, +0.008° at res 40 | PASS |
| Plasma slab | thickness 0.20 and 0.40; 3 frequencies | transfer matrix | h=0.02 | not in this package | T −0.0015 dB (0.40) | +0.025° | absorption proxy +3.2e-4 | — | R+T=0.99968 | — | — | PASS |
| Quartz/plasma/quartz | 0.10/0.40/0.10; 3 frequencies | transfer matrix | h=0.02 | not run | T −0.0015 dB | +0.015° | — | — | — | — | h=0.04 is an invalid grid, not a solver fail | PASS |
| Bare disk, r=0.230, plasma | 3.85 GHz; forward, back, sides, near | Mie | h=0.04, 0.02, 0.012 | prior plasma, not this package | worst probe −0.0015 dB | −0.031° | analytic absorption +1.34e-5 | old cut rel L2 0.052% | — | — | 1.68e6 DOFs | PASS FEM; FAIL Meep plasma |
| Bare r=0.12 and r=0.40 | same probes | Mie | h=0.04, 0.02 | not run | −0.0052 dB and −0.0022 dB | −0.130° and −0.062° | — | — | — | — | — | PASS |
| Dielectric cylinder ε=3.8 | r=0.230 | Mie | h=0.02 | res 24 and 40 | FEM −0.0003 dB | −0.007° | analytic contour ~1e-12 | — | ~0 | — | Meep forward +8.6e-5 dB | PASS |
| ε=−1.5 cylinder | r=0.230, near the ε=−1 feature | Mie | h=0.04 fails; h=0.02 passes | not run | h=0.02 worst −0.0195 dB | −0.420° | — | — | — | — | h=0.04 was −0.119 dB, −1.60° | PASS at h=0.02 |
| ε=−8 cylinder | r=0.230 | Mie | h=0.04 only | not run | −0.0083 dB | −0.194° | — | — | — | — | one mesh | PASS at that mesh |
| Mie frequency sweep | 3.2–6.5 GHz, 102 samples | Mie | FEM only at 3.85 GHz | not run | peak max\|T\|=0.999 at 5.30 GHz | — | — | — | — | — | — | UNRESOLVED for FEM |
| Coated bulb | r=0.230 / 0.325 / 0.375, ε_q=3.8 | multilayer T-matrix | h=0.04, 0.02 | not run | worst −0.0037 dB | −0.081° | analytic absorption +1.81e-5 | — | — | 5.5e-14 on 6 pairs | 604838 DOFs | PASS |
| Coated frequencies | 3.50, 4.20, 5.30 GHz | T-matrix | not run | not run | analytic only | — | — | — | — | — | — | UNRESOLVED |
| Two bare cylinders | pitch 1.0 | T-matrix | h=0.04, 0.02 | not run | worst −0.0128 dB | −0.083° | — | — | — | 3.1e-15 on 6 pairs | — | PASS FEM; Meep UNRESOLVED |
| Two cylinders, 30° | one mesh | T-matrix of the rotated pair | h=0.04 | not run | −0.0316 dB | −0.232° | — | — | — | — | one mesh | UNRESOLVED |
| Coated pair | pitch 1.0 | T-matrix | h=0.04 | not run | −0.0477 dB | −0.318° | — | — | — | — | one mesh | UNRESOLVED as a mesh sequence |
| Three coated bulbs | nearest-neighbor triangle | T-matrix | h=0.02 and older C/F | not run | new-domain worst −0.0088 dB; old FEM-F −0.00017 dB | −0.178°; −0.110° | — | old cut rel L2 0.18% | — | — | side +60 and −60 differ because the triangle is not y-symmetric | PASS FEM; Meep UNRESOLVED |
| Seven bare | hex pitch 1.0 | T-matrix | h=0.04, 0.02 | not run | h=0.02 worst −0.0160 dB | −0.224° | — | — | — | — | shell phase step 0.54° | PASS FEM; Meep UNRESOLVED |
| Seven coated, forward, old domain | source (−6,0) | multilayer T-matrix, m_max=12 | FEM-C/F/X | prior, not re-run | FEM-X −0.0023 dB, abs 0.00313 | −0.100° | — | exterior rel L2 0.21% | — | — | last step 0.0030 dB, 0.13° | PASS |
| Seven coated, multi-probe, new domain | forward, back, ±60°, gap, outside shell | same T-matrix | h=0.04, 0.02 | not run | h=0.02 forward −0.0106 dB; shell −0.0248 dB | −0.301°; −0.412° | — | full-field map not made | — | — | shell phase step 0.77°, not yet ≪ 1° | UNRESOLVED |
| Seven coated Meep | original grid, late DFT | the FEM-X / T-matrix value | — | prior registrations, not averaged | about 0.6 dB | about 15° | — | — | — | — | registration spread kept separate | FAIL |
| T-matrix truncation | m_max=2…14; 1, 3, and 7 cylinders | previous m_max | — | — | coated-7 forward change < 3e-5 by m_max=8; < 1e-8 by m_max=14 | — | — | — | — | — | cond ~2e10 at m_max=8; 1e19–1e25 at m_max=12–14 | PASS probe; conditioning UNRESOLVED |
| One-cylinder reduction | nmax=12 | Mie | — | — | 6.8e-21 | — | — | — | — | — | Graf worst ~1e-11 | PASS |
| PML slab | dpml 0.6/1.4, sigma ×0.5/×2, pad shift −0.4 | transfer matrix | one mesh | not run | ≤ 0.0013 dB | ≤ 0.03° | — | — | — | — | cylinder and oblique sweeps not run | PASS slab; other PML rows UNRESOLVED |
| Quartz-inclusion reciprocity | 15 pairs | G(a,b)=G(b,a) | h=0.025 | not run | — | — | — | — | — | max rel 1.7e-15 | — | PASS |
| Linear residual | SuperLU | — | scatter solves | — | — | — | — | — | — | — | \|\|Ax−b\|\|/\|\|b\|\| ~ 1e-14 | PASS |
| Geometry | production disk and quartz annulus | πr² | h=0.012 and h_edge=0.008 | — | — | — | — | — | — | — | plasma area rel 7.8e-5; annulus ~8e-7; ~6 elements across the wall | PASS |

## Solve counts in this package

| Kind | Count | What is counted |
|---|---|---|
| Analytic evaluations stored | 275 | 54 truncation, 15 Mie radius×material, 102 Mie frequencies, 4 coated frequencies, 9 cluster geometries, 4 contours, 12 Fresnel coefficients, 24 Fresnel fields, 30 slabs, 8 homogeneous, 9 guide families, 4 MMS tensors |
| FEM solves | 187 | 16 MMS, 64 planar, 45 canonical (guide, homogeneous power, PML), 7 refined guides, 52 scatterer solves (26 cases × vacuum and object), 3 reciprocity factorizations |
| Meep simulations | 8 | quartz slab and dielectric cylinder, two resolutions each, vacuum reference plus object |

The earlier plasma Meep campaign is not included in the Meep count. It is the evidence for the FAIL rows, and those registrations were not averaged.

Prior production-domain FEM-C, FEM-F, and FEM-X factorizations for the seven-bulb cluster are additional to the 187. FEM-X is 3,028,995 nodes.

## Closeout addition

| Kind | Added | Running total |
|---|---|---|
| Analytic evaluations | 92 | 367 |
| FEM solves | 77 | 264 |
| Meep simulations | 36 | 44 |

The 36 Meep simulations are 2 guide runs and 17 scatter cases × (vacuum + object). They are not averaged across grid registrations. Plasma Meep remains outside the frozen tolerances. The FEM failures that keep the package at `B0_FEM_PARTIALLY_VALIDATED` are the guide gradient at h=0.0008 (L2 0.105%), the coated cylinder at 5.30 GHz (−4.73°), and the oblique PML harness.
