# Gyrotropic validation matrix

Scope is the 2D Hz finite-element operator with a cold-plasma tensor, `B` along `z`. Thresholds are the ones frozen in `GYROTROPIC_PASS_CRITERIA.md` before these numbers were computed. The 91-bulb circulator was not run.

| Item | Result | Status |
|---|---|---|
| Production core radius | `Add_Bulb` profile 0 uses `4.6/6.5 * 0.325 a = 0.230 a`. `sixport_common.r_plasma = 0.250 a` is unused rod metadata | PASS |
| Tensor vs scalar Drude at `B = 0` | worst absolute difference `1.78e-15`, including the stored production value | PASS |
| `ε_xx = ε_yy`, `ε_xy = -ε_yx`, diagonal even in `B`, off-diagonal odd in `B` | exact on 75 cases | PASS |
| `ρ ε = I` and analytic inverse vs `numpy.linalg.inv` | worst residuals `7.5e-16` and `2.3e-15` | PASS |
| Onsager `ε(B) = ε(-B)^T` | residual `0` | PASS |
| Passive dissipation eigenvalues | minimum eigenvalue `0` when lossless and `≥ 1.2e-5` on the lossy grid | PASS |
| Manufactured solution, physical `ρ`, P1 | L2 order `1.999`, H1 order `0.999` at `n = 64` | PASS |
| Homogeneous Voigt wave | Hz relative L2 `1.84e-5` at `h = 0.01`. `β(+B) = β(-B)`. `E_x` reverses with `B` and matches the analytic polarization | PASS |
| Faraday `κ` along `B` | Lorentz `κ(+B) = -0.078460/a`. Stored Meep res 64 is `-0.078479/a`. Difference `1.87e-5`. `κ(0) = 0`. `κ(-B) = -κ(+B)` | PASS |
| Hz FEM as a Faraday propagator | `∂/∂z = 0`, so `k` is perpendicular to `B`. This solver is not asked to rotate a wave along `z` | PASS |
| Interior FEM matrix Onsager | interior-block residual `0` at `B = 0` and at `fc = ±0.05`. Green function `2.5e-15`. Same-`B` nonreciprocity on that mesh `1.3e-3` | PASS |
| Oblique magnetoplasma interface | FEM field relative error `≤ 4.3e-4`. `R(k_y,+B) = R(-k_y,-B)` to `0`. `R(k_y,+B)` differs from `R(-k_y,+B)` | PASS |
| Normal slab, propagating | field relative error `1.17e-4`. Analytic `T(+B) = T(-B)`. Stack residual `≤ 2.7e-16` over four thicknesses | PASS |
| Normal slab, production 3.85 GHz | field relative error `1.08e-5`. `R` matches the B=0 stack routine to `2e-16`. `±B` transmissions agree | PASS |
| Plane-wave dissipation identity | `dS_x/dx + p_abs` residual `1.6e-9` | PASS |
| FEM absorbed power | positive, linear in the collision rate, and equal for `+B` and `-B` | PASS |
| FEM contour at collision rate `0.01` | outward flux plus volume absorption within `0.54%` of the absorbed power at `h = 0.01` | PASS |
| One gyrotropic cylinder | FEM versus independent `T_n`: worst about `1.0e-4 dB` and `7.6e-4 deg`. `B = 0` reduces to isotropic Mie within `6e-17`. Side probes swap under `B → -B` | PASS |
| Two and three cylinders | FEM versus Graf T-matrix: worst `2.6e-4 dB` and `2.5e-3 deg`. Finer mesh `7.2e-5 dB`. Analytic Onsager `≤ 7e-16`. Flipped pair matches | PASS |
| Three nodal ports | `S_ij(+B)` versus `S_ji(-B)` at `6e-15`. Same-`B` nonreciprocity `0.63%`, stable from `h = 0.025` to `0.0125`. At `B = 0` the residual is `9e-16`. Repeat solve residual `0` | PASS |
| Straight-guide even-mode projection | `±B` difference `2e-5`, comparable to the mesh residual. The full-width integral cancels the odd scattered field. Not used as evidence of reciprocity | PASS |
| Meep cylinder | not used as a reference. The B=0 campaign already separated finite-grid Meep from continuum plasma cylinders | PASS |
| Full 91-bulb magnetized circulator | not run | NOT RUN |
