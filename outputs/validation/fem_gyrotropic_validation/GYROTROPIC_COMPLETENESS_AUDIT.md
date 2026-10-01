# Gyrotropic completeness audit

Final status: **GYROTROPIC_FEM_VERIFIED_AND_VALIDATED**

B=0 remains frozen. This audit does not reopen it. Production source was not modified. The acceptance file `GYROTROPIC_PASS_CRITERIA.md` is unchanged since it was written, before the comparison runs.

| Required check | Evidence | Status |
|---|---|---|
| Where `0.250 a` and `0.230 a` come from | `GYROTROPIC_GEOMETRY_PRECHECK.md` | PASS |
| Production material radius | `0.230 a` in `Add_Bulb` profile 0 | PASS |
| Independent tensor derivation, including collisions and `ρ` | `GYROTROPIC_TENSOR_DERIVATION.md` | PASS |
| Identity grid: frequency, `|B|`, sign, collisions, plasma frequency | `tensor_identities.json`, 75 cases, residuals `≤ 2e-15` | PASS |
| Physical-`ρ` manufactured solution, both `B` signs | `gyrotropic_mms.json`, orders 1.999 and 0.999 | PASS |
| Homogeneous analytic `β`, polarization, `E`, `B → 0` | `gyrotropic_homo.json` | PASS |
| Faraday analytic sweep and one Meep comparison | `gyrotropic_faraday.json`, stored `faraday_res64.json` | PASS |
| Statement that the Hz FEM is Voigt, not `k ∥ B` | derivation and report | PASS |
| `ε(B) = ε(-B)^T` and the FEM interior operator | identities and `gyrotropic_onsager.json` | PASS |
| Oblique interface, complex `R`, several `k_y` and `B` signs | `gyrotropic_oblique.json` | PASS |
| Slab, several thicknesses, `±B`, production frequency, `B → 0` | `gyrotropic_slab.json`, `gyrotropic_slab_sweep.json`, `gyrotropic_prod_slab.json` | PASS |
| Dissipation formula, sign of `B`, collision rate `→ 0` in the volume | `gyrotropic_passivity.json` and the plane-wave residual `1.6e-9` | PASS |
| Port relation `S_ij(+B) = S_ji(-B)`, and same-`B` nonreciprocity stated separately | `gyrotropic_ports.json` | PASS |
| One magnetized cylinder, near-field probes, `B` reversal, `B → 0` | `gyrotropic_cylinder.json` | PASS |
| Two and three cylinders, mesh, orientation, Onsager | `gyrotropic_cluster.json` | PASS |
| Meep classified apart from the analytic reference | Faraday cell only; cylinder Meep not treated as truth | PASS |
| Small multiport before any 91-bulb model | three nodal ports plus the straight-guide diagnostic | PASS |
| 91-bulb magnetized forward run | not run | NOT RUN |
| Adjoints and optimization | not run | NOT RUN |

No required FEM-versus-analytic row is FAIL. The straight-guide full-width integral is recorded so it is not mistaken for a reciprocal tensor. Thresholds were not changed after the results.
