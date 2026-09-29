# Horn localization

**Status: FAILED**

B≠0, adjoints, and optimization were not run. Quartz shells and the 91-bulb device were not run, because the six-horn test failed the 0.25 dB gate.

Artifacts are in `horn_localization/`. The polygons FEM solved are `horn_localization/meep_horn_polygons.json`, written from `horn_for_port` → the four quads that `_mount_horn_walls` passes to `Add_Prism` when `horn_walls="prism"`. There is one vertex list.

---

## Discrepancies found

| # | What | What was done |
|---|---|---|
| 1 | FEM metal mask was a tube of radius `0.55 * wall_thickness` around prism edges. Area Jaccard vs the filled prisms was 0.32 in the earlier full-device raster and 0.444 on a 0.015a raster here. | FEM horns-only solves mark triangles whose centroids lie inside the exported quads. The two polygon definitions are the same arrays. Shoelace area of those quads is **0.007152 m²**, equal to the expected prism area. |
| 2 | A one-horn FEM box placed entirely outside the PML bands produced a closed cavity. Time-average flux was exactly 0. | The local box was mapped onto `[0, L]` so the outer 1.0a is a PML. Inward/monitor ratio then matched Meep. |
| 3 | Adding polygon vertices onto a grid produced duplicate nodes and an exactly singular matrix. | Coordinates were uniquified; empty rows are pinned. |
| 4 | The first six-horn Meep job deadlocked after the timestep loop because only rank 0 called `get_dft_array`. | Every rank calls it. Powers below are from the completed rerun. |
| 5 | Six-horn port powers still differ by up to 1.25 dB after the geometry fix. | Stopped. No bulbs added. No scale factor applied. |

---

## 1. Was the previous FEM horn geometry equivalent to Meep?

No. The vertex list in `full_horns` was already the Meep prism list, but the metal region FEM actually used was not those solids. `wall_triangle_mask` painted a neighborhood of the boundary polyline. Filling the true quads on the earlier full-device FEM-M run did not recover the Meep port distribution; that run still contained the 91 bulbs and a non-converged mesh.

## 2. Why was the Jaccard index only ~0.32?

The comparison was filled prism versus edge tube, on the same quads.

- Capture radius 0.55 × 4 mm = 2.2 mm. Meep half-thickness is 2.0 mm (`stage2/geometry_compare.json`).
- Earlier raster: polyline-tube area 0.00814 m², polygon area 0.00747 m², intersection small enough that Jaccard = **0.323** (`stage5/polygon_walls_FEM-M.json`).
- This raster (pitch 0.015a): filled area 17.82 a², tube area 40.11 a², Jaccard **0.444**.
- Identical filled-vs-filled Jaccard on the exported vertices is **1.0**.
- Maximum coordinate discrepancy between the Meep `Add_Prism` vertices and the vertices FEM consumed is **0**. The Meep log for the six-horn run prints the same numbers (aperture corner `(4.70513, 2.6)`, throat `(9.15513, 1.2)`).

The pictures can look alike because the tube follows the wall centerline neighborhood. The areas are not the same.

Overlay: `horn_localization/meep_horn_polygons.png`.

## 3. What PEC condition does Meep’s prism correspond to?

`Add_Prism(..., PEC=True)` appends `mp.Prism` with `mp.perfect_electric_conductor`. The six-horn log shows a finite-area prism, height infinite, epsilon diagonal `(-1e20, -1e20, -1e20)`. Meep staircases that material onto the Yee grid. It is a volume of enormous permittivity, not a zero-thickness boundary and not a forced `Hz = 0` condition.

For this TE polarization (only `Hz`, `e^{-iωt}`, `Ez = 0`):

- `Ex = (i /(ω ε)) ∂Hz/∂y`, `Ey = -(i /(ω ε)) ∂Hz/∂x`.
- PEC requires `n × E = 0`. On a wall with in-plane normal, that sets the tangential `E` to zero, which forces `∂Hz/∂n = 0` on the air side.
- `Hz = 0` would instead be the PMC condition `n × H = 0`, because `Hz` is tangential to a z-extruded wall.

## 4. Was FEM enforcing that condition?

In the solves reported here, yes, as a volume condition. Triangles inside the exported prisms have `ρ = 0` (`ε → ∞`). The weak form has no explicit boundary integral, so the natural condition on the air side of that interface is `∂Hz/∂n = 0`. No Dirichlet `Hz = 0` row replacement was applied.

The straight-guide test (transmission difference **6.9×10⁻⁵ dB**) and the one-horn flare ratio below are the checks that this matches Meep’s PEC blocks/prisms for a single guide. The old full-device runs that zeroed `Hz` on wall nodes were a different, incorrect condition; those numbers are not used here.

## 5. Does one-horn FEM match Meep?

The flare ratio does. One port-1 prism, no bulbs, cosine `Hz` line, resolution 40 in Meep, FEM `h = 0.04a`.

| Quantity | Meep | FEM | Δ |
|---|---:|---:|---:|
| `\|P(throat − 1.2a)\| / \|P(monitor)\|` | 0.7507 | 0.7562 | **0.032 dB** |
| `\|P(source + 1.2a)\| / \|P(monitor)\|` | 0.8679 | 1.0978 | **1.02 dB** |

Both monitors on the horn side of the source agree. The outward monitor, which sits toward the PML behind the source, does not. That 1.02 dB is a source-termination difference (FEM PML thickness 1.0a and a quadratic conductivity, Meep PML 0.8a). It is not a flare-shape difference. Signs agree: power at the monitor flows inward (negative along the outward normal) in both solvers.

## 6. Does six-horn, no-bulb FEM match Meep?

No. Stopped here.

Meep: `device_mode="horns_only"`, `num_mode_guide_normal`, 25 ppc (res 50), `run_time=20`, production straight-feed `P_inc = 2909.721`. FEM: same polygons, same numerical-mode samples, guide-normal `∫ S·n̂`, each mesh normalized by its own straight-feed `P_inc`. Source P1. FEM-H is the finest of the three existing meshes.

| Port | Meep power | FEM-H power | Meep dB | FEM-H dB | FEM−Meep dB |
|---|---:|---:|---:|---:|---:|
| P1 | −0.9860 | −0.9237 | — | — | — |
| P2 | 0.037632 | 0.050140 | −14.244 | −12.998 | **+1.246** |
| P3 | 0.085539 | 0.098278 | −10.678 | −10.075 | **+0.603** |
| P4 | 0.696684 | 0.569300 | −1.570 | −2.447 | **−0.877** |
| P5 | 0.085539 | 0.093173 | −10.678 | −10.307 | **+0.371** |
| P6 | 0.037632 | 0.049075 | −14.244 | −13.091 | **+1.153** |

P1 is negative at the source monitor in both solvers. The cached Meep `FluxData` is still empty, so `load_minus_flux_data` did not produce a reflection coefficient. P1 is omitted from the dB gate. The transmission ports are the comparison.

Mirror symmetry across the P1–P4 axis requires P2 = P6 and P3 = P5. Meep matches to the printed digits. FEM-H splits P2/P6 by 0.09 dB and P3/P5 by 0.23 dB. The mesh is not symmetric.

The errors are not one scale factor: P4 is low while P2 and P6 are high.

## 7. Does horn-only FEM converge with mesh refinement?

Not on these three meshes, and the opposite port moves away from Meep.

| Grade | DOFs | min metal edge | factor time | P1→P2 | P1→P3 | P1→P4 | P1→P5 | P1→P6 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| FEM-L | 89062 | 6.4 µm | 1.46 s | 0.04743 | 0.08106 | 0.6563 | 0.08905 | 0.04080 |
| FEM-M | 234307 | 2.7 µm | 2.97 s | 0.05395 | 0.09470 | 0.6142 | 0.09090 | 0.05644 |
| FEM-H | 661890 | 0.66 µm | 11.73 s | 0.05014 | 0.09828 | 0.5693 | 0.09317 | 0.04907 |
| Meep res 50 | — | Yee 0.4 mm | 182 s wall | 0.03763 | 0.08554 | 0.6967 | 0.08554 | 0.03763 |

The minimum edge is the shortest edge of a metal triangle, including slivers, so it is smaller than the mesh’s design `h_crit` (0.8 / 0.4 / 0.2 mm). P14 changes by 0.33 dB from FEM-M to FEM-H (0.614 → 0.569) while Meep is 0.697. Refinement of this family is not approaching the Meep column. FEM-L is closer on P4 (−0.26 dB) but still misses P2 by 1.00 dB, so no grade passes 0.25 dB on every meaningful port.

## 8. Where do the complex fields first diverge?

`horn_localization/field_lines_horns_only.png`, `horn_localization/field_line_compare.json`. FEM-M versus Meep DFT, one complex scale removed so absolute source strength cancels. Inner 80% of each line, to drop samples that snap onto the metal:

| Line | rel L2 of Hz | phase std |
|---|---:|---:|
| source | 0.025 | 1.37° |
| monitor | 0.013 | 0.05° |
| throat | 0.009 | 0.29° |
| center, between horns | **0.112** | **4.9°** |

The feed, the monitor, and the throat agree. The first line that does not is the open segment through the middle of the six-horn cavity. The lobe count and the phase slope still match there; the lobe amplitudes do not. That is the same pattern as the port table: power is leaving the opposite horn low and the side horns high. The horn cross-section itself is not where the solutions separate.

## 9. Do quartz-only structures match?

Not run. Six-horn transmission already exceeds 0.25 dB.

## 10. Does B=0 plasma then match?

Not run.

## 11. What is the corrected full-device B=0 port error?

There is no corrected full-device number. The earlier FEM-H column (max |Δ| 9.32 dB) used a different wall mask and included the bulbs. It is not reused.

## 12. Status

**FAILED**

What is established:

- Ex/Ey and `∫ S·n̂` on a straight guide.
- The Meep prism vertices, now shared.
- The TE PEC condition: Neumann `∂Hz/∂n = 0` via `ρ → 0`, consistent with Meep’s `ε = −1e20` prisms.
- One horn’s flare ratio, to 0.03 dB.

What failed:

- Six-horn coupling across the open cavity, by up to 1.25 dB, with the mesh trend moving away from Meep.

Next physical question, still at B=0 and still without bulbs: why the open region between six correct horns does not propagate the same way (PML profile, numerical dispersion on the coarse cavity mesh, nearest-node source inside a multi-port box). That is not a license to retune amplitudes until the ports match.
