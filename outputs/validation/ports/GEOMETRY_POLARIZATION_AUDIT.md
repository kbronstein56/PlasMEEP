# Geometry, units, and polarization audit — agent/eigenmode-ports

**Date:** 2026-08-27  
**Sources:** `scripts/validation/sixport_common.py`, `Sketchbook_PMMCirculator.ipynb`, beam-steering paper (Rodriguez et al., ACS Photonics 2026), PlasMEEP `Add_Bulb` / `Add_*_Source`.

---

## 1. What is `a`?

| Quantity | Value |
|---|---|
| **`a` in code** | `0.028` **metres** |
| **Physical length of 1 Meep unit** | **2.8 cm** |
| **Is `a` the lattice spacing?** | **No** |

Throughout PlasMEEP/Meep for this project, **`a` is the characteristic length used to nondimensionalize Maxwell’s equations** (Meep “a-units”). Frequencies enter as \(f_a = f\,a/c\).

The triangular-lattice pitch is set separately:

```text
d_exp = 0.020 / a   # 20 mm center-to-center → 0.714286 a-units
```

So **`a = 28 mm = 1.4 × (20 mm lattice pitch)`**. The choice is inherited from the Fullfields / circulator notebook (`a = 0.028  # 1 Meep length unit = 2.8 cm`); it is a unit convention, not a hardware dimension.

**Answer: `a` is 2.8 cm (0.028 m).**

---

## 2. Geometry vs beam-steering paper

Paper (METHODS): 91 discharges, triangular lattice, six bulbs per side, **20 mm** pitch, quartz OD **15 mm**, wall **1 mm**, ID **13 mm**, quartz **ε ≈ 3.8**.

| Feature | Paper | Current Meep harness | Match? |
|---|---|---|---|
| Element count | 91 | 91 (`side_dim=6`) | yes |
| Lattice | triangular / hex | `Rod_Array_Hexagon` triangular | yes |
| Pitch | 20 mm | `0.020/a` → 20 mm | yes |
| Quartz OD | 15 mm | `r_bulb_outer=0.0075/a` → 15 mm | yes |
| Quartz ID | 13 mm | `r_bulb_inner=0.0065/a` → 13 mm | yes |
| Quartz ε | ≈ 3.8 | `Get_Med(3.8)` in `Add_Bulb` | yes |
| Plasma fill radius | (model-dependent) | `r_plasma=5 mm` (~4.6/6.5×ID) | model choice |
| Horns | commercial MW horns (3-port steering) | synthetic PEC flares (6-port); wall 4 mm, aperture 104 mm, throat 48 mm, depth 89 mm + 60 mm feed | **not paper geometry** |
| Ports | 1 in / 2 out | 6-port circulator | different device |

**Discrepancies to keep visible:** six-port PEC horn model ≠ paper’s three commercial horns; plasma radial profile is a PlasMEEP model (not a paper-stated radius); `a` is a unit scale, not a hardware length.

---

## 3. Polarization audit (components, not TE/TM labels)

### Published beam-steering / Ceviche (paper)

> “The source mode E-field is polarized along the length of the plasma discharges (**Ez** with z out of the side of the page).”

Nonzero fields (2D): **Ez, Hx, Hy**. Plots show **|Ez|²**.

### Current Meep baseline and `te1_hz_line`

- Notebook / harness: `component=mp.Hz`, or `Add_Cont_Source(..., E=False)` with `Pol=[0,0,1]`.
- TE1 launch: same **Hz** points with cos envelope on the true tangent.
- Flux / DFT receivers use **Ex, Ey, Hz** (Poynting in-plane).

Nonzero fields (2D): **Hz, Ex, Ey**.

### PlasMEEP naming trap

`Add_*_Source(..., E=True/False)` maps `Pol` to **E** or **H** Cartesian components. Notebook horn block is labeled “TM orientation” (legacy PMMInverse horn sizes) while the **source is Hz** — **do not trust TE/TM labels**.

### Gyrotropy with **B ∥ z** (circulator)

Cold-plasma tensor with **B = Bẑ** couples **Ex ↔ Ey** via ε_g. **Ez** only sees ε_∥ and does **not** experience that Hall/Faraday coupling in 2D.

| Question | Answer |
|---|---|
| Same polarization as published/experimental beam-steering PMM? | **No** — paper/Ceviche/experiment: **Ez**; current circulator model: **Hz**. |
| Wrong for a **magnetized** 2D circulator with **B ∥ discharge axis**? | **No** — **Hz / in-plane E** is the gyrotropy-active polarization. Switching to Ez would remove the intended B≠0 nonreciprocity mechanism in this 2D model. |

**Resolution campaign decision:** keep **`te1_hz_line` (Hz)** fixed. Treat Ez-matching as a **separate** experimental-fidelity task, not as a prerequisite to grid-convergence of the present circulator model.

---

## 4. Measurement-formulation result (already completed)

No alternative receiver **materially** beat `te1_hz_line`+axis flux toward the 0.2 dB gate (cheap res32): best TE1+DFT S·n ≈ **0.32 dB** vs TE1+axis **0.38 dB** vs baseline+axis **0.18 dB**. **Not continuing measurement A/B.**

Completed TE1 P1↔P2: res32 0.38 → res48 1.28 → res64 1.03 → **res96 0.682 dB** (np=4, 7494 s; incident 32.72 vs 33.82, ~3.3%).
