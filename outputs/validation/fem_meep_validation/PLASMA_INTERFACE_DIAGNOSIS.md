# Plasma interface diagnosis

B=0 only. No 91-bulb refinement, no gyrotropic run, no adjoints.

The first geometry on which FEM and Meep diverge under refinement is one bare circular plasma disk. They agree on that disk at 25 and 35 cells per wavelength, and the FEM field matches a point-source Mie series. At 50 cells per wavelength Meep keeps a field ringing on the curved boundary after the source has died, and the scattered-field DFT moves away from FEM. An equal-area square of the same plasma does not ring. One quartz-coated bulb at 25 and 35 cells per wavelength agrees with FEM to about the same tolerance as the bare disk. Neither one disk nor one bulb reproduces the multi-dB full-device failure already present at 35 cells per wavelength.

Plots: `plasma_interface/disk_diameter.png`, `plasma_interface/disk_convergence.png`.

## Setup

Shared frequency `fs = 3.85 GHz` (`fs_a = 0.256844`). Drude parameters are the production values `fp = 8 GHz`, `γ = 1 MHz`. Both codes use

```
ε(fs) = −3.31775987 + 0.00112150 i
```

Meep 1.30.0 `Medium.epsilon(fs)` matches that value to `1e−8` before any time step. The disk radius is the production plasma radius `R = 0.230 a` (`a = 20 mm`). The domain is `18 a × 14 a` with a quadratic PML of thickness `2 a` and `σ_max = 25.904`. The source is a point `Hz` Gaussian at `x = −6 a`, `fwidth = 0.20 fs`. Vacuum and inclusion runs share that source, so the reported quantities are ratios.

The quartz bulb keeps the production stack: plasma core `r = 0.230 a`, vacuum gap out to the inner quartz radius `0.325 a`, quartz shell `ε = 3.8` out to `0.375 a` (1.00 mm wall). The square has the same area, half-side `0.2038 a`, edges parallel to the axes.

FEM is the canonical P1 operator: element-wise `ρ = 1/ε`, PEC excised, Meep PML profile. Two body-fitted meshes, FEM-C (`h = 0.030`) and FEM-F (`h = 0.018`). Meep resolutions: 25, 35, and 50 cells per wavelength (`res = 50, 70, 100`).

## Answers

1. **One bare plasma disk matches at 25 and 35 cells per wavelength, and FEM matches the continuum solution. Meep at 50 cells per wavelength does not.**

   FEM-F has 1,343,200 DOFs. Against FEM-C the forward ratio changes by `+0.0019 dB`, the interior ratio by `+0.0016 dB`, and the field at `x = +0.30 a` by `+0.0042 dB`. A point-source Mie series (Graf expansion of the Hankel wave from `x = −6`, then the usual cylinder coefficients) reproduces the FEM-F ratio `Hz_disk / Hz_vac` with a maximum absolute error `0.0024` on the diameter and `0.0004` at `x = +3`. The air Mie check gives `max |a_n| = 2×10−17`.

   Forward ratio `T = Hz(x=+3, with disk) / Hz(x=+3, vacuum)`:

   | solver | T | |T| vs FEM-F |
   | --- | --- | --- |
   | FEM-F | `1.08364 + 0.18542 i` | 0 |
   | Mie | `1.0995` at phase `+9.73°` | `0.0004` in the complex ratio |
   | Meep 25 | `1.07110 + 0.15695 i` | `−0.134 dB` |
   | Meep 35 | `1.08105 + 0.17909 i` | `−0.028 dB` |
   | Meep 35, run extended to 40 | `1.08321 + 0.17996 i` | `−0.010 dB` |
   | Meep 50 | `1.03931 + 0.32038 i` | `−0.094 dB`, complex distance `0.142` |
   | Meep 50, run extended to 100 | `1.12391 + 0.31939 i` | `+0.529 dB`, complex distance `0.140` |

2. **From 25 to 35 cells per wavelength Meep moves onto FEM. The 50-cell run leaves that trend.**

   Complex distance of `T` from FEM-F: `0.031` (25), `0.0068` (35), `0.142` (50). The same pattern holds just outside the disk (`x = +0.30 a`): `−0.626 dB`, `+0.018 dB`, then `+0.853 dB`. Extending the 35-cell run from 20 to 40 time units after the source tightens `T` to `−0.010 dB`. Extending the 50-cell run from 20 to 100 makes the forward error larger (`+0.529 dB`) while the interior ratio improves to `+0.053 dB`.

   The cause is visible in the time series. After the source is off, vacuum `|Hz|` is below `10−9`. On the disk at 35 cells per wavelength the boundary probes fall through `3×10−3` at `t = 210` to below `10−3` by `t = 230`. At 50 cells per wavelength the same probes are still `0.023` and `0.031` at `t = 290`, while the center is `0.0015`. The hotspot in the DFT sits on the boundary (`r ≈ 0.226` versus `R = 0.230`) and reaches `|Hz| = 4.75`, against `1.89` at 35 cells per wavelength. Waiting longer does not remove it from the DFT.

3. **The square does not do this.** At 50 cells per wavelength the square’s late-time `|Hz|` at `t = 230` is `~10−7`, the same floor as vacuum. There is no boundary ring. The interior ratio matches FEM-F to `−0.026 dB`. Forward `|T|` is `−0.119 dB` at 35 cells per wavelength and `−0.094 dB` at 50, so the amplitude is moving toward FEM, with a remaining phase gap of about `2.4°` (`14.04°` versus FEM-F `16.5°`). The sample just outside a face is still `−0.57 dB` away, and the FEM square itself moved `0.145 dB` there between meshes, so that one station is not a settled comparison. Corner singularities of a negative-ε square are a plausible reason the face value converges slowly. The sharp circle-versus-square difference is the ring: it appears on the curved boundary and is absent on the axis-aligned boundary.

4. **Meep cannot time-step this fixed-frequency ε without a dispersive susceptibility.**

   Inspected API: Meep 1.30.0. `Medium(epsilon=complex)` stores the value, and `Medium.epsilon(fs)` echoes `−3.31776 + 0.0011215 i`, but `Simulation` initialization raises `TypeError: '<' not supported between instances of 'complex' and 'int'`. A real ε plus `D_conductivity` can be chosen so that the documented conductivity formula

   ```
   ε_eff = ε (1 + i σ / (2π f))
   ```

   (`meep/geom.py`, `_get_epsmu`) hits the Drude value exactly: `σ = Im(ε) · 2π fs / Re(ε) = −5.455×10−4`, and `Medium.epsilon(fs) − ε_Drude = 0`. A 4×4 cell at resolution 10 then diverges: after `until = 5`, `|Hz| ~ 2×10^18`. The same run in air stays at `|Hz| ~ 7×10−7`. Meep warns that `ε < 1` may require a Courant change. No other nondispersive fixed-frequency object in this installation runs an overdense plasma. Drude is the only stable representation that realizes `ε(fs)`.

5. **There is no Meep fixed-ε result to compare with Drude.** The conductivity material matches `ε(fs)` algebraically and is unstable in the time stepper. It was not used for the disk.

6. **FEM fixed-ε matches the Mie series and matches Meep Drude at 25 and 35 cells per wavelength.** The unstable conductivity run is not a third solver.

7. **One quartz-coated bulb matches at 35 cells per wavelength.** FEM-F (1,342,729 DOFs) versus FEM-C: forward ratio `+0.0019 dB`, interior `−0.0005 dB`, just outside `+0.006 dB`. Meep minus FEM-F on the forward ratio: `−0.154 dB` at 25 cells per wavelength, `−0.025 dB` at 35. Just outside the shell: `−0.946 dB` then `−0.001 dB`. Interior: `−0.087 dB` then `−0.040 dB`. A 50-cell bulb was not run. The bare disk at that resolution is already contaminated by the boundary ring, so a 50-cell bulb DFT would not be a clean nested-interface test.

8. **On the clean runs the residual is smallest at the center and in transmission, and larger in the small reflected wave. On the 50-cell disk the disagreement starts at the boundary.**

   Along the diameter, `|Hz_disk / Hz_vac − Mie|` for FEM-F is `0.0009` at `x = −0.5`, `0.0023` at the upstream boundary, `0` at the center, `0.0024` at the downstream boundary, and `0.0005` at `x = +2`. Meep 35 stays within `0.02` of that ratio everywhere on `|x| < 1.2`. Meep 50 stays close at the center (`|Δratio| = 0.036`) and reaches `0.71–0.75` on the boundary, then `0.14` at `x = +3`. The error is a boundary field plus the wave it radiates. It is not a phase offset accumulated while crossing the plasma: the interior ratio moves back toward FEM when the run is extended, and the boundary field does not.

9. **The FEM interface condition is the one implied by the weak form, applied element-wise.**

   For `div(ρ ∇Hz) + k0² Hz = 0` with `ρ = 1/ε` finite, a pillbox across the material jump gives continuity of the normal flux `ρ ∂Hz/∂n`. The trial space is continuous, so `Hz` is continuous. Together those are continuity of `Hz` and of tangential `E` for this polarization (`E = (i/ω) ρ (∂y Hz, −∂x Hz)`). `ε` itself jumps, and `∂Hz/∂n` jumps with it.

   `assemble_anisotropic` multiplies the element stiffness by the element’s own `ρ_xx, ρ_xy, ρ_yx, ρ_yy`. Those values are constant on each triangle. `rho_on` sets them from the centroid: inside the constrained circle, `ρ = 1/ε`; outside, `ρ = 1`. Circles and the square are constrained edges, so a triangle is not cut by the interface and there is nothing to average. There is no nodal, arithmetic, or harmonic average of `ε`. The natural boundary term of the weak form is what enforces `ρ ∂n Hz` continuity. The flat-slab check from the previous campaign (analytic `|T|²` within `0.09 dB`) is the same operator on a flat jump.

10. **Meep staircases the Drude susceptibility. Subpixel smoothing does not act on it.**

    In Meep 1.30.0 `src/meepgeom.cpp` (the source matching this install):

    - `geom_epsilon::eff_chi1inv_matrix` subpixel-averages the instantaneous permittivity with the Kottke split (harmonic average on the normal component, arithmetic average on the tangential components).
    - `geom_epsilon::sigma_row` returns the susceptibility amplitude by a point test, `material_of_unshifted_point_in_tree_inobject`. Whichever object contains the Yee location contributes its full Drude `σ`. There is no fill fraction.
    - `eff_chi1inv_row_disp` builds the dispersive tensor at a single point. The adjoint comment in that file states that the subpixel routine neglects susceptibilities, and that the dispersive routine does not apply subpixel smoothing.

    The disk is built as `epsilon=1` plus `DrudeSusceptibility`. Air is also `epsilon=1`. The instantaneous-ε average therefore has nothing to change, which is why `eps_averaging` on and off at 25 cells per wavelength produced forward ratios differing by `1.4×10−15`. The curved Drude boundary is a staircase of which `Ex` and `Ey` locations fall inside the cylinder. Those two components sit on offset Yee nodes, so they see different staircases. Ordinary nondispersive subpixel behavior does not describe this interface.

11. **What the disagreement is.**

    - The FEM plasma operator on a body-fitted curve is the continuum solution. The Mie series says so.
    - A nondispersive Meep stand-in for this ε does not exist in a form the time stepper can run.
    - Meep’s curved Drude boundary, at 50 cells per wavelength, supports a long-lived interface ring that the DFT folds into the `fs` bin. The axis-aligned square of the same material does not ring. That points at the curved dispersive staircase.
    - The nested quartz/vacuum/plasma bulb agrees with FEM where the ring has died (25 and 35 cells per wavelength). The nested interface is not the first failure.
    - Accumulation across many bulbs is still unmeasured. It has to be, because one bulb at 35 cells per wavelength agrees and the 91-bulb device at 35 cells per wavelength does not.

12. **Smallest geometry that reproduces the full-device failure: not found.**

    One bare disk and one coated bulb, at the 35-cell resolution where the full device is already far from FEM (P4 about `2.5 dB`, P2 about `13 dB` in the overnight record), sit within `0.03 dB` of FEM on the forward ratio. The 50-cell disk does diverge, and in the same direction as the device (further refinement makes Meep worse), but the forward-power error is `0.1–0.5 dB`, not the device’s many dB. The local boundary field is where that disk failure is large. A single inclusion does not reproduce the port-power failure.

13. **Next test.**

    Three neighboring production bulbs, B=0, same point-source ratios, FEM on a body-fitted mesh and Meep at 25 and 35 cells per wavelength only. Leave 50 cells per wavelength out until the interface ring is out of the DFT. If three bulbs still agree, repeat at seven before touching 91. In parallel, one Harminv at the 50-cell disk boundary, after the source is off, to record the ring frequency and Q. Do not run `B ≠ 0`.

## Field and power notes

Closed-contour outward Poynting on the square `|x|,|y| ≤ 1`, source outside, is a few times `10−6` in the FEM normalization for both vacuum and the disk. With `Im ε ~ 10−3` the absorbed power is too small, relative to contour cancellation, to rank the solvers. The diameter ratios are the comparison.

`eps_averaging=False` was tested on the single disk at 25 cells per wavelength only. It did not move the result. No 91-bulb averaging sweep was run.
