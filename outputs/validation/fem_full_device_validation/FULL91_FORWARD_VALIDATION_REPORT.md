# Full 91-bulb gyrotropic forward model

Label: **FULL91_GYROTROPIC_FORWARD_VALIDATED**

Reference: mesh FULL91-Q for the final continuum estimate (6,052,164 nodes);
FULL91-F remains the archived six-port magnetized matrix. Operating point
`f = 3.85 GHz`, `a = 0.020 m`, `fp = 8.00 GHz`, `gamma = 1.00 MHz`,
production `|B| = 0.05 T`. Time convention `e^{-iωt}`. Tensor convention is
the validated Lorentz tensor: `ω_c = e B_z / m > 0` for `B_z > 0`,
`ε_xx = ε_yy = ε_⊥`, `ε_xy = +iη`, `ε_yx = -iη`.
P1 in the port list below is code port 0: the horn whose outward normal
is nearest `+x`. Ports then increase counterclockwise.

## Closeout addendum (air, interface, and closed balance)

The earlier M→F B = 0 step (`+0.226 dB` through) is retained as a failed
pair. Holding plasma/quartz at the F size and refining cavity air produced
meshes A then Q. A→Q receiving-port changes are at most `0.060 dB` and
`0.35°`, inside the frozen `0.10 dB / 1°` gate. Magnetized F→Q is inside
`0.047 dB / 0.53°`.

Crossed refinement (adversarial check of the VALIDATED label): mesh **X2**
keeps cavity air at the A size (`h_air = 0.022 a`, `h_pml = 0.06 a`) and
refines plasma/quartz beyond F (`h_iface = 0.0042 a`, 11.9 elements across
the quartz wall; F used `0.005 a`). Preferred X at `h_iface = 0.004 a`
reached 6.64e6 nodes and aborted in SuperLU; X2 at 6,233,120 nodes
factored. At B = 0, source P1:

| Comparison | max \|ΔdB\| | max \|Δphase\| |
|------------|------------:|---------------:|
| X2 vs A (iface only) | 0.022 | 0.12° |
| X2 vs Q (crossed) | 0.038 | 0.22° |
| A vs Q (air only) | 0.060 | 0.35° |

All three are inside `0.10 dB / 1°`. Through power: F `0.4237`, A `0.4327`,
X2 `0.4341`, Q `0.4366`. At `+0.05 T`, X2 versus Q is `0.010 dB` /
`0.12°` over all six sources. Observables are therefore independently
converged under cavity-air refinement and under plasma/quartz-interface
refinement, not only within a one-parameter family.

Closed control volume (PML-inner rectangle with the source slot excised):
on mesh F, max relative residual `|F_rect + F_slot + P_abs| /
(|F_rect|+|F_slot|+|P_abs|) = 8.36e-3` against the frozen `0.02` at
B = 0 and both production signs, all six sources.

Port-station flux uncertainty at the production span is `≤ 0.023 dB`.
The Hamid ellipse paper reproduction remains FAIL and is not used as
continuum evidence.

Machine-readable solves:

`outputs/validation/fem_full_device_validation/full91_<grade>_B<sign><tesla>_f<GHz>.json`

Column `j` of every matrix is the response to source port `j`.
The reaction matrix is `S_ij = b_i · u_j` (non-conjugated). That is the
pairing for which Onsager–Casimir is an identity of the complex-symmetric
operator. Horn power is the outward Poynting flux through the monitor
line and is a different functional.

## Geometry

Audit: `FULL91_GEOMETRY_AUDIT.md`. Status PASS.

- 91 bulbs, pointy-top triangular lattice, basis `[[0, 1], [√3/2, 1/2]]`, pitch `1.0`
- Plasma radius `0.230 a` (4.60 mm). `sixport_common.r_plasma = 0.250 a` is unused metadata
- Quartz inner `0.325 a`, outer `0.375 a`, vacuum gap `0.095 a`, quartz `ε = 3.8`
- Domain `30 a × 28 a`, PML `2.0 a`
- PEC is natural Neumann on Hz (`ρ = 0` on the horn-wall triangles)

The first mesh batch had an air seed inside a bulb, so the connected air
region was unconstrained and edges of length `~1–3 a` appeared. Those
solves were discarded. The meshes used here seed the air at the port-0
source center and constrain the source and monitor transects. A mesh is
rejected if its longest edge exceeds eight times the air/PML size target.

## B = 0 bridge

On FULL91-C the gyrotropic assembly at `B = 0` and the scalar assembly
differ by a relative data norm of **0**. Reaction matrices on C, M, and F
are symmetric to `2–3 × 10^{-14}`.

Poynting power for source P1, normalized by the absolute driven-port flux:

| Port | FULL91-C | FULL91-M | FULL91-F | M→F (dB) | M→F phase (deg) |
|------|----------|----------|----------|----------|-----------------|
| P2 adjacent | 0.165228 | 0.153875 | 0.148919 | -0.142 | +1.87 |
| P3 next | 0.086118 | 0.073954 | 0.068863 | -0.310 | +1.24 |
| P4 through | 0.352037 | 0.402207 | 0.423715 | +0.226 | +1.03 |
| P5 next | 0.086113 | 0.073939 | 0.068894 | -0.307 | +1.23 |
| P6 adjacent | 0.165155 | 0.153906 | 0.148906 | -0.143 | +1.89 |

Absolute driven-port flux: C `0.061364`, M `0.062837`, F `0.063355`.
Absorption, source P1: C `1.903e-3`, M `1.880e-3`, F `1.865e-3`.
Sixfold symmetry of the power matrix is at the `4e-4` relative level on F.
The M→F step is outside `0.10 dB / 1 deg`. C→M→F is monotone on every
channel above. That pair is superseded by the air-refined A→Q sequence and
by the independent crossed interface mesh X2 (below). FULL91-S
(`h_iface = 0.0035 a`, 7,154,769 nodes) and the preferred crossed mesh X
(`h_iface = 0.004 a` with A air, 6,641,826 nodes) aborted in SuperLU and
were not repeated.

## Production field, mesh F

`B = +0.05 T`, `f_c a/c = 0.09337289556608207`

`ε_xx = -3.9752973044933535 + 0.0016858808147405032 i`

`ε_xy = +0.0010826788464825367 + 1.8087134773285605 i`

At `B = -0.05 T`, `ε_xy` changes sign and `ε_xx` is unchanged.

### Source P1, complex monitor amplitude

Weights are the trapezoid weights on 41 samples across `0.92` of the feed.
Phase is in degrees.

| Port | +B magnitude | +B phase | +B power | +B power dB vs \|P1\| | −B magnitude | −B phase | −B power | −B power dB vs \|P1\| |
|------|-------------:|---------:|---------:|----------------------:|-------------:|---------:|---------:|----------------------:|
| P1 | 0.68596649 | 166.17246 | -0.01430256 | 0 | 0.68596346 | 166.17218 | -0.01432006 | 0 |
| P2 | 0.04371073 | -35.28977 | 0.00060221 | -13.756 | 0.02315675 | 8.02760 | 0.00023625 | -17.824 |
| P3 | 0.06010985 | -163.27976 | 0.00083112 | -12.357 | 0.05062577 | -158.02363 | 0.00076446 | -12.724 |
| P4 | 0.18564183 | -163.87951 | 0.00783449 | -2.614 | 0.18564300 | -163.88242 | 0.00783355 | -2.620 |
| P5 | 0.05061690 | -158.02299 | 0.00076440 | -12.720 | 0.06010271 | -163.26857 | 0.00083099 | -12.361 |
| P6 | 0.02317267 | 7.97719 | 0.00023619 | -17.821 | 0.04372220 | -35.28443 | 0.00060261 | -13.757 |

Normalized outgoing power at `+B`, source P1:
P2 `0.042105`, P3 `0.058110`, P4 `0.547769`, P5 `0.053445`, P6 `0.016514`.
At `−B` the adjacent pair swaps: P2 `0.016498`, P6 `0.042081`.

M→F at `+0.05 T` for those same normalized powers: `-0.013`, `-0.062`,
`-0.017`, `-0.062`, `-0.042` dB, with monitor phase changes
`-0.006`, `+0.811`, `+0.446`, `+0.146`, `+0.756` deg.
At `−0.05 T` the largest M→F change is `0.083 dB` and `0.87 deg`.

### Poynting matrix at +0.05 T (mesh F)

```
        P1          P2          P3          P4          P5          P6
P1  -0.01430256  0.00023652  0.00076429  0.00783462  0.00083176  0.00060277
P2   0.00060221 -0.01431292  0.00023592  0.00076443  0.00783342  0.00083048
P3   0.00083112  0.00060165 -0.01431680  0.00023616  0.00076407  0.00783295
P4   0.00783449  0.00083130  0.00060281 -0.01430589  0.00023638  0.00076455
P5   0.00076440  0.00783463  0.00083045  0.00060251 -0.01430436  0.00023605
P6   0.00023619  0.00076425  0.00783314  0.00083150  0.00060151 -0.01431493
```

### Poynting matrix at −0.05 T (mesh F)

```
        P1          P2          P3          P4          P5          P6
P1  -0.01432006  0.00060202  0.00083102  0.00783404  0.00076540  0.00023613
P2   0.00023625 -0.01430886  0.00060198  0.00083085  0.00783442  0.00076396
P3   0.00076446  0.00023622 -0.01430205  0.00060301  0.00083061  0.00783360
P4   0.00783355  0.00076472  0.00023600 -0.01431924  0.00060195  0.00083121
P5   0.00083099  0.00783405  0.00076389  0.00023597 -0.01431556  0.00060156
P6   0.00060261  0.00083071  0.00783305  0.00076502  0.00023594 -0.01430509
```

Power Onsager residual `max |P_ij(+B) - P_ji(−B)| = 1.75e-5`
(RMS `5.13e-6`). The circulation contrast on the adjacent ports is
`3.7e-4`, about twenty times that residual.

The reaction matrices are stored as `reaction_reim` in
`full91_F_Bp0.0500_f3.8500.json` and `full91_F_Bm0.0500_f3.8500.json`.
`max |S_ij(+B) - S_ji(−B)| = 1.172918e-13`. Significant-channel phase
residual is `2.6e-11` deg. Same-B nonsymmetry `max |S - S^T| = 0.03167`
(`0.108` of the largest entry). The same identity on FULL91-M holds to
`7.4e-14` or better at every sampled `|B|`.

## Device characterization at the current geometry

These numbers describe the unoptimized production layout. They are not a
circulator grade.

Mesh F, source P1, powers relative to the absolute driven-port flux:

| Quantity | B = 0 | +0.05 T | −0.05 T |
|----------|------:|--------:|--------:|
| Through, P4 / \|P1\| | 0.423715 | 0.547769 | 0.547033 |
| Through, dB | -3.729 | -2.614 | -2.620 |
| CCW adjacent P2 / \|P1\| | 0.148919 | 0.042105 | 0.016498 |
| CW adjacent P6 / \|P1\| | 0.148906 | 0.016514 | 0.042081 |
| P2 versus P6, dB | 0.000 | +4.065 | -4.067 |
| Next ports P3 and P5 | 0.06886 / 0.06889 | 0.05811 / 0.05345 | 0.05338 / 0.05803 |

At `+B`, power prefers the counterclockwise adjacent horn over the
clockwise horn by `4.07 dB`. At `−B` the preference reverses. The through
horn still carries more power than either adjacent horn (`0.55` versus
`0.042`). The other five sources repeat this pattern by rotation:
through-port insertion is `-2.614` to `-2.620 dB` at `+B` and the
CCW/CW ratio is `+4.055` to `+4.074 dB`.

## Passivity

Absorption uses `(ω/2) E† ((ε − ε†)/(2i)) E` on plasma elements only.

| Case | Absorption, all six sources | Horn residual `(−P_jj − Σ_{i≠j} P_ij − P_abs) / (−P_jj)` |
|------|-----------------------------:|----------------------------------------------------------:|
| F, B = 0 | `1.8644e-3` to `1.8648e-3` | 0.111 |
| F, +0.05 T | `6.9962e-4` to `6.9974e-4` | 0.233 |
| F, −0.05 T | `6.9961e-4` to `6.9973e-4` | 0.234 |

Absorption is positive and the same for both B signs. The horn residual
is the same for every source to three digits, which follows the sixfold
layout. It is not closed because each horn monitor is an open line inside
the feeds. The inner-PML rectangle flux is a different surface and is
larger (`0.225` at B = 0, `0.043` at `±0.05 T` for source P1) because the
source lies inside that rectangle. The six-port powers are not forced to
sum to one.

## B sweep (FULL91-M, 3.85 GHz, source P1)

| B (T) | Driven flux | P2 | P3 | P4 through | P5 | P6 | P2 − P6 | Absorption |
|------:|------------:|---:|---:|-----------:|---:|---:|--------:|-----------:|
| 0 | -0.062837 | 0.009669 | 0.004647 | 0.025273 | 0.004646 | 0.009671 | -2.0e-6 | 1.880e-3 |
| +0.0125 | -0.066077 | 0.013083 | 0.002369 | 0.032337 | 0.005623 | 0.004317 | +8.77e-3 | 1.787e-3 |
| +0.0250 | -0.057283 | 0.005973 | 0.003295 | 0.037549 | 0.001172 | 0.002221 | +3.75e-3 | 1.453e-3 |
| +0.0375 | -0.038951 | 0.002650 | 0.004994 | 0.013804 | 0.003252 | 0.001908 | +7.42e-4 | 1.671e-3 |
| +0.0500 | -0.014289 | 0.000603 | 0.000842 | 0.007859 | 0.000775 | 0.000238 | +3.65e-4 | 6.949e-4 |

Negative B reproduces the row with P2 and P6 exchanged, and likewise P3
and P5. `Im(ε_xy)` moves smoothly through `0`, `0.396`, `0.812`, `1.272`,
`1.809`. Through power peaks near `0.025 T` and then falls. That track is
mirrored at `−B`, absorption stays positive, and the reaction Onsager
residual stays below `8e-14`, so the peak is kept as device behavior.

## Frequency sweep (FULL91-M, B = +0.05 T, source P1)

The tensor is rebuilt at each frequency. `f_c` is held at the `0.05 T` value.

| f (GHz) | Re ε_xx | Driven flux | P2 | P4 through | P6 | Absorption |
|--------:|--------:|------------:|---:|-----------:|---:|-----------:|
| 3.6575 | -4.60501 | -0.030394 | 0.003487 | 0.016031 | 0.001832 | 5.940e-4 |
| 3.7538 | -4.27544 | -0.080410 | 0.013156 | 0.041902 | 0.000459 | 1.746e-3 |
| 3.8500 | -3.97530 | -0.014289 | 0.000603 | 0.007859 | 0.000238 | 6.949e-4 |
| 3.9463 | -3.70106 | -0.033249 | 0.001255 | 0.005192 | 0.000844 | 4.053e-3 |
| 4.0425 | -3.44974 | -0.018067 | 0.002842 | 0.002880 | 0.001480 | 1.118e-3 |

`Re ε_xx` is monotone. Absorption is positive at every sample. The
operating point sits between a high-transmission sample at `3.754 GHz`
and a much smaller through power at `3.850 GHz`. Neighboring samples are
not a flat baseline.

## Comparison with existing Meep data

Meep is not the continuum reference. No saved full-device magnetized Meep
result was found, and no new Meep run was launched.

Existing B = 0 case `outputs/validation/fem_meep_validation/overnight_b0/meep_plasma_ppc50.json`:

- resolution 100 (`dx = 0.01 a`, 0.2 mm), 50 points per cm, runtime 20, wall time 4338 s
- formulation `num_mode_guide_normal`, `eps_averaging` on
- `P_inc = 2910.05055`
- normalized flux `[−0.672167, 0.002130, 0.002142, 0.031266, 0.002142, 0.002130]`

Dividing the outgoing Meep fluxes by the absolute driven-port flux gives
through `0.0465` and adjacent `0.00317`. FULL91-F at B = 0 gives through
`0.424` and adjacent `0.149`. The two side classes in Meep differ by
`0.024 dB`; in the FEM they differ by `3.35 dB` and that split is stable
from C to F and matches the two crystallographic port classes of the hex
layout.

The earlier body-fitted B = 0 snapshot `fem_plasma.json` level FEM-VH
(916,818 nodes) has outgoing/driven ratios through `0.341`, adjacent
`0.167`, next `0.089`. Those sit next to FULL91-C (`0.352`, `0.165`,
`0.086`). FULL91-Q through/driven is `0.437`. The Meep side remains the
finite-grid dispersive-interface disagreement already separated in the
frozen B = 0 campaign. This FEM was not adjusted toward the Meep numbers.

## Field plot

`full91_hz_port0.png` is `log10 |Hz|` on mesh C for source P1 at
`B = 0`, `+0.05 T`, and `−0.05 T`. It is a visualization of the coarse
mesh, not the production field.

## Open items

- FULL91-S cannot be factored on this SuperLU build. Continuity evidence
  uses A→Q instead.
- The 3.85 GHz, 0.05 T operating point remains steep; dense samples are
  recorded and are physical structure, not a solver failure.
- Hamid & Cooray ellipse echo-width reproduction is FAIL under the frozen
  criteria and is not used as continuum evidence.
- No full-device magnetized Meep comparison exists in the saved results.
