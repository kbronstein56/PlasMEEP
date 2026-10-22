# Experimental mode resolution feasibility

**Scope:** feasibility study, not hardware control.  
**Claim level:** TE0/TE1 decomposition in the **numerical port basis is not
currently available** from ordinary VNA/horn coupling alone.

## What a VNA/horn measurement gives today

Standard six-port S-parameters from horn coupling report **port-to-port power
and phase** in the receive-horn basis. That is **not** automatically a
projection onto simulated flat-port TE0/TE1 eigenmodes.

## Options

### A. Near-field scan across the output aperture (most direct)

- Measure complex amplitude and phase vs transverse position (probe + VNA).
- Project measured \(H_z(y)\) (or \(E\)) onto simulated TE0/TE1 basis with the
  same transverse weight/measure used in FEM.
- **Needs:** XY scanner or rail, calibrated probe, phase-stable cables,
  positioning accuracy ≪ guide width, careful de-embedding of probe coupling.
- **Feasible** in a microwave lab with modest hardware; not free.

### B. Two+ spatial probe locations (even/odd)

- Place probes at even/odd sampling points to separate TE0 (even) vs TE1 (odd).
- Cheaper than a full scan; more sensitive to placement error and higher-order
  content.
- Useful as a **qualitative** TE1 detector; quantitative power needs calibration.

### C. Mode-selective waveguide adapter / mode filter

- Overmoded guide + mode filter to reject TE1 before the horn.
- Realistic only if we redesign the feed; not a drop-in for current horns.
- Would change the **experimental objective** (force TE0), which is a different
  experiment than “what mode is present.”

### D. Measurement-aware objective (scientifically valuable)

- Simulate or measure the **actual receive-horn mode**.
- Optimize FEM overlap with that measured/estimated receive pattern rather than
  abstract flat-port TE0.
- **Needs:** calibrated horn pattern or a one-time near-field characterization of
  the receive chain; then freeze that pattern into the FoM.
- Aligns simulation with what the VNA actually sees.

## Recommendation

| Goal | Approach |
|------|----------|
| Answer “which mode carries power?” | Near-field scan (A) or dual-probe (B) |
| Match experiment without mode resolution | Measurement-aware objective (D) |
| Force TE0 experimentally | Mode filter (C) — changes the device target |

**Do not claim** experimental TE0/TE1 resolution is available until A or B is
demonstrated on the bench with a published calibration residual.
