# Specialist panel review — B=0 six-port reciprocity residual

Independent findings before new Meep runs. Evidence from notebook canonical
matrix (res=64, run_time=80, uniform fp=8 GHz, B=0).

Canonical residual: mean |Pij−Pji| ≈ 0.75 dB, max ≈ 1.64 dB.

## Physics Validator

**Assessment:** At B=0 the gyrotropic term is off; Faraday benchmark already
PASS for PlasMEEP material path. Residual is unlikely to be magnetized-plasma
physics. Passive plasma+quartz+PEC horns should be reciprocal up to numerical
error for consistent port conventions.

### Confirmed
- Faraday material path validated independently (PASS).
- B=0 circulation preference ≈ 0 in notebook (consistent with reciprocity of
  handedness, not of matrix elements).

### Physics blockers
1. None unique to plasma physics for the B=0 residual.
2. Residual must be attributed to numerics / ports / measurement definition.

## Numerics Validator

**Assessment:** Pair structure strongly indicates Cartesian-grid rotated-feed
discretization error.

### Evidence
From canonical dB matrix, unique pairs:

| |diff| (dB) | pairs | geometry class |
|---:|---|---|
| 1.64 | P1↔P2, P1↔P6, P3↔P4, P4↔P5 | axis-aligned ↔ ±60° diagonal |
| 1.19 | P1↔P3, P1↔P5, P2↔P4, P4↔P6 | axis ↔ other diagonals |
| 0.00 | P2↔P3, P2↔P5, P2↔P6, P3↔P5, P3↔P6, P5↔P6, P1↔P4 | diagonal↔diagonal or opposite-axis |

Ports: P1/P4 along ±x; P2/P3/P5/P6 at 60° increments — poorly aligned with the
Cartesian Yee grid relative to P1/P4.

### Highest-value tests
1. Resolution sweep on P1↔P2 (worst) — error should fall if discretization.
2. Contrast P1↔P2 vs P2↔P3 at same cost — orientation dependence check.
3. Runtime sweep — if flat, not field-decay truncation.

## Ports & Measurements Validator

**Assessment:** Measurement definition (signed flux / incident subtraction /
per-port normalization) can contribute, but the *pattern* of which pairs fail
points first at rotated source/monitor sampling, not a global normalization
bug (a global scale error would not zero some pairs while failing others).

### Confirmed
- Standardized `make_port_source` used for all ports in canonical path.
- Incident subtraction applied at source port.

### Risks
1. Flux monitor normals/signs on rotated ports.
2. Source point sampling on odd/even pixel alignments differing by orientation.
3. Power matrix is normalized flux ratios, not complex S — reciprocity of
   power ratios is approximate when reflections/modes mix; still expect much
   smaller error than 1.6 dB if ports are equivalent.

## Code Integration Validator

**Assessment:** Notebook `simulate_circulator` rebuilds materials per B; B=0
path uses ordinary Drude. Production `PMMCirculatorInverse` still lacks
ported S-parameter methods — investigation must use notebook-faithful harness
(`scripts/validation/sixport_common.py`), not unfinished production API.

### Confirmed
- Extracted harness asserts 91 elements; geometry matches notebook ports.

### Risks
1. Cache/reuse of incident flux data across incompatible resolutions.
2. Port-subset matrix indexing bugs in diagnostics (guard with tests).

## Joint ranking (pre-run)

1. **Most likely:** Cartesian discretization / rotated-port inequivalence
2. **Plausible contributor:** monitor/source sampling details on diagonals
3. **Less likely alone:** insufficient run_time (to be checked)
4. **Unlikely:** PlasMEEP B=0 material bug (Faraday/B=0 Drude path OK)
5. **Out of scope for residual pattern:** true nonreciprocal physics at B=0
