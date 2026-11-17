# Circulator physical-limits context (literature)

**Purpose:** qualitative context for the PlasMEEP 91-element feasibility study.  
**Stage-6:** NOT STARTED.  
**Does not change** the production FEM / constitutive model.

## Production target vs related experiments

| Item | This campaign | Closest logged experiment (LIT-003) |
|------|---------------|-------------------------------------|
| Band | 4–6 GHz | 4–6 GHz (WR187 Y-junction) |
| Ports | **6** | **3** |
| Geometry | 91-element PMM array + horns | Single plasma column in Y-junction |
| B (production) | **0.05 T** | ~**0.104 T** center (NdFeB / coil setup) |
| Reported 10 dB isolation BW | *not achieved on continuum audit* | ~**1.6 GHz** (40% BW) at peak iso ~18.6 dB |
| Loss | FEM includes plasma collisions | Experiment: large absorbed fraction (~50%+) |

Source: Li et al., *Phys. Rev. Applied* **18**, 064066 (2022); arXiv:2206.12507 — already in `docs/circulator_literature/LITERATURE_LOG.md`.

## What the literature suggests (qualitative)

1. **Broadband microwave plasma circulation is experimentally real** near 4–6 GHz — but demonstrated on a **3-port Y-junction**, not a 6-port 91-element array.
2. Isolation bandwidth of order **~1 GHz** (not a flat 2 GHz continuum without notches) appears in that Y-junction family when B ~ 0.1 T.
3. Raising **B** increased peak isolation with little frequency shift in Li et al. simulations — consistent with gyrotropy / Voigt-parameter strength as a bottleneck.
4. Plasma **absorption** and mismatch are first-class IL mechanisms; deep FEM notches with tiny T_des are compatible with resonant absorption / transmission zeros.

## Bounds that do **not** directly apply

- Classical ferrite circulator Bode-type or matching bounds generally assume different materials and 3-port junction topology.
- Passivity/causality bounds on nonreciprocal response exist in abstract multiport theory, but translating them into a **quantitative A2 floor for this lossy gyrotropic FEM device** would require assumptions we have **not** validated. We therefore treat them as **context only**, not as a proof of impossibility.

## Implication for PlasMEEP

Asking the **current** 91-element / **0.05 T** / **C6-16** package for a **uniformly good 2 GHz** continuum (4–6 GHz) is a **harder** ask than the successful ~1–1.6 GHz isolation windows reported for related 3-port plasma circulators at higher B. That does **not** prove the geometry cannot do it; it **does** motivate:

- a **bandwidth ladder** (can we achieve ~0.2–1.0 GHz convincingly?);
- a **B sensitivity** probe only at documented values (0.05 T production; ~0.1 T as literature-motivated extrapolation, clearly labeled);
- a **DOF** probe (C6 vs C3) before blaming geometry alone.

```
STAGE-6: NOT STARTED / AWAITING AUTHORIZATION
```
