# Meep PEC / Cartesian grid handling (horns_only study)

## What the harness uses today

- Horn walls are `mp.Prism` extrusions with `material=PEC` (via `PlasMEEP.Add_Prism(..., PEC=True)`).
- Meep rasterizes **all** geometry—including PEC—onto the **Cartesian Yee grid**.
- Angled polygon edges become **staircased** voxels; there is **no dielectric subpixel averaging** for perfect conductors.

## Implications for rotated horns

| Effect | Axis-aligned port (P1) | Rotated port (P2, ~60°) |
|---|---|---|
| Wall staircasing | minimal (edges ∥ grid) | strong (oblique edges) |
| Effective aperture | stable vs resolution | shifts with sub-cell placement |
| Flux monitor (axis-aligned `FluxRegion`) | samples along grid line | samples along grid line, not along physical mouth |

This matches the observed pattern:

- P2↔P3 (both diagonal): reciprocity ≈ 0 dB
- P1↔P2 (axis↔diagonal): ~0.6 dB at res32, **non-monotonic** with resolution
- Error persists in `horns_only` (no PMM)

## Legitimate alternatives to test

1. **Fractional-cell geometry registration** — translate all horns together; if reciprocity varies strongly, staircasing/registration is confirmed.
2. **Finite-ε metal** (`high_eps` wall mode) — ε≈10⁶ smooths interfaces differently from true PEC; compare on horns_only.
3. **Guide-normal flux integration** — already tested; does not fix axis error (~0.60 dB).
4. **Subpixel smoothing** — Meep's `eps_averaging` applies to **dielectric** interfaces, not PEC prisms in this build.
5. **Finer local resolution** — only after a representation shows coherent convergence; not brute-force res160+.

## Working conclusion (pre-sweep)

The reciprocity artifact is **consistent with rotated PEC horn staircasing + axis-aligned flux sampling on a Cartesian grid**, not plasma physics.
