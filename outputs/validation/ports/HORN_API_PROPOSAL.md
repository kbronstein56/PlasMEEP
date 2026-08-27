# Proposed reusable PlasMEEP horn API (deferred implementation)

**Do not implement until resolution/polarization gates are settled.** Validation should keep using `sixport_common` / `port_formulations` until then.

## Goal

Replace duplicated horn construction in validation scripts with one PlasMEEP-facing builder that supports:

- arbitrary in-plane orientation (outward normal)
- PEC or finite-ε walls
- roles: `source` | `probe` | `passive`
- explicit field component: `Ez` | `Hz` | `Ex`/`Ey` (avoid bare TE/TM names in the public API)
- optional amplitude profile along the aperture (uniform, cos TE1, custom)
- monitor: axis-aligned flux and/or guide-normal / DFT hooks

## Sketch

```python
horn = plasmeep.ports.Horn(
    open_center=...,
    outward_dir=...,
    throat_width_m=0.048,
    aperture_width_m=0.104,
    flare_depth_m=0.089,
    wall_thickness_m=0.004,
    feed_length_m=0.060,
    wall="pec",                 # or epsilon=1e6
)
sources = horn.make_sources(
    frequency=fs_a,
    fwidth=df,
    component="Hz",             # explicit
    profile="te1_cos",
)
monitors = horn.make_flux_monitors(component_frame="axis_aligned")
```

## Where code should live

| Today | Target |
|---|---|
| `sixport_common.get_full_horn` / port dirs | `plasmeep/ports/horn.py` (or under `inverse_pmm`) |
| `port_formulations.make_*_sources` | thin wrappers over horn API |
| validation scripts | call library; keep experiment loops only |

Production notebook / `PMMCirculatorInverse.py` stay untouched until a validated library API exists.
