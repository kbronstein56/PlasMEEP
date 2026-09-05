# Phase 1 — S21/S12 calculation audit

Generated: 2026-09-05T05:18:04.606128+00:00
Resolution: 50 points/cm, Meep res=100, a=0.02 m, dx=0.2 mm

## How S21 / S12 are formed

| Quantity | Formula | Uses `load_minus_flux_data`? |
|---|---|---|
| Incident P1 | `\|flux\|` on straight-feed reference at P1 | no |
| Incident P2 | `\|flux\|` on straight-feed reference at P2 | no |
| S21-like | `flux(P2\|drive P1) / incident_P1` | **no** |
| S12-like | `flux(P1\|drive P2) / incident_P2` | **no** |
| S11-like | `flux(P1\|drive P1 after subtract) / incident_P1` | **yes** |

**Critical finding:** load_minus_flux_data cannot explain S21!=S12 transmission mismatch: it only modifies the source-port (reflection) monitor. S21 and S12 are raw receive fluxes divided by per-port incident powers.

## Reference roles

- P1 device run uses **P1** reference only (`incident_flux_data_by_port[0]`).
- P2 device run uses **P2** reference only (`incident_flux_data_by_port[1]`).
- References are **not** shared across incompatible ports.
- Cache keys include res, run_time, formulation, grid/monitor offsets, horn walls, rotation.

## Port Yee registration (50 points/cm)

| Port | angle° | src frac→Hz | mon frac→Hz | src–mon mm |
|---|---:|---:|---:|---:|
| P1 | -0.0 | (+0.013,-0.500) | (+0.013,-0.500) | 24.00 |
| P2 | 60.0 | (+0.256,+0.223) | (+0.256,+0.300) | 24.00 |
| P3 | 120.0 | (-0.256,+0.223) | (-0.256,+0.300) | 24.00 |

Full JSON: `/home/daq_user/PlasMEEP/outputs/validation/sparam_audit/phase1_sparam_path_audit.json`
