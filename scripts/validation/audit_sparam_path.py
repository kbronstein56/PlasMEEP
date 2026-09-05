#!/usr/bin/env python3
"""
Static audit of S21/S12 calculation path + Yee-grid port registration.

No FDTD. Documents references, cache keys, and grid DOFs for P1/P2/P3.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from physical_units import resolve_simulation_resolution, physical_resolution_report
from sixport_common import (
    a,
    clear_width,
    effective_port_dir,
    fs_a,
    horn_for_port,
    monitor_center_for_port,
    set_geometry_context,
    wall_thickness,
)


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def yee_hz_index(xy: np.ndarray, res: int) -> Dict[str, Any]:
    """Hz is sampled at integer cell centers in Meep's Yee grid."""
    xy = np.asarray(xy, dtype=float)
    # Meep Hz lives at (i+0.5, j+0.5)/res in some conventions; for point
    # sources Meep snaps to nearest Hz location. Report both requested and
    # nearest half-integer cell index.
    ix = xy[0] * res
    iy = xy[1] * res
    snap_x = (np.floor(ix) + 0.5) / res
    snap_y = (np.floor(iy) + 0.5) / res
    return {
        "requested_xy_a": [float(xy[0]), float(xy[1])],
        "requested_xy_mm": [float(xy[0] * a * 1e3), float(xy[1] * a * 1e3)],
        "continuous_cell": [float(ix), float(iy)],
        "frac_offset_from_hz_center": [
            float(ix - (np.floor(ix) + 0.5)),
            float(iy - (np.floor(iy) + 0.5)),
        ],
        "nearest_hz_xy_a": [float(snap_x), float(snap_y)],
        "nearest_hz_cell_index": [float(np.floor(ix)), float(np.floor(iy))],
    }


def port_audit(port_index: int, res: int) -> Dict[str, Any]:
    u = np.asarray(effective_port_dir(port_index), dtype=float)
    u = u / np.linalg.norm(u)
    tangent = np.array([-u[1], u[0]])
    horn = horn_for_port(port_index, res)
    src = np.asarray(horn["source_center"], dtype=float)
    mon = np.asarray(monitor_center_for_port(port_index, res), dtype=float)
    # guide-normal sample count matches formulation default
    n_points = 31
    span = 0.96 * clear_width
    offsets = np.linspace(-span / 2, span / 2, n_points)
    samples = [mon + s * tangent for s in offsets]
    return {
        "port": f"P{port_index + 1}",
        "port_index": port_index,
        "outward_normal": u.tolist(),
        "tangent": tangent.tolist(),
        "angle_deg": float(np.degrees(np.arctan2(u[1], u[0]))),
        "source": yee_hz_index(src, res),
        "monitor": yee_hz_index(mon, res),
        "source_to_monitor_distance_a": float(np.linalg.norm(mon - src)),
        "source_to_monitor_distance_mm": float(np.linalg.norm(mon - src) * a * 1e3),
        "guide_normal_n_points": n_points,
        "guide_normal_span_a": float(span),
        "guide_normal_span_mm": float(span * a * 1e3),
        "first_sample": yee_hz_index(samples[0], res),
        "mid_sample": yee_hz_index(samples[n_points // 2], res),
        "last_sample": yee_hz_index(samples[-1], res),
        "clear_width_mm": float(clear_width * a * 1e3),
        "wall_thickness_mm": float(wall_thickness * a * 1e3),
    }


def calculation_path_audit() -> Dict[str, Any]:
    return {
        "quantity_name": "power_matrix[receive, source] ≈ |S_rs|^2 (power ratio, not complex S)",
        "S21_like": "power_matrix[P2, P1] = flux_out_P2(drive_P1) / incident_P1",
        "S12_like": "power_matrix[P1, P2] = flux_out_P1(drive_P2) / incident_P2",
        "dB": "10*log10(power_matrix) when positive",
        "reference_simulation": {
            "function": "normalize_port() in sixport_common.py",
            "geometry": "isolated straight PEC/high-eps parallel-plate feed at that port's source_center orientation",
            "source": "same formulation source as device (e.g. numerical mode Hz line)",
            "saved": [
                "incident_power = |signed outward flux| at monitor plane",
                "incident_flux_data = Meep FluxData DFT for load_minus_flux_data (flux formulations only)",
            ],
            "P1_reference": "built with port_index=0 orientation only",
            "P2_reference": "built with port_index=1 orientation only",
            "shared_between_ports": False,
        },
        "device_simulation": {
            "function": "simulate_circulator() in sixport_common.py",
            "per_source": "one FDTD run; monitors on all requested ports",
            "load_minus_flux_data": {
                "called": "only on SOURCE port flux monitor",
                "subtracts": "stored reference incident DFT flux for THAT same source port",
                "affects_reflection_diagonal": True,
                "affects_transmission_S21_S12": False,
                "code": "sim_i.load_minus_flux_data(monitors_i[src_local], incident_flux_data_by_port[source_port])",
            },
            "transmitted_power": "signed guide-normal flux at receive port (NO subtraction)",
            "reflected_power": "source-port flux AFTER load_minus_flux_data / incident",
            "normalization": "column_i = flux_i / incident_power_by_port[source_port]",
        },
        "cache_keys": {
            "norm_pickle": [
                "res",
                "run_time",
                "formulation",
                "grid_offset_cells",
                "monitor_offset_cells",
                "horn_walls",
                "coord_rotation_deg",
                "port subset in filename",
            ],
            "mode_json": [
                "res directory (res{N}/numerical_mode_P{k}.json)",
                "frequency_a checked at load",
                "per-port file (P1..P6 not shared)",
            ],
            "missing_from_norm_key_explicitly": [
                "points_per_cm (implied by res + a)",
                "a_m (implicit via module constant; changing a invalidates physics of old caches)",
            ],
        },
        "orientation_dependence": [
            "source amplitudes follow per-port numerical mode + local tangent",
            "flux regions use outward normal (n_x, n_y) weights along true tangent",
            "P1 is axis-aligned (horizontal); P2/P3 are ±60° — FluxRegion axis-alignment requires guide-normal point decomposition",
        ],
        "critical_finding": (
            "load_minus_flux_data cannot explain S21!=S12 transmission mismatch: "
            "it only modifies the source-port (reflection) monitor. "
            "S21 and S12 are raw receive fluxes divided by per-port incident powers."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points-per-cm", type=float, default=50.0)
    parser.add_argument("--ports", type=str, default="0,1,2")
    parser.add_argument(
        "--json-out",
        type=str,
        default="/home/daq_user/PlasMEEP/outputs/validation/sparam_audit/phase1_sparam_path_audit.json",
    )
    args = parser.parse_args()

    class _A:
        points_per_cm = args.points_per_cm
        res = None

    sim_res, ppc, report = resolve_simulation_resolution(_A(), a_m=a)
    set_geometry_context(res=sim_res)
    ports = [int(x) for x in args.ports.split(",")]

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "physical_resolution": report,
        "calculation_path": calculation_path_audit(),
        "ports": [port_audit(p, sim_res) for p in ports],
        "reference_vs_device_monitor_identity": {
            "same_functions": [
                "monitor_center_for_port(port, res)",
                "make_flux_region_for_formulation(...)",
                "horn_for_port(port, res)['source_center'] for launch",
            ],
            "reference_differs_by": "geometry content (straight feed vs full device), NOT monitor coordinate formulas",
            "yee_dofs_match_if": "same res, same grid_offset, same monitor_offset, same a",
        },
    }

    out = Path(args.json_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(_json_safe(payload), indent=2) + "\n", encoding="utf-8")

    md = out.with_suffix(".md")
    lines = [
        "# Phase 1 — S21/S12 calculation audit",
        "",
        f"Generated: {payload['timestamp_utc']}",
        f"Resolution: {ppc:g} points/cm, Meep res={sim_res}, a={a} m, dx={report['dx_mm']} mm",
        "",
        "## How S21 / S12 are formed",
        "",
        "| Quantity | Formula | Uses `load_minus_flux_data`? |",
        "|---|---|---|",
        r"| Incident P1 | `|flux|` on straight-feed reference at P1 | no |",
        r"| Incident P2 | `|flux|` on straight-feed reference at P2 | no |",
        r"| S21-like | `flux(P2|drive P1) / incident_P1` | **no** |",
        r"| S12-like | `flux(P1|drive P2) / incident_P2` | **no** |",
        r"| S11-like | `flux(P1|drive P1 after subtract) / incident_P1` | **yes** |",
        "",
        f"**Critical finding:** {payload['calculation_path']['critical_finding']}",
        "",
        "## Reference roles",
        "",
        "- P1 device run uses **P1** reference only (`incident_flux_data_by_port[0]`).",
        "- P2 device run uses **P2** reference only (`incident_flux_data_by_port[1]`).",
        "- References are **not** shared across incompatible ports.",
        "- Cache keys include res, run_time, formulation, grid/monitor offsets, horn walls, rotation.",
        "",
        "## Port Yee registration (50 points/cm)",
        "",
        "| Port | angle° | src frac→Hz | mon frac→Hz | src–mon mm |",
        "|---|---:|---:|---:|---:|",
    ]
    for p in payload["ports"]:
        sf = p["source"]["frac_offset_from_hz_center"]
        mf = p["monitor"]["frac_offset_from_hz_center"]
        lines.append(
            f"| {p['port']} | {p['angle_deg']:.1f} | "
            f"({sf[0]:+.3f},{sf[1]:+.3f}) | ({mf[0]:+.3f},{mf[1]:+.3f}) | "
            f"{p['source_to_monitor_distance_mm']:.2f} |"
        )
    lines += [
        "",
        "Full JSON: `" + str(out) + "`",
        "",
    ]
    md.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out}")
    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
