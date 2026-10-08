#!/usr/bin/env python3
"""
Phase 2 ContinuousSource FDTD screen (benchmark-only).

solve_cw cannot handle Drude; this tests single-frequency ContinuousSource
time-domain settling as an alternative.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import meep as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
VAL = ROOT / "scripts" / "validation"
for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import sixport_common as sc  # noqa: E402
from physical_units import meep_resolution_from_points_per_cm, physical_resolution_report  # noqa: E402
from port_formulations import (  # noqa: E402
    add_flux_monitor,
    extract_flux_powers,
    get_formulation,
    make_flux_region_for_formulation,
    port_measure_center,
)
from plasmeep.ports.mode_registry import get_numerical_mode  # noqa: E402
from plasmeep.ports.numerical_launch import port_tangent  # noqa: E402
from sixport_common import (  # noqa: E402
    build_circulator_device,
    default_uniform_rho,
    set_geometry_context,
    normalize_port,
)


FORMULATION = "num_mode_guide_normal"
RECEIVE_PORTS = list(range(6))
SOURCE_PORT = 0


def _mpi_rank() -> int:
    try:
        from mpi4py import MPI
        return int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        return int(os.environ.get("PMI_RANK", "0"))


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if obj is None or isinstance(obj, (str, bool, int, float)):
        return obj
    return str(obj)


def make_cw_numerical_sources(port_index: int, frequency: float) -> List[mp.Source]:
    """Same spatial mode as num_mode_guide_normal, ContinuousSource envelope."""
    res = sc.current_res()
    u = sc.effective_port_dir(port_index)
    tangent = port_tangent(u)
    mode = get_numerical_mode(
        port_index,
        res=res,
        frequency_a=frequency,
        horn_walls=sc.get_horn_walls(),
        grid_offset_cells=sc.get_grid_offset_cells(),
        coord_rotation_deg=sc.get_coord_rotation_deg(),
        validate_alignment=(u, tangent),
        log_audit=True,
    )
    center = np.asarray(sc.horn_for_port(port_index, res)["source_center"], dtype=float)
    weights = mode.source_amplitudes(target_power=1.0)
    offsets = np.asarray(mode.offsets_a, dtype=float)
    sources: List[mp.Source] = []
    for s, amp in zip(offsets, weights):
        xy = center + s * tangent
        sources.append(
            mp.Source(
                src=mp.ContinuousSource(frequency=frequency, is_integrated=True),
                component=mp.Hz,
                center=mp.Vector3(xy[0], xy[1], 0),
                amplitude=complex(amp),
            )
        )
    return sources


def run_device_cw(*, res: int, until: float, incident_power: float, incident_flux) -> Dict[str, Any]:
    form = get_formulation(FORMULATION)
    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, P_device, _wp = build_circulator_device(rho, B, res=res, device_mode="full", wall_pec=True)
    # Monkeypatch sources for this object only
    P_device.sources = make_cw_numerical_sources(SOURCE_PORT, sc.fs_a)
    sim = P_device.Get_Sim()
    monitors, signs = [], []
    for out_p in RECEIVE_PORTS:
        measure_xy = port_measure_center(form.name, out_p)
        regions, sign = make_flux_region_for_formulation(
            form.name, measure_xy, sc.effective_port_dir(out_p), out_p
        )
        monitors.append(add_flux_monitor(sim, regions))
        signs.append(sign)
    # For CW, flux subtraction of Gaussian-norm data is wrong. Skip subtraction;
    # report raw normalized by a CW reference power measured separately.
    # If incident_flux is None, no subtraction.
    if incident_flux is not None:
        sim.load_minus_flux_data(monitors[0], incident_flux)

    print(f"Running ContinuousSource device until={until} ...")
    t0 = time.perf_counter()
    sim.run(until=until)
    t_dev = time.perf_counter() - t0
    flux = extract_flux_powers(sim, monitors, signs)
    col = flux / incident_power
    by_port = {int(p): float(col[k]) for k, p in enumerate(RECEIVE_PORTS)}
    return {"device_wall_s": t_dev, "by_port": by_port, "raw_flux": {int(p): float(flux[k]) for k,p in enumerate(RECEIVE_PORTS)}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cycles", type=float, required=True, help="Run duration in periods of fs")
    ap.add_argument("--points-per-cm", type=float, default=25.0)
    ap.add_argument("--json-out", type=str, required=True)
    args = ap.parse_args()

    res = meep_resolution_from_points_per_cm(args.points_per_cm, a_m=sc.a)
    set_geometry_context(res=res, horn_walls="prism", grid_offset_cells=(0, 0), coord_rotation_deg=0.0)
    period = 1.0 / float(sc.fs_a)
    until = float(args.cycles) * period

    print(sc.geometry_summary())
    print(f"CW-FDTD cycles={args.cycles} period={period:.6f} until={until:.6f} res={res}")

    form = get_formulation(FORMULATION)
    from sixport_common import (
        PMMI,
        rotated_wall,
        clear_width,
        wall_thickness,
        nx_ports,
        ny_ports,
        dpml_ports,
        effective_port_dir,
        horn_for_port,
        a,
    )

    # CW reference: same geometry as normalize_port, ContinuousSource (no Formulation monkeypatch).
    t0 = time.perf_counter()
    u = np.asarray(effective_port_dir(SOURCE_PORT), dtype=float)
    u = u / np.linalg.norm(u)
    n = np.array([-u[1], u[0]])
    pmm_ref = PMMI(a=a, res=res, nx=nx_ports, ny=ny_ports, dpml=dpml_ports, B=np.zeros(3))
    P_ref = pmm_ref.Build_Sim()
    reference_center = np.asarray(horn_for_port(SOURCE_PORT, res)["source_center"], dtype=float)
    guide_length = 2.0 * np.hypot(nx_ports, ny_ports)
    wall_center_offset = clear_width / 2 + wall_thickness / 2
    for sgn in (+1, -1):
        wall = rotated_wall(
            centerline=reference_center,
            axis=u,
            normal=n,
            normal_offset=sgn * wall_center_offset,
            length=guide_length,
            thickness=wall_thickness,
        )
        P_ref.Add_Prism(vertices=wall, axis=np.array([0, 0, 1]), PEC=True)
    P_ref.sources = make_cw_numerical_sources(SOURCE_PORT, sc.fs_a)
    ref_sim = P_ref.Get_Sim()
    measure_xy = port_measure_center(form.name, SOURCE_PORT)
    regions, sign = make_flux_region_for_formulation(
        form.name, measure_xy, effective_port_dir(SOURCE_PORT), SOURCE_PORT
    )
    mon = add_flux_monitor(ref_sim, regions)
    print("Running CW reference...")
    ref_sim.run(until=until)
    raw = sign * mp.get_fluxes(mon)[0]
    incident_power = abs(raw)
    incident_flux = ref_sim.get_flux_data(mon)
    print("CW incident_power", incident_power)
    t_norm = time.perf_counter() - t0

    device = run_device_cw(res=res, until=until, incident_power=incident_power, incident_flux=incident_flux)
    by = device["by_port"]
    payload = {
        "label": f"cw_fdtd_cyc{args.cycles:g}_ppc{args.points_per_cm:g}",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "phase": 2,
        "method": "ContinuousSource_FDTD",
        "cycles": float(args.cycles),
        "period": period,
        "until": until,
        "timings_s": {"normalization": t_norm, "device": device["device_wall_s"], "total": t_norm + device["device_wall_s"]},
        "incident_power": incident_power,
        "normalized_power_by_port": by,
        "P11": by.get(0),
        "P1_to_P2": by.get(1),
        "settings": {"points_per_cm": args.points_per_cm, "res": res, "formulation": FORMULATION},
        "physical_resolution": physical_resolution_report(points_per_cm=args.points_per_cm),
        "note": "solve_cw unsupported for Drude; ContinuousSource FDTD with matched CW reference + flux subtraction",
    }
    if _mpi_rank() == 0:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(_json_safe(payload), indent=2) + "\n")
        print("Wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
