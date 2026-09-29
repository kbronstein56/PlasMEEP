#!/usr/bin/env python3
"""
Trusted Meep P1→P1..P6 reference for FEM validation.

Uses production num_mode_guide_normal + source_df=0.10*fs_a + production norm cache.
Does NOT modify production defaults permanently (process-local only).
"""
from __future__ import annotations

import json
import os
import pickle
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
VAL = ROOT / "scripts" / "validation"
sys.path[:0] = [str(VAL), str(ROOT / "scripts")]

import sixport_common as sc  # noqa: E402
from physical_units import meep_resolution_from_points_per_cm, physical_resolution_report  # noqa: E402
from port_formulations import (  # noqa: E402
    add_flux_monitor,
    extract_flux_powers,
    get_formulation,
    make_flux_region_for_formulation,
    port_measure_center,
)
from sixport_common import (  # noqa: E402
    build_circulator_device,
    default_uniform_rho,
    ensure_normalizations,
    set_geometry_context,
)

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "phase2"
FORMULATION = "num_mode_guide_normal"
SOURCE_PORT = 0
RECEIVE_PORTS = list(range(6))


def _rank() -> int:
    try:
        from mpi4py import MPI
        return int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        return 0


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    points_per_cm = 25.0
    run_time = 20.0
    # Ensure production 0.10 definition
    assert abs(sc.source_df - 0.10 * sc.fs_a) < 1e-12, (sc.source_df, 0.10 * sc.fs_a)

    res = int(meep_resolution_from_points_per_cm(points_per_cm, a_m=sc.a))
    set_geometry_context(
        res=res,
        horn_walls="prism",
        grid_offset_cells=(0.0, 0.0),
        monitor_offset_cells=(0.0, 0.0),
        coord_rotation_deg=0.0,
    )

    # Production norm cache for port 0 only (compatible with P1 column)
    prod_norm = (
        VAL / ".cache" /
        f"norm_{FORMULATION}_full_prism_rot0_pec_g0_0_m0_0_"
        f"ppc{points_per_cm:g}_res{res}_rt{run_time:g}_p0.pkl"
    )
    print(f"Using norm cache path: {prod_norm} exists={prod_norm.is_file()}")
    t0 = time.perf_counter()
    incident_cache = ensure_normalizations(
        res=res,
        run_time=run_time,
        force=False,
        ports=[SOURCE_PORT],
        cache_path=str(prod_norm),
        verbose=True,
        formulation=FORMULATION,
    )
    t_norm = time.perf_counter() - t0

    form = get_formulation(FORMULATION)
    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, P_device, wp = build_circulator_device(rho, B, res=res, device_mode="full", wall_pec=True)

    incident_power = float(incident_cache["incident_power_by_port"][SOURCE_PORT])
    incident_flux = incident_cache["incident_flux_data_by_port"][SOURCE_PORT]

    P_device.sources = form.make_sources(SOURCE_PORT)
    sim = P_device.Get_Sim()

    monitors = []
    signs = []
    measure_centers = {}
    for out_p in RECEIVE_PORTS:
        measure_xy = port_measure_center(form.name, out_p)
        measure_centers[out_p] = np.asarray(measure_xy, dtype=float).tolist()
        regions, sign = make_flux_region_for_formulation(
            form.name, measure_xy, sc.effective_port_dir(out_p), out_p
        )
        monitors.append(add_flux_monitor(sim, regions))
        signs.append(sign)

    sim.load_minus_flux_data(monitors[RECEIVE_PORTS.index(SOURCE_PORT)], incident_flux)

    print("Running Meep P1 → P1..P6 ...")
    t1 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    t_dev = time.perf_counter() - t1

    flux = extract_flux_powers(sim, monitors, signs)
    column = flux / incident_power
    by_port = {int(p): float(column[k]) for k, p in enumerate(RECEIVE_PORTS)}
    raw = {int(p): float(flux[k]) for k, p in enumerate(RECEIVE_PORTS)}
    db = {
        str(p): (10.0 * np.log10(v) if v > 0 else float("nan"))
        for p, v in by_port.items()
    }

    if _rank() == 0:
        for p, v in by_port.items():
            print(f"  P1->P{p+1}: {v:.8e}  ({db[str(p)]:.3f} dB)")

    # Export numerical mode for FEM
    from plasmeep.ports.mode_registry import get_numerical_mode
    from plasmeep.ports.numerical_launch import port_tangent

    mode = get_numerical_mode(
        SOURCE_PORT,
        res=res,
        frequency_a=sc.fs_a,
        horn_walls="prism",
        grid_offset_cells=(0.0, 0.0),
        coord_rotation_deg=0.0,
    )
    mode_path = OUT / "numerical_mode_P1_res50.json"
    if _rank() == 0:
        mode.save(mode_path)

    src_center = np.asarray(sc.horn_for_port(SOURCE_PORT, res)["source_center"], dtype=float)
    result = {
        "label": "fem_validation_meep_P1_all6_ppc25_df0.10",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {
            "points_per_cm": points_per_cm,
            "res": res,
            "run_time": run_time,
            "source_df": float(sc.source_df),
            "df_frac": 0.10,
            "formulation": FORMULATION,
            "device_mode": "full",
            "B": [0.0, 0.0, 0.0],
            "horn_walls": "prism",
            "grid_offset_cells": [0.0, 0.0],
            "coord_rotation_deg": 0.0,
            "rho": "uniform_fp_8GHz",
        },
        "physical_resolution": physical_resolution_report(points_per_cm=points_per_cm, a_m=sc.a),
        "norm_cache_path": str(prod_norm),
        "timings_s": {"normalization_or_load": t_norm, "device": t_dev, "total": t_norm + t_dev},
        "incident_power": incident_power,
        "source_center_cell": src_center.tolist(),
        "source_tangent": port_tangent(sc.effective_port_dir(SOURCE_PORT)).tolist(),
        "source_outward": np.asarray(sc.effective_port_dir(SOURCE_PORT), dtype=float).tolist(),
        "monitor_centers_cell": measure_centers,
        "clear_width": float(sc.clear_width),
        "span_factor": 0.96,
        "n_flux_points": 31,
        "normalized_power_by_port": by_port,
        "raw_signed_flux_by_port": raw,
        "power_dB_by_port": db,
        "P11": by_port[0],
        "P1_to_P2": by_port[1],
        "P1_to_P3": by_port[2],
        "P1_to_P4": by_port[3],
        "P1_to_P5": by_port[4],
        "P1_to_P6": by_port[5],
        "mode_json": str(mode_path),
        "trusted_crosscheck": {
            "note": "compare to prior P1P2 study",
            "prior_P11": 0.16826699269846146,
            "prior_P1_to_P2": 0.013729681958625567,
            "dP11": by_port[0] - 0.16826699269846146,
            "dP12": by_port[1] - 0.013729681958625567,
        },
    }
    out_path = OUT / "meep_P1_all6_ppc25_df0.10.json"
    if _rank() == 0:
        out_path.write_text(json.dumps(result, indent=2) + "\n")
        print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
