#!/usr/bin/env python3
"""Instrumented trusted 50 ppc one-source profile (benchmark-only). Production defaults."""
from __future__ import annotations

import json
import os
import pickle
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import meep as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
VAL = ROOT / "scripts" / "validation"
for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import sixport_common as sc  # noqa: E402
from physical_units import physical_resolution_report  # noqa: E402
from port_formulations import (  # noqa: E402
    add_flux_monitor,
    extract_flux_powers,
    get_formulation,
    make_flux_region_for_formulation,
    port_measure_center,
)
from sixport_common import build_circulator_device, default_uniform_rho, set_geometry_context  # noqa: E402

FORMULATION = "num_mode_guide_normal"
RECEIVE = list(range(6))
SRC = 0


def main() -> int:
    res = 100
    run_time = 20.0
    assert abs(sc.source_df / sc.fs_a - 0.10) < 1e-12, "must use production 0.10 fs"
    set_geometry_context(res=res, horn_walls="prism", grid_offset_cells=(0, 0), coord_rotation_deg=0.0)

    # Reuse existing full-device port-0,1 norm; take port 0 only
    norm_path = VAL / ".cache" / (
        "norm_num_mode_guide_normal_full_prism_rot0_pec_g0_0_m0_0_ppc50_res100_rt20_p01.pkl"
    )
    with open(norm_path, "rb") as f:
        norm = pickle.load(f)
    assert 0 in norm["incident_power_by_port"]
    print(f"Loaded production-compatible norm {norm_path}")
    print(f"source_df={sc.source_df} (=0.10*fs)")

    form = get_formulation(FORMULATION)
    rho = default_uniform_rho()
    B = np.zeros(3)

    t_setup0 = time.perf_counter()
    _pmm, P_device, _wp = build_circulator_device(rho, B, res=res, device_mode="full", wall_pec=True)
    P_device.sources = form.make_sources(SRC)
    sim = P_device.Get_Sim()
    monitors, signs = [], []
    for out_p in RECEIVE:
        xy = port_measure_center(form.name, out_p)
        regions, sign = make_flux_region_for_formulation(
            form.name, xy, sc.effective_port_dir(out_p), out_p
        )
        monitors.append(add_flux_monitor(sim, regions))
        signs.append(sign)
    sim.load_minus_flux_data(monitors[0], norm["incident_flux_data_by_port"][SRC])
    t_setup = time.perf_counter() - t_setup0

    print("Running profiled Meep device (50 ppc, P1 source, P1–P6 monitors)...")
    t0 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    t_dev = time.perf_counter() - t0
    flux = extract_flux_powers(sim, monitors, signs)
    inc = float(norm["incident_power_by_port"][SRC])
    col = flux / inc
    by = {int(p): float(col[k]) for k, p in enumerate(RECEIVE)}

    T = 2 * 5 / sc.source_df + run_time
    dt = 0.5 / res
    nsteps = int(round(T / dt))
    ncells = 30 * res * 28 * res
    payload = {
        "label": "profile_trusted_ppc50_P1_multirec",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "phase": 5,
        "df_frac": 0.10,
        "run_time": run_time,
        "T_end_expected": T,
        "nsteps_expected": nsteps,
        "ncells": ncells,
        "timings_s": {
            "python_setup_build": t_setup,
            "device": t_dev,
            "seconds_per_timestep": t_dev / nsteps,
            "cell_updates_per_sec": (ncells * nsteps) / t_dev,
        },
        "incident_power": inc,
        "normalized_power_by_port": by,
        "P11": by[0],
        "P1_to_P2": by[1],
        "settings": {"points_per_cm": 50, "res": 100, "formulation": FORMULATION, "np": 32},
        "physical_resolution": physical_resolution_report(points_per_cm=50),
        "norm_cache_path": str(norm_path),
    }
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
    except Exception:
        rank = 0
    if rank == 0:
        out = Path(__file__).resolve().parent / "trusted_ppc50_profile.json"
        out.write_text(json.dumps(payload, indent=2) + "\n")
        print("Wrote", out)
        print(json.dumps(payload["timings_s"], indent=2))
        print("P1 column", by)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
