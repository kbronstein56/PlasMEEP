#!/usr/bin/env python3
"""
Performance campaign Phase 1 driver (benchmark-only).

df_frac fixed (default 0.20); vary post-source run_time.
One source P1, monitors P1–P6. Fwidth+rt tagged norm caches.
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import meep as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[3]  # performance_campaign -> validation -> outputs -> wait
# file at outputs/validation/performance_campaign/phase1/run_ringdown.py
# parents[0]=phase1, [1]=performance_campaign, [2]=validation, [3]=outputs, [4]=PlasMEEP
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
from sixport_common import (  # noqa: E402
    build_circulator_device,
    default_uniform_rho,
    ensure_normalizations,
    set_geometry_context,
)

GAUSSIAN_CUTOFF = 5.0
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--df-frac", type=float, default=0.20)
    ap.add_argument("--run-time", type=float, required=True)
    ap.add_argument("--points-per-cm", type=float, default=25.0)
    ap.add_argument("--json-out", type=str, required=True)
    ap.add_argument("--norm-dir", type=str, default="")
    ap.add_argument("--force-norm", action="store_true")
    args = ap.parse_args()

    df_frac = float(args.df_frac)
    run_time = float(args.run_time)
    source_df = df_frac * float(sc.fs_a)
    sc.source_df = source_df

    res = meep_resolution_from_points_per_cm(args.points_per_cm, a_m=sc.a)
    set_geometry_context(
        grid_offset_cells=(0.0, 0.0),
        monitor_offset_cells=(0.0, 0.0),
        res=res,
        horn_walls="prism",
        coord_rotation_deg=0.0,
    )

    width = 1.0 / source_df
    t_src = 2.0 * GAUSSIAN_CUTOFF * width
    t_end = t_src + run_time
    dt = 0.5 / res

    norm_dir = Path(args.norm_dir) if args.norm_dir else Path(__file__).resolve().parent / "norm_cache"
    norm_path = norm_dir / (
        f"norm_{FORMULATION}_full_prism_ppc{args.points_per_cm:g}_res{res}"
        f"_rt{run_time:g}_p0_df{df_frac:.2f}.pkl"
    )

    print(sc.geometry_summary())
    print(f"PHASE1 df_frac={df_frac} run_time={run_time} source_df={source_df:.8f}")
    print(f"T_src={t_src:.4f} T_end={t_end:.4f} steps~{t_end/dt:.1f}")
    print(f"norm_path={norm_path}")

    # Build or load tagged norm
    if norm_path.is_file() and not args.force_norm:
        with open(norm_path, "rb") as f:
            norm = pickle.load(f)
        if abs(float(norm.get("df_frac", -1)) - df_frac) > 1e-12 or abs(float(norm.get("source_df", -1)) - source_df) > 1e-12:
            raise RuntimeError("norm tag mismatch")
        if abs(float(norm.get("run_time", -1)) - run_time) > 1e-12:
            raise RuntimeError("norm run_time mismatch")
        print(f"Loaded tagged norm {norm_path}")
        t_norm = 0.0
    else:
        print("Building NEW normalization (matched df+rt)")
        t0 = time.perf_counter()
        norm = ensure_normalizations(
            res=res, run_time=run_time, force=True, ports=[SOURCE_PORT],
            cache_path=str(norm_path), verbose=True, formulation=FORMULATION,
        )
        t_norm = time.perf_counter() - t0
        norm = dict(norm)
        norm["df_frac"] = df_frac
        norm["source_df"] = source_df
        norm["run_time"] = run_time
        if _mpi_rank() == 0:
            norm_path.parent.mkdir(parents=True, exist_ok=True)
            with open(norm_path, "wb") as f:
                pickle.dump(norm, f, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"Wrote tagged norm {norm_path}")

    assert abs(float(norm["df_frac"]) - df_frac) < 1e-12
    assert abs(float(norm["source_df"]) - source_df) < 1e-12
    print(f"PRE-DEVICE ASSERT OK df={norm['df_frac']} rt={norm.get('run_time')}")

    form = get_formulation(FORMULATION)
    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, P_device, _wp = build_circulator_device(rho, B, res=res, device_mode="full", wall_pec=True)
    inc_p = float(norm["incident_power_by_port"][SOURCE_PORT])
    inc_flux = norm["incident_flux_data_by_port"][SOURCE_PORT]

    P_device.sources = form.make_sources(SOURCE_PORT)
    sim = P_device.Get_Sim()
    monitors, signs = [], []
    for out_p in RECEIVE_PORTS:
        measure_xy = port_measure_center(form.name, out_p)
        regions, sign = make_flux_region_for_formulation(
            form.name, measure_xy, sc.effective_port_dir(out_p), out_p
        )
        monitors.append(add_flux_monitor(sim, regions))
        signs.append(sign)
    sim.load_minus_flux_data(monitors[RECEIVE_PORTS.index(SOURCE_PORT)], inc_flux)

    print("Running Meep device...")
    t0 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    t_dev = time.perf_counter() - t0
    flux = extract_flux_powers(sim, monitors, signs)
    col = flux / inc_p
    by_port = {int(p): float(col[k]) for k, p in enumerate(RECEIVE_PORTS)}
    print("P1 column:", by_port)

    payload = {
        "label": f"phase1_df{df_frac:.2f}_rt{run_time:g}_ppc{args.points_per_cm:g}",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "phase": 1,
        "df_frac": df_frac,
        "source_df": source_df,
        "run_time": run_time,
        "T_src": t_src,
        "T_end_expected": t_end,
        "timesteps_expected": t_end / dt,
        "timings_s": {"normalization": t_norm, "device": t_dev, "total": t_norm + t_dev},
        "incident_power": inc_p,
        "normalized_power_by_port": by_port,
        "P11": by_port.get(0),
        "P1_to_P2": by_port.get(1),
        "P1_to_P3": by_port.get(2),
        "P1_to_P4": by_port.get(3),
        "P1_to_P5": by_port.get(4),
        "P1_to_P6": by_port.get(5),
        "settings": {
            "points_per_cm": args.points_per_cm,
            "res": res,
            "formulation": FORMULATION,
            "device_mode": "full",
            "horn_walls": "prism",
        },
        "physical_resolution": physical_resolution_report(points_per_cm=args.points_per_cm),
        "norm_cache_path": str(norm_path),
    }
    if _mpi_rank() == 0:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(_json_safe(payload), indent=2) + "\n")
        print("Wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
