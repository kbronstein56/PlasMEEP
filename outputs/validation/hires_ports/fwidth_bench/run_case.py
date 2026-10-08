#!/usr/bin/env python3
"""
Benchmark-only fwidth sweep driver (NOT production).

One source (port 0 / P1), monitors on receivers 0..5 in the SAME Meep solve.
Overrides sixport_common.source_df for this process only; writes fwidth-tagged
normalization pickles so production 0.10 caches cannot be reused.
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

ROOT = Path(__file__).resolve().parents[4]
VAL = ROOT / "scripts" / "validation"
for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import sixport_common as sc  # noqa: E402
from physical_units import (  # noqa: E402
    meep_resolution_from_points_per_cm,
    physical_resolution_report,
)
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
        return int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))


def _mpi_info() -> Dict[str, Any]:
    try:
        from mpi_runner import mpi_info

        return mpi_info()
    except Exception:
        return {
            "ranks": int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1"))),
            "rank": _mpi_rank(),
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
        }


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, complex):
        return {"re": float(obj.real), "im": float(obj.imag)}
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if obj is None or isinstance(obj, (str, bool, int, float)):
        return obj
    return str(obj)


def source_timing(source_df: float, run_time: float) -> Dict[str, float]:
    width = 1.0 / float(source_df)
    t_src = 2.0 * GAUSSIAN_CUTOFF * width
    t_end = t_src + float(run_time)
    return {
        "gaussian_cutoff": GAUSSIAN_CUTOFF,
        "gaussian_width": width,
        "T_src": t_src,
        "post_source_run_time": float(run_time),
        "T_end_expected": t_end,
    }


def tagged_norm_path(cache_dir: Path, *, points_per_cm: float, res: int, run_time: float, df_frac: float) -> Path:
    return cache_dir / (
        f"norm_{FORMULATION}_full_prism_rot0_pec_g0_0_m0_0_"
        f"ppc{points_per_cm:g}_res{res}_rt{run_time:g}_p0_df{df_frac:.2f}.pkl"
    )


def assert_norm_matches_df(norm: Dict[str, Any], df_frac: float, source_df: float) -> None:
    got_frac = norm.get("df_frac")
    got_df = norm.get("source_df")
    if got_frac is None or got_df is None:
        raise RuntimeError(
            "Normalization cache missing df_frac/source_df tags — refusing to use "
            "an untagged (likely 0.10) production pickle."
        )
    if abs(float(got_frac) - float(df_frac)) > 1e-12:
        raise RuntimeError(f"Norm df_frac={got_frac} != requested {df_frac}")
    if abs(float(got_df) - float(source_df)) > 1e-12:
        raise RuntimeError(f"Norm source_df={got_df} != requested {source_df}")
    if abs(float(sc.source_df) - float(source_df)) > 1e-12:
        raise RuntimeError("sixport_common.source_df was mutated unexpectedly")


def build_tagged_normalization(
    *,
    res: int,
    run_time: float,
    df_frac: float,
    source_df: float,
    cache_path: Path,
    force: bool,
) -> Dict[str, Any]:
    """Always write/read only fwidth-tagged pickles; never production untagged paths."""
    cache_path = Path(cache_path)
    if cache_path.is_file() and not force:
        with open(cache_path, "rb") as f:
            loaded = pickle.load(f)
        assert_norm_matches_df(loaded, df_frac, source_df)
        print(f"Loaded fwidth-tagged normalizations from {cache_path}")
        return loaded

    # Ensure production untagged caches cannot be selected by ensure_normalizations.
    # We pass an explicit tagged path and force a fresh normalize_port via force=True
    # when the tagged file is absent.
    print(f"Building NEW normalization for df_frac={df_frac:.2f} source_df={source_df:.8f}")
    print(f"Tagged norm path: {cache_path}")
    t0 = time.perf_counter()
    norm = ensure_normalizations(
        res=res,
        run_time=run_time,
        force=True,
        ports=[SOURCE_PORT],
        cache_path=str(cache_path),
        verbose=True,
        formulation=FORMULATION,
    )
    t_norm = time.perf_counter() - t0
    norm = dict(norm)
    norm["df_frac"] = float(df_frac)
    norm["source_df"] = float(source_df)
    norm["fs_a"] = float(sc.fs_a)
    norm["benchmark_norm_wall_s"] = float(t_norm)
    if _mpi_rank() == 0:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump(norm, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Wrote tagged normalization to {cache_path}")
    assert_norm_matches_df(norm, df_frac, source_df)
    return norm


def simulate_one_source_multi_receiver(
    *,
    res: int,
    run_time: float,
    incident_cache: Dict[str, Any],
    receive_ports: List[int],
) -> Dict[str, Any]:
    """One full-device Meep solve: source P1, flux monitors on receive_ports."""
    form = get_formulation(FORMULATION)
    if form.measurement != "flux":
        raise RuntimeError(f"expected flux measurement, got {form.measurement}")

    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, P_device, wp = build_circulator_device(
        rho, B, res=res, device_mode="full", wall_pec=True
    )

    incident_power = float(incident_cache["incident_power_by_port"][SOURCE_PORT])
    incident_flux = incident_cache["incident_flux_data_by_port"][SOURCE_PORT]

    P_device.sources = form.make_sources(SOURCE_PORT)
    sim = P_device.Get_Sim()

    monitors = []
    signs = []
    for out_p in receive_ports:
        measure_xy = port_measure_center(form.name, out_p)
        regions, sign = make_flux_region_for_formulation(
            form.name, measure_xy, sc.effective_port_dir(out_p), out_p
        )
        monitors.append(add_flux_monitor(sim, regions))
        signs.append(sign)

    src_mon_idx = receive_ports.index(SOURCE_PORT)
    sim.load_minus_flux_data(monitors[src_mon_idx], incident_flux)

    print("Running Meep (one source P1, multi-receiver)...")
    print(f"  sc.source_df={sc.source_df:.8f}  receive_ports={receive_ports}")
    t0 = time.perf_counter()
    sim.run(until_after_sources=run_time)
    t_device = time.perf_counter() - t0

    flux = extract_flux_powers(sim, monitors, signs)
    column = flux / incident_power

    # Meep end time / steps from simulation state when available
    actual_t = float(getattr(sim, "meep_time", lambda: float("nan"))())
    try:
        # sim.fields.t is timestep index in some Meep versions
        nsteps = int(sim.fields.t) if hasattr(sim, "fields") and sim.fields is not None else None
    except Exception:
        nsteps = None

    by_port = {int(p): float(column[k]) for k, p in enumerate(receive_ports)}
    flux_by_port = {int(p): float(flux[k]) for k, p in enumerate(receive_ports)}

    print("P1 column (normalized power):")
    for p, val in by_port.items():
        print(f"  -> P{p + 1}: {val:.6e}")

    return {
        "device_wall_s": float(t_device),
        "meep_time_actual": actual_t,
        "nsteps_guess": nsteps,
        "incident_power": incident_power,
        "normalized_power_by_port": by_port,
        "raw_signed_flux_by_port": flux_by_port,
        "receive_ports": list(receive_ports),
        "wp": np.asarray(wp, dtype=float),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--df-frac", type=float, required=True)
    parser.add_argument("--points-per-cm", type=float, default=25.0)
    parser.add_argument("--run-time", type=float, default=20.0)
    parser.add_argument("--json-out", type=str, required=True)
    parser.add_argument(
        "--norm-dir",
        type=str,
        default=str(Path(__file__).resolve().parent / "norm_cache"),
    )
    parser.add_argument(
        "--force-norm",
        action="store_true",
        help="Rebuild tagged norm even if present",
    )
    args = parser.parse_args()

    df_frac = float(args.df_frac)
    if df_frac <= 0:
        raise SystemExit("df-frac must be positive")

    # Process-local override — does not edit production source files.
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

    timing = source_timing(source_df, args.run_time)
    # dt = Courant/res with default Courant 0.5
    dt = 0.5 / float(res)
    timing["dt_expected"] = dt
    timing["timesteps_expected"] = timing["T_end_expected"] / dt

    print(sc.geometry_summary())
    print()
    print(f"BENCHMARK fwidth case df_frac={df_frac:.2f}")
    print(f"  fs_a={sc.fs_a:.12f}  source_df={source_df:.12f}")
    print(f"  T_src={timing['T_src']:.6f}  T_end_expected={timing['T_end_expected']:.6f}")
    print(f"  points_per_cm={args.points_per_cm:g} res={res} run_time={args.run_time}")
    print(f"  formulation={FORMULATION}  source_port={SOURCE_PORT}  receivers={RECEIVE_PORTS}")

    # Refuse accidental use of production untagged pickle names.
    prod_untagged = VAL / ".cache" / (
        f"norm_{FORMULATION}_full_prism_rot0_pec_g0_0_m0_0_"
        f"ppc{args.points_per_cm:g}_res{res}_rt{args.run_time:g}_p0.pkl"
    )
    print(f"  production_untagged_path (must not load)={prod_untagged}")

    norm_path = tagged_norm_path(
        Path(args.norm_dir),
        points_per_cm=args.points_per_cm,
        res=res,
        run_time=args.run_time,
        df_frac=df_frac,
    )

    t_all0 = time.perf_counter()
    t_norm0 = time.perf_counter()
    norm = build_tagged_normalization(
        res=res,
        run_time=args.run_time,
        df_frac=df_frac,
        source_df=source_df,
        cache_path=norm_path,
        force=bool(args.force_norm) or (not norm_path.is_file()),
    )
    t_norm = time.perf_counter() - t_norm0
    if "benchmark_norm_wall_s" in norm and not args.force_norm and norm_path.is_file():
        # If loaded from disk, wall is ~0; keep stored build time separately.
        norm_build_s = float(norm.get("benchmark_norm_wall_s", t_norm))
    else:
        norm_build_s = float(norm.get("benchmark_norm_wall_s", t_norm))

    assert_norm_matches_df(norm, df_frac, source_df)
    print(
        f"PRE-DEVICE ASSERT OK: using norm df_frac={norm['df_frac']} "
        f"source_df={norm['source_df']}"
    )

    device = simulate_one_source_multi_receiver(
        res=res,
        run_time=args.run_time,
        incident_cache=norm,
        receive_ports=RECEIVE_PORTS,
    )
    t_total = time.perf_counter() - t_all0

    P = device["normalized_power_by_port"]
    payload = {
        "label": f"fwidth_ppc{args.points_per_cm:g}_df{df_frac:.2f}",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "benchmark": "fwidth_stage1",
        "mpi": _mpi_info(),
        "df_frac": df_frac,
        "source_df": source_df,
        "fs_a": float(sc.fs_a),
        "fs_Hz": float(sc.fs_Hz),
        "source_timing": timing,
        "settings": {
            "points_per_cm": float(args.points_per_cm),
            "res": int(res),
            "run_time": float(args.run_time),
            "a": float(sc.a),
            "B": [0.0, 0.0, 0.0],
            "rho": "uniform_fp_8GHz",
            "port_formulation": FORMULATION,
            "device_mode": "full",
            "horn_walls": "prism",
            "grid_offset_cells": [0.0, 0.0],
            "coord_rotation_deg": 0.0,
            "source_port": SOURCE_PORT,
            "receive_ports": RECEIVE_PORTS,
            "measurement_note": (
                "One Meep solve: source P1 only; flux monitors on P1..P6. "
                "Production simulate_circulator uses the same list for sources "
                "and receivers; this harness separates them without changing production."
            ),
        },
        "physical_resolution": physical_resolution_report(
            points_per_cm=float(args.points_per_cm)
        ),
        "norm_cache_path": str(norm_path),
        "timings_s": {
            "normalization_build_or_load": t_norm,
            "normalization_build_stored": norm_build_s,
            "device": device["device_wall_s"],
            "total": t_total,
        },
        "incident_power": device["incident_power"],
        "normalized_power_by_port": device["normalized_power_by_port"],
        "raw_signed_flux_by_port": device["raw_signed_flux_by_port"],
        "P11": P.get(0),
        "P1_to_P2": P.get(1),
        "P1_to_P3": P.get(2),
        "P1_to_P4": P.get(3),
        "P1_to_P5": P.get(4),
        "P1_to_P6": P.get(5),
        "meep_time_actual": device["meep_time_actual"],
        "nsteps_guess": device["nsteps_guess"],
    }

    if _mpi_rank() == 0:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(_json_safe(payload), f, indent=2)
            f.write("\n")
        print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    # Avoid meep importing issues with unused mp in some paths
    _ = mp
    raise SystemExit(main())
