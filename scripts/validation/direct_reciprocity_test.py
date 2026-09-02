#!/usr/bin/env python3
"""
Direct Lorentz reciprocity test without port normalization.

Supports point-probe and matched discrete-overlap observables for B=0 diagnostics.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
for p in (VAL, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from plasmeep.ports.lorentz_probe import (  # noqa: E402
    ProbeSites,
    build_hz_grid_patch,
    discrete_reciprocity_specification,
    evaluate_discrete_reciprocity_pair,
    evaluate_lorentz_pair,
    lorentz_test_specification,
    probe_sites_grid_audit,
)
from sixport_common import (  # noqa: E402
    a,
    build_circulator_device,
    default_uniform_rho,
    effective_port_dir,
    fs_a,
    horn_for_port,
    monitor_center_for_port,
    physical_resolution_report,
    set_geometry_context,
    source_df,
)


OUT = os.path.join(ROOT, "outputs", "validation", "lorentz_direct")


def _sites(port: int, res: int) -> ProbeSites:
    horn = horn_for_port(port, res)
    return ProbeSites(
        port_index=port,
        source_xy=np.asarray(horn["source_center"], dtype=float),
        monitor_xy=np.asarray(monitor_center_for_port(port, res), dtype=float),
        outward_dir=np.asarray(effective_port_dir(port), dtype=float),
    )


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _mpi_info() -> Dict[str, Any]:
    try:
        from mpi_runner import mpi_info

        return mpi_info()
    except Exception:
        size = int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1")))
        rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
        return {
            "ranks": size,
            "rank": rank,
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
            "fi_provider": os.environ.get("FI_PROVIDER", ""),
            "mpich_ch4_netmod": os.environ.get("MPICH_CH4_NETMOD", ""),
        }


def _rank0() -> bool:
    try:
        from mpi4py import MPI

        return int(MPI.COMM_WORLD.Get_rank()) == 0
    except Exception:
        return int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0"))) == 0


def _run_point_pair(
    p_device,
    pa: int,
    pb: int,
    *,
    res: int,
    run_time: float,
) -> Dict[str, Any]:
    pair = evaluate_lorentz_pair(
        p_device,
        _sites(pa, res),
        _sites(pb, res),
        res=res,
        frequency=fs_a,
        fwidth=source_df,
        run_time=run_time,
    )
    d = pair.as_dict()
    d["observable"] = "point"
    d["fields_large_enough"] = d["abs_a_to_b"] > 1e-3 and d["abs_b_to_a"] > 1e-3
    if _rank0():
        print(f"  [point] H_{pa}->{pb} = {pair.a_to_b.hz_complex:.6g}", flush=True)
        print(f"  [point] H_{pb}->{pa} = {pair.b_to_a.hz_complex:.6g}", flush=True)
        print(
            f"  [point] amp_err={pair.amp_err_dB:.4f} dB  "
            f"phase={pair.phase_diff_deg:.2f}°  "
            f"|dH|/mean|H|={pair.complex_sym_err:.4e}",
            flush=True,
        )
    return d


def _run_discrete_pair(
    p_device,
    pa: int,
    pb: int,
    *,
    res: int,
    run_time: float,
    half_width_cells: int,
    weight_kind: str,
) -> Dict[str, Any]:
    sites_a = _sites(pa, res)
    sites_b = _sites(pb, res)
    patch_a = build_hz_grid_patch(
        sites_a.source_xy,
        res=res,
        half_width_cells=half_width_cells,
        weight_kind=weight_kind,  # type: ignore[arg-type]
        label=f"P{pa+1}_source",
    )
    patch_b = build_hz_grid_patch(
        sites_b.source_xy,
        res=res,
        half_width_cells=half_width_cells,
        weight_kind=weight_kind,  # type: ignore[arg-type]
        label=f"P{pb+1}_source",
    )
    pair = evaluate_discrete_reciprocity_pair(
        p_device,
        patch_a,
        patch_b,
        frequency=fs_a,
        fwidth=source_df,
        run_time=run_time,
    )
    d = pair.as_dict()
    d["observable"] = "discrete_overlap"
    if _rank0():
        print(f"  [discrete] G_{pa}->{pb} = {pair.a_to_b.response:.6g}", flush=True)
        print(f"  [discrete] G_{pb}->{pa} = {pair.b_to_a.response:.6g}", flush=True)
        print(
            f"  [discrete] amp_err={pair.amp_err_dB:.4f} dB  "
            f"phase={pair.phase_diff_deg:.2f}°  "
            f"|dG|/mean|G|={pair.complex_sym_err:.4e}",
            flush=True,
        )
    return d


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=20.0)
    parser.add_argument("--pairs", type=str, default="0,1")
    parser.add_argument("--horn-walls", type=str, default="prism")
    parser.add_argument("--device-mode", type=str, default="horns_only", choices=["full", "horns_only"])
    parser.add_argument(
        "--plasma-fill",
        type=str,
        default="active",
        choices=["active", "geometry_only", "dielectric_fill"],
    )
    parser.add_argument("--n-bulbs", type=int, default=None)
    parser.add_argument(
        "--susceptibility-mode",
        type=str,
        default="auto",
        choices=["auto", "force_gyrotropic_b0"],
    )
    parser.add_argument(
        "--observable",
        type=str,
        default="both",
        choices=["point", "discrete", "both"],
        help="point probe, matched discrete overlap, or both",
    )
    parser.add_argument("--patch-half-width", type=int, default=1)
    parser.add_argument("--patch-weights", type=str, default="uniform", choices=["uniform", "gaussian"])
    parser.add_argument("--label", type=str, default="")
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    set_geometry_context(res=args.res, horn_walls=args.horn_walls)

    pairs: List[tuple[int, int]] = []
    for chunk in args.pairs.split(";"):
        a_p, b_p = [int(x) for x in chunk.split(",")]
        pairs.append((a_p, b_p))

    os.makedirs(OUT, exist_ok=True)

    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, p_device, wp_values = build_circulator_device(
        rho,
        B,
        res=args.res,
        device_mode=args.device_mode,
        plasma_fill=args.plasma_fill,
        n_bulbs=args.n_bulbs,
        susceptibility_mode=args.susceptibility_mode,
    )

    probe_audit = []
    grid_audit = []
    for port in sorted({p for pair in pairs for p in pair}):
        s = _sites(port, args.res)
        probe_audit.append(
            {
                "port": port,
                "source_xy": s.source_xy.tolist(),
                "monitor_xy": s.monitor_xy.tolist(),
                "outward_dir": s.outward_dir.tolist(),
                "tangent": s.tangent.tolist(),
            }
        )
        grid_audit.append(probe_sites_grid_audit(s, res=args.res))

    results: List[Dict[str, Any]] = []
    t0 = time.perf_counter()
    for pa, pb in pairs:
        if _rank0():
            print(f"\n=== reciprocity pair P{pa+1}<->P{pb+1} ===", flush=True)
        pair_result: Dict[str, Any] = {"port_a": pa, "port_b": pb}
        if args.observable in ("point", "both"):
            pair_result["point"] = _run_point_pair(
                p_device, pa, pb, res=args.res, run_time=args.run_time
            )
        if args.observable in ("discrete", "both"):
            pair_result["discrete"] = _run_discrete_pair(
                p_device,
                pa,
                pb,
                res=args.res,
                run_time=args.run_time,
                half_width_cells=args.patch_half_width,
                weight_kind=args.patch_weights,
            )
        results.append(pair_result)

    payload = {
        "label": args.label,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "test_specification": {
            "point": lorentz_test_specification(),
            "discrete": discrete_reciprocity_specification(),
        },
        "settings": {
            "res": args.res,
            "run_time": args.run_time,
            "pairs": pairs,
            "horn_walls": args.horn_walls,
            "device_mode": args.device_mode,
            "plasma_fill": args.plasma_fill,
            "n_bulbs": args.n_bulbs,
            "susceptibility_mode": args.susceptibility_mode,
            "observable": args.observable,
            "patch_half_width_cells": args.patch_half_width,
            "patch_weights": args.patch_weights,
            "a_m": a,
            "fs_a": fs_a,
            "fwidth_a": source_df,
            "rho": "uniform_fp_8GHz",
            "B": [0.0, 0.0, 0.0],
            "wp_mean": float(np.mean(wp_values)) if len(wp_values) else 0.0,
        },
        "probe_sites": probe_audit,
        "grid_audit": grid_audit,
        "mpi": _mpi_info(),
        "physical_resolution": physical_resolution_report(args.res),
        "timings_s": {"total": time.perf_counter() - t0},
        "pairs": results,
    }

    out_path = args.json_out or os.path.join(
        OUT, f"lorentz_{args.label or 'run'}_res{args.res}_rt{args.run_time:g}.json"
    )
    if _rank0():
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(_json_safe(payload), f, indent=2)
            f.write("\n")
        print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
