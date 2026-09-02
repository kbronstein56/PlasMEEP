#!/usr/bin/env python3
"""
Direct Lorentz reciprocity test (horns_only, B=0) without port normalization.

Compares complex Hz transfer for swapped localized point sources at two ports.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, List

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
for p in (VAL, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from plasmeep.ports.lorentz_probe import (  # noqa: E402
    ProbeSites,
    evaluate_lorentz_pair,
)
from sixport_common import (  # noqa: E402
    a,
    build_circulator_device,
    default_uniform_rho,
    fs_a,
    horn_for_port,
    monitor_center_for_port,
    physical_resolution_report,
    set_geometry_context,
    source_df,
    effective_port_dir,
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=40.0)
    parser.add_argument("--pairs", type=str, default="0,1;1,2")
    parser.add_argument("--horn-walls", type=str, default="prism")
    parser.add_argument(
        "--device-mode",
        type=str,
        default="horns_only",
        choices=["full", "horns_only"],
    )
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
    _pmm, p_device, _ = build_circulator_device(
        rho, B, res=args.res, device_mode=args.device_mode
    )

    results: List[Dict[str, Any]] = []
    t0 = time.perf_counter()
    for pa, pb in pairs:
        print(f"\n=== Lorentz pair P{pa+1}<->P{pb+1} ===", flush=True)
        pair = evaluate_lorentz_pair(
            p_device,
            _sites(pa, args.res),
            _sites(pb, args.res),
            res=args.res,
            frequency=fs_a,
            fwidth=source_df,
            run_time=args.run_time,
        )
        d = pair.as_dict()
        results.append(d)
        print(
            f"  H_{pa}->{pb} = {pair.a_to_b.hz_complex:.6g}",
            flush=True,
        )
        print(
            f"  H_{pb}->{pa} = {pair.b_to_a.hz_complex:.6g}",
            flush=True,
        )
        print(
            f"  amp_err={pair.amp_err_dB:.4f} dB  "
            f"phase_diff={pair.phase_diff_deg:.3f} deg  "
            f"|dH|/mean|H|={pair.complex_sym_err:.4e}",
            flush=True,
        )

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {
            "res": args.res,
            "run_time": args.run_time,
            "pairs": pairs,
            "horn_walls": args.horn_walls,
            "device_mode": args.device_mode,
            "a_m": a,
            "fs_a": fs_a,
            "source": "localized Hz point (Gaussian), no flux normalization",
        },
        "mpi": _mpi_info(),
        "physical_resolution": physical_resolution_report(args.res),
        "timings_s": {"total": time.perf_counter() - t0},
        "pairs": results,
    }

    out_path = args.json_out or os.path.join(
        OUT, f"lorentz_direct_{args.device_mode}_res{args.res}.json"
    )
    try:
        from mpi4py import MPI

        rank = int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        rank = 0

    if rank == 0:
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(_json_safe(payload), f, indent=2)
            f.write("\n")
        print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
