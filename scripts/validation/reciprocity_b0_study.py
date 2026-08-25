#!/usr/bin/env python3
"""
B=0 reciprocity study for the six-port PMM circulator forward model.

Runs uniform-rho simulate_circulator at B=0, prints reciprocity diagnostics,
and writes a JSON summary. Does not modify production PlasMEEP / PMMI code.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from datetime import datetime, timezone
from typing import Any, Dict, List

import numpy as np

from sixport_common import (
    ensure_normalizations,
    default_uniform_rho,
    reciprocity_metrics,
    simulate_circulator,
    geometry_summary,
    nx_ports,
    ny_ports,
    dpml_ports,
    a,
    fs_Hz,
    fp_Hz,
)


def _parse_ports(s: str | None) -> List[int] | None:
    if s is None or s.strip() == "":
        return None
    ports = [int(p.strip()) for p in s.split(",") if p.strip() != ""]
    for p in ports:
        if p < 0 or p > 5:
            raise argparse.ArgumentTypeError(
                f"port indices must be in 0..5, got {p}"
            )
    if len(ports) != len(set(ports)):
        raise argparse.ArgumentTypeError("duplicate port indices")
    return ports


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="B=0 six-port circulator reciprocity study (uniform rho)."
    )
    parser.add_argument("--res", type=int, default=64, help="Meep resolution")
    parser.add_argument(
        "--run-time",
        type=float,
        default=80.0,
        help="until_after_sources run time",
    )
    parser.add_argument(
        "--mpi-note",
        type=str,
        default="",
        help="Free-text note about MPI / machine / ranks",
    )
    parser.add_argument(
        "--json-out",
        type=str,
        default="",
        help="Path for JSON results (default: reciprocity_b0_<label>.json)",
    )
    parser.add_argument(
        "--label",
        type=str,
        default="b0",
        help="Label embedded in default output filename and JSON",
    )
    parser.add_argument(
        "--ports",
        type=str,
        default="",
        help="Optional subset, e.g. '0,1' for a faster 2x2 diagnostic",
    )
    parser.add_argument(
        "--skip-norm-if-cached",
        action="store_true",
        help="Reuse pickled incident normalizations when present for this res/run-time",
    )
    args = parser.parse_args()

    ports = _parse_ports(args.ports if args.ports else None)

    here = os.path.dirname(os.path.abspath(__file__))
    cache_dir = os.path.join(here, ".cache")
    os.makedirs(cache_dir, exist_ok=True)
    port_tag = (
        "p" + "".join(str(p) for p in ports) if ports is not None else "pall"
    )
    cache_path = os.path.join(
        cache_dir, f"norm_res{args.res}_rt{args.run_time:g}_{port_tag}.pkl"
    )

    json_out = args.json_out
    if not json_out:
        port_tag = (
            "ports" + "".join(str(p) for p in ports) if ports is not None else "all"
        )
        json_out = os.path.join(
            here, f"reciprocity_b0_{args.label}_r{args.res}_{port_tag}.json"
        )

    print(geometry_summary())
    print()
    print(f"label={args.label}  res={args.res}  run_time={args.run_time}")
    print(f"ports={ports if ports is not None else list(range(6))}")
    print(f"mpi_note={args.mpi_note!r}")
    print(f"norm cache={cache_path}")
    print(f"json_out={json_out}")

    t0 = time.perf_counter()
    norm = ensure_normalizations(
        res=args.res,
        run_time=args.run_time,
        force=not args.skip_norm_if_cached,
        ports=ports,
        cache_path=cache_path,
        verbose=True,
    )
    t_norm = time.perf_counter() - t0

    rho = default_uniform_rho()
    B = np.array([0.0, 0.0, 0.0])

    t1 = time.perf_counter()
    result = simulate_circulator(
        rho,
        B,
        res=args.res,
        run_time=args.run_time,
        direction="CCW",
        verbose=True,
        incident_cache=norm,
        ports=ports,
    )
    t_device = time.perf_counter() - t1
    t_total = time.perf_counter() - t0

    metrics = reciprocity_metrics(result["power_matrix_dB"])

    print()
    print("=" * 54)
    print("B=0 RECIPROCITY DIAGNOSTICS")
    print("=" * 54)
    print(f"mean |Pij-Pji| = {metrics['mean_abs_diff_dB']} dB")
    print(f"max  |Pij-Pji| = {metrics['max_abs_diff_dB']} dB")
    for pair in metrics["all_pairs"]:
        diff = pair["abs_diff_dB"]
        print(
            f"  {pair['label']}: "
            f"{'nan' if diff is None else f'{diff:.3f} dB'}"
        )
    if metrics["worst_pairs"]:
        print("worst pairs:")
        for pair in metrics["worst_pairs"]:
            print(f"  {pair['label']}: {pair['abs_diff_dB']:.3f} dB")

    payload: Dict[str, Any] = {
        "label": args.label,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "mpi_note": args.mpi_note,
        "settings": {
            "res": args.res,
            "run_time": args.run_time,
            "B": [0.0, 0.0, 0.0],
            "ports": result["ports"],
            "nx_ports": nx_ports,
            "ny_ports": ny_ports,
            "dpml_ports": dpml_ports,
            "a": a,
            "fs_Hz": fs_Hz,
            "fp_Hz": fp_Hz,
            "rho": "uniform_fp_8GHz",
            "skip_norm_if_cached": bool(args.skip_norm_if_cached),
            "norm_cache_path": cache_path,
        },
        "timings_s": {
            "normalization": t_norm,
            "device": t_device,
            "total": t_total,
        },
        "incident_power_by_port": {
            str(k): float(v) for k, v in norm["incident_power_by_port"].items()
        },
        "power_matrix": result["power_matrix"],
        "power_matrix_dB": result["power_matrix_dB"],
        "reciprocity": {
            "mean_abs_diff_dB": metrics["mean_abs_diff_dB"],
            "max_abs_diff_dB": metrics["max_abs_diff_dB"],
            "worst_pairs": metrics["worst_pairs"],
            "all_pairs": metrics["all_pairs"],
            "error_matrix_dB": metrics["error_matrix_dB"],
        },
        "negative_entries": result["negative_entries"].tolist(),
        "objective": result["objective"],
    }

    with open(json_out, "w", encoding="utf-8") as f:
        json.dump(_json_safe(payload), f, indent=2)
        f.write("\n")

    print()
    print(f"Wrote {json_out}")


if __name__ == "__main__":
    main()
