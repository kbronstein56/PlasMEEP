#!/usr/bin/env python3
"""
A/B decomposition of port launch and receiver formulations on horns_only.

Compares te1 vs numerical launch crossed with flux vs modal-overlap receiver.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
for p in (VAL, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from plasmeep.ports.modal_receiver import modal_coefficient_metrics  # noqa: E402
from port_formulations import uses_modal_measurement  # noqa: E402
from sixport_common import (  # noqa: E402
    default_uniform_rho,
    ensure_normalizations,
    physical_resolution_report,
    reciprocity_metrics,
    set_geometry_context,
    simulate_circulator,
)

OUT = os.path.join(ROOT, "outputs", "validation", "modal_ports")
CACHE_DIR = os.path.join(VAL, ".cache")


def norm_cache_path(
    formulation: str,
    device_mode: str,
    res: int,
    run_time: float,
    ports: List[int],
) -> str:
    port_tag = "p" + "".join(str(p) for p in ports)
    return os.path.join(
        CACHE_DIR,
        f"norm_{formulation}_{device_mode}_prism_rot0_pec_g0_0_m0_0"
        f"_res{res}_rt{run_time:g}_{port_tag}.pkl",
    )

FORMULATIONS = (
    "te1_hz_line",
    "num_mode_hz_line",
    "te1_hz_modal",
    "num_mode_modal",
)


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _pair_modal_reciprocity(
    result: Dict[str, Any], pa: int, pb: int
) -> Optional[Dict[str, float]]:
    coeffs = result.get("modal_coefficients")
    if not coeffs:
        return None
    cab = coeffs.get((pb, pa))
    cba = coeffs.get((pa, pb))
    if cab is None or cba is None:
        return None
    metrics = modal_coefficient_metrics(complex(cab), complex(cba))
    metrics["coeff_ab"] = complex(cab)
    metrics["coeff_ba"] = complex(cba)
    return metrics


def run_case(
    formulation: str,
    ports: Tuple[int, int],
    *,
    res: int,
    run_time: float,
    device_mode: str,
    skip_norm_if_cached: bool,
) -> Dict[str, Any]:
    set_geometry_context(res=res, horn_walls="prism")
    port_list = list(ports)
    cache_path = norm_cache_path(
        formulation, device_mode, res, run_time, port_list
    )

    t0 = time.perf_counter()
    norm = ensure_normalizations(
        res,
        run_time=run_time,
        force=not skip_norm_if_cached,
        ports=port_list,
        cache_path=cache_path,
        verbose=False,
        formulation=formulation,
    )
    t_norm = time.perf_counter() - t0

    rho = default_uniform_rho()
    B = np.array([0.0, 0.0, 0.0])
    t1 = time.perf_counter()
    result = simulate_circulator(
        rho,
        B,
        res=res,
        run_time=run_time,
        verbose=False,
        incident_cache=norm,
        ports=port_list,
        formulation=formulation,
        device_mode=device_mode,
    )
    t_sim = time.perf_counter() - t1

    metrics = reciprocity_metrics(result["power_matrix_dB"], port_ids=result["ports"])
    pa, pb = ports
    modal_pair = _pair_modal_reciprocity(result, pa, pb)

    return {
        "formulation": formulation,
        "ports": list(ports),
        "reciprocity": metrics,
        "max_abs_diff_dB": metrics["max_abs_diff_dB"],
        "modal_reciprocity": modal_pair,
        "power_matrix_dB": result["power_matrix_dB"].tolist(),
        "timings_s": {
            "normalization": t_norm,
            "simulate": t_sim,
            "total": t_norm + t_sim,
        },
        "norm_cache_path": cache_path,
        "norm_cache_hit": skip_norm_if_cached and os.path.isfile(cache_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=5.0)
    parser.add_argument("--device-mode", type=str, default="horns_only")
    parser.add_argument("--pairs", type=str, default="0,1;1,2")
    parser.add_argument(
        "--formulations",
        type=str,
        default=",".join(FORMULATIONS),
    )
    parser.add_argument("--skip-norm-if-cached", action="store_true", default=True)
    parser.add_argument("--force-norm", action="store_true")
    args = parser.parse_args()

    pairs: List[Tuple[int, int]] = []
    for chunk in args.pairs.split(";"):
        a_p, b_p = [int(x) for x in chunk.split(",")]
        pairs.append((a_p, b_p))

    formulations = [f.strip() for f in args.formulations.split(",") if f.strip()]
    os.makedirs(OUT, exist_ok=True)

    payload: Dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {
            "res": args.res,
            "run_time": args.run_time,
            "device_mode": args.device_mode,
            "pairs": pairs,
            "formulations": formulations,
        },
        "physical_resolution": physical_resolution_report(args.res),
        "cases": [],
    }

    skip = args.skip_norm_if_cached and not args.force_norm
    for form in formulations:
        for ports in pairs:
            label = f"{form}_P{ports[0]+1}P{ports[1]+1}"
            print(f"\n=== {label} ===", flush=True)
            case = run_case(
                form,
                ports,
                res=args.res,
                run_time=args.run_time,
                device_mode=args.device_mode,
                skip_norm_if_cached=skip,
            )
            payload["cases"].append(case)
            print(
                f"  |dP|_max={case['max_abs_diff_dB']:.4f} dB  "
                f"wall={case['timings_s']['total']:.1f}s",
                flush=True,
            )
            if case.get("modal_reciprocity"):
                mr = case["modal_reciprocity"]
                print(
                    f"  modal amp_err={mr['amp_err_dB']:.4f} dB  "
                    f"phase={mr['phase_diff_deg']:.2f}°",
                    flush=True,
                )

    out_path = os.path.join(
        OUT, f"modal_ablation_res{args.res}_rt{args.run_time:g}.json"
    )
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(_json_safe(payload), f, indent=2)
        f.write("\n")
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
