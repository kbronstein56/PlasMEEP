#!/usr/bin/env python3
"""
Horns_only port-formulation gate campaign.

Compares launch/receiver variants, grid-offset robustness, and optional
matched discrete FDTD reciprocity control (np=32 MPI).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "port_gates"
STUDY = VAL / "reciprocity_b0_study.py"

for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from mpi_runner import DEFAULT_PYTHON, DEFAULT_RANKS, run_mpi_script, verify_mpi_ranks  # noqa: E402


FORMULATIONS = (
    "te1_hz_line",
    "num_mode_hz_line",
    "num_mode_guide_normal",
    "num_mode_yee_sdotn",
    "num_mode_yee_guide_normal",
)

GRID_OFFSETS: Tuple[Tuple[float, float], ...] = (
    (0.0, 0.0),
    (-0.5, 0.0),
    (0.5, 0.0),
    (0.0, -0.5),
    (0.0, 0.5),
    (-0.5, -0.5),
    (0.5, 0.5),
)


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _pair_err(data: Dict[str, Any], pa: int, pb: int) -> Optional[float]:
    for pair in data.get("reciprocity", {}).get("all_pairs", []):
        if pair.get("i") == pa and pair.get("j") == pb:
            return pair.get("abs_diff_dB")
    return None


def _run_case(
    *,
    label: str,
    formulation: str,
    ports: str,
    res: int,
    run_time: float,
    grid_offset: Tuple[float, float],
    device_mode: str,
    ranks: int,
    discrete_control: bool,
    skip_existing: bool,
) -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    gox, goy = grid_offset
    json_out = OUT / f"{label}.json"
    log_out = OUT / f"{label}.log"
    if skip_existing and json_out.is_file():
        data = json.loads(json_out.read_text(encoding="utf-8"))
        pa, pb = [int(x) for x in ports.split(",")]
        return {
            "label": label,
            "skipped": True,
            "json_out": str(json_out),
            "port_err_dB": _pair_err(data, pa, pb),
            "discrete_amp_err_dB": data.get("discrete_control", {}).get("amp_err_dB"),
        }

    args = [
        str(STUDY),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--ports",
        ports,
        "--port-formulation",
        formulation,
        "--device-mode",
        device_mode,
        f"--grid-offset-cells={gox},{goy}",
        "--label",
        label,
        "--json-out",
        str(json_out),
        "--skip-norm-if-cached",
    ]
    if discrete_control:
        args.append("--discrete-control")

    t0 = time.perf_counter()
    proc = run_mpi_script(args, ranks=ranks, log_path=log_out)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"case {label} failed exit={proc.returncode} log={log_out}")

    data = json.loads(json_out.read_text(encoding="utf-8"))
    pa, pb = [int(x) for x in ports.split(",")]
    return {
        "label": label,
        "skipped": False,
        "json_out": str(json_out),
        "formulation": formulation,
        "ports": ports,
        "grid_offset": list(grid_offset),
        "port_err_dB": _pair_err(data, pa, pb),
        "discrete_amp_err_dB": data.get("discrete_control", {}).get("amp_err_dB"),
        "elapsed_s": elapsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=5.0)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--device-mode", type=str, default="horns_only")
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--no-skip-existing", dest="skip_existing", action="store_false")
    parser.add_argument("--grid-offset-sweep", action="store_true", default=False)
    parser.add_argument("--full-device", action="store_true", default=False)
    parser.add_argument("--discrete-control", action="store_true", default=True)
    parser.add_argument(
        "--formulations",
        type=str,
        default=",".join(FORMULATIONS),
        help="Comma-separated subset of formulations",
    )
    args = parser.parse_args()

    verify_mpi_ranks(ranks=args.ranks, python=DEFAULT_PYTHON, cwd=VAL)
    device_mode = "full" if args.full_device else args.device_mode
    formulations = [f.strip() for f in args.formulations.split(",") if f.strip()]

    cases: List[Dict[str, Any]] = []
    offsets = GRID_OFFSETS if args.grid_offset_sweep else ((0.0, 0.0),)
    for form in formulations:
        for ports in ("0,1", "1,2"):
            for gox, goy in offsets:
                tag = f"{form}_P{ports.replace(',', 'P')}_g{gox:g}_{goy:g}"
                if args.full_device:
                    tag = f"full_{tag}"
                cases.append(
                    {
                        "label": f"{tag}_res{args.res}_rt{args.run_time:g}",
                        "formulation": form,
                        "ports": ports,
                        "grid_offset": (gox, goy),
                    }
                )

    rows: List[Dict[str, Any]] = []
    t_campaign = time.perf_counter()
    for spec in cases:
        print(f"\n=== {spec['label']} ===", flush=True)
        row = _run_case(
            label=spec["label"],
            formulation=spec["formulation"],
            ports=spec["ports"],
            res=args.res,
            run_time=args.run_time,
            grid_offset=spec["grid_offset"],
            device_mode=device_mode,
            ranks=args.ranks,
            discrete_control=args.discrete_control,
            skip_existing=args.skip_existing,
        )
        rows.append(row)
        print(
            f"  port_err={row.get('port_err_dB')} dB  "
            f"discrete={row.get('discrete_amp_err_dB')} dB  "
            f"skipped={row.get('skipped')}",
            flush=True,
        )

    ablation = [
        {
            "formulation": r["formulation"],
            "ports": r["ports"],
            "grid_offset": r.get("grid_offset"),
            "port_err_dB": r.get("port_err_dB"),
            "discrete_amp_err_dB": r.get("discrete_amp_err_dB"),
        }
        for r in rows
        if r.get("grid_offset") == [0.0, 0.0] or r.get("grid_offset") == (0.0, 0.0)
    ]

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {
            "res": args.res,
            "run_time": args.run_time,
            "ranks": args.ranks,
            "device_mode": device_mode,
            "grid_offset_sweep": args.grid_offset_sweep,
            "gates": {
                "P1P2_le_0.2_dB": 0.2,
                "P2P3_le_0.05_dB": 0.05,
            },
        },
        "ablation_zero_offset": ablation,
        "rows": rows,
        "elapsed_s": time.perf_counter() - t_campaign,
    }
    summary_path = OUT / "campaign_summary.json"
    summary_path.write_text(json.dumps(_json_safe(summary), indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
