#!/usr/bin/env python3
"""
Full-device P1↔P2 reciprocity vs run_time for num_mode_hz_line.

Runs increasing run_time values until reciprocity stabilizes or the list ends.
Reuses cached normalizations per (res, run_time) when available.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
OUT = Path(ROOT) / "outputs" / "validation" / "modal_ports"


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def run_one(
    python: str,
    *,
    res: int,
    run_time: float,
    formulation: str,
    skip_norm_cached: bool,
) -> Dict[str, Any]:
    label = f"full_{formulation}_P1P2_res{res}_rt{run_time:g}"
    json_out = OUT / f"{label}.json"
    cmd = [
        python,
        os.path.join(VAL, "reciprocity_b0_study.py"),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--ports",
        "0,1",
        "--port-formulation",
        formulation,
        "--device-mode",
        "full",
        "--label",
        label,
        "--json-out",
        str(json_out),
    ]
    if skip_norm_cached:
        cmd.append("--skip-norm-if-cached")
    print(f"\n>>> {' '.join(cmd)}", flush=True)
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, check=False)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        return {"run_time": run_time, "error": f"exit {proc.returncode}", "wall_s": wall}
    if not json_out.is_file():
        return {"run_time": run_time, "error": "missing json", "wall_s": wall}
    data = json.loads(json_out.read_text(encoding="utf-8"))
    rec = data.get("reciprocity", {})
    pm = data.get("power_matrix_dB", [])
    t12 = t21 = None
    if len(pm) == 2:
        t12 = pm[0][1]
        t21 = pm[1][0]
    return {
        "run_time": run_time,
        "json_out": str(json_out),
        "reciprocity_dB": rec.get("max_abs_diff_dB"),
        "T_P1P2_dB": t12,
        "T_P2P1_dB": t21,
        "incident_mismatch": data.get("incident_power_mismatch"),
        "timings_s": data.get("timings_s"),
        "wall_s": wall,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=64)
    parser.add_argument(
        "--run-times",
        type=str,
        default="10,20,40,60,80",
        help="Comma-separated run_time values to test in order",
    )
    parser.add_argument(
        "--formulation",
        type=str,
        default="num_mode_hz_line",
    )
    parser.add_argument("--tol-dB", type=float, default=0.05)
    parser.add_argument("--skip-norm-if-cached", action="store_true", default=True)
    parser.add_argument("--python", type=str, default=sys.executable)
    args = parser.parse_args()

    run_times = [float(x) for x in args.run_times.split(",") if x.strip()]
    OUT.mkdir(parents=True, exist_ok=True)

    # Audit mode cache first
    audit_cmd = [
        args.python,
        os.path.join(VAL, "audit_mode_cache.py"),
        "--res",
        str(args.res),
        "--json-out",
        str(OUT / f"mode_cache_audit_res{args.res}.json"),
    ]
    subprocess.run(audit_cmd, cwd=ROOT, check=True)

    results: List[Dict[str, Any]] = []
    prev_rec = None
    for rt in run_times:
        row = run_one(
            args.python,
            res=args.res,
            run_time=rt,
            formulation=args.formulation,
            skip_norm_cached=args.skip_norm_if_cached,
        )
        results.append(row)
        rec = row.get("reciprocity_dB")
        print(
            f"  rt={rt:g}: |dP|={rec} dB  T12={row.get('T_P1P2_dB')}  "
            f"T21={row.get('T_P2P1_dB')}  wall={row.get('wall_s', 0):.0f}s",
            flush=True,
        )
        if rec is not None and prev_rec is not None:
            delta = abs(rec - prev_rec)
            row["delta_vs_prev_dB"] = delta
            if delta < args.tol_dB and rec < 0.2:
                print(
                    f"  early stop: |dP| change {delta:.4f} dB < {args.tol_dB}",
                    flush=True,
                )
                break
        prev_rec = rec

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "res": args.res,
        "formulation": args.formulation,
        "run_times_requested": run_times,
        "results": results,
    }
    summary_path = OUT / f"full_device_rt_convergence_res{args.res}.json"
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(_json_safe(payload), f, indent=2)
        f.write("\n")
    print(f"\nWrote {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
