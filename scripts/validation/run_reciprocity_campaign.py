#!/usr/bin/env python3
"""
Overnight B=0 reciprocity investigation orchestrator.

Runs discriminating tests in priority order:
  1) orientation pair contrast (P1-P2 vs P2-P3) at modest res
  2) resolution sweep on worst pair
  3) run-time sweep at fixed res
  4) optional full 6x6 at best practical res

Writes a rolling JSON index under outputs/validation/reciprocity/.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
PY = "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python"
MPIRUN = "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/mpirun"


def run_study(
    *,
    label: str,
    res: int,
    run_time: float,
    ports: str,
    ranks: int,
    out_dir: str,
    skip_norm_if_cached: bool = True,
) -> dict:
    json_out = os.path.join(out_dir, f"{label}.json")
    log_out = os.path.join(out_dir, f"{label}.log")
    cmd = [
        MPIRUN,
        "-np",
        str(ranks),
        PY,
        os.path.join(VAL, "reciprocity_b0_study.py"),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--ports",
        ports,
        "--label",
        label,
        "--json-out",
        json_out,
        "--mpi-note",
        f"np={ranks}",
    ]
    if skip_norm_if_cached:
        cmd.append("--skip-norm-if-cached")

    env = os.environ.copy()
    env["PYTHONPATH"] = VAL + os.pathsep + os.path.join(ROOT, "scripts") + os.pathsep + env.get("PYTHONPATH", "")

    print(f"\n[{datetime.now().isoformat(timespec='seconds')}] START {label}", flush=True)
    print(" ".join(cmd), flush=True)
    t0 = time.time()
    with open(log_out, "w", encoding="utf-8") as log:
        proc = subprocess.run(
            cmd,
            cwd=VAL,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    wall = time.time() - t0
    entry = {
        "label": label,
        "res": res,
        "run_time": run_time,
        "ports": ports,
        "ranks": ranks,
        "wall_s": wall,
        "returncode": proc.returncode,
        "json": json_out,
        "log": log_out,
    }
    if proc.returncode == 0 and os.path.exists(json_out):
        with open(json_out, encoding="utf-8") as handle:
            data = json.load(handle)
        metrics = data.get("reciprocity", data.get("metrics", {}))
        entry["metrics"] = metrics
        # also pull common keys
        for key in ("mean_abs_dB", "max_abs_dB", "worst_pairs"):
            if key in data:
                entry[key] = data[key]
        if "reciprocity_metrics" in data:
            entry["reciprocity_metrics"] = data["reciprocity_metrics"]
    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] DONE {label} "
        f"wall={wall/60:.1f} min rc={proc.returncode}",
        flush=True,
    )
    return entry


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument(
        "--phase",
        choices=["all", "pairs", "res-sweep", "runtime-sweep", "full6"],
        default="all",
    )
    parser.add_argument("--budget-hours", type=float, default=10.0)
    args = parser.parse_args()

    out_dir = os.path.join(ROOT, "outputs", "validation", "reciprocity")
    os.makedirs(out_dir, exist_ok=True)
    index_path = os.path.join(out_dir, "index.json")
    t_start = time.time()
    budget_s = args.budget_hours * 3600.0
    entries = []

    def within_budget(reserve_min: float = 30.0) -> bool:
        return (time.time() - t_start) < (budget_s - reserve_min * 60.0)

    def save_index() -> None:
        payload = {
            "started": datetime.fromtimestamp(t_start, tz=timezone.utc).isoformat(),
            "updated": datetime.now(timezone.utc).isoformat(),
            "ranks": args.ranks,
            "entries": entries,
        }
        with open(index_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    # Phase A: orientation contrast at modest cost
    if args.phase in ("all", "pairs") and within_budget():
        entries.append(
            run_study(
                label="pair_P1P2_res32_rt40",
                res=32,
                run_time=40.0,
                ports="0,1",
                ranks=args.ranks,
                out_dir=out_dir,
            )
        )
        save_index()
        entries.append(
            run_study(
                label="pair_P2P3_res32_rt40",
                res=32,
                run_time=40.0,
                ports="1,2",
                ranks=args.ranks,
                out_dir=out_dir,
            )
        )
        save_index()

    # Phase B: resolution sweep on worst pair (axis-diagonal)
    if args.phase in ("all", "res-sweep") and within_budget():
        for res in (32, 48, 64, 96):
            if not within_budget(reserve_min=45.0):
                break
            entries.append(
                run_study(
                    label=f"pair_P1P2_res{res}_rt80",
                    res=res,
                    run_time=80.0,
                    ports="0,1",
                    ranks=args.ranks,
                    out_dir=out_dir,
                )
            )
            save_index()

    # Phase C: runtime sweep at fixed res
    if args.phase in ("all", "runtime-sweep") and within_budget():
        for rt in (40.0, 80.0, 160.0):
            if not within_budget(reserve_min=45.0):
                break
            entries.append(
                run_study(
                    label=f"pair_P1P2_res64_rt{int(rt)}",
                    res=64,
                    run_time=rt,
                    ports="0,1",
                    ranks=args.ranks,
                    out_dir=out_dir,
                )
            )
            save_index()

    # Phase D: full 6x6 at practical res if budget remains
    if args.phase in ("all", "full6") and within_budget(reserve_min=20.0):
        entries.append(
            run_study(
                label="full6_res48_rt80",
                res=48,
                run_time=80.0,
                ports="0,1,2,3,4,5",
                ranks=args.ranks,
                out_dir=out_dir,
            )
        )
        save_index()

    print(f"\nIndex written to {index_path}")
    print(f"Total wall time: {(time.time()-t_start)/3600:.2f} h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
