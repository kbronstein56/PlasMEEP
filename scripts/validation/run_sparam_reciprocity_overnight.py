#!/usr/bin/env python3
"""
Overnight S-parameter reciprocity decision-tree campaign.

Runs dual-receiver diagnostics (flux + modal) through horns-only then full
device at 50 points/cm, with optional monitor-depth and subtraction-bypass
variants. One np=32 job at a time.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "sparam_audit"
DUAL = VAL / "dual_receiver_diagnostic.py"

for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from mpi_runner import DEFAULT_RANKS, run_mpi_script, verify_mpi_ranks  # noqa: E402


def _run(
    *,
    label: str,
    ports: str,
    device_mode: str,
    points_per_cm: float,
    run_time: float,
    skip_sub: bool,
    monitor_cells: float,
    ranks: int,
    skip_existing: bool,
) -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    json_out = OUT / f"{label}.json"
    log_out = OUT / f"{label}.log"
    if skip_existing and json_out.is_file():
        data = json.loads(json_out.read_text(encoding="utf-8"))
        pairs = data.get("results", {}).get("pair_summaries", [])
        return {
            "label": label,
            "skipped": True,
            "json_out": str(json_out),
            "pairs": pairs,
        }

    args = [
        str(DUAL),
        "--points-per-cm",
        str(points_per_cm),
        "--run-time",
        str(run_time),
        "--ports",
        ports,
        "--device-mode",
        device_mode,
        "--label",
        label,
        "--json-out",
        str(json_out),
        f"--monitor-outward-cells={monitor_cells}",
        "--discrete-control",
    ]
    if skip_sub:
        args.append("--skip-flux-subtraction")

    t0 = time.perf_counter()
    proc = run_mpi_script(args, ranks=ranks, log_path=log_out)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"{label} failed exit={proc.returncode} log={log_out}")
    data = json.loads(json_out.read_text(encoding="utf-8"))
    return {
        "label": label,
        "skipped": False,
        "json_out": str(json_out),
        "wall_s": elapsed,
        "pairs": data.get("results", {}).get("pair_summaries", []),
        "discrete_amp_err_dB": (data.get("results", {}) or {})
        .get("discrete_control", {})
        .get("amp_err_dB"),
        "flux_subtraction_healthy": (data.get("results", {}) or {})
        .get("flux_subtraction_health", {})
        .get("healthy"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points-per-cm", type=float, default=50.0)
    parser.add_argument("--run-time", type=float, default=20.0)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument(
        "--stages",
        type=str,
        default="horns,full,monitor",
        help="Comma list: horns,full,monitor",
    )
    args = parser.parse_args()

    verify_mpi_ranks(ranks=args.ranks)
    stages = {s.strip() for s in args.stages.split(",") if s.strip()}
    results: List[Dict[str, Any]] = []

    # Phase 3: horns-only P1P2 and P2P3
    if "horns" in stages:
        for ports, tag in (("0,1", "P1P2"), ("1,2", "P2P3")):
            for skip_sub in (False, True):
                sub = "nosub" if skip_sub else "sub"
                label = f"horns_dual_{tag}_ppc{args.points_per_cm:g}_rt{args.run_time:g}_{sub}"
                print(f"\n>>> {label}", flush=True)
                row = _run(
                    label=label,
                    ports=ports,
                    device_mode="horns_only",
                    points_per_cm=args.points_per_cm,
                    run_time=args.run_time,
                    skip_sub=skip_sub,
                    monitor_cells=0.0,
                    ranks=args.ranks,
                    skip_existing=args.skip_existing,
                )
                print(row)
                results.append(row)

    # Phase 2/9: full device — only with subtraction on (matches production);
    # nosub transmission should match sub (Phase 1 finding) — verify once on P1P2
    if "full" in stages:
        for ports, tag in (("0,1", "P1P2"), ("1,2", "P2P3")):
            label = f"full_dual_{tag}_ppc{args.points_per_cm:g}_rt{args.run_time:g}_sub"
            print(f"\n>>> {label}", flush=True)
            row = _run(
                label=label,
                ports=ports,
                device_mode="full",
                points_per_cm=args.points_per_cm,
                run_time=args.run_time,
                skip_sub=False,
                monitor_cells=0.0,
                ranks=args.ranks,
                skip_existing=args.skip_existing,
            )
            print(row)
            results.append(row)

        # One nosub confirmation on full P1P2
        label = f"full_dual_P1P2_ppc{args.points_per_cm:g}_rt{args.run_time:g}_nosub"
        print(f"\n>>> {label}", flush=True)
        row = _run(
            label=label,
            ports="0,1",
            device_mode="full",
            points_per_cm=args.points_per_cm,
            run_time=args.run_time,
            skip_sub=True,
            monitor_cells=0.0,
            ranks=args.ranks,
            skip_existing=args.skip_existing,
        )
        print(row)
        results.append(row)

    # Phase 6: monitor farther into guide (horns-only, cheap)
    if "monitor" in stages:
        for cells in (2.0, 5.0):
            label = f"horns_dual_P1P2_ppc{args.points_per_cm:g}_rt{args.run_time:g}_mon{cells:g}"
            print(f"\n>>> {label}", flush=True)
            row = _run(
                label=label,
                ports="0,1",
                device_mode="horns_only",
                points_per_cm=args.points_per_cm,
                run_time=args.run_time,
                skip_sub=False,
                monitor_cells=cells,
                ranks=args.ranks,
                skip_existing=args.skip_existing,
            )
            print(row)
            results.append(row)

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "points_per_cm": args.points_per_cm,
        "run_time": args.run_time,
        "cases": results,
    }
    path = OUT / "overnight_sparam_summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
