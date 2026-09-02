#!/usr/bin/env python3
"""
Six-hour autonomous port-development block campaign.

Phases (skip completed JSON):
  1. horns_only formulation gates (incl. num_mode_yee_guide_normal)
  2. grid-offset robustness on best candidates
  3. full-device res32 rt convergence
  4. res64 confirmation (if res32 P1↔P2 ≤0.2 dB)
  5. forward-eval timing benchmark
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "port_block"

for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from mpi_runner import DEFAULT_PYTHON, DEFAULT_RANKS, verify_mpi_ranks  # noqa: E402
from run_port_gates_campaign import (  # noqa: E402
    GRID_OFFSETS,
    _json_safe,
    _run_case,
)


HORNS_FORMULATIONS = (
    "num_mode_yee_guide_normal",
    "num_mode_guide_normal",
)

GRID_FORMULATIONS = (
    "num_mode_yee_guide_normal",
    "num_mode_guide_normal",
)

RT_VALUES = (5.0, 10.0, 20.0, 40.0)


def _phase_horns(res: int, run_time: float, ranks: int) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for form in HORNS_FORMULATIONS:
        for ports in ("0,1", "1,2"):
            label = f"horns_{form}_P{ports.replace(',', 'P')}_g0_0_res{res}_rt{run_time:g}"
            print(f"\n=== phase1 {label} ===", flush=True)
            row = _run_case(
                label=label,
                formulation=form,
                ports=ports,
                res=res,
                run_time=run_time,
                grid_offset=(0.0, 0.0),
                device_mode="horns_only",
                ranks=ranks,
                discrete_control=True,
                skip_existing=True,
            )
            rows.append(row)
            print(f"  port_err={row.get('port_err_dB')} dB", flush=True)
    return rows


def _phase_grid(res: int, run_time: float, ranks: int) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for form in GRID_FORMULATIONS:
        for gox, goy in GRID_OFFSETS:
            label = f"grid_{form}_P0P1_g{gox:g}_{goy:g}_res{res}_rt{run_time:g}"
            print(f"\n=== phase2 {label} ===", flush=True)
            row = _run_case(
                label=label,
                formulation=form,
                ports="0,1",
                res=res,
                run_time=run_time,
                grid_offset=(gox, goy),
                device_mode="horns_only",
                ranks=ranks,
                discrete_control=True,
                skip_existing=True,
            )
            rows.append(row)
            print(f"  port_err={row.get('port_err_dB')} dB", flush=True)
    return rows


def _phase_full_rt(
    formulation: str, res: int, ranks: int
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    prev_err: float | None = None
    for rt in RT_VALUES:
        label = f"full_{formulation}_P0P1_g0_0_res{res}_rt{rt:g}"
        print(f"\n=== phase3 {label} ===", flush=True)
        row = _run_case(
            label=label,
            formulation=formulation,
            ports="0,1",
            res=res,
            run_time=rt,
            grid_offset=(0.0, 0.0),
            device_mode="full",
            ranks=ranks,
            discrete_control=True,
            skip_existing=True,
        )
        row["run_time"] = rt
        rows.append(row)
        err = row.get("port_err_dB")
        print(f"  port_err={err} dB  discrete={row.get('discrete_amp_err_dB')}", flush=True)
        if err is not None and prev_err is not None:
            if abs(float(err) - float(prev_err)) < 0.01 and float(err) <= 0.2:
                print("  rt convergence: stable below gate — stopping sweep", flush=True)
                break
        prev_err = float(err) if err is not None else None
    return rows


def _phase_res64(
    formulation: str, run_time: float, ranks: int
) -> Dict[str, Any]:
    """res64 full-device confirmation; yee_sdotn fallback if flux subtraction fails."""
    rows: List[Dict[str, Any]] = []
    for form in (formulation, "num_mode_yee_sdotn"):
        if form != formulation and rows:
            prev = rows[-1]
            if prev.get("flux_subtraction_healthy", True):
                break
        label = f"full_{form}_P0P1_g0_0_res64_rt{run_time:g}"
        print(f"\n=== phase4 {label} ===", flush=True)
        row = _run_case(
            label=label,
            formulation=form,
            ports="0,1",
            res=64,
            run_time=run_time,
            grid_offset=(0.0, 0.0),
            device_mode="full",
            ranks=ranks,
            discrete_control=True,
            skip_existing=True,
        )
        json_path = Path(row.get("json_out", ""))
        if json_path.is_file():
            data = json.loads(json_path.read_text(encoding="utf-8"))
            fs = data.get("flux_subtraction", {})
            row["flux_subtraction_healthy"] = bool(fs.get("healthy", True))
        rows.append(row)
        if row.get("port_err_dB") is not None and float(row["port_err_dB"]) <= 0.2:
            break
    return rows[-1] if len(rows) == 1 else {"primary": rows[0], "fallback": rows[-1]}


def _phase_benchmark(formulation: str, res: int, run_time: float, ranks: int) -> Dict[str, Any]:
    """Time one full P1↔P2 reciprocity study (norm + device + discrete)."""
    label = f"bench_{formulation}_P0P1_res{res}_rt{run_time:g}"
    json_path = OUT / f"{label}.json"
    if json_path.is_file():
        data = json.loads(json_path.read_text(encoding="utf-8"))
        return data
    row = _run_case(
        label=label,
        formulation=formulation,
        ports="0,1",
        res=res,
        run_time=run_time,
        grid_offset=(0.0, 0.0),
        device_mode="full",
        ranks=ranks,
        discrete_control=True,
        skip_existing=False,
    )
    bench = {
        "label": label,
        "wall_s": row.get("elapsed_s"),
        "mpi_ranks": ranks,
        "res": res,
        "run_time": run_time,
        "formulation": formulation,
        "note": "includes norm cache miss path if first run",
    }
    (OUT / f"{label}_bench.json").write_text(
        json.dumps(_json_safe(bench), indent=2) + "\n", encoding="utf-8"
    )
    return bench


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--horns-rt", type=float, default=5.0)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--skip-res64", action="store_true", default=False)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    verify_mpi_ranks(ranks=args.ranks, python=DEFAULT_PYTHON, cwd=VAL)

    t0 = time.perf_counter()
    report: Dict[str, Any] = {"phases": {}}

    report["phases"]["horns_gates"] = _phase_horns(32, args.horns_rt, args.ranks)
    report["phases"]["grid_offset"] = _phase_grid(32, args.horns_rt, args.ranks)

    horns_p12 = next(
        (
            r
            for r in report["phases"]["horns_gates"]
            if r.get("ports") == "0,1"
        ),
        None,
    )
    best = "num_mode_guide_normal"
    if horns_p12:
        candidates = [r for r in report["phases"]["horns_gates"] if r.get("ports") == "0,1"]
        best_row = min(
            candidates,
            key=lambda r: float(r.get("port_err_dB") or 999),
        )
        best = str(best_row.get("formulation", best))
    if horns_p12 and horns_p12.get("port_err_dB") is not None:
        if float(horns_p12.get("port_err_dB", 999)) <= 0.2:
            report["phases"]["full_rt"] = _phase_full_rt(best, 32, args.ranks)
            last = report["phases"]["full_rt"][-1]
            min_rt = float(last.get("run_time", 20.0)) if last else 20.0
            if (
                not args.skip_res64
                and last.get("port_err_dB") is not None
                and float(last["port_err_dB"]) <= 0.2
            ):
                report["phases"]["res64"] = _phase_res64(best, min_rt, args.ranks)
            report["phases"]["benchmark"] = _phase_benchmark(
                best, 32, min_rt, args.ranks
            )

    report["timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    report["elapsed_s"] = time.perf_counter() - t0
    report["recommended_formulation"] = best

    out_path = OUT / "block_summary.json"
    out_path.write_text(json.dumps(_json_safe(report), indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
