#!/usr/bin/env python3
"""
High-resolution B=0 port reciprocity campaign (>= 50 points/cm).

Runs full 91-bulb device with the best current numerical-mode formulation
(num_mode_guide_normal) plus matched discrete Yee-grid reciprocity control.
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
OUT = ROOT / "outputs" / "validation" / "hires_ports"
MODE_STUDY = VAL / "mode_profile_study.py"
STUDY = VAL / "reciprocity_b0_study.py"

for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from mpi_runner import DEFAULT_RANKS, run_mpi_script, verify_mpi_ranks  # noqa: E402
from physical_units import (  # noqa: E402
    meep_resolution_from_points_per_cm,
    physical_resolution_report,
)
from sixport_common import a  # noqa: E402

PAIRS = (
    ("P1P2", "0,1"),
    ("P2P3", "1,2"),
)
MODE_PORTS = "0,1,2"


def _ensure_mode_cache(
    *,
    points_per_cm: float,
    run_time: float,
    ranks: int,
) -> None:
    res = meep_resolution_from_points_per_cm(points_per_cm, a_m=a)
    cache_dir = ROOT / "outputs" / "validation" / "mode_profiles" / f"res{res}"
    port_ids = [int(x) for x in MODE_PORTS.split(",")]
    needed = [cache_dir / f"numerical_mode_P{i + 1}.json" for i in port_ids]
    if all(p.is_file() for p in needed):
        print(f"  mode cache res{res}: ports {MODE_PORTS} present")
        return
    log_out = OUT / f"mode_profiles_ppc{points_per_cm:g}_res{res}.log"
    args = [
        str(MODE_STUDY),
        "--points-per-cm",
        str(points_per_cm),
        "--run-time",
        str(run_time),
        "--ports",
        MODE_PORTS,
    ]
    proc = run_mpi_script(args, ranks=ranks, log_path=log_out)
    if proc.returncode != 0:
        raise RuntimeError(
            f"mode extraction failed ppc={points_per_cm} res={res} log={log_out}"
        )


def _pair_err(data: Dict[str, Any], pa: int, pb: int) -> Optional[float]:
    for pair in data.get("reciprocity", {}).get("all_pairs", []):
        if pair.get("i") == pa and pair.get("j") == pb:
            return pair.get("abs_diff_dB")
    return None


def _run_case(
    *,
    label: str,
    points_per_cm: float,
    run_time: float,
    ports: str,
    ranks: int,
    skip_existing: bool,
) -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    json_out = OUT / f"{label}.json"
    log_out = OUT / f"{label}.log"
    if skip_existing and json_out.is_file():
        data = json.loads(json_out.read_text(encoding="utf-8"))
        pa, pb = [int(x) for x in ports.split(",")]
        return {
            "label": label,
            "skipped": True,
            "json_out": str(json_out),
            "points_per_cm": points_per_cm,
            "port_err_dB": _pair_err(data, pa, pb),
            "discrete_amp_err_dB": (data.get("discrete_control") or {}).get(
                "amp_err_dB"
            ),
            "wall_s": data.get("timings_s", {}).get("total"),
        }

    args = [
        str(STUDY),
        "--points-per-cm",
        str(points_per_cm),
        "--run-time",
        str(run_time),
        "--ports",
        ports,
        "--port-formulation",
        "num_mode_guide_normal",
        "--device-mode",
        "full",
        "--label",
        label,
        "--json-out",
        str(json_out),
        "--discrete-control",
    ]

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
        "log_out": str(log_out),
        "points_per_cm": points_per_cm,
        "meep_resolution": data.get("settings", {}).get("res"),
        "port_err_dB": _pair_err(data, pa, pb),
        "discrete_amp_err_dB": (data.get("discrete_control") or {}).get("amp_err_dB"),
        "flux_subtraction_healthy": data.get("flux_subtraction", {}).get("healthy"),
        "incident_power_mismatch": data.get("incident_power_mismatch"),
        "wall_s": elapsed,
        "timings_s": data.get("timings_s"),
    }


def _estimate_100ppc(cases: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Rough wall-time / cell-count scaling from completed lower-ppc runs."""
    refs = [c for c in cases if not c.get("skipped") and c.get("wall_s")]
    if not refs:
        return {"feasible": False, "reason": "no completed timing references"}
    # Use highest ppc reference available.
    ref = max(refs, key=lambda c: float(c["points_per_cm"]))
    ppc_ref = float(ref["points_per_cm"])
    wall_ref = float(ref["wall_s"])
    res_ref = int(ref.get("meep_resolution") or 0)
    res_100 = meep_resolution_from_points_per_cm(100.0, a_m=a)
    if res_ref <= 0:
        return {"feasible": False, "reason": "missing meep resolution in reference"}
    scale = (100.0 / ppc_ref) ** 2
    est_wall = wall_ref * scale
    est_cells = (res_100 / res_ref) ** 2
    return {
        "feasible": True,
        "reference_label": ref["label"],
        "reference_ppc": ppc_ref,
        "reference_wall_s": wall_ref,
        "estimated_100ppc_wall_s": est_wall,
        "estimated_100ppc_cell_scale": est_cells,
        "meep_res_100": res_100,
        "recommend_run_100ppc": est_wall < 6 * 3600,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="High-res port reciprocity campaign.")
    parser.add_argument(
        "--points-per-cm-list",
        type=str,
        default="50,75",
        help="Comma-separated physical grid densities to run sequentially.",
    )
    parser.add_argument("--run-time", type=float, default=20.0)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument(
        "--include-100",
        action="store_true",
        help="After 50/75 ppc, run 100 ppc if timing estimate is reasonable.",
    )
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    verify_mpi_ranks(ranks=args.ranks)
    ppc_list = [float(x.strip()) for x in args.points_per_cm_list.split(",") if x.strip()]

    results: List[Dict[str, Any]] = []
    for ppc in ppc_list:
        res_report = physical_resolution_report(points_per_cm=ppc, a_m=a)
        print(
            f"\n=== ppc={ppc:g}  res={int(res_report['meep_resolution'])}  "
            f"dx_mm={res_report['dx_mm']:.4f}  a={res_report['a_m']} m ==="
        )
        _ensure_mode_cache(
            points_per_cm=ppc,
            run_time=max(args.run_time, 20.0),
            ranks=args.ranks,
        )
        for pair_label, ports in PAIRS:
            label = f"full_num_mode_guide_normal_{pair_label}_ppc{ppc:g}_rt{args.run_time:g}"
            row = _run_case(
                label=label,
                points_per_cm=ppc,
                run_time=args.run_time,
                ports=ports,
                ranks=args.ranks,
                skip_existing=args.skip_existing,
            )
            print(
                f"  {pair_label}: port_err={row.get('port_err_dB')} dB  "
                f"discrete={row.get('discrete_amp_err_dB')} dB  "
                f"wall={row.get('wall_s', 0):.1f}s"
            )
            results.append(row)

    estimate = _estimate_100ppc(results)
    if args.include_100 and estimate.get("recommend_run_100ppc"):
        ppc = 100.0
        print(
            f"\n=== estimated 100 ppc wall ~{estimate['estimated_100ppc_wall_s']:.0f}s; running ==="
        )
        for pair_label, ports in PAIRS:
            label = f"full_num_mode_guide_normal_{pair_label}_ppc100_rt{args.run_time:g}"
            row = _run_case(
                label=label,
                points_per_cm=ppc,
                run_time=args.run_time,
                ports=ports,
                ranks=args.ranks,
                skip_existing=args.skip_existing,
            )
            print(
                f"  {pair_label}: port_err={row.get('port_err_dB')} dB  "
                f"discrete={row.get('discrete_amp_err_dB')} dB  "
                f"wall={row.get('wall_s', 0):.1f}s"
            )
            results.append(row)
    elif args.include_100:
        print(
            f"\nSkipping 100 ppc: estimate={estimate.get('estimated_100ppc_wall_s')} s "
            f"({estimate.get('reason', 'over budget')})"
        )

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "a_m": a,
        "a_m_previous": 0.028,
        "run_time": args.run_time,
        "formulation": "num_mode_guide_normal",
        "device_mode": "full",
        "cases": results,
        "estimate_100ppc": estimate,
    }
    summary_path = OUT / "hires_port_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {summary_path}")


if __name__ == "__main__":
    main()
