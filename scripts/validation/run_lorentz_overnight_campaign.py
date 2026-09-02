#!/usr/bin/env python3
"""
Overnight direct-Lorentz reciprocity diagnostic campaign.

Phases:
  1. Test specification audit (embedded in each JSON)
  2. Runtime convergence (full PMM, res32)
  3. Material-complexity ladder (res32, shortest validated rt)
  4. Optional res64 confirmation (one case, if ladder identifies failure)

All Meep runs use mpi_runner (np=32).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "lorentz_overnight"
SCRIPT = VAL / "direct_reciprocity_test.py"

for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from mpi_runner import DEFAULT_PYTHON, DEFAULT_RANKS, run_mpi_script, verify_mpi_ranks  # noqa: E402


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _run_case(
    *,
    label: str,
    res: int,
    run_time: float,
    device_mode: str = "full",
    plasma_fill: str = "active",
    n_bulbs: Optional[int] = None,
    susceptibility_mode: str = "auto",
    ranks: int = DEFAULT_RANKS,
    skip_existing: bool = True,
) -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    json_out = OUT / f"{label}.json"
    log_out = OUT / f"{label}.log"
    if skip_existing and json_out.is_file():
        data = json.loads(json_out.read_text(encoding="utf-8"))
        pair = data.get("pairs", [{}])[0]
        return {
            "label": label,
            "skipped": True,
            "json_out": str(json_out),
            "amp_err_dB": pair.get("amp_err_dB"),
            "abs_a_to_b": pair.get("abs_a_to_b"),
            "abs_b_to_a": pair.get("abs_b_to_a"),
            "mpi_ranks": data.get("mpi", {}).get("ranks"),
        }

    args = [
        str(SCRIPT),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--pairs",
        "0,1",
        "--device-mode",
        device_mode,
        "--plasma-fill",
        plasma_fill,
        "--susceptibility-mode",
        susceptibility_mode,
        "--label",
        label,
        "--json-out",
        str(json_out),
    ]
    if n_bulbs is not None:
        args.extend(["--n-bulbs", str(n_bulbs)])

    t0 = time.perf_counter()
    proc = run_mpi_script(
        args,
        ranks=ranks,
        python=DEFAULT_PYTHON,
        cwd=ROOT,
        verify_ranks=True,
        log_path=log_out,
    )
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        return {"label": label, "error": proc.returncode, "wall_s": wall, "log": str(log_out)}
    if not json_out.is_file():
        return {"label": label, "error": "missing_json", "wall_s": wall}
    data = json.loads(json_out.read_text(encoding="utf-8"))
    pair = data.get("pairs", [{}])[0]
    return {
        "label": label,
        "json_out": str(json_out),
        "run_time": run_time,
        "res": res,
        "device_mode": device_mode,
        "plasma_fill": plasma_fill,
        "n_bulbs": n_bulbs,
        "susceptibility_mode": susceptibility_mode,
        "abs_a_to_b": pair.get("abs_a_to_b"),
        "abs_b_to_a": pair.get("abs_b_to_a"),
        "amp_err_dB": pair.get("amp_err_dB"),
        "phase_diff_deg": pair.get("phase_diff_deg"),
        "complex_sym_err": pair.get("complex_sym_err"),
        "fields_large_enough": pair.get("fields_large_enough"),
        "mpi_ranks": data.get("mpi", {}).get("ranks"),
        "wall_s": wall,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--phase", type=str, default="all", help="all|rt|ladder|res64")
    args = parser.parse_args()

    print(f"Verifying MPI np={args.ranks} ...", flush=True)
    verify_mpi_ranks(ranks=args.ranks, python=DEFAULT_PYTHON, cwd=ROOT)

    summary: Dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "mpi_ranks": args.ranks,
        "phases": {},
    }
    ladder_rows: List[Dict[str, Any]] = []

    # Phase 2: runtime convergence (full PMM, res32)
    if args.phase in ("all", "rt"):
        print("\n=== PHASE 2: direct Lorentz runtime convergence (res32, full) ===", flush=True)
        rt_results = []
        for rt in (5.0, 10.0, 20.0, 40.0):
            row = _run_case(
                label=f"rt_conv_full_res32_rt{rt:g}",
                res=32,
                run_time=rt,
                device_mode="full",
                plasma_fill="active",
                ranks=args.ranks,
                skip_existing=args.skip_existing,
            )
            rt_results.append(row)
            print(
                f"  rt={rt:g}: amp_err={row.get('amp_err_dB')} dB  "
                f"|H|={row.get('abs_a_to_b')},{row.get('abs_b_to_a')}  "
                f"wall={row.get('wall_s', 0):.0f}s",
                flush=True,
            )
        summary["phases"]["runtime_convergence_res32"] = rt_results

    # Phase 3: material ladder at rt=20 (validated shortest for full at res32)
    if args.phase in ("all", "ladder"):
        print("\n=== PHASE 3: material-complexity ladder (res32, rt=20) ===", flush=True)
        ladder_cases = [
            ("A_horns_only", dict(device_mode="horns_only")),
            ("B_geometry_91bulbs", dict(plasma_fill="geometry_only", n_bulbs=None)),
            ("C_dielectric_91bulbs", dict(plasma_fill="dielectric_fill", n_bulbs=None)),
            ("D_drude_active_1bulb", dict(plasma_fill="active", n_bulbs=1)),
            ("D_drude_active_7bulbs", dict(plasma_fill="active", n_bulbs=7)),
            ("D_drude_active_91bulbs", dict(plasma_fill="active", n_bulbs=None)),
            ("E_gyrotropic_b0_1bulb", dict(plasma_fill="active", n_bulbs=1, susceptibility_mode="force_gyrotropic_b0")),
            ("E_gyrotropic_b0_91bulbs", dict(plasma_fill="active", n_bulbs=None, susceptibility_mode="force_gyrotropic_b0")),
        ]
        for label, kwargs in ladder_cases:
            dm = kwargs.pop("device_mode", "full")
            row = _run_case(
                label=f"ladder_{label}_res32_rt20",
                res=32,
                run_time=20.0,
                device_mode=dm,
                ranks=args.ranks,
                skip_existing=args.skip_existing,
                **kwargs,
            )
            ladder_rows.append(row)
            print(
                f"  {label}: amp_err={row.get('amp_err_dB')} dB  "
                f"sym_err={row.get('complex_sym_err')}",
                flush=True,
            )
        summary["phases"]["material_ladder_res32"] = ladder_rows

    # Phase 9: conditional res64 (one confirmation of full active plasma)
    if args.phase in ("all", "res64"):
        print("\n=== PHASE 9: res64 confirmation (full active, rt=20) ===", flush=True)
        row = _run_case(
            label="confirm_full_active_res64_rt20",
            res=64,
            run_time=20.0,
            device_mode="full",
            plasma_fill="active",
            ranks=args.ranks,
            skip_existing=args.skip_existing,
        )
        summary["phases"]["res64_confirmation"] = row
        ladder_rows.append(row)

    summary["comparison_table"] = ladder_rows
    out_path = OUT / "overnight_summary.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(_json_safe(summary), f, indent=2)
        f.write("\n")
    print(f"\nWrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
