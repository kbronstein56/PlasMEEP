#!/usr/bin/env python3
"""Grid-offset robustness for selected port formulations."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "port_gates"

for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from run_port_gates_campaign import GRID_OFFSETS, _run_case  # noqa: E402
from mpi_runner import DEFAULT_PYTHON, DEFAULT_RANKS, verify_mpi_ranks  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=5.0)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument(
        "--formulations",
        type=str,
        default="num_mode_guide_normal,num_mode_hz_line",
    )
    parser.add_argument("--ports", type=str, default="0,1")
    args = parser.parse_args()

    verify_mpi_ranks(ranks=args.ranks, python=DEFAULT_PYTHON, cwd=VAL)
    forms = [f.strip() for f in args.formulations.split(",") if f.strip()]
    rows: List[Dict[str, Any]] = []
    for form in forms:
        for gox, goy in GRID_OFFSETS:
            label = f"grid_{form}_P{args.ports.replace(',', 'P')}_g{gox:g}_{goy:g}_res{args.res}_rt{args.run_time:g}"
            print(f"\n=== {label} ===", flush=True)
            row = _run_case(
                label=label,
                formulation=form,
                ports=args.ports,
                res=args.res,
                run_time=args.run_time,
                grid_offset=(gox, goy),
                device_mode="horns_only",
                ranks=args.ranks,
                discrete_control=True,
                skip_existing=True,
            )
            rows.append(row)
            print(f"  port_err={row.get('port_err_dB')} dB", flush=True)

    out = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": vars(args),
        "rows": rows,
    }
    path = OUT / f"grid_offset_P{args.ports.replace(',', 'P')}_res{args.res}.json"
    path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
