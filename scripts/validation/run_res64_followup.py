#!/usr/bin/env python3
"""
Follow-up tests after res64 full-device axis failure diagnosis.

Runs sequentially (one np=32 job at a time):
  1. full-device res64 P2↔P3 (diagonal pair control)
  2. full-device res32 num_mode_yee_sdotn P1↔P2 baseline
  3. horns_only grid-offset sweep for num_mode_yee_sdotn (res32, rt=5)
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
from run_port_gates_campaign import GRID_OFFSETS, _json_safe, _run_case  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--no-skip-existing", dest="skip_existing", action="store_false")
    args = parser.parse_args()

    verify_mpi_ranks(ranks=args.ranks, python=DEFAULT_PYTHON, cwd=VAL)
    OUT.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    rows: List[Dict[str, Any]] = []

    cases = [
        dict(
            label="full_num_mode_guide_normal_P1P2_g0_0_res64_rt20",
            formulation="num_mode_guide_normal",
            ports="1,2",
            res=64,
            run_time=20.0,
            device_mode="full",
            grid_offset=(0.0, 0.0),
        ),
        dict(
            label="full_num_mode_yee_sdotn_P0P1_g0_0_res32_rt20",
            formulation="num_mode_yee_sdotn",
            ports="0,1",
            res=32,
            run_time=20.0,
            device_mode="full",
            grid_offset=(0.0, 0.0),
        ),
    ]
    for gox, goy in GRID_OFFSETS:
        cases.append(
            dict(
                label=f"grid_num_mode_yee_sdotn_P0P1_g{gox:g}_{goy:g}_res32_rt5",
                formulation="num_mode_yee_sdotn",
                ports="0,1",
                res=32,
                run_time=5.0,
                device_mode="horns_only",
                grid_offset=(gox, goy),
            )
        )

    for spec in cases:
        print(f"\n=== {spec['label']} ===", flush=True)
        row = _run_case(
            label=spec["label"],
            formulation=spec["formulation"],
            ports=spec["ports"],
            res=spec["res"],
            run_time=spec["run_time"],
            grid_offset=spec["grid_offset"],
            device_mode=spec["device_mode"],
            ranks=args.ranks,
            discrete_control=True,
            skip_existing=args.skip_existing,
        )
        rows.append(row)
        print(
            f"  port_err={row.get('port_err_dB')} dB  "
            f"discrete={row.get('discrete_amp_err_dB')} dB",
            flush=True,
        )

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_s": time.perf_counter() - t0,
        "cases": rows,
    }
    out_path = OUT / "res64_followup_summary.json"
    out_path.write_text(json.dumps(_json_safe(summary), indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
