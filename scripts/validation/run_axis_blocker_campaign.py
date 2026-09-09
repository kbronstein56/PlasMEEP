#!/usr/bin/env python3
"""
Axis-blocker campaign: P1/P2 direction split, mode reextract, grid-index ports.

Runs sequentially (one np=32 job at a time). Skips existing JSON artifacts.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

VAL = Path(__file__).resolve().parent
ROOT = VAL.parents[1]
OUT = ROOT / "outputs" / "validation" / "port_axis"

for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from mpi_runner import DEFAULT_PYTHON, DEFAULT_RANKS, run_mpi_script, verify_mpi_ranks  # noqa: E402
from run_port_gates_campaign import GRID_OFFSETS, _json_safe, _run_case  # noqa: E402


def _run_py(script: str, args: List[str], *, mpi: bool, log: Path) -> int:
    if mpi:
        proc = run_mpi_script([script, *args], ranks=DEFAULT_RANKS, log_path=log, cwd=VAL)
        return int(proc.returncode)
    proc = subprocess.run(
        [DEFAULT_PYTHON, str(VAL / script), *args],
        cwd=VAL,
        capture_output=False,
    )
    return int(proc.returncode)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--phase", type=str, default="all", help="audit|diag|reextract|grid|all")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    verify_mpi_ranks(ranks=DEFAULT_RANKS, python=DEFAULT_PYTHON, cwd=VAL)
    t0 = time.perf_counter()
    report: Dict[str, Any] = {"phases": {}}

    if args.phase in ("audit", "all"):
        audit_out = OUT / "port_coord_audit_res32.json"
        if not (args.skip_existing and audit_out.is_file()):
            rc = _run_py(
                "audit_port_coordinates.py",
                ["--res", "32", "--json-out", str(audit_out)],
                mpi=False,
                log=OUT / "audit_port_coordinates.log",
            )
            if rc != 0:
                raise RuntimeError(f"audit_port_coordinates failed rc={rc}")
        report["phases"]["coord_audit"] = str(audit_out)

    if args.phase in ("reextract", "all"):
        for res in (32, 64):
            out = OUT / f"mode_reextract_res{res}.json"
            if args.skip_existing and out.is_file():
                report["phases"][f"reextract_res{res}"] = {"skipped": True, "path": str(out)}
                continue
            rc = _run_py(
                "compare_mode_reextract.py",
                [
                    "--res",
                    str(res),
                    "--run-time",
                    "40" if res >= 64 else "20",
                    "--ports",
                    "0,1",
                    "--json-out",
                    str(out),
                ],
                mpi=True,
                log=OUT / f"mode_reextract_res{res}.log",
            )
            if rc != 0:
                raise RuntimeError(f"mode reextract res{res} failed rc={rc}")
            report["phases"][f"reextract_res{res}"] = {"path": str(out)}

    diag_cases = [
        ("horns_res32", "horns_only", 32, 5.0),
        ("full_res32", "full", 32, 20.0),
        ("horns_res64", "horns_only", 64, 5.0),
        ("full_res64", "full", 64, 20.0),
    ]
    if args.phase in ("diag", "all"):
        rows: List[Dict[str, Any]] = []
        for tag, dmode, res, rt in diag_cases:
            label = f"axis_diag_{tag}_guide_normal_rt{rt:g}"
            out = OUT / f"{label}.json"
            if args.skip_existing and out.is_file():
                rows.append({"label": label, "skipped": True, "json": str(out)})
                continue
            rc = _run_py(
                "port_axis_diagnostic.py",
                [
                    "--res",
                    str(res),
                    "--run-time",
                    str(rt),
                    "--device-mode",
                    dmode,
                    "--formulation",
                    "num_mode_guide_normal",
                    "--label",
                    label,
                    "--json-out",
                    str(out),
                ],
                mpi=True,
                log=OUT / f"{label}.log",
            )
            if rc != 0:
                raise RuntimeError(f"axis diag {label} failed rc={rc}")
            data = json.loads(out.read_text(encoding="utf-8"))
            rows.append(
                {
                    "label": label,
                    "flux_recip_dB": data["reciprocity"]["flux_recip_err_dB"],
                    "modal_recip_dB": data["reciprocity"]["modal_recip_err_dB"],
                    "json": str(out),
                }
            )
        report["phases"]["direction_diag"] = rows

    if args.phase in ("grid", "all"):
        grid_rows: List[Dict[str, Any]] = []
        for form in (
            "num_mode_guide_normal",
            "num_mode_yee_guide_normal",
            "num_mode_grid_index",
        ):
            for gox, goy in GRID_OFFSETS:
                label = f"grid_{form}_P0P1_g{gox:g}_{goy:g}_res32_rt5"
                row = _run_case(
                    label=label,
                    formulation=form,
                    ports="0,1",
                    res=32,
                    run_time=5.0,
                    grid_offset=(gox, goy),
                    device_mode="horns_only",
                    ranks=DEFAULT_RANKS,
                    discrete_control=True,
                    skip_existing=args.skip_existing,
                )
                grid_rows.append(row)
        report["phases"]["grid_index_offsets"] = grid_rows

    report["timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    report["elapsed_s"] = time.perf_counter() - t0
    summary = OUT / "axis_blocker_summary.json"
    summary.write_text(json.dumps(_json_safe(report), indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
