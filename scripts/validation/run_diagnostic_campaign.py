#!/usr/bin/env python3
"""
Cheap discriminating diagnostics for B=0 reciprocity (Phase 4).

Runs only cases not already present in outputs/validation/diagnostics/.
Reuses reciprocity_b0_study.py via MPI — does not duplicate simulation logic.
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

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "scripts" / "validation"
OUT = ROOT / "outputs" / "validation" / "diagnostics"
PY = os.environ.get(
    "PLASMEEP_PYTHON",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python",
)
MPIRUN = os.environ.get(
    "PLASMEEP_MPIRUN",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/mpirun",
)

# Cases: (label_suffix, formulation, device_mode, ports, res, rt)
CASES = [
    ("horns_te1hz_P1P2", "te1_hz_line", "horns_only", "0,1", 32, 40),
    ("horns_te1hz_P2P3", "te1_hz_line", "horns_only", "1,2", 32, 40),
    ("horns_te1gn_P1P2", "te1_guide_normal", "horns_only", "0,1", 32, 40),
    ("horns_te1ez_P1P2", "te1_ez_line", "horns_only", "0,1", 32, 40),
    ("horns_te1ez_P2P3", "te1_ez_line", "horns_only", "1,2", 32, 40),
    ("full_te1ez_P1P2", "te1_ez_line", "full", "0,1", 32, 40),
    ("full_te1ez_P2P3", "te1_ez_line", "full", "1,2", 32, 40),
]


def _mpi_env() -> dict:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("FI_PROVIDER", "tcp")
    env.setdefault("MPICH_CH4_NETMOD", "ofi")
    env.setdefault("UCX_TLS", "tcp,self")
    env["PYTHONPATH"] = f"{VAL}:{ROOT / 'scripts'}:" + env.get("PYTHONPATH", "")
    return env


def run_case(
    label: str,
    formulation: str,
    device_mode: str,
    ports: str,
    res: int,
    rt: float,
    ranks: int,
) -> Dict[str, Any]:
    json_out = OUT / f"pair_{label}_res{res}_rt{int(rt)}.json"
    log_out = OUT / f"pair_{label}_res{res}_rt{int(rt)}.log"
    if json_out.is_file():
        with open(json_out, encoding="utf-8") as f:
            data = json.load(f)
        return {
            "label": label,
            "skipped": True,
            "path": str(json_out),
            "P12_dB": data.get("reciprocity", {}).get("max_abs_diff_dB"),
        }

    OUT.mkdir(parents=True, exist_ok=True)
    cmd = [
        MPIRUN,
        "-np",
        str(ranks),
        PY,
        str(VAL / "reciprocity_b0_study.py"),
        "--res",
        str(res),
        "--run-time",
        str(rt),
        "--ports",
        ports,
        "--port-formulation",
        formulation,
        "--device-mode",
        device_mode,
        "--label",
        label,
        "--json-out",
        str(json_out),
        "--mpi-note",
        f"np={ranks}",
        "--skip-norm-if-cached",
    ]
    t0 = time.time()
    with open(log_out, "w", encoding="utf-8") as log:
        log.write(f"# {' '.join(cmd)}\n")
        log.flush()
        proc = subprocess.run(
            cmd, cwd=str(ROOT), env=_mpi_env(), stdout=log, stderr=subprocess.STDOUT
        )
    wall = time.time() - t0
    row: Dict[str, Any] = {
        "label": label,
        "formulation": formulation,
        "device_mode": device_mode,
        "ports": ports,
        "res": res,
        "run_time": rt,
        "ranks": ranks,
        "wall_s": wall,
        "returncode": proc.returncode,
        "path": str(json_out),
        "log": str(log_out),
    }
    if json_out.is_file():
        with open(json_out, encoding="utf-8") as f:
            data = json.load(f)
        row["P12_dB"] = data.get("reciprocity", {}).get("max_abs_diff_dB")
        row["incident_rel"] = (data.get("incident_power_mismatch") or {}).get(
            "rel_spread"
        )
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--cases", type=str, default="", help="Comma filter on labels")
    args = parser.parse_args()

    selected = CASES
    if args.cases.strip():
        want = {x.strip() for x in args.cases.split(",") if x.strip()}
        selected = [c for c in CASES if c[0] in want]

    index: Dict[str, Any] = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "ranks": args.ranks,
        "runs": [],
    }
    for label, formul, mode, ports, res, rt in selected:
        print(f"\n=== {label} ===", flush=True)
        row = run_case(label, formul, mode, ports, res, rt, args.ranks)
        index["runs"].append(row)
        print(
            f"  rc={row.get('returncode')} skipped={row.get('skipped')} "
            f"P12={row.get('P12_dB')} wall={row.get('wall_s', 0):.1f}s",
            flush=True,
        )
        with open(OUT / "index.json", "w", encoding="utf-8") as f:
            json.dump(index, f, indent=2)

    # Refresh master registry
    subprocess.run([PY, str(VAL / "validation_registry.py")], cwd=str(ROOT), check=False)
    print(f"\nWrote {OUT / 'index.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
