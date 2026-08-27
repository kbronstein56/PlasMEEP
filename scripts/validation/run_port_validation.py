#!/usr/bin/env python3
"""
Discriminating A/B validation for six-port formulations.

Runs (by default):
  1) incident-power uniformity at res=32 and res=48 (all six ports)
  2) P1↔P2 (axis–diag) and P2↔P3 (diag–diag) at res=32, rt=40
  3) confirmation of the best candidate at res=64, rt=80 (P1↔P2 + incident)

Reuses overnight baseline numbers from outputs/validation/reciprocity/ when
comparing; does not re-run expensive baseline resolution sweeps.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs", "validation", "ports")
PY = os.environ.get(
    "PLASMEEP_PYTHON",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python",
)
MPIRUN = os.environ.get(
    "PLASMEEP_MPIRUN",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/mpirun",
)


def _rank0_print(msg: str) -> None:
    print(msg, flush=True)


def run_mpi(cmd: List[str], ranks: int, log_path: str) -> int:
    full = [MPIRUN, "-np", str(ranks)] + cmd
    _rank0_print(f"$ {' '.join(full)}")
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(f"# {' '.join(full)}\n")
        log.flush()
        proc = subprocess.run(
            full,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={
                **os.environ,
                "PYTHONPATH": f"{VAL}:{os.path.join(ROOT, 'scripts')}:"
                + os.environ.get("PYTHONPATH", ""),
            },
        )
    return int(proc.returncode)


def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def overnight_baseline_row() -> Dict[str, Any]:
    """Known baseline from completed campaign (do not re-simulate)."""
    recip = os.path.join(ROOT, "outputs", "validation", "reciprocity")
    p12_32 = load_json(os.path.join(recip, "pair_P1P2_res32_rt40.json"))
    p23_32 = load_json(os.path.join(recip, "pair_P2P3_res32_rt40.json"))
    p12_64 = load_json(os.path.join(recip, "pair_P1P2_res64_rt80.json"))
    incid = load_json(os.path.join(recip, "incident_power_vs_res.json"))
    incid_by_res = {int(r["res"]): r for r in incid["rows"]}
    return {
        "formulation": "baseline_hz_line",
        "source": "overnight_campaign",
        "P1P2_res32_rt40_dB": p12_32["reciprocity"]["max_abs_diff_dB"],
        "P2P3_res32_rt40_dB": p23_32["reciprocity"]["max_abs_diff_dB"],
        "P1P2_res64_rt80_dB": p12_64["reciprocity"]["max_abs_diff_dB"],
        "incident": {
            str(res): {
                "axis_over_diag": incid_by_res[res]["axis_over_diag"],
                "max_over_min": incid_by_res[res]["max_over_min"],
            }
            for res in (32, 48, 64, 96)
            if res in incid_by_res
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument(
        "--formulations",
        type=str,
        default="aligned_hz_line,te1_hz_line,eigenmode",
        help="Comma list of NEW formulations to run (baseline reused from overnight)",
    )
    parser.add_argument(
        "--phase",
        type=str,
        default="all",
        choices=["smoke", "cheap", "confirm", "all"],
    )
    parser.add_argument("--skip-confirm", action="store_true")
    args = parser.parse_args()

    os.makedirs(OUT, exist_ok=True)
    formulations = [f.strip() for f in args.formulations.split(",") if f.strip()]
    index: Dict[str, Any] = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "ranks": args.ranks,
        "formulations": formulations,
        "baseline_overnight": overnight_baseline_row(),
        "runs": [],
    }

    def record(label: str, path: str, rc: int, wall_s: float) -> None:
        entry: Dict[str, Any] = {
            "label": label,
            "path": path,
            "returncode": rc,
            "wall_s": wall_s,
        }
        if rc == 0 and path.endswith(".json") and os.path.isfile(path):
            try:
                entry["data"] = load_json(path)
            except Exception as exc:  # noqa: BLE001
                entry["load_error"] = str(exc)
        index["runs"].append(entry)
        with open(os.path.join(OUT, "index.json"), "w", encoding="utf-8") as f:
            json.dump(index, f, indent=2)

    # --- smoke: one-port eigenmode/aligned norm ---
    if args.phase in ("smoke", "cheap", "all"):
        for formul in formulations:
            json_out = os.path.join(OUT, f"smoke_incident_{formul}_res32.json")
            log = os.path.join(OUT, f"smoke_incident_{formul}_res32.log")
            cmd = [
                PY,
                os.path.join(VAL, "incident_power_orientation.py"),
                "--res-list",
                "32",
                "--run-time",
                "40",
                "--port-formulation",
                formul,
                "--json-out",
                json_out,
            ]
            t0 = time.time()
            rc = run_mpi(cmd, args.ranks, log)
            record(f"smoke_incident_{formul}_res32", json_out, rc, time.time() - t0)
            if rc != 0:
                _rank0_print(f"SMOKE FAILED for {formul}; see {log}")
                # Continue other formulations rather than aborting whole cycle.

    def smoke_ok(formul: str) -> bool:
        """Reject formulations with collapsed/unphysical incident power."""
        smoke = os.path.join(OUT, f"smoke_incident_{formul}_res32.json")
        if not os.path.isfile(smoke):
            return False
        try:
            row = load_json(smoke)["rows"][0]
            ratio = float(row["axis_over_diag"])
            mom = float(row["max_over_min"])
            powers = [float(v) for v in row["powers"].values()]
            if min(powers) <= 0.0:
                _rank0_print(
                    f"SMOKE REJECT {formul}: non-positive port power {powers}"
                )
                return False
            if not (0.5 <= ratio <= 1.5):
                _rank0_print(
                    f"SMOKE REJECT {formul}: axis/diag={ratio:.4f} outside [0.5,1.5]"
                )
                return False
            if mom > 5.0:
                _rank0_print(
                    f"SMOKE REJECT {formul}: max/min={mom:.3f} > 5 (orientation collapse)"
                )
                return False
            return True
        except Exception as exc:  # noqa: BLE001
            _rank0_print(f"SMOKE REJECT {formul}: parse error {exc}")
            return False

    if args.phase in ("cheap", "all"):
        for formul in formulations:
            # Skip if smoke failed hard (no json) or failed quality gate
            if not smoke_ok(formul):
                _rank0_print(f"Skipping cheap pair tests for {formul}")
                continue

            # Incident at res=48 as second resolution
            json_out = os.path.join(OUT, f"incident_{formul}_res32_48.json")
            log = os.path.join(OUT, f"incident_{formul}_res32_48.log")
            cmd = [
                PY,
                os.path.join(VAL, "incident_power_orientation.py"),
                "--res-list",
                "32,48",
                "--run-time",
                "40",
                "--port-formulation",
                formul,
                "--json-out",
                json_out,
            ]
            t0 = time.time()
            rc = run_mpi(cmd, args.ranks, log)
            record(f"incident_{formul}_res32_48", json_out, rc, time.time() - t0)

            for ports, tag in (("0,1", "P1P2"), ("1,2", "P2P3")):
                json_out = os.path.join(
                    OUT, f"pair_{tag}_{formul}_res32_rt40.json"
                )
                log = os.path.join(OUT, f"pair_{tag}_{formul}_res32_rt40.log")
                cmd = [
                    PY,
                    os.path.join(VAL, "reciprocity_b0_study.py"),
                    "--res",
                    "32",
                    "--run-time",
                    "40",
                    "--ports",
                    ports,
                    "--port-formulation",
                    formul,
                    "--label",
                    f"{tag}_{formul}_res32_rt40",
                    "--json-out",
                    json_out,
                    "--mpi-note",
                    f"np={args.ranks}",
                    "--skip-norm-if-cached",
                ]
                t0 = time.time()
                rc = run_mpi(cmd, args.ranks, log)
                record(
                    f"pair_{tag}_{formul}_res32_rt40",
                    json_out,
                    rc,
                    time.time() - t0,
                )

    # Choose best candidate by (1) lower P1P2 residual, (2) closer axis/diag to 1
    best: Optional[str] = None
    best_score = None
    for formul in formulations:
        p12 = os.path.join(OUT, f"pair_P1P2_{formul}_res32_rt40.json")
        inc = os.path.join(OUT, f"smoke_incident_{formul}_res32.json")
        if not (os.path.isfile(p12) and os.path.isfile(inc)):
            continue
        d12 = load_json(p12)["reciprocity"]["max_abs_diff_dB"]
        ratio = load_json(inc)["rows"][0]["axis_over_diag"]
        # score: residual + 2*|log10(ratio)|*10 as soft penalty for incident bias
        score = float(d12) + 20.0 * abs(__import__("math").log10(ratio))
        if best_score is None or score < best_score:
            best_score = score
            best = formul
    index["best_cheap"] = {"formulation": best, "score": best_score}

    if (
        args.phase in ("confirm", "all")
        and not args.skip_confirm
        and best is not None
    ):
        formul = best
        json_out = os.path.join(OUT, f"incident_{formul}_res64.json")
        log = os.path.join(OUT, f"incident_{formul}_res64.log")
        cmd = [
            PY,
            os.path.join(VAL, "incident_power_orientation.py"),
            "--res-list",
            "64",
            "--run-time",
            "80",
            "--port-formulation",
            formul,
            "--json-out",
            json_out,
        ]
        t0 = time.time()
        rc = run_mpi(cmd, args.ranks, log)
        record(f"incident_{formul}_res64", json_out, rc, time.time() - t0)

        json_out = os.path.join(OUT, f"pair_P1P2_{formul}_res64_rt80.json")
        log = os.path.join(OUT, f"pair_P1P2_{formul}_res64_rt80.log")
        cmd = [
            PY,
            os.path.join(VAL, "reciprocity_b0_study.py"),
            "--res",
            "64",
            "--run-time",
            "80",
            "--ports",
            "0,1",
            "--port-formulation",
            formul,
            "--label",
            f"P1P2_{formul}_res64_rt80",
            "--json-out",
            json_out,
            "--mpi-note",
            f"np={args.ranks}",
            "--skip-norm-if-cached",
        ]
        t0 = time.time()
        rc = run_mpi(cmd, args.ranks, log)
        record(f"pair_P1P2_{formul}_res64_rt80", json_out, rc, time.time() - t0)

        # Also confirm at res=48 for resolution trend of the candidate
        json_out = os.path.join(OUT, f"pair_P1P2_{formul}_res48_rt80.json")
        log = os.path.join(OUT, f"pair_P1P2_{formul}_res48_rt80.log")
        cmd = [
            PY,
            os.path.join(VAL, "reciprocity_b0_study.py"),
            "--res",
            "48",
            "--run-time",
            "80",
            "--ports",
            "0,1",
            "--port-formulation",
            formul,
            "--label",
            f"P1P2_{formul}_res48_rt80",
            "--json-out",
            json_out,
            "--mpi-note",
            f"np={args.ranks}",
            "--skip-norm-if-cached",
        ]
        t0 = time.time()
        rc = run_mpi(cmd, args.ranks, log)
        record(f"pair_P1P2_{formul}_res48_rt80", json_out, rc, time.time() - t0)

    index["finished_utc"] = datetime.now(timezone.utc).isoformat()
    with open(os.path.join(OUT, "index.json"), "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2)
    _rank0_print(f"Wrote {os.path.join(OUT, 'index.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
