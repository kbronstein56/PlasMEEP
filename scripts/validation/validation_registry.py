#!/usr/bin/env python3
"""
Scan validation JSON artifacts and emit a consolidated registry.

Does not re-run simulations. Safe to call after any campaign completes.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "validation"

SCRIPT_ROLES = {
    "sixport_common.py": "shared geometry, device build, simulate, norm cache",
    "port_formulations.py": "source/measurement formulation registry",
    "reciprocity_b0_study.py": "B=0 reciprocity pair/subset runs → JSON",
    "incident_power_orientation.py": "six-port incident uniformity",
    "run_port_validation.py": "A/B formulation orchestrator",
    "run_reciprocity_campaign.py": "overnight baseline reciprocity batch",
    "mpi_scale_oneport.py": "MPI scaling + throughput benchmarks",
    "faraday_benchmark.py": "isolated gyrotropy / Faraday rotation",
    "make_advisor_figures.py": "plots from existing JSON only",
    "validation_registry.py": "this consolidation tool",
    "run_diagnostic_campaign.py": "cheap orientation / device-mode diagnostics",
}


def _load_json(path: Path) -> Dict[str, Any] | None:
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def _classify_reciprocity(path: Path, data: Dict[str, Any]) -> Dict[str, Any]:
    settings = data.get("settings", {})
    recip = data.get("reciprocity", {})
    pr = data.get("physical_resolution", {})
    inc = data.get("incident_power_mismatch") or {}
    if not inc and "incident_power_by_port" in data:
        vals = [float(v) for v in data["incident_power_by_port"].values()]
        if vals:
            mean = sum(vals) / len(vals)
            inc = {"rel_spread": (max(vals) - min(vals)) / mean if mean else None}
    return {
        "path": str(path.relative_to(ROOT)),
        "label": data.get("label"),
        "formulation": settings.get("port_formulation"),
        "device_mode": settings.get("device_mode", "full"),
        "res": settings.get("res"),
        "run_time": settings.get("run_time"),
        "ports": settings.get("ports"),
        "P12_dB": recip.get("max_abs_diff_dB"),
        "incident_rel_spread": inc.get("rel_spread"),
        "dx_mm": pr.get("grid_spacing_mm"),
        "px_per_cm": pr.get("pixels_per_cm"),
        "mpi_note": data.get("mpi_note"),
        "total_s": (data.get("timings_s") or {}).get("total"),
        "timestamp_utc": data.get("timestamp_utc"),
    }


def build_registry() -> Dict[str, Any]:
    reciprocity_rows: List[Dict[str, Any]] = []
    mpi_rows: List[Dict[str, Any]] = []
    faraday_rows: List[Dict[str, Any]] = []

    for path in sorted(OUT.rglob("*.json")):
        if path.name in ("index.json", "MASTER_VALIDATION_REGISTRY.json"):
            continue
        data = _load_json(path)
        if not data or not isinstance(data, dict):
            continue
        if "reciprocity" in data and "settings" in data:
            reciprocity_rows.append(_classify_reciprocity(path, data))
        if "results" in data and isinstance(data["results"], list):
            if any("t_device_s" in r for r in data["results"]):
                mpi_rows.append(
                    {
                        "path": str(path.relative_to(ROOT)),
                        "res": data.get("res"),
                        "run_time": data.get("run_time"),
                        "formulation": data.get("formulation"),
                        "results": data["results"],
                    }
                )
        if "B=0" in str(data.get("status", "")) or "kappa_rad_per_a" in str(data):
            if data.get("resolution") or data.get("res"):
                faraday_rows.append(
                    {
                        "path": str(path.relative_to(ROOT)),
                        "res": data.get("res") or data.get("resolution"),
                        "status": data.get("status"),
                    }
                )

    # De-duplicate reciprocity by label keeping newest timestamp
    by_label: Dict[str, Dict[str, Any]] = {}
    for row in reciprocity_rows:
        key = row.get("label") or row["path"]
        prev = by_label.get(key)
        if prev is None or (row.get("timestamp_utc") or "") >= (
            prev.get("timestamp_utc") or ""
        ):
            by_label[key] = row
    reciprocity_rows = sorted(
        by_label.values(),
        key=lambda r: (r.get("formulation") or "", r.get("res") or 0),
    )

    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "branch_hint": "agent/eigenmode-ports",
        "scripts": SCRIPT_ROLES,
        "reciprocity_runs": reciprocity_rows,
        "mpi_benchmarks": mpi_rows,
        "faraday_runs": faraday_rows,
        "completed_te1_hz_line_P1P2": [
            r
            for r in reciprocity_rows
            if r.get("formulation") == "te1_hz_line"
            and r.get("ports") == [0, 1]
            and r.get("device_mode", "full") == "full"
        ],
        "notes": {
            "do_not_rerun": [
                "te1_hz_line P1P2 res32-128 full device (completed)",
                "measurement formulation cheap screen res32 (completed)",
                "faraday res32/64 PASS (completed)",
            ],
            "obsolete_or_superseded": [
                "te1_axis_at_source (broken source-plane flux)",
                "eigenmode on PEC horns (MPB failure)",
            ],
        },
    }


def main() -> int:
    reg = build_registry()
    out_path = OUT / "MASTER_VALIDATION_REGISTRY.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(reg, f, indent=2)
        f.write("\n")
    print(f"Wrote {out_path}")
    print(f"  reciprocity runs: {len(reg['reciprocity_runs'])}")
    print(f"  mpi benchmarks: {len(reg['mpi_benchmarks'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
