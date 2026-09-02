#!/usr/bin/env python3
"""
Matched discrete-overlap vs point-probe reciprocity campaign.

Compares observables across material-complexity ladder and bulb-count scaling.
All Meep runs use mpi_runner (np=32).
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
OUT = ROOT / "outputs" / "validation" / "discrete_reciprocity"
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


def _extract_pair_metrics(data: Dict[str, Any]) -> Dict[str, Any]:
    pair = data.get("pairs", [{}])[0]
    out: Dict[str, Any] = {
        "label": data.get("label"),
        "mpi_ranks": data.get("mpi", {}).get("ranks"),
        "n_bulbs": data.get("settings", {}).get("n_bulbs"),
        "plasma_fill": data.get("settings", {}).get("plasma_fill"),
        "device_mode": data.get("settings", {}).get("device_mode"),
    }
    if "point" in pair:
        p = pair["point"]
        out["point_amp_err_dB"] = p.get("amp_err_dB")
        out["point_abs_a_to_b"] = p.get("abs_a_to_b")
        out["point_abs_b_to_a"] = p.get("abs_b_to_a")
        out["point_complex_sym_err"] = p.get("complex_sym_err")
    if "discrete" in pair:
        d = pair["discrete"]
        out["discrete_amp_err_dB"] = d.get("amp_err_dB")
        out["discrete_abs_g_a_to_b"] = d.get("abs_g_a_to_b")
        out["discrete_abs_g_b_to_a"] = d.get("abs_g_b_to_a")
        out["discrete_complex_sym_err"] = d.get("complex_sym_err")
    return out


def _run_case(
    *,
    label: str,
    res: int,
    run_time: float,
    device_mode: str = "full",
    plasma_fill: str = "active",
    n_bulbs: Optional[int] = None,
    observable: str = "both",
    ranks: int = DEFAULT_RANKS,
    skip_existing: bool = True,
) -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    json_out = OUT / f"{label}.json"
    log_out = OUT / f"{label}.log"
    if skip_existing and json_out.is_file():
        data = json.loads(json_out.read_text(encoding="utf-8"))
        row = _extract_pair_metrics(data)
        row["skipped"] = True
        row["json_out"] = str(json_out)
        return row

    args = [
        str(SCRIPT),
        "--res",
        str(res),
        "--run-time",
        str(run_time),
        "--pairs",
        "0,1",
        "--horn-walls",
        "prism",
        "--device-mode",
        device_mode,
        "--plasma-fill",
        plasma_fill,
        "--observable",
        observable,
        "--label",
        label,
        "--json-out",
        str(json_out),
    ]
    if n_bulbs is not None:
        args.extend(["--n-bulbs", str(n_bulbs)])

    t0 = time.perf_counter()
    proc = run_mpi_script(args, ranks=ranks, log_path=log_out)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"case {label} failed: exit={proc.returncode} log={log_out}")

    data = json.loads(json_out.read_text(encoding="utf-8"))
    row = _extract_pair_metrics(data)
    row["skipped"] = False
    row["json_out"] = str(json_out)
    row["elapsed_s"] = elapsed
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--run-time", type=float, default=20.0)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--no-skip-existing", dest="skip_existing", action="store_false")
    parser.add_argument("--include-res64", action="store_true", default=False)
    args = parser.parse_args()

    verify_mpi_ranks(ranks=args.ranks, python=DEFAULT_PYTHON, cwd=VAL)

    cases: List[Dict[str, Any]] = [
        {
            "label": f"compare_horns_only_res{args.res}_rt{args.run_time:g}",
            "device_mode": "horns_only",
            "plasma_fill": "active",
            "n_bulbs": None,
        },
        {
            "label": f"compare_quartz91_res{args.res}_rt{args.run_time:g}",
            "device_mode": "full",
            "plasma_fill": "geometry_only",
            "n_bulbs": None,
        },
        {
            "label": f"compare_drude7_res{args.res}_rt{args.run_time:g}",
            "device_mode": "full",
            "plasma_fill": "active",
            "n_bulbs": 7,
        },
        {
            "label": f"compare_drude91_res{args.res}_rt{args.run_time:g}",
            "device_mode": "full",
            "plasma_fill": "active",
            "n_bulbs": None,
        },
    ]

    bulb_counts = [0, 1, 7, 19, 37, 91]
    for n in bulb_counts:
        if n == 0:
            cases.append(
                {
                    "label": f"scale_n0_horns_res{args.res}_rt{args.run_time:g}",
                    "device_mode": "horns_only",
                    "plasma_fill": "active",
                    "n_bulbs": None,
                }
            )
        else:
            cases.append(
                {
                    "label": f"scale_n{n}_drude_res{args.res}_rt{args.run_time:g}",
                    "device_mode": "full",
                    "plasma_fill": "active",
                    "n_bulbs": n,
                }
            )

    if args.include_res64:
        cases.append(
            {
                "label": f"compare_drude91_res64_rt{args.run_time:g}",
                "device_mode": "full",
                "plasma_fill": "active",
                "n_bulbs": None,
                "res": 64,
            }
        )

    seen = set()
    rows: List[Dict[str, Any]] = []
    t_campaign = time.perf_counter()
    for spec in cases:
        label = spec["label"]
        if label in seen:
            continue
        seen.add(label)
        res = int(spec.get("res", args.res))
        print(f"\n=== {label} ===", flush=True)
        row = _run_case(
            label=label,
            res=res,
            run_time=args.run_time,
            device_mode=spec["device_mode"],
            plasma_fill=spec["plasma_fill"],
            n_bulbs=spec.get("n_bulbs"),
            observable="both",
            ranks=args.ranks,
            skip_existing=args.skip_existing,
        )
        rows.append(row)
        print(
            f"  point={row.get('point_amp_err_dB')} dB  "
            f"discrete={row.get('discrete_amp_err_dB')} dB  "
            f"skipped={row.get('skipped')}",
            flush=True,
        )

    comparison = []
    for row in rows:
        if row["label"].startswith("compare_"):
            comparison.append(
                {
                    "case": row["label"].replace(f"_res{args.res}_rt{args.run_time:g}", "").replace("_res64", ""),
                    "geometry_material": row["label"].split("_res")[0].replace("compare_", ""),
                    "point_amp_err_dB": row.get("point_amp_err_dB"),
                    "discrete_amp_err_dB": row.get("discrete_amp_err_dB"),
                    "point_complex_sym_err": row.get("point_complex_sym_err"),
                    "discrete_complex_sym_err": row.get("discrete_complex_sym_err"),
                    "mpi_ranks": row.get("mpi_ranks"),
                }
            )

    scaling = [
        {
            "n_bulbs": (0 if r["label"].startswith("scale_n0") else int(r["label"].split("_n")[1].split("_")[0])),
            "point_amp_err_dB": r.get("point_amp_err_dB"),
            "discrete_amp_err_dB": r.get("discrete_amp_err_dB"),
        }
        for r in rows
        if r["label"].startswith("scale_")
    ]

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {
            "res": args.res,
            "run_time": args.run_time,
            "ranks": args.ranks,
            "observable": "point vs discrete_overlap",
            "patch_half_width_cells": 1,
            "patch_weights": "uniform",
        },
        "comparison_table": comparison,
        "bulb_scaling": scaling,
        "rows": rows,
        "elapsed_s": time.perf_counter() - t_campaign,
    }
    summary_path = OUT / "campaign_summary.json"
    table_path = OUT / "comparison_table.json"
    summary_path.write_text(json.dumps(_json_safe(summary), indent=2) + "\n", encoding="utf-8")
    table_path.write_text(json.dumps(_json_safe(comparison), indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {summary_path}")
    print(f"Wrote {table_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
