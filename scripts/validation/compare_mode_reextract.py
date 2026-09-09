#!/usr/bin/env python3
"""Re-extract P1/P2 numerical modes and compare to cached JSON profiles."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

import numpy as np

VAL = Path(__file__).resolve().parent
ROOT = VAL.parents[1]
for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from mode_profile_study import extract_hz_line_profile  # noqa: E402
from plasmeep.ports.mode_registry import audit_mode_cache, get_numerical_mode  # noqa: E402
from plasmeep.ports.numerical_mode import NumericalPortMode  # noqa: E402
from sixport_common import (  # noqa: E402
    build_circulator_device,
    default_uniform_rho,
    effective_port_dir,
    fs_a,
    set_geometry_context,
    source_df,
)

OUT = ROOT / "outputs" / "validation" / "port_axis"


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def compare_port_mode(
    port_index: int,
    *,
    res: int,
    run_time: float,
) -> Dict[str, Any]:
    set_geometry_context(res=res)
    audit = audit_mode_cache(port_index, res=res, frequency_a=fs_a)
    cached = get_numerical_mode(port_index, res=res, frequency_a=fs_a)

    rho = default_uniform_rho()
    B = np.zeros(3)
    _pmm, p_device, _ = build_circulator_device(
        rho, B, res=res, device_mode="horns_only"
    )
    profile = extract_hz_line_profile(
        p_device,
        port_index,
        res,
        frequency=fs_a,
        fwidth=source_df,
        run_time=run_time,
    )
    fresh = NumericalPortMode.from_profile(
        profile,
        label=f"P{port_index + 1}_reextract_res{res}",
        frequency_a=fs_a,
        res=res,
        outward_dir=effective_port_dir(port_index),
    )

    ov = cached.overlap_with(fresh)
    power_ov = float(np.abs(ov) ** 2)
    max_amp_diff = float(
        np.max(np.abs(cached.normalized_field() - fresh.normalized_field()))
    )
    phase_diff_deg = float(np.degrees(np.angle(ov)))

    return {
        "port_index": port_index,
        "cached_label": cached.label,
        "cached_path": audit.cache_path,
        "cached_res": cached.res,
        "fresh_res": fresh.res,
        "complex_overlap": complex(ov),
        "power_overlap": power_ov,
        "phase_diff_deg": phase_diff_deg,
        "max_normalized_amp_diff": max_amp_diff,
        "monitor_xy_cached": list(cached.monitor_xy),
        "monitor_xy_fresh": list(fresh.monitor_xy),
        "n_samples_cached": len(cached.offsets_a),
        "n_samples_fresh": len(fresh.offsets_a),
        "identical_within_1e-6": bool(power_ov > 1.0 - 1e-6),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=64)
    parser.add_argument("--run-time", type=float, default=40.0)
    parser.add_argument("--ports", type=str, default="0,1")
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    ports = [int(x) for x in args.ports.split(",") if x.strip()]
    rows = [compare_port_mode(p, res=args.res, run_time=args.run_time) for p in ports]
    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "res": args.res,
        "run_time": args.run_time,
        "comparisons": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    out = Path(args.json_out) if args.json_out else OUT / f"mode_reextract_res{args.res}.json"
    out.write_text(json.dumps(_json_safe(payload), indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out}")
    for row in rows:
        print(
            f"P{row['port_index']+1}: power_overlap={row['power_overlap']:.6f} "
            f"max_amp_diff={row['max_normalized_amp_diff']:.3e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
