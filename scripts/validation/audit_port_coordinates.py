#!/usr/bin/env python3
"""Audit port source/monitor coordinates and grid-index cross-sections."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

VAL = Path(__file__).resolve().parent
ROOT = VAL.parents[1]
for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from plasmeep.ports.grid_index_port import build_grid_index_port_line  # noqa: E402
from plasmeep.ports.lorentz_probe import hz_yee_site  # noqa: E402
from plasmeep.ports.mode_registry import get_numerical_mode  # noqa: E402
from plasmeep.ports.numerical_launch import port_tangent  # noqa: E402
from port_formulations import make_numerical_mode_sources  # noqa: E402
from sixport_common import (  # noqa: E402
    clear_width,
    effective_port_dir,
    fs_a,
    horn_for_port,
    monitor_center_for_port,
    set_geometry_context,
)


OUT = ROOT / "outputs" / "validation" / "port_axis"


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def audit_port(
    port_index: int,
    *,
    res: int,
    grid_offset: Tuple[float, float],
) -> Dict[str, Any]:
    set_geometry_context(res=res, grid_offset_cells=grid_offset)
    horn = horn_for_port(port_index, res)
    src_req = np.asarray(horn["source_center"], dtype=float)
    mon_req = np.asarray(monitor_center_for_port(port_index, res), dtype=float)
    u = effective_port_dir(port_index)
    tangent = port_tangent(u)
    span = 0.96 * clear_width

    mode = get_numerical_mode(port_index, res=res, frequency_a=fs_a)
    launch_coords = []
    for s, amp in zip(mode.offsets_a, mode.source_amplitudes(target_power=1.0)):
        xy = src_req + float(s) * tangent
        site = hz_yee_site(xy, res=res)
        launch_coords.append(
            {
                "offset_a": float(s),
                "requested_xy": [float(xy[0]), float(xy[1])],
                "yee_xy": list(site.grid_xy),
                "cell_index": list(site.cell_index),
                "frac_offset_cells": list(site.fractional_offset_cells),
                "amplitude_abs": float(abs(amp)),
            }
        )

    grid_line = build_grid_index_port_line(
        port_index,
        res=res,
        center_xy=mon_req,
        outward_dir=u,
        span_a=span,
    )

    return {
        "port": port_index,
        "label": f"P{port_index + 1}",
        "res": res,
        "grid_offset_cells": list(grid_offset),
        "outward_dir": u.tolist(),
        "tangent": tangent.tolist(),
        "source_center_requested": src_req.tolist(),
        "source_center_yee": hz_yee_site(src_req, res=res).as_dict(),
        "monitor_center_requested": mon_req.tolist(),
        "monitor_center_yee": hz_yee_site(mon_req, res=res).as_dict(),
        "launch_samples": launch_coords,
        "grid_index_monitor_line": grid_line.as_dict(),
        "n_launch_dofs": len(launch_coords),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, default=32)
    parser.add_argument("--ports", type=str, default="0,1,2")
    parser.add_argument(
        "--offsets",
        type=str,
        default="0,0;0,0.5;0,-0.5;0.5,0;-0.5,0",
        help="semicolon-separated ox,oy pairs",
    )
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    ports = [int(x) for x in args.ports.split(",") if x.strip()]
    offset_list: List[Tuple[float, float]] = []
    for chunk in args.offsets.split(";"):
        ox, oy = [float(v) for v in chunk.split(",")]
        offset_list.append((ox, oy))

    rows = []
    for gox, goy in offset_list:
        for p in ports:
            rows.append(audit_port(p, res=args.res, grid_offset=(gox, goy)))

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "res": args.res,
        "audits": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    out = (
        Path(args.json_out)
        if args.json_out
        else OUT / f"port_coord_audit_res{args.res}.json"
    )
    out.write_text(json.dumps(_json_safe(payload), indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
