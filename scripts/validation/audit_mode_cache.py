#!/usr/bin/env python3
"""Print numerical-mode cache identity for every six-port device port."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.dirname(os.path.abspath(__file__))
for p in (VAL, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from plasmeep.ports.mode_registry import audit_mode_cache  # noqa: E402
from plasmeep.ports.numerical_launch import port_tangent  # noqa: E402
from sixport_common import (  # noqa: E402
    effective_port_dir,
    fs_a,
    set_geometry_context,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", type=int, required=True)
    parser.add_argument("--horn-walls", type=str, default="prism")
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    set_geometry_context(res=args.res, horn_walls=args.horn_walls)
    audits = []
    print(f"\n=== mode cache audit res{args.res} ===", flush=True)
    for port in range(6):
        audit = audit_mode_cache(
            port,
            res=args.res,
            frequency_a=fs_a,
            horn_walls=args.horn_walls,
        )
        u = effective_port_dir(port)
        t = port_tangent(u)
        import numpy as np

        mo = np.array(audit.outward_dir)
        mt = np.array(audit.tangent)
        outward_dot = float(np.dot(u / np.linalg.norm(u), mo))
        tangent_dot = float(np.dot(t, mt / np.linalg.norm(mt)))
        row = audit.as_dict()
        row["device_outward_dot"] = outward_dot
        row["device_tangent_dot"] = tangent_dot
        audits.append(row)
        status = "OK" if audit.res_match and tangent_dot > 0.999 else "MISMATCH"
        print(
            f"  P{port + 1}: {status} file={Path(audit.cache_path).name} "
            f"stored_res={audit.mode_res} tangent_dot={tangent_dot:.4f}",
            flush=True,
        )

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "res": args.res,
        "frequency_a": fs_a,
        "horn_walls": args.horn_walls,
        "ports": audits,
    }
    if args.json_out:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
            f.write("\n")
        print(f"\nWrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
