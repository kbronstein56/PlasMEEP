#!/usr/bin/env python3
"""Compare FEM vs Meep P1→all6 normalized powers; write comparison table JSON."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict


def db(x: float) -> float:
    return 10.0 * math.log10(x) if x > 0 else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--meep", required=True)
    ap.add_argument("--fem", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    meep = json.loads(Path(args.meep).read_text())
    fem = json.loads(Path(args.fem).read_text())
    mpow = {int(k): float(v) for k, v in meep["normalized_power_by_port"].items()}
    fpow = {int(k): float(v) for k, v in fem["normalized_power_by_port"].items()}

    rows = []
    dbs = []
    for p in range(6):
        mp, fp = mpow[p], fpow[p]
        mdb, fdb = db(mp), db(fp)
        ddb = fdb - mdb if (mp > 0 and fp > 0) else float("nan")
        rel = (fp - mp) / abs(mp) if abs(mp) > 0 else float("nan")
        if mp > 0 and fp > 0:
            dbs.append(ddb)
        rows.append(
            {
                "port": f"P1->P{p+1}",
                "meep_power": mp,
                "fem_power": fp,
                "meep_dB": mdb,
                "fem_dB": fdb,
                "FEM_minus_Meep_dB": ddb,
                "relative_power_error": rel,
            }
        )

    abs_db = [abs(x) for x in dbs]
    summary = {
        "meep_file": args.meep,
        "fem_file": args.fem,
        "rows": rows,
        "max_abs_dB_error": max(abs_db) if abs_db else float("nan"),
        "rms_dB_error": (sum(x * x for x in abs_db) / len(abs_db)) ** 0.5 if abs_db else float("nan"),
        "transmitted_ports_P2_P6_max_abs_dB": max(abs(rows[i]["FEM_minus_Meep_dB"]) for i in range(1, 6)
                                                   if rows[i]["meep_power"] > 0 and rows[i]["fem_power"] > 0),
        "verdict_band": None,
    }
    m = summary["max_abs_dB_error"]
    if m < 0.05:
        band = "excellent (<0.05 dB)"
    elif m < 0.10:
        band = "very good (0.05–0.10 dB)"
    elif m < 0.25:
        band = "possibly acceptable (0.10–0.25 dB)"
    elif m < 1.0:
        band = "not yet validated (0.25–1 dB)"
    else:
        band = "inconsistent (>1 dB)"
    summary["verdict_band"] = band

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
