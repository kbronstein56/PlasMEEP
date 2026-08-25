#!/usr/bin/env python3
"""
Cheap diagnostic: incident-power inequivalence across port orientations.

For each requested resolution, normalize all six ports (or a subset) and
report launched incident power. Axis ports (P1,P4) vs diagonal ports should
converge toward equality under grid refinement if the residual is
discretization of rotated sources.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

VAL = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(VAL, "..", ".."))
sys.path.insert(0, VAL)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from sixport_common import ensure_normalizations, port_dirs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res-list", type=str, default="32,48,64")
    parser.add_argument("--run-time", type=float, default=80.0)
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    resolutions = [int(x) for x in args.res_list.split(",") if x.strip()]
    rows = []
    for res in resolutions:
        t0 = time.time()
        cache = ensure_normalizations(
            res, run_time=args.run_time, force=True, ports=list(range(6)), verbose=True
        )
        powers = {int(k): float(v) for k, v in cache["incident_power_by_port"].items()}
        axis = np.mean([powers[0], powers[3]])
        diag = np.mean([powers[1], powers[2], powers[4], powers[5]])
        row = {
            "res": res,
            "run_time": args.run_time,
            "wall_s": time.time() - t0,
            "powers": powers,
            "axis_mean": float(axis),
            "diag_mean": float(diag),
            "axis_over_diag": float(axis / diag),
            "max_over_min": float(max(powers.values()) / min(powers.values())),
            "port_dirs": {i: port_dirs[i].tolist() for i in range(6)},
        }
        rows.append(row)
        print(
            f"res={res}: axis/diag={row['axis_over_diag']:.6f} "
            f"max/min={row['max_over_min']:.6f} wall={row['wall_s']:.1f}s"
        )

    out = {"rows": rows}
    json_out = args.json_out or os.path.join(
        ROOT, "outputs", "validation", "reciprocity", "incident_power_vs_res.json"
    )
    os.makedirs(os.path.dirname(json_out), exist_ok=True)
    # rank-0 write
    try:
        from mpi4py import MPI

        rank = MPI.COMM_WORLD.Get_rank()
    except Exception:
        rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", "0"))
    if rank == 0:
        with open(json_out, "w", encoding="utf-8") as handle:
            json.dump(out, handle, indent=2)
        print(f"Wrote {json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
