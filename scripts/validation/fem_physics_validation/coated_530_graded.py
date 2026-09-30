#!/usr/bin/env python3
"""Graded meshes for the coated cylinder at 5.30 GHz.

The uniform h=0.006 mesh (6.71e6 nodes) is built, then SuperLU aborts in
gstrf with "malloc fails for local dworkptr[]". These meshes keep that
resolution only inside a disk that covers the scatterer and the probes, and
use a coarser size in the outer air.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_sweeps import plasma  # noqa: E402
from fem_scatterers import solve_case  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def main():
    saved = float(sc.fs_a)
    # h_far, h_wall, h_near, r_near. The first two levels are already in the JSON.
    levels = [(0.02, 0.0022, 0.0032, 3.2)]
    path = OUT / "coated_530_graded.json"
    rows = json.loads(path.read_text()) if path.exists() else []
    try:
        sc.fs_a = saved * 5.30 / 3.85
        eps = plasma(sc.fs_a)
        for spec in levels:
            t0 = time.perf_counter()
            got = solve_case("coated_5.30GHz_graded", np.zeros((1, 2)), 0.230, eps, True, [spec])
            rec = got[0]
            rec["seconds"] = time.perf_counter() - t0
            rec["ru_maxrss_kb"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            rows.append(rec)
            (OUT / "coated_530_graded.json").write_text(json.dumps(rows, indent=2) + "\n")
            print("SAVED", spec, "forward", rec["probes"]["forward"], "ring", rec.get("ring_l2"), flush=True)
            fwd = rec["probes"]["forward"]
            near = rec["probes"]["near_quartz"]
            if abs(fwd["db"]) <= 0.05 and abs(fwd["phase_deg"]) <= 0.5 and abs(near["db"]) <= 0.05 and rec["ring_l2"] <= 0.01:
                print("GRADED_BARS_MET", flush=True)
                break
    finally:
        sc.fs_a = saved
    print("COATED_530_GRADED_DONE", flush=True)


if __name__ == "__main__":
    main()
