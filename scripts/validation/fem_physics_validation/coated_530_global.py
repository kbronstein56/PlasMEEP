#!/usr/bin/env python3
"""Finer global meshes for the coated cylinder at 5.30 GHz.

h_edge is half of h, so the quartz wall, the plasma, and the air are refined
together. The analytic reference is the same multilayer solution as the
coarser rows.
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
    levels = [(0.008, 0.004), (0.006, 0.003)]
    rows = []
    try:
        sc.fs_a = saved * 5.30 / 3.85
        eps = plasma(sc.fs_a)
        for h, he in levels:
            t0 = time.perf_counter()
            got = solve_case("coated_5.30GHz_global", np.zeros((1, 2)), 0.230, eps, True, [(h, he)])
            rec = got[0]
            rec["seconds"] = time.perf_counter() - t0
            rec["ru_maxrss_kb"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            rows.append(rec)
            (OUT / "coated_530_global.json").write_text(json.dumps(rows, indent=2) + "\n")
            print("SAVED", h, "forward", rec["probes"]["forward"], "ring", rec.get("ring_l2"), flush=True)
    finally:
        sc.fs_a = saved
    print("COATED_530_GLOBAL_DONE", flush=True)


if __name__ == "__main__":
    main()
