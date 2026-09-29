#!/usr/bin/env python3
"""Time the quartz FEM-VB mesh that passed the Meep resolution gate.

One factorization, then six lu.solve calls. P1 and P2 use the cached numerical
modes. Ports P3-P6 reuse the P2 load so the extra solves measure LU backsolves
rather than a new mode computation.
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
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "overnight_b0"


def main() -> int:
    import sixport_common as sc
    from b0_campaign import build_symmetric_mesh
    from fem_validated_solver import (
        ElementSampler,
        assemble_anisotropic,
        guide_normal_flux,
        inject_mode_consistent,
        load_numerical_mode,
    )
    from scipy.sparse.linalg import splu
    from sixport_common import set_geometry_context

    set_geometry_context(res=50, horn_walls="prism", grid_offset_cells=(0, 0), coord_rotation_deg=0.0)
    t0 = time.perf_counter()
    points, tris = build_symmetric_mesh(0.040, 0.012, 0.022, 0.028, h_quartz=0.010)
    t_mesh = time.perf_counter() - t0
    from overnight_b0 import material_rho

    t0 = time.perf_counter()
    rho, *_ = material_rho(points, tris, "quartz")
    t_mat = time.perf_counter() - t0
    k0 = 2 * np.pi * float(sc.fs_a)
    pec = np.zeros(len(points), dtype=bool)
    sampler = ElementSampler(points, tris)
    t0 = time.perf_counter()
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    t_asm = time.perf_counter() - t0
    src1 = load_numerical_mode()
    b1 = inject_mode_consistent(points, tris, src1, sampler)
    t0 = time.perf_counter()
    lu = splu(A.tocsc())
    t_fac = time.perf_counter() - t0
    t0 = time.perf_counter()
    x1 = lu.solve(b1)
    t_p1 = time.perf_counter() - t0
    p2_path = ROOT / "outputs" / "validation" / "mode_profiles" / "res50" / "numerical_mode_P2.json"
    src2 = load_numerical_mode(p2_path)
    b2 = inject_mode_consistent(points, tris, src2, sampler)
    extra = []
    t0 = time.perf_counter()
    x2 = lu.solve(b2)
    extra.append(time.perf_counter() - t0)
    for _ in range(4):
        t0 = time.perf_counter()
        lu.solve(b2)
        extra.append(time.perf_counter() - t0)
    t0 = time.perf_counter()
    powers = [
        guide_normal_flux(points, x1, None, None, p, sampler=sampler, rho=rho, omega=k0) for p in range(6)
    ]
    t_flux = time.perf_counter() - t0
    rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rec = {
        "mesh": "quartz FEM-VB",
        "dofs": int(len(points)),
        "mesh_s": t_mesh,
        "material_s": t_mat,
        "assemble_s": t_asm,
        "factor_s": t_fac,
        "p1_solve_s": t_p1,
        "p2_to_p6_solve_s": extra,
        "p2_to_p6_sum_s": float(sum(extra)),
        "port_integration_s": t_flux,
        "six_port_after_factor_s": float(t_p1 + sum(extra) + t_flux),
        "complete_s": float(t_mesh + t_mat + t_asm + t_fac + t_p1 + sum(extra) + t_flux),
        "ru_maxrss_gb": rss_kb / (1024.0 * 1024.0),
        "p1_raw": [complex(v).real for v in powers],
        "note": "P3-P6 backsolves reuse the P2 load. They measure solve time, not new port physics.",
    }
    print(json.dumps(rec, indent=2))
    (OUT / "quartz_vb_timing.json").write_text(json.dumps(rec, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
