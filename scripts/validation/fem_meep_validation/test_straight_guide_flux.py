#!/usr/bin/env python3
"""
Cheap straight-guide FEM test: numerical-mode source + guide-normal flux sign/conservation.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "scripts" / "validation" / "fem_meep_validation"), str(ROOT / "scripts" / "validation")]

from fem_validated_solver import (  # noqa: E402
    OUT,
    assemble_anisotropic,
    fields_from_hz,
    guide_normal_flux,
    inject_numerical_mode_rhs,
    load_mesh,
    load_numerical_mode,
    log,
    solve_system,
)
import sixport_common as sc  # noqa: E402


def main() -> int:
    grade = "FEM-M"  # cheaper than H for sign check
    src = load_numerical_mode(OUT / "phase2" / "numerical_mode_P1_res50.json")
    points, tris = load_mesh(grade)
    T = len(tris)
    rho = (
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    )
    # Straight feed walls only
    pec = np.zeros(len(points), dtype=bool)
    u, n, center = src["outward"], src["tangent"], src["center_mesh"]
    wall_off = sc.clear_width / 2 + sc.wall_thickness / 2
    axis_len = float(2.0 * np.hypot(sc.nx_ports, sc.ny_ports))
    rad = 0.55 * sc.wall_thickness
    for sign in (+1.0, -1.0):
        c0 = center + sign * wall_off * n - 0.5 * axis_len * u
        c1 = center + sign * wall_off * n + 0.5 * axis_len * u
        L = float(np.linalg.norm(c1 - c0))
        ns = max(2, int(L / max(0.5 * rad, 1e-6)))
        for t in np.linspace(0, 1, ns):
            c = c0 * (1 - t) + c1 * t
            pec |= (points[:, 0] - c[0]) ** 2 + (points[:, 1] - c[1]) ** 2 <= rad * rad

    k0 = 2 * np.pi * sc.fs_a
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    b = inject_numerical_mode_rhs(points, src, pec)
    x, _, t_fac, t_sol, resid = solve_system(A, b, pec)
    Ex, Ey = fields_from_hz(points, tris, x, *rho, 2 * np.pi * sc.fs_a)
    P0 = guide_normal_flux(points, x, Ex, Ey, 0)

    # Sample flux slightly inward (−n̂) and outward (+n̂) of source along guide
    # by shifting monitor idea: compare |Hz| energy left/right of source along axis
    # For port0, outward = +x in origin frame → +x in mesh. Source launches into device = −outward.
    center = src["center_mesh"]
    # Probe nodes ±1.0 along outward from source
    i_plus = int(np.argmin(np.sum((points - (center + 1.0 * u)) ** 2, axis=1)))
    i_minus = int(np.argmin(np.sum((points - (center - 1.0 * u)) ** 2, axis=1)))
    e_plus = float(np.abs(x[i_plus]) ** 2)
    e_minus = float(np.abs(x[i_minus]) ** 2)

    out = {
        "grade": grade,
        "P_raw_port0_monitor": float(P0),
        "P_inc_abs": float(abs(P0)),
        "true_residual": resid,
        "factor_s": t_fac,
        "solve_s": t_sol,
        "energy_plus_outward": e_plus,
        "energy_minus_outward": e_minus,
        "launch_into_device_expected": "energy_minus > energy_plus for P1 (source launches −n̂)",
        "launch_ok": bool(e_minus > e_plus),
        "flux_sign_note": (
            "Monitor at throat with +n̂ outward. Incident traveling −n̂ into device "
            "gives S·n̂ < 0, so P_raw_ref typically negative; Meep uses P_inc=|P_raw|."
        ),
    }
    path = OUT / "phase4" / "straight_guide_flux_sign.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    log(json.dumps(out, indent=2))
    return 0 if out["launch_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
