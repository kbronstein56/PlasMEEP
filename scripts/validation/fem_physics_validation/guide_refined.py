#!/usr/bin/env python3
"""Finer PEC-guide meshes and a centered-element power integral."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_maxwell import kz_of  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from planar_fem import assign_rho, dirichlet, l2_error, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def run(width, k0, h):
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    x0, y0, length = 2.2, 2.2, 1.6

    def ufn(x, y):
        return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

    pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
    bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0)
    values = ufn(pts[:, 0], pts[:, 1])
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    b = np.zeros(len(pts), dtype=np.complex128)
    A, b = dirichlet(A, b, bound, values)
    uh = splu(A.tocsc()).solve(b)
    rel = l2_error(pts, tris, uh, ufn)[1]
    # Element power on triangles whose centroid is within h of x = x0+0.8, average Sx * width.
    xcut = x0 + 0.8
    acc = 0.0
    wsum = 0.0
    for tri in tris:
        xy = pts[tri]
        mid = xy.mean(0)
        if abs(mid[0] - xcut) > 0.5 * h:
            continue
        twice = (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
        area = 0.5 * abs(twice)
        bx = np.array([xy[1, 1] - xy[2, 1], xy[2, 1] - xy[0, 1], xy[0, 1] - xy[1, 1]]) / twice
        dhz = bx @ uh[tri]
        hz = np.mean(uh[tri])
        ey = (1j / k0) * (-dhz)
        sx = 0.5 * np.real(ey * np.conj(hz))
        acc += sx * area
        wsum += area
    sx_mean = acc / wsum
    p_fem = sx_mean * width
    p_exact = 0.5 * np.real(beta) * (width / 2.0) / k0
    rec = {
        "width": width,
        "k0": k0,
        "h": h,
        "rel_L2": rel,
        "power_rel": float(abs(p_fem - p_exact) / abs(p_exact)),
        "beta": [beta.real, beta.imag],
        "dofs": int(len(pts)),
    }
    print(rec, flush=True)
    return rec


def main():
    rows = []
    for width, k0 in ((1.0, 5.0), (0.8, 7.5), (1.4, 7.5)):
        for h in (0.005, 0.0025):
            rows.append(run(width, k0, h))
    (OUT / "guide_refined.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("GUIDE_REFINED", flush=True)


if __name__ == "__main__":
    main()
