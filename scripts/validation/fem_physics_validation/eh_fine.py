#!/usr/bin/env python3
"""One finer mesh for the E/H samples that were still above 0.1%."""
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

from eh_orders import element_grad, rel_l2  # noqa: E402
from fem_canonical import plasma  # noqa: E402
from analytic_maxwell import kz_of  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from planar_fem import assign_rho, dirichlet, l2_error, rect_mesh, solve_dirichlet  # noqa: E402
import sixport_common as sc  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def plasma_one(h):
    saved = float(sc.fs_a)
    k0 = 2 * np.pi * saved
    eps = plasma(saved)
    kx = kz_of(k0, eps)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    pts, tris, _, _ = rect_mesh(1.5, 2.5, 1.5, 2.2, h)

    def ufn(x, y, kx=kx):
        return np.exp(1j * kx * (x - 1.5))

    bound = (
        (np.abs(pts[:, 0] - 1.5) < 1e-10) | (np.abs(pts[:, 0] - 2.5) < 1e-10)
        | (np.abs(pts[:, 1] - 1.5) < 1e-10) | (np.abs(pts[:, 1] - 2.2) < 1e-10)
    )
    rho = assign_rho(pts[tris].mean(1), lambda x, y, eps=eps: eps)
    uh = solve_dirichlet(pts, tris, rho, k0, ufn(pts[:, 0], pts[:, 1]), bound)
    _, rel = l2_error(pts, tris, uh, ufn)
    dux, duy, area, cents = element_grad(pts, tris, uh)
    u_c = ufn(cents[:, 0], cents[:, 1])
    ey = (1j / k0) * (1.0 / eps) * (-dux)
    exact = kx / (k0 * eps)
    interior = (cents[:, 0] > 1.5 + 2 * h) & (cents[:, 0] < 2.5 - 2 * h) & (cents[:, 1] > 1.5 + 2 * h) & (cents[:, 1] < 2.2 - 2 * h)
    point = np.abs(ey[interior] / u_c[interior] - exact) / np.abs(exact)
    grad = rel_l2(np.column_stack([dux, duy]), np.column_stack([1j * kx * u_c, np.zeros_like(u_c)]), area)
    rec = {"case": "plasma_fs", "h": h, "dofs": int(len(pts)), "hz_l2": rel, "grad_l2": grad, "element_EH_median": float(np.median(point)), "element_EH_max": float(np.max(point))}
    print("fine", rec, flush=True)
    return rec


def guide_one(h):
    width, k0 = 1.0, 5.0
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    x0, y0, length = 2.2, 2.2, 1.6

    def ufn(x, y):
        return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

    pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
    bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    b = np.zeros(len(pts), np.complex128)
    A, b = dirichlet(A, b, bound, ufn(pts[:, 0], pts[:, 1]))
    uh = splu(A.tocsc()).solve(b)
    _, rel = l2_error(pts, tris, uh, ufn)
    dux, duy, area, cents = element_grad(pts, tris, uh)
    u_c = ufn(cents[:, 0], cents[:, 1])
    gex = 1j * beta * u_c
    gey = -ky * np.sin(ky * (cents[:, 1] - y0)) * np.exp(1j * beta * (cents[:, 0] - x0))
    grad = rel_l2(np.column_stack([dux, duy]), np.column_stack([gex, gey]), area)
    ey = (1j / k0) * (-dux)
    exact = beta / k0
    band = (np.abs(cents[:, 1] - (y0 + 0.25 * width)) < 1.5 * h) & (cents[:, 0] > x0 + 3 * h) & (cents[:, 0] < x0 + length - 3 * h)
    point = np.abs(ey[band] / u_c[band] - exact) / np.abs(exact)
    rec = {
        "case": "guide_w1_k0_5",
        "h": h,
        "dofs": int(len(pts)),
        "hz_l2": rel,
        "grad_l2": grad,
        "element_EH_median": float(np.median(point)),
        "element_EH_max": float(np.max(point)),
    }
    print("fine", rec, flush=True)
    return rec


def main():
    rows = [plasma_one(0.0015625), guide_one(0.00125)]
    if rows[-1]["element_EH_median"] > 1e-3:
        rows.append(guide_one(0.0008))
    prev = json.loads((OUT / "eh_orders.json").read_text())
    prev["fine"] = rows
    (OUT / "eh_orders.json").write_text(json.dumps(prev, indent=2) + "\n")
    print("EH_FINE_DONE", flush=True)


if __name__ == "__main__":
    main()
