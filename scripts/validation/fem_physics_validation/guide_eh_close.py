#!/usr/bin/env python3
"""Finer guide mesh for the frozen 0.1% direct P1 E/H threshold.

Pass quantities are the element gradient L2 and the element-center Ey/Hz
sample. Nodal recovery is not used. The threshold is not changed.
"""
from __future__ import annotations

import json
import sys
import time
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
from eh_orders import element_grad, rel_l2  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from planar_fem import assign_rho, dirichlet, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def guide_close(h: float) -> dict:
    width, k0 = 1.0, 5.0
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    x0, y0, length = 2.2, 2.2, 1.6

    def ufn(x, y):
        return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

    t0 = time.perf_counter()
    pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
    print(f"mesh dofs {len(pts)} tris {len(tris)}", flush=True)
    bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    b = np.zeros(len(pts), np.complex128)
    A, b = dirichlet(A, b, bound, ufn(pts[:, 0], pts[:, 1]))
    t_fac = time.perf_counter()
    uh = splu(A.tocsc()).solve(b)
    fac_s = time.perf_counter() - t_fac
    dux, duy, area, cents = element_grad(pts, tris, uh)
    u_c = ufn(cents[:, 0], cents[:, 1])
    hz_c = uh[tris].mean(1)
    rel = float(np.sqrt(np.sum(np.abs(hz_c - u_c) ** 2 * area) / np.sum(np.abs(u_c) ** 2 * area)))
    gex = 1j * beta * u_c
    gey = -ky * np.sin(ky * (cents[:, 1] - y0)) * np.exp(1j * beta * (cents[:, 0] - x0))
    grad = rel_l2(np.column_stack([dux, duy]), np.column_stack([gex, gey]), area)
    ey = (1j / k0) * (-dux)
    exact = beta / k0
    band = (np.abs(cents[:, 1] - (y0 + 0.25 * width)) < 1.5 * h) & (cents[:, 0] > x0 + 3 * h) & (cents[:, 0] < x0 + length - 3 * h)
    point = np.abs(ey[band] / u_c[band] - exact) / np.abs(exact)
    # Area-weighted mean Sx on the strip x = x0+0.8, times the guide width.
    xcut = x0 + 0.8
    strip = np.abs(cents[:, 0] - xcut) <= 0.5 * h
    hz = uh[tris].mean(1)
    sx = 0.5 * np.real(ey * np.conj(hz))
    sx_mean = float(np.sum(sx[strip] * area[strip]) / np.sum(area[strip]))
    p_fem = sx_mean * width
    p_exact = 0.5 * float(np.real(beta)) * (width / 2.0) / k0
    rec = {
        "case": "guide_w1_k0_5",
        "h": h,
        "h_realized_x": float((length) / round(length / h)),
        "dofs": int(len(pts)),
        "hz_l2": rel,
        "grad_l2": grad,
        "element_EH_median": float(np.median(point)),
        "element_EH_max": float(np.max(point)),
        "n_element_samples": int(np.count_nonzero(band)),
        "power_fem": p_fem,
        "power_exact": p_exact,
        "power_rel": float(abs(p_fem - p_exact) / abs(p_exact)),
        "factor_seconds": fac_s,
        "seconds": time.perf_counter() - t0,
    }
    print("guide_close", rec, flush=True)
    return rec


def main():
    rec = guide_close(0.00064)
    prev = json.loads((OUT / "eh_orders.json").read_text())
    fine = [r for r in prev.get("fine", []) if not (r.get("case") == rec["case"] and abs(r.get("h", 9) - rec["h"]) < 1e-8)]
    fine.append(rec)
    prev["fine"] = fine
    (OUT / "eh_orders.json").write_text(json.dumps(prev, indent=2) + "\n")
    (OUT / "guide_eh_close.json").write_text(json.dumps(rec, indent=2) + "\n")
    print("GUIDE_EH_CLOSE_DONE", flush=True)


if __name__ == "__main__":
    main()
