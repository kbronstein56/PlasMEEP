#!/usr/bin/env python3
"""One overdense plasma disk in the validated PML box. FEM now, Meep with argv meep."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "overnight_b0"


def fem():
    from scipy.spatial import Delaunay

    import sixport_common as sc
    from fem_validated_solver import (
        ElementSampler,
        assemble_anisotropic,
        eps_tensor_at_bias,
        rho_from_eps,
        solve_system,
    )

    nx, ny, h = 16.0, 16.0, 0.025
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, 2.0
    xs = np.arange(0.0, nx + 1e-9, h)
    ys = np.arange(0.0, ny + 1e-9, h)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    tris = Delaunay(pts).simplices.astype(int)
    k0 = 2 * np.pi * float(sc.fs_a)
    rxx, rxy, ryx, ryy = rho_from_eps(*eps_tensor_at_bias(0.0))
    cents = pts[tris].mean(1)
    center = np.array([8.0, 8.0])
    r_p = 4.6 * float(sc.r_bulb_inner) / 6.5
    inside = np.sum((cents - center) ** 2, axis=1) <= r_p**2
    T = len(tris)
    rho = [np.ones(T, dtype=np.complex128), np.zeros(T, dtype=np.complex128), np.zeros(T, dtype=np.complex128), np.ones(T, dtype=np.complex128)]
    rho[0][inside] = rxx
    rho[3][inside] = ryy
    pec = np.zeros(len(pts), dtype=bool)
    A = assemble_anisotropic(pts, tris, *rho, k0, pec)
    sampler = ElementSampler(pts, tris)
    src = np.array([[4.0, 8.0]])
    b = np.zeros(len(pts), dtype=np.complex128)
    t = int(sampler.locate(src)[0])
    w = sampler._bary(t, src[0])
    for k in range(3):
        b[int(tris[t, k])] += complex(w[k])
    x, *_ = solve_system(A, b, pec)
    radii = [0.05, 0.12, 0.20, 0.28, 0.40, 0.80]
    ang = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    rows = []
    for r in radii:
        xy = center + r * np.column_stack([np.cos(ang), np.sin(ang)])
        hz, _, _ = sampler.fields(x, *rho, k0, xy)
        rows.append({"r_a": r, "mean_abs_Hz": float(np.mean(np.abs(hz))), "inside_disk": bool(r < r_p)})
    rec = {"r_p": r_p, "eps": [complex(eps_tensor_at_bias(0.0)[0]).real, complex(eps_tensor_at_bias(0.0)[0]).imag], "rings": rows, "n_inside": int(inside.sum())}
    print(json.dumps(rec, indent=2))
    (OUT / "plasma_disk_fem.json").write_text(json.dumps(rec, indent=2) + "\n")


if __name__ == "__main__":
    fem()
