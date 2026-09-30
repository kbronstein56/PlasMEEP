#!/usr/bin/env python3
"""Oblique PML reflection from a Bloch waveguide, replacing the broken 13b box.

The old box imposed the unstretched incident wave on three walls, so its
interior error did not move with the PML. This test launches one oblique
mode with a Robin condition on the left, identifies the top and bottom by
the Bloch rule u_top = exp(i ky Ly) u_bottom, and terminates only the right
end in the PML. The test-function factor is the conjugate of that phase.

On a horizontal line in the homogeneous region the field is fit to

    A exp(+i kx x) + B exp(-i kx x)

and R_PML = B/A. The stretch in assemble_anisotropic is the standard
complex-coordinate form: xx weight sy/sx, yy weight sx/sy, mass sx*sy.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
import fem_validated_solver as fv  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from planar_fem import assign_rho, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
FS0 = float(sc.fs_a)


def solve_one(angle_deg: float, dp: float, sigma_scale: float, length: float, ly: float, h: float) -> dict:
    k0 = 2 * np.pi * FS0
    ang = np.deg2rad(angle_deg)
    kx = k0 * np.cos(ang)
    ky = k0 * np.sin(ang)
    fv.PML_SIGMA_SCALE = sigma_scale
    nx = length + 2 * dp
    ny = ly + 2 * dp
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, dp
    sc.fs_a = FS0
    x0, x1 = dp, dp + length + dp
    y0, y1 = dp, dp + ly
    pts, tris, _, _ = rect_mesh(x0, x1, y0, y1, h)
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool)).tolil()
    b = np.zeros(len(pts), np.complex128)
    left = np.abs(pts[:, 0] - x0) < 1e-8
    top = np.abs(pts[:, 1] - y1) < 1e-8
    bot = np.abs(pts[:, 1] - y0) < 1e-8
    # Full left edge, including the top corner. That corner is a Bloch slave;
    # its Robin residual is carried into the bottom test function below.
    order = np.flatnonzero(left)
    order = order[np.argsort(pts[order, 1])]
    for i, j in zip(order[:-1], order[1:]):
        hy = float(pts[j, 1] - pts[i, 1])
        ui = np.exp(1j * (kx * pts[i, 0] + ky * pts[i, 1]))
        uj = np.exp(1j * (kx * pts[j, 0] + ky * pts[j, 1]))
        for p, q, w in ((i, i, hy / 3), (j, j, hy / 3), (i, j, hy / 6), (j, i, hy / 6)):
            A[p, q] = A[p, q] - 1j * kx * w
        b[i] += -2j * kx * (hy / 3 * ui + hy / 6 * uj)
        b[j] += -2j * kx * (hy / 6 * ui + hy / 3 * uj)
    # Nodal trial values obey u_top = exp(i ky Ly) u_bottom. The form is
    # sesquilinear on real hats, so the test function contributes conj(phase).
    # Using phase on the test side leaves a net top/bottom flux and a decaying mode.
    phase = np.exp(1j * ky * ly)
    top_at = {round(float(pts[i, 0]), 8): int(i) for i in np.flatnonzero(top)}
    pairs = [(int(i), top_at[round(float(pts[i, 0]), 8)]) for i in np.flatnonzero(bot)]
    top_to_bottom = {t: bottom for bottom, t in pairs}
    co = A.tocoo()
    rows, cols, data = [], [], []
    for r, c, val in zip(co.row, co.col, co.data):
        r, c = int(r), int(c)
        fac = 1.0 + 0j
        if r in top_to_bottom:
            r = top_to_bottom[r]
            fac *= np.conj(phase)
        if c in top_to_bottom:
            c = top_to_bottom[c]
            fac *= phase
        rows.append(r)
        cols.append(c)
        data.append(fac * val)
    A = coo_matrix((data, (rows, cols)), shape=co.shape).tolil()
    b_new = np.zeros_like(b)
    for i, val in enumerate(b):
        if i in top_to_bottom:
            b_new[top_to_bottom[i]] += np.conj(phase) * val
        else:
            b_new[i] += val
    b = b_new
    for bottom, top_i in pairs:
        A[top_i, :] = 0.0
        A[:, top_i] = 0.0
        A[top_i, top_i] = 1.0
        A[top_i, bottom] = -phase
        b[top_i] = 0.0
    t0 = time.perf_counter()
    uh = splu(A.tocsc()).solve(b)
    fac_s = time.perf_counter() - t0
    # Propagating modal amplitude. Integration against exp(-i ky y) removes
    # evanescent Bloch harmonics that a single-y sample would mix into the fit.
    x_pml = dp + length
    x_nodes = np.unique(np.round(pts[:, 0], 8))
    x_nodes = x_nodes[(x_nodes > x0 + 0.5) & (x_nodes < x_pml - 0.35)]
    modal = []
    for x in x_nodes:
        col = np.flatnonzero(np.abs(pts[:, 0] - x) < 1e-8)
        col = col[np.argsort(pts[col, 1])]
        yy = pts[col, 1]
        uu = uh[col] * np.exp(-1j * ky * yy)
        trap = np.sum(0.5 * (uu[:-1] + uu[1:]) * np.diff(yy))
        modal.append((float(x), trap / ly))
    xx = np.array([m[0] for m in modal])
    vv = np.array([m[1] for m in modal])
    design = np.column_stack([np.exp(1j * kx * xx), np.exp(-1j * kx * xx)])
    coef, residual, _, _ = np.linalg.lstsq(design, vv, rcond=None)
    amp, refl = coef
    fit = float(np.linalg.norm(design @ coef - vv) / (np.linalg.norm(vv) + 1e-30))
    fv.PML_SIGMA_SCALE = 1.0
    return {
        "angle_deg": angle_deg,
        "dpml": dp,
        "sigma_scale": sigma_scale,
        "length": length,
        "ly": ly,
        "h": h,
        "dofs": int(len(pts)),
        "kx": float(kx),
        "ky": float(ky),
        "R_abs": float(abs(refl / amp)),
        "R_phase_deg": float(np.angle(refl / amp) * 180 / np.pi),
        "fit_rel": fit,
        "factor_seconds": fac_s,
    }


def main():
    rows = []
    cases = []
    for angle in (0.0, 25.0, 60.0):
        for dp in (0.4, 0.8, 1.2):
            cases.append((angle, dp, 1.0, 2.4, 1.0, 0.02))
    cases += [
        (0.0, 1.2, 0.5, 2.4, 1.0, 0.02),
        (0.0, 1.2, 2.0, 2.4, 1.0, 0.02),
        (25.0, 1.2, 0.5, 2.4, 1.0, 0.02),
        (25.0, 1.2, 2.0, 2.4, 1.0, 0.02),
        (60.0, 1.6, 1.0, 2.4, 1.0, 0.02),
        (60.0, 1.2, 2.0, 2.4, 1.0, 0.02),
        (0.0, 1.2, 1.0, 3.6, 1.0, 0.02),
        (25.0, 1.2, 1.0, 3.6, 1.0, 0.02),
        (25.0, 1.2, 1.0, 2.4, 1.0, 0.01),
    ]
    try:
        for case in cases:
            rec = solve_one(*case)
            rows.append(rec)
            print(rec, flush=True)
    finally:
        fv.PML_SIGMA_SCALE = 1.0
        sc.fs_a = FS0
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    (OUT / "oblique_pml_fit.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("OBLIQUE_PML_FIT_DONE", flush=True)


if __name__ == "__main__":
    main()
