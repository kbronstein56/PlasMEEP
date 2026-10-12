#!/usr/bin/env python3
"""Straight parallel-plate guide with modal Robin (DtN) ends — no end PML."""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from feed_constants import clear_width, fs_a, wall_thickness  # noqa: E402
from fem_local import ElementSampler, assemble_anisotropic  # noqa: E402
from dtn_assembly import add_modal_incident_rhs, add_modal_robin  # noqa: E402
from modes import mode_constants, normalize_on_line, power_scale, transverse_mode  # noqa: E402


def build_mesh(L: float, W: float, wall: float, h: float) -> tuple[np.ndarray, np.ndarray]:
    y0, y1 = -0.5 * W - wall, 0.5 * W + wall
    nx = max(10, int(np.ceil(L / h)))
    ny = max(8, int(np.ceil((y1 - y0) / h)))
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    tri = Delaunay(pts)
    return pts.astype(float), np.asarray(tri.simplices, dtype=int)


def pec_wall_mask(points, tris, W: float, wall: float) -> np.ndarray:
    cents = points[tris].mean(axis=1)
    return np.abs(cents[:, 1]) > 0.5 * W


def port_nodes(points, x_port: float, W: float, tol: float = 1e-8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nodes on clear aperture at x=x_port, sorted by y."""
    mask = (np.abs(points[:, 0] - x_port) < tol) & (np.abs(points[:, 1]) <= 0.5 * W + tol)
    idx = np.flatnonzero(mask)
    order = idx[np.argsort(points[idx, 1])]
    y = points[order, 1]
    span = float(y[-1] - y[0]) if len(y) > 1 else W
    if len(y) == 1:
        w = np.array([span])
    else:
        ds = np.diff(y)
        w = np.zeros(len(y))
        w[0] = 0.5 * ds[0]
        w[-1] = 0.5 * ds[-1]
        w[1:-1] = 0.5 * (ds[:-1] + ds[1:])
    return order, y, w


def solve_guide(
    *,
    L: float = 10.0,
    h: float = 0.08,
    k0: Optional[float] = None,
    drive_mode: int = 0,
    drive_side: str = "left",
    n_modes_robin: int = 2,
    n_evanescent: int = 2,
) -> dict:
    W = float(clear_width)
    wall = float(wall_thickness)
    k0 = float(2 * np.pi * fs_a if k0 is None else k0)
    points, tris = build_mesh(L, W, wall, h)
    T = len(tris)
    rho = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    wmask = pec_wall_mask(points, tris, W, wall)
    for arr in rho:
        arr[wmask] = 0.0
    A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    b = np.zeros(len(points), dtype=np.complex128)

    modes = mode_constants(k0, W, n_modes=n_modes_robin + n_evanescent)
    left_idx, y_l, w_l = port_nodes(points, 0.0, W)
    right_idx, y_r, w_r = port_nodes(points, L, W)

    for md in modes:
        m = md["m"]
        phi_l, _ = normalize_on_line(y_l, w_l, m, W)
        phi_r, _ = normalize_on_line(y_r, w_r, m, W)
        beta = md["beta"]
        if md["propagating"] or md["m"] < n_modes_robin + n_evanescent:
            A = add_modal_robin(A, left_idx, y_l, w_l, phi_l, beta)
            A = add_modal_robin(A, right_idx, y_r, w_r, phi_r, beta)

    # Unit power incident on drive_mode
    md = modes[drive_mode]
    phi_drv, _ = normalize_on_line(y_l if drive_side == "left" else y_r,
                                   w_l if drive_side == "left" else w_r, drive_mode, W)
    # amplitude a such that |a|^2 = 1: use field amplitude a_f = 1/sqrt(scale^2) ...
    scale = power_scale(md["beta"], k0)
    a_field = 1.0 / max(abs(scale), 1e-30)
    if drive_side == "left":
        b = add_modal_incident_rhs(b, left_idx, y_l, w_l, phi_drv, md["beta"], a_field, side="left")
    else:
        b = add_modal_incident_rhs(b, right_idx, y_r, w_r, phi_drv, md["beta"], a_field, side="right")

    lu = splu(A.tocsc())
    hz = lu.solve(b)
    sampler = ElementSampler(points, tris)

    def modal_ab(port_idx, y, w, m):
        phi, _ = normalize_on_line(y, w, m, W)
        beta = modes[m]["beta"]
        xy = np.column_stack([np.full(len(y), points[port_idx[0], 0]), y])
        # sample Hz on port
        vals = []
        dn = []
        n_hat = np.array([-1.0, 0.0]) if port_idx is left_idx else np.array([1.0, 0.0])
        for j in range(len(y)):
            p = np.array([points[port_idx[j], 0], y[j]])
            tid = sampler.locate(p[None, :])[0]
            if tid < 0:
                i = int(port_idx[j])
                vals.append(hz[i])
                dn.append(0.0)
                continue
            nodes = tris[tid]
            ell = sampler._bary(int(tid), p)
            un = hz[nodes]
            hzv = complex(np.dot(ell, un))
            gx = complex(np.dot(sampler.bx[tid], un))
            gy = complex(np.dot(sampler.by[tid], un))
            vals.append(hzv)
            dn.append(n_hat[0] * gx + n_hat[1] * gy)
        vals = np.asarray(vals, np.complex128)
        dn = np.asarray(dn, np.complex128)
        cH = np.sum(w * phi * vals)
        cD = np.sum(w * phi * dn)
        bprop = complex(beta)
        if abs(bprop.imag) > 0:
            bprop = 1j * modes[m]["alpha"]
        sum_ab = cH
        diff = cD / (1j * bprop) if abs(bprop) > 1e-30 else 0.0
        a_f = 0.5 * (sum_ab - diff)
        b_f = 0.5 * (sum_ab + diff)
        scf = power_scale(bprop, k0)
        return a_f * scf, b_f * scf

    results = {"ports": {}}
    for name, idx, y, w in (("left", left_idx, y_l, w_l), ("right", right_idx, y_r, w_r)):
        results["ports"][name] = {}
        for m in range(n_modes_robin):
            a, bb = modal_ab(idx, y, w, m)
            results["ports"][name][f"m{m}"] = {
                "a": [a.real, a.imag],
                "b": [bb.real, bb.imag],
                "P_inc": float(abs(a) ** 2),
                "P_out": float(abs(bb) ** 2),
            }

    # S: column = drive mode at drive_side
    a_in = results["ports"]["left" if drive_side == "left" else "right"][f"m{drive_mode}"]["P_inc"]
    a_in = complex(*results["ports"]["left" if drive_side == "left" else "right"][f"m{drive_mode}"]["a"])
    S = {}
    for pname, idx, y, w in (("left", left_idx, y_l, w_l), ("right", right_idx, y_r, w_r)):
        for m in range(n_modes_robin):
            bout = complex(*results["ports"][pname][f"m{m}"]["b"])
            S[(pname, m)] = bout / a_in if abs(a_in) > 1e-30 else np.nan

    results["S"] = {f"{k[0]}_m{k[1]}": [v.real, v.imag] for k, v in S.items()}
    results["drive"] = {"mode": drive_mode, "side": drive_side, "a_in": [a_in.real, a_in.imag]}
    results["k0"] = k0
    results["L"] = L
    results["n_nodes"] = len(points)
    return results


if __name__ == "__main__":
    for dm in (0, 1):
        r = solve_guide(drive_mode=dm)
        print("drive m=", dm, "S", r["S"], flush=True)
