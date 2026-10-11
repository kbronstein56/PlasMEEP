#!/usr/bin/env python3
"""Cheap 2-port straight parallel-plate FEM for modal-S validation.

Geometry (c=1 units, clear_width from sixport_common):
  x ∈ [0, L], y ∈ [-W/2 − wall, +W/2 + wall]
  PEC sidewalls via ρ→0 strips; x-PML at both ends.
  Flat Hz line source near left; monitors at ports 0 (left) and 1 (right).

Outward normals: port0 +x̂? No — for S-parameters we use *port outward from device*:
  Here the "device" is the guide interior between monitors.
  Port L (0): outward = −x̂ (toward left PML); incident into guide = +x̂.
  Port R (1): outward = +x̂; incident into guide = −x̂.

This matches the six-port convention (outward away from junction).
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from fem_validated_solver import ElementSampler  # noqa: E402
from modal_s_core import extract_ab, flat_mode_weights  # noqa: E402


def build_rect_mesh(L: float, W: float, h: float, wall: float) -> tuple[np.ndarray, np.ndarray]:
    """Structured-ish Delaunay mesh of the guide rectangle including wall strips."""
    y0, y1 = -0.5 * W - wall, 0.5 * W + wall
    nx = max(8, int(np.ceil(L / h)))
    ny = max(6, int(np.ceil((y1 - y0) / h)))
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    # Jitter interior slightly to avoid perfect grid slivers
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    rng = np.random.default_rng(0)
    interior = (pts[:, 0] > 1e-9) & (pts[:, 0] < L - 1e-9) & (pts[:, 1] > y0 + 1e-9) & (pts[:, 1] < y1 - 1e-9)
    pts[interior, 0] += rng.normal(0, 0.05 * h, interior.sum())
    pts[interior, 1] += rng.normal(0, 0.05 * h, interior.sum())
    tri = Delaunay(pts)
    return pts.astype(float), np.asarray(tri.simplices, dtype=int)


def assemble_guide(
    points,
    tris,
    rho_xx,
    rho_xy,
    rho_yx,
    rho_yy,
    k0: float,
    L: float,
    dpml: float,
    sigma_max: float,
) -> sparse.csr_matrix:
    """Anisotropic weak form with x-directed PML only (guide ends)."""
    n = len(points)
    cents = points[tris].mean(axis=1)
    omega = float(k0)
    x = cents[:, 0]
    sigx = np.zeros_like(x)
    m = x < dpml
    sigx[m] = sigma_max * ((dpml - x[m]) / dpml) ** 2
    m = x > L - dpml
    sigx[m] = sigma_max * ((x[m] - (L - dpml)) / dpml) ** 2
    sx = 1.0 + 1j * sigx / omega
    sy = np.ones_like(sx)

    px = points[tris, 0]
    py = points[tris, 1]
    twice = (px[:, 1] - px[:, 0]) * (py[:, 2] - py[:, 0]) - (px[:, 2] - px[:, 0]) * (py[:, 1] - py[:, 0])
    area = 0.5 * np.abs(twice)
    bx = np.stack([py[:, 1] - py[:, 2], py[:, 2] - py[:, 0], py[:, 0] - py[:, 1]], axis=1) / twice[:, None]
    by = np.stack([px[:, 2] - px[:, 1], px[:, 0] - px[:, 2], px[:, 1] - px[:, 0]], axis=1) / twice[:, None]
    mass = (k0**2) * sx * sy
    ax = sy / sx
    ay = sx / sy
    pec_el = (
        (np.abs(rho_xx) < 1e-10)
        & (np.abs(rho_xy) < 1e-10)
        & (np.abs(rho_yx) < 1e-10)
        & (np.abs(rho_yy) < 1e-10)
    )
    I = []
    J = []
    V = []
    for e in range(len(tris)):
        if area[e] < 1e-18 or pec_el[e]:
            continue
        be = bx[e]
        ce = by[e]
        Ae = area[e]
        rxx, rxy, ryx, ryy = rho_xx[e], rho_xy[e], rho_yx[e], rho_yy[e]
        Ke = np.zeros((3, 3), dtype=np.complex128)
        for a in range(3):
            for b in range(3):
                Ke[a, b] = Ae * (
                    rxx * ax[e] * be[a] * be[b]
                    + rxy * ce[a] * be[b]
                    + ryx * be[a] * ce[b]
                    + ryy * ay[e] * ce[a] * ce[b]
                    - mass[e] * (1.0 / 12.0) * (1.0 + (a == b))
                )
        nodes = tris[e]
        for a in range(3):
            for b in range(3):
                I.append(nodes[a])
                J.append(nodes[b])
                V.append(Ke[a, b])
    A = sparse.coo_matrix((V, (I, J)), shape=(n, n)).tocsr()
    # Pin disconnected metal / empty nodes
    deg = np.asarray(A.getnnz(axis=1)).ravel()
    empty = deg == 0
    if np.any(empty):
        A = A.tolil()
        for i in np.flatnonzero(empty):
            A[i, i] = 1.0
        A = A.tocsr()
    return A


def inject_flat(points, tris, sampler: ElementSampler, xy: np.ndarray, weights: np.ndarray) -> np.ndarray:
    b = np.zeros(len(points), dtype=np.complex128)
    tid = sampler.locate(xy)
    for s, (p, w) in enumerate(zip(xy, weights)):
        t = int(tid[s])
        if t < 0:
            i = int(np.argmin(np.sum((points - p) ** 2, axis=1)))
            b[i] += complex(w)
            continue
        nodes = tris[t]
        ell = sampler._bary(t, p)
        for k in range(3):
            b[int(nodes[k])] += complex(w) * float(ell[k])
    return b


def make_port_line(x_mon: float, outward: np.ndarray, W: float, n_points: int = 41) -> dict:
    span = 0.92 * W
    offsets, weights, _ = flat_mode_weights(n_points, span)
    tang = np.array([-outward[1], outward[0]], float)
    mon = np.array([x_mon, 0.0])
    return {
        "outward": np.asarray(outward, float) / np.linalg.norm(outward),
        "tangent": tang,
        "monitor_xy": np.vstack([mon + float(s) * tang for s in offsets]),
        "weights": weights,
        "offsets": offsets,
    }


def solve_straight_guide(
    *,
    k0: Optional[float] = None,
    L: float = 8.0,
    h: float = 0.12,
    dpml: float = 1.0,
    eps_slab: Optional[complex] = None,
    slab_x0: float = 3.5,
    slab_x1: float = 4.5,
    drive_port: int = 0,
) -> dict:
    """Solve 2-port guide; optional dielectric slab for known reflection."""
    W = float(sc.clear_width)
    wall = float(sc.wall_thickness)
    k0 = float(2 * np.pi * sc.fs_a if k0 is None else k0)
    points, tris = build_rect_mesh(L, W, h, wall)
    T = len(tris)
    cents = points[tris].mean(axis=1)
    rho = [
        np.ones(T, np.complex128),
        np.zeros(T, np.complex128),
        np.zeros(T, np.complex128),
        np.ones(T, np.complex128),
    ]
    # PEC walls outside clear aperture
    yabs = np.abs(cents[:, 1])
    wall_mask = yabs > 0.5 * W
    for arr in rho:
        arr[wall_mask] = 0.0
    if eps_slab is not None:
        slab = (cents[:, 0] >= slab_x0) & (cents[:, 0] <= slab_x1) & (~wall_mask)
        inv = 1.0 / complex(eps_slab)
        rho[0][slab] = inv
        rho[3][slab] = inv

    sigma_max = -np.log(1e-15) * 3.0 / (2.0 * dpml)
    A = assemble_guide(points, tris, *rho, k0, L, dpml, sigma_max)
    sampler = ElementSampler(points, tris)

    # Source / monitors inside non-PML region
    x_src_L = dpml + 0.6
    x_mon_L = dpml + 1.4
    x_mon_R = L - dpml - 1.4
    x_src_R = L - dpml - 0.6
    line_L = make_port_line(x_mon_L, np.array([-1.0, 0.0]), W)
    line_R = make_port_line(x_mon_R, np.array([+1.0, 0.0]), W)
    lines = [line_L, line_R]

    span = 0.92 * W
    offsets, weights, _ = flat_mode_weights(41, span)
    if drive_port == 0:
        src_c = np.array([x_src_L, 0.0])
        tang = np.array([0.0, 1.0])
    else:
        src_c = np.array([x_src_R, 0.0])
        tang = np.array([0.0, 1.0])
    src_xy = np.vstack([src_c + float(s) * tang for s in offsets])
    bvec = inject_flat(points, tris, sampler, src_xy, weights)

    lu = splu(A.tocsc())
    x = lu.solve(bvec)

    amps = [extract_ab(sampler, x, rho, k0, line) for line in lines]
    # Driven port index
    j = int(drive_port)
    a_j = amps[j].a
    S = np.array([amps[i].b / a_j for i in range(2)], dtype=np.complex128)

    # Analytic slab reflection for normal-incidence TEM (β=k0 in air)
    analytic = None
    if eps_slab is not None:
        # Transfer-matrix for Hz ~ E_y TEM: Z ∝ 1/n, n=√ε
        n = np.sqrt(complex(eps_slab))
        d = slab_x1 - slab_x0
        k1 = k0 * n
        # Interface air(0)->slab(1): r01=(Z1-Z0)/(Z1+Z0), Z∝1/n so Z1/Z0=1/n
        # For E-parallel / TEM in this 2D Hz formulation, power wave impedance
        # scales as 1/β_group; use Fresnel for magnetic-wall consistent TEM:
        # r = (1-n)/(1+n) for E-field; for our Hz-flat mode S· uses β/ω same form
        # as optical TEM with n=√ε.
        r01 = (1 - n) / (1 + n)
        t01 = 2 / (1 + n)
        r10 = -r01
        t10 = 2 * n / (1 + n)
        phase = np.exp(1j * k1 * d)
        # Multiple reflection
        denom = 1 - r10 * r10 * phase * phase
        r = r01 + t01 * t10 * r10 * phase * phase / denom
        t = t01 * t10 * phase / denom
        # Power transmission |t|² * (n?); for same media in/out, T_power = |t|²
        analytic = {
            "r": [complex(r).real, complex(r).imag],
            "t": [complex(t).real, complex(t).imag],
            "R": float(abs(r) ** 2),
            "T": float(abs(t) ** 2),
            "R_plus_T": float(abs(r) ** 2 + abs(t) ** 2),
            "eps": [complex(eps_slab).real, complex(eps_slab).imag],
            "d": float(d),
            "n": [complex(n).real, complex(n).imag],
        }

    return {
        "k0": k0,
        "L": L,
        "W": W,
        "n_nodes": int(len(points)),
        "n_tris": int(len(tris)),
        "drive_port": j,
        "a": [[amps[i].a.real, amps[i].a.imag] for i in range(2)],
        "b": [[amps[i].b.real, amps[i].b.imag] for i in range(2)],
        "P_inc": [amps[i].P_inc for i in range(2)],
        "P_out": [amps[i].P_out for i in range(2)],
        "P_net_modal": [amps[i].P_net_modal for i in range(2)],
        "P_net_poynting": [amps[i].P_net_poynting for i in range(2)],
        "overlap_quality": [amps[i].overlap_quality for i in range(2)],
        "S_col": [[S[i].real, S[i].imag] for i in range(2)],
        "S11": [S[0].real, S[0].imag] if j == 0 else None,
        "S21": [S[1].real, S[1].imag] if j == 0 else None,
        "|S11|^2": float(abs(S[0]) ** 2) if j == 0 else float(abs(S[1]) ** 2),
        "|S21|^2": float(abs(S[1]) ** 2) if j == 0 else float(abs(S[0]) ** 2),
        "analytic": analytic,
        "wall_frac": float(wall_mask.mean()),
    }
