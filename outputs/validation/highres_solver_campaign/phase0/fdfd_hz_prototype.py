#!/usr/bin/env python3
"""
Benchmark-only scalar-Hz FDFD prototype for PlasMEEP feasibility (NOT production).

TE-to-z (Hz, Ex, Ey), μ=μ0, d/dz=0.
B=0:  ∇·((1/ε)∇Hz) + k0² Hz = 0  (+ CFS-like complex stretch PML)
B≠0:  ∇·(ρ ∇Hz) + k0² Hz = 0 with ρ = ε_2D^{-1} from Faraday/Meep tensor.

Stages: level1 straight guide | level2 horns_only | level3 full 91-bulb PMM
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import gmres, bicgstab, spilu, LinearOperator, spsolve

ROOT = Path(__file__).resolve().parents[3]
VAL = ROOT / "scripts" / "validation"
for p in (str(VAL), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import sixport_common as sc  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402
from physical_units import meep_resolution_from_points_per_cm  # noqa: E402
from sixport_common import (  # noqa: E402
    build_circulator_device,
    default_uniform_rho,
    set_geometry_context,
)

OUT = Path(__file__).resolve().parent


def _rss_gib() -> float:
    # Linux: ru_maxrss is KiB
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 ** 2)


def log(msg: str) -> None:
    print(msg, flush=True)


@dataclass
class Grid:
    res: int
    Nx: int
    Ny: int
    dx: float
    x: np.ndarray
    y: np.ndarray


def make_grid(res: int) -> Grid:
    Nx, Ny = int(sc.nx_ports * res), int(sc.ny_ports * res)
    dx = 1.0 / res
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dx
    return Grid(res=res, Nx=Nx, Ny=Ny, dx=dx, x=x, y=y)


def eps_tensor_components(fs_a: float, fp_a: float, gamma_a: float, bias_a: float):
    """Meep/Faraday e^{-iωt} in-plane tensor: [[eps_perp, -1j*eta],[+1j*eta, eps_perp]]."""
    eps_perp, eta = gyrotropic_drude_eps_eta(fs_a, fp_a, gamma_a, bias_a)
    eps_xx = eps_perp
    eps_yy = eps_perp
    eps_xy = -1j * eta
    eps_yx = +1j * eta
    det = eps_xx * eps_yy - eps_xy * eps_yx
    # ρ = ε^{-1}
    rho_xx = eps_yy / det
    rho_yy = eps_xx / det
    rho_xy = -eps_xy / det
    rho_yx = -eps_yx / det
    return {
        "eps_perp": eps_perp,
        "eta": eta,
        "eps_xx": eps_xx,
        "eps_xy": eps_xy,
        "eps_yx": eps_yx,
        "eps_yy": eps_yy,
        "rho_xx": rho_xx,
        "rho_xy": rho_xy,
        "rho_yx": rho_yx,
        "rho_yy": rho_yy,
        "det": det,
    }


def material_maps(
    level: str, grid: Grid, bias_a: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns rho_xx, rho_xy, rho_yx, rho_yy, pec  (all Ny,Nx; rho complex).
    Outside plasma: isotropic ε=1 => ρ=I. Quartz ε=3.8. PEC mask for horn/guide walls.
    """
    Nx, Ny, dx = grid.Nx, grid.Ny, grid.dx
    X, Y = np.meshgrid(grid.x, grid.y)
    xc = X - sc.nx_ports / 2
    yc = Y - sc.ny_ports / 2

    rho_xx = np.ones((Ny, Nx), dtype=np.complex128)
    rho_yy = np.ones((Ny, Nx), dtype=np.complex128)
    rho_xy = np.zeros((Ny, Nx), dtype=np.complex128)
    rho_yx = np.zeros((Ny, Nx), dtype=np.complex128)
    pec = np.zeros((Ny, Nx), dtype=bool)

    ten = eps_tensor_components(sc.fs_a, sc.fp_a, sc.gamma_a, bias_a)
    eps_q = 3.8 + 0j

    if level == "level1":
        half = 0.5 * sc.clear_width
        wall = (
            (np.abs(yc) > half)
            & (np.abs(yc) < half + sc.wall_thickness)
            & (np.abs(xc) < sc.nx_ports / 2 - sc.dpml_ports)
        )
        pec |= wall
        return rho_xx, rho_xy, rho_yx, rho_yy, pec

    set_geometry_context(res=grid.res, horn_walls="prism", grid_offset_cells=(0, 0), coord_rotation_deg=0.0)
    rho_opt = default_uniform_rho()
    B = np.array([0.0, 0.0, 0.0])
    device_mode = "horns_only" if level == "level2" else "full"
    pmm, _P, _ = build_circulator_device(rho_opt, B, res=grid.res, device_mode=device_mode, wall_pec=True)

    if level == "level3":
        locs = np.asarray(pmm.train_elem_locs, dtype=float)
        r_q = sc.r_bulb_outer
        r_vac = sc.r_bulb_inner
        r_p = 4.6 * sc.r_bulb_inner / 6.5
        for i in range(len(locs)):
            cx0, cy0 = float(locs[i, 0]), float(locs[i, 1])
            r2 = (X - cx0) ** 2 + (Y - cy0) ** 2
            mask_q = (r2 <= r_q**2) & (r2 > r_vac**2)
            rho_xx[mask_q] = 1.0 / eps_q
            rho_yy[mask_q] = 1.0 / eps_q
            rho_xy[mask_q] = 0.0
            rho_yx[mask_q] = 0.0
            mask_p = r2 <= r_p**2
            rho_xx[mask_p] = ten["rho_xx"]
            rho_yy[mask_p] = ten["rho_yy"]
            rho_xy[mask_p] = ten["rho_xy"]
            rho_yx[mask_p] = ten["rho_yx"]

    for horn in sc.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], dtype=float)
            for k in range(len(poly) - 1):
                p0, p1 = poly[k], poly[k + 1]
                nseg = max(int(np.linalg.norm(p1 - p0) / (0.5 * dx)), 2)
                ts = np.linspace(0, 1, nseg)
                pxs = p0[0] * (1 - ts) + p1[0] * ts
                pys = p0[1] * (1 - ts) + p1[1] * ts
                ixs = ((pxs + sc.nx_ports / 2) / dx).astype(int)
                iys = ((pys + sc.ny_ports / 2) / dx).astype(int)
                rad = max(1, int(0.5 * sc.wall_thickness / dx))
                for jx, jy in zip(ixs, iys):
                    y0, y1 = max(0, jy - rad), min(Ny, jy + rad + 1)
                    x0, x1 = max(0, jx - rad), min(Nx, jx + rad + 1)
                    pec[y0:y1, x0:x1] = True

    return rho_xx, rho_xy, rho_yx, rho_yy, pec


def pml_sigma_1d(n: int, dpml_cells: int, sigma_max: float = 2.0) -> np.ndarray:
    s = np.zeros(n, dtype=np.float64)
    if dpml_cells <= 0:
        return s
    for i in range(dpml_cells):
        u = (dpml_cells - i) / dpml_cells
        s[i] = sigma_max * u**2
        s[n - 1 - i] = sigma_max * u**2
    return s


def build_operator(
    grid: Grid,
    rho_xx: np.ndarray,
    rho_xy: np.ndarray,
    rho_yx: np.ndarray,
    rho_yy: np.ndarray,
    pec: np.ndarray,
    k0: float,
) -> sparse.csr_matrix:
    """
    Vectorized 5-point (isotropic) / 9-point (gyrotropic) discretization of
        ∂x(sx^{-1}[ρxx ∂x + ρxy ∂y]Hz) + ∂y(sy^{-1}[ρyx ∂x + ρyy ∂y]Hz) + k0² sx sy Hz
    PEC nodes: Dirichlet row (unused DOFs); air adjacent to PEC gets zero normal flux (TE Neumann).
    """
    Nx, Ny, dx = grid.Nx, grid.Ny, grid.dx
    N = Nx * Ny
    omega = 2 * np.pi * sc.fs_a
    dpml_cells = int(round(sc.dpml_ports * grid.res))
    sigx = pml_sigma_1d(Nx, dpml_cells)
    sigy = pml_sigma_1d(Ny, dpml_cells)
    sx = 1.0 + 1j * sigx / max(omega, 1e-30)
    sy = 1.0 + 1j * sigy / max(omega, 1e-30)
    Sx = np.broadcast_to(sx[None, :], (Ny, Nx)).copy()
    Sy = np.broadcast_to(sy[:, None], (Ny, Nx)).copy()

    # Zero ρ inside PEC so flux into metal vanishes
    pec = pec.astype(bool)
    for arr in (rho_xx, rho_xy, rho_yx, rho_yy):
        arr[pec] = 0.0

    inv_dx2 = 1.0 / (dx * dx)
    inv_dx2_4 = 1.0 / (4.0 * dx * dx)

    def harm(a, b):
        return 2 * a * b / (a + b + 1e-30)

    # Face-averaged ρxx / ρyy for principal axes
    rxx_xp = np.zeros_like(rho_xx)
    rxx_xm = np.zeros_like(rho_xx)
    ryy_yp = np.zeros_like(rho_yy)
    ryy_ym = np.zeros_like(rho_yy)
    rxx_xp[:, :-1] = harm(rho_xx[:, :-1], rho_xx[:, 1:])
    rxx_xm[:, 1:] = harm(rho_xx[:, 1:], rho_xx[:, :-1])
    ryy_yp[:-1, :] = harm(rho_yy[:-1, :], rho_yy[1:, :])
    ryy_ym[1:, :] = harm(rho_yy[1:, :], rho_yy[:-1, :])

    # Cross terms: cell-centered ρxy, ρyx (simple; sufficient for feasibility)
    # ∂x(ρxy ∂y u) ≈ (ρxy_{i+1/2} (u_{i+1,j+1}-u_{i+1,j-1}) - ρxy_{i-1/2} (u_{i-1,j+1}-u_{i-1,j-1})) / (4 dx²)
    # Use arithmetic face averages of ρxy / ρyx
    rxy_xp = np.zeros_like(rho_xy)
    rxy_xm = np.zeros_like(rho_xy)
    ryx_yp = np.zeros_like(rho_yx)
    ryx_ym = np.zeros_like(rho_yx)
    rxy_xp[:, :-1] = 0.5 * (rho_xy[:, :-1] + rho_xy[:, 1:])
    rxy_xm[:, 1:] = 0.5 * (rho_xy[:, 1:] + rho_xy[:, :-1])
    ryx_yp[:-1, :] = 0.5 * (rho_yx[:-1, :] + rho_yx[1:, :])
    ryx_ym[1:, :] = 0.5 * (rho_yx[1:, :] + rho_yx[:-1, :])

    cx = inv_dx2 / Sx
    cy = inv_dx2 / Sy

    diag = -cx * (rxx_xp + rxx_xm) - cy * (ryy_yp + ryy_ym) + (k0**2) * Sx * Sy
    # Cross-term diagonal contributions cancel at lowest order for centered diffs; omit.

    def idx(iy, ix):
        return iy * Nx + ix

    # Build diagonals via shifts (principal part)
    data = []
    rows = []
    cols = []

    def add_diag_array(coeff: np.ndarray, diy: int, dix: int):
        iy0 = max(0, -diy)
        iy1 = Ny - max(0, diy)
        ix0 = max(0, -dix)
        ix1 = Nx - max(0, dix)
        c = coeff[iy0:iy1, ix0:ix1].ravel()
        if c.size == 0:
            return
        IY, IX = np.mgrid[iy0:iy1, ix0:ix1]
        r = (IY * Nx + IX).ravel()
        cc = ((IY + diy) * Nx + (IX + dix)).ravel()
        mask = np.abs(c) > 0
        # skip pec rows — handled separately
        pec_flat = pec[iy0:iy1, ix0:ix1].ravel()
        mask &= ~pec_flat
        if not np.any(mask):
            return
        rows.append(r[mask])
        cols.append(cc[mask])
        data.append(c[mask])

    add_diag_array(diag, 0, 0)
    add_diag_array(cx * rxx_xp, 0, +1)
    add_diag_array(cx * rxx_xm, 0, -1)
    add_diag_array(cy * ryy_yp, +1, 0)
    add_diag_array(cy * ryy_ym, -1, 0)

    # Gyrotropic cross terms (only if any nonzero)
    if np.max(np.abs(rho_xy)) + np.max(np.abs(rho_yx)) > 0:
        cxy = inv_dx2_4 / Sx
        cyx = inv_dx2_4 / Sy
        # ∂x(ρxy ∂y u): +rxy_xp * (u_{+1,+1} - u_{+1,-1}) - rxy_xm * (u_{-1,+1} - u_{-1,-1})
        add_diag_array(+cxy * rxy_xp, +1, +1)
        add_diag_array(-cxy * rxy_xp, -1, +1)
        add_diag_array(-cxy * rxy_xm, +1, -1)
        add_diag_array(+cxy * rxy_xm, -1, -1)
        # ∂y(ρyx ∂x u): +ryx_yp * (u_{+1,+1} - u_{-1,+1}) - ryx_ym * (u_{+1,-1} - u_{-1,-1})
        add_diag_array(+cyx * ryx_yp, +1, +1)
        add_diag_array(-cyx * ryx_yp, +1, -1)
        add_diag_array(-cyx * ryx_ym, -1, +1)
        add_diag_array(+cyx * ryx_ym, -1, -1)

    # PEC Dirichlet rows
    pec_idx = np.flatnonzero(pec.ravel())
    if pec_idx.size:
        rows.append(pec_idx)
        cols.append(pec_idx)
        data.append(np.ones(pec_idx.size, dtype=np.complex128))

    r = np.concatenate(rows) if rows else np.array([], dtype=np.int64)
    c = np.concatenate(cols) if cols else np.array([], dtype=np.int64)
    v = np.concatenate(data) if data else np.array([], dtype=np.complex128)
    A = sparse.coo_matrix((v, (r, c)), shape=(N, N)).tocsr()
    return A


def _to_cell(xy_array_centered: np.ndarray) -> np.ndarray:
    """Array-centered coords → Meep cell coords (array at nx/2, ny/2)."""
    return np.asarray(xy_array_centered, dtype=float) + np.array(
        [sc.nx_ports / 2, sc.ny_ports / 2]
    )


def inject_mode_source(grid: Grid, port_index: int = 0, level: str = "level2") -> np.ndarray:
    Nx, Ny, dx = grid.Nx, grid.Ny, grid.dx
    b = np.zeros(Nx * Ny, dtype=np.complex128)
    if level == "level1":
        cx = sc.nx_ports / 2
        cy = sc.ny_ports / 2
        span = 0.96 * sc.clear_width
        for s in np.linspace(-span / 2, span / 2, 31):
            xy = np.array([cx, cy + s])
            ix = int(xy[0] / dx)
            iy = int(xy[1] / dx)
            if 0 <= ix < Nx and 0 <= iy < Ny:
                amp = max(np.cos(np.pi * s / span), 0.0)
                b[iy * Nx + ix] += amp
        nrm = np.linalg.norm(b)
        if nrm > 0:
            b /= nrm
        return b
    set_geometry_context(res=grid.res, horn_walls="prism")
    horn = sc.horn_for_port(port_index, grid.res)
    center = _to_cell(horn["source_center"])
    u = np.asarray(sc.effective_port_dir(port_index), dtype=float)
    u /= np.linalg.norm(u)
    tangent = np.array([-u[1], u[0]])
    span = 0.96 * sc.clear_width
    for s in np.linspace(-span / 2, span / 2, 31):
        xy = center + s * tangent
        ix = int(xy[0] / dx)
        iy = int(xy[1] / dx)
        if 0 <= ix < Nx and 0 <= iy < Ny:
            amp = max(np.cos(np.pi * s / span), 0.0)
            b[iy * Nx + ix] += amp
    nrm = np.linalg.norm(b)
    if nrm > 0:
        b /= nrm
    return b


def flux_proxy_through_port(grid: Grid, u: np.ndarray, port_index: int, pec: Optional[np.ndarray] = None) -> complex:
    """
    Rough power proxy: integrate |Hz|^2 on monitor line (skip PEC samples).
    Not exact Poynting; relative GO/NO-GO only after normalizing to port0.
    """
    Nx, Ny, dx = grid.Nx, grid.Ny, grid.dx
    U = u.reshape(Ny, Nx)
    set_geometry_context(res=grid.res, horn_walls="prism")
    center = _to_cell(sc.monitor_center_for_port(port_index, grid.res))
    uhat = np.asarray(sc.effective_port_dir(port_index), dtype=float)
    uhat /= np.linalg.norm(uhat)
    tangent = np.array([-uhat[1], uhat[0]])
    span = 0.96 * sc.clear_width
    acc = 0.0
    wsum = 0.0
    for s in np.linspace(-span / 2, span / 2, 63):
        xy = center + s * tangent
        ix = int(xy[0] / dx)
        iy = int(xy[1] / dx)
        if 0 <= ix < Nx and 0 <= iy < Ny:
            if pec is not None and pec[iy, ix]:
                continue
            w = max(np.cos(np.pi * s / span), 0.0)
            y0, y1 = max(0, iy - 1), min(Ny, iy + 2)
            x0, x1 = max(0, ix - 1), min(Nx, ix + 2)
            patch = U[y0:y1, x0:x1]
            acc += w * float(np.mean(np.abs(patch) ** 2))
            wsum += w
    return acc / (wsum + 1e-30)


def make_preconditioner(A: sparse.csr_matrix, kind: str) -> Tuple[LinearOperator, float, Dict[str, Any]]:
    t0 = time.perf_counter()
    if kind == "jacobi":
        diag = A.diagonal()
        diag = np.where(np.abs(diag) > 1e-30, diag, 1.0 + 0j)
        inv = 1.0 / diag

        def mv(x):
            return inv * x

        M = LinearOperator(A.shape, matvec=mv, dtype=np.complex128)
        meta: Dict[str, Any] = {"kind": "jacobi"}
    elif kind == "shifted":
        # Poor-man complex-shifted Jacobi: diag(A - i σ)^{-1}, σ ~ k0-ish scale from mean |diag|
        diag = A.diagonal().astype(np.complex128)
        scale = float(np.median(np.abs(diag))) + 1e-30
        sigma = 0.5j * scale
        d = diag - sigma
        d = np.where(np.abs(d) > 1e-30, d, 1.0 + 0j)
        inv = 1.0 / d

        def mv(x):
            return inv * x

        M = LinearOperator(A.shape, matvec=mv, dtype=np.complex128)
        meta = {"kind": "shifted", "sigma": complex(sigma)}
    elif kind == "ilu0":
        ilu = spilu(A.tocsc(), drop_tol=1e-2, fill_factor=1.0)
        M = LinearOperator(A.shape, matvec=ilu.solve, dtype=np.complex128)
        meta = {"kind": "ilu0", "ilu": ilu}
    else:
        raise ValueError(kind)
    return M, time.perf_counter() - t0, meta


def solve_system(
    A: sparse.csr_matrix,
    b: np.ndarray,
    prec: Optional[Any] = None,
    precond_kind: str = "jacobi",
    tol: float = 1e-5,
    maxiter: int = 500,
    restart: int = 40,
    solver: str = "bicgstab",
) -> Dict[str, Any]:
    from scipy.sparse.linalg import splu

    residuals: List[float] = []
    niter = [0]

    if solver == "direct":
        if prec is not None and "lu" in prec:
            t_prec = 0.0
            lu = prec["lu"]
            reusable = prec
            t0 = time.perf_counter()
            x = lu.solve(b)
            t_solve = time.perf_counter() - t0
        else:
            t0 = time.perf_counter()
            lu = splu(A.tocsc())
            t_prec = time.perf_counter() - t0
            t1 = time.perf_counter()
            x = lu.solve(b)
            t_solve = time.perf_counter() - t1
            reusable = {"lu": lu, "kind": "splu", "M": None}
        info = 0
        niter[0] = -1
    else:
        if prec is None:
            M, t_prec, meta = make_preconditioner(A, precond_kind)
            reusable = {"M": M, **meta}
        else:
            M = prec["M"]
            t_prec = 0.0
            reusable = prec

        def cb(arg):
            niter[0] += 1
            if niter[0] % 25 == 0:
                log(f"    {solver} iter {niter[0]}")
                if np.isscalar(arg) or isinstance(arg, (float, np.floating)):
                    residuals.append(float(arg))
                    log(f"      pr_norm={float(arg):.3e}")
                else:
                    residuals.append(
                        float(np.linalg.norm(A @ arg - b) / (np.linalg.norm(b) + 1e-30))
                    )

        t1 = time.perf_counter()
        if solver == "gmres":
            x, info = gmres(
                A, b, M=M, rtol=tol, atol=0, maxiter=maxiter, restart=restart,
                callback=cb, callback_type="pr_norm",
            )
        else:
            x, info = bicgstab(A, b, M=M, rtol=tol, atol=0, maxiter=maxiter, callback=cb)
        t_solve = time.perf_counter() - t1

    resid = float(np.linalg.norm(A @ x - b) / (np.linalg.norm(b) + 1e-30))
    return {
        "x": x,
        "info": int(info) if info is not None else 0,
        "residual": resid,
        "iterations": int(niter[0]),
        "residual_history": residuals,
        "timings_s": {
            "preconditioner": t_prec,
            "solve": t_solve,
            "total_solver": t_prec + t_solve,
        },
        "prec": reusable,
        "precond_kind": reusable.get("kind", precond_kind),
        "solver": solver,
    }


def run_case(
    level: str,
    points_per_cm: float,
    bias_a: float = 0.0,
    second_rhs: bool = False,
    precond: str = "jacobi",
    maxiter: int = 500,
    restart: int = 40,
    tol: float = 1e-5,
    solver: str = "bicgstab",
) -> Dict[str, Any]:
    res = meep_resolution_from_points_per_cm(points_per_cm, a_m=sc.a)
    grid = make_grid(res)
    k0 = 2 * np.pi * sc.fs_a
    ten = eps_tensor_components(sc.fs_a, sc.fp_a, sc.gamma_a, bias_a)
    log(
        f"=== {level} ppc={points_per_cm} res={res} N={grid.Nx}x{grid.Ny}={grid.Nx*grid.Ny} "
        f"bias={bias_a} solver={solver}/{precond} ==="
    )

    t_mat0 = time.perf_counter()
    rho_xx, rho_xy, rho_yx, rho_yy, pec = material_maps(level, grid, bias_a=bias_a)
    t_mat = time.perf_counter() - t_mat0
    log(f"  material maps: {t_mat:.2f}s  pec_frac={pec.mean():.4f}  rss={_rss_gib():.2f} GiB")

    t_op0 = time.perf_counter()
    A = build_operator(grid, rho_xx, rho_xy, rho_yx, rho_yy, pec, k0)
    t_op = time.perf_counter() - t_op0
    log(f"  operator: {t_op:.2f}s  nnz={A.nnz}  rss={_rss_gib():.2f} GiB")

    # Matvec throughput for projection
    v = np.random.randn(A.shape[0]) + 1j * np.random.randn(A.shape[0])
    v = v.astype(np.complex128)
    for _ in range(3):
        _ = A @ v
    tmv0 = time.perf_counter()
    n_mv = 20
    for _ in range(n_mv):
        _ = A @ v
    matvec_s = (time.perf_counter() - tmv0) / n_mv
    log(f"  matvec: {matvec_s*1e3:.2f} ms  ({1.0/matvec_s:.1f}/s)")

    b0 = inject_mode_source(grid, 0, level=level)
    sol = solve_system(
        A, b0, precond_kind=precond, maxiter=maxiter, restart=restart, tol=tol, solver=solver
    )
    log(
        f"  solve1: prec={sol['timings_s']['preconditioner']:.2f}s "
        f"krylov={sol['timings_s']['solve']:.2f}s iters={sol['iterations']} "
        f"info={sol['info']} resid={sol['residual']:.2e} rss={_rss_gib():.2f} GiB"
    )
    u = sol["x"]

    n_ports = 1 if level == "level1" else 6
    proxies = {p: float(np.real(flux_proxy_through_port(grid, u, p, pec=pec))) for p in range(n_ports)}
    p0 = proxies.get(0, 1.0) + 1e-30
    norms = {str(p): float(proxies[p] / p0) for p in proxies}

    out: Dict[str, Any] = {
        "level": level,
        "points_per_cm": points_per_cm,
        "res": int(res),
        "Nx": grid.Nx,
        "Ny": grid.Ny,
        "unknowns": grid.Nx * grid.Ny,
        "nnz": int(A.nnz),
        "bias_a": bias_a,
        "precond": precond,
        "krylov_solver": solver,
        "matvec_s": matvec_s,
        "tensor": {
            "eps_perp_re": float(np.real(ten["eps_perp"])),
            "eps_perp_im": float(np.imag(ten["eps_perp"])),
            "eta_re": float(np.real(ten["eta"])),
            "eta_im": float(np.imag(ten["eta"])),
            "eps_xy": "-1j*eta",
            "eps_yx": "+1j*eta",
        },
        "timings_s": {
            "material": t_mat,
            "operator": t_op,
            "preconditioner": sol["timings_s"]["preconditioner"],
            "solve": sol["timings_s"]["solve"],
            "wall_total_first_rhs": t_mat + t_op + sol["timings_s"]["total_solver"],
        },
        "solver": {
            "info": sol["info"],
            "residual": sol["residual"],
            "iterations": sol["iterations"],
            "residual_history": sol["residual_history"],
        },
        "normalized_power_proxy_by_port": norms,
        "peak_rss_gib": _rss_gib(),
        "pec_fraction": float(np.mean(pec)),
    }

    if second_rhs and level != "level1":
        b1 = inject_mode_source(grid, 1, level=level)
        t2 = time.perf_counter()
        sol2 = solve_system(
            A, b1, prec=sol["prec"], maxiter=maxiter, restart=restart, tol=tol, solver=solver
        )
        t_second = time.perf_counter() - t2
        log(
            f"  solve2 (reuse prec): krylov={sol2['timings_s']['solve']:.2f}s "
            f"iters={sol2['iterations']} info={sol2['info']} resid={sol2['residual']:.2e}"
        )
        out["second_rhs"] = {
            "port": 1,
            "timings_s": sol2["timings_s"],
            "wall_s": t_second,
            "iterations": sol2["iterations"],
            "info": sol2["info"],
            "residual": sol2["residual"],
        }
        out["timings_s"]["estimated_six_source_s"] = (
            out["timings_s"]["wall_total_first_rhs"] + 5.0 * sol2["timings_s"]["solve"]
        )
        # Clarify: for direct, preconditioner==factorization; solve==substitution
        out["timings_s"]["factorization_s"] = sol["timings_s"]["preconditioner"]
        out["timings_s"]["substitution_per_rhs_s"] = sol2["timings_s"]["solve"]

    # Level1 qualitative: mid-guide |Hz| energy left vs right of source
    if level == "level1":
        U = np.abs(u.reshape(grid.Ny, grid.Nx)) ** 2
        mid_y = slice(grid.Ny // 2 - grid.res, grid.Ny // 2 + grid.res)
        left = float(U[mid_y, : grid.Nx // 2].sum())
        right = float(U[mid_y, grid.Nx // 2 :].sum())
        out["level1_energy_left"] = left
        out["level1_energy_right"] = right
        out["level1_right_over_left"] = right / (left + 1e-30)
        log(f"  level1 energy R/L={out['level1_right_over_left']:.3f}")

    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", choices=["level1", "level2", "level3"], required=True)
    ap.add_argument("--points-per-cm", type=float, default=25.0)
    ap.add_argument("--bias-a", type=float, default=0.0)
    ap.add_argument("--second-rhs", action="store_true")
    ap.add_argument("--precond", choices=["jacobi", "ilu0", "shifted"], default="shifted")
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--restart", type=int, default=40)
    ap.add_argument("--tol", type=float, default=1e-5)
    ap.add_argument("--solver", choices=["bicgstab", "gmres", "direct"], default="bicgstab")
    ap.add_argument("--json-out", type=str, required=True)
    args = ap.parse_args()
    result = run_case(
        args.level,
        args.points_per_cm,
        bias_a=args.bias_a,
        second_rhs=args.second_rhs,
        precond=args.precond,
        maxiter=args.maxiter,
        restart=args.restart,
        tol=args.tol,
        solver=args.solver,
    )
    Path(args.json_out).write_text(json.dumps(result, indent=2) + "\n")
    # print summary without huge history duplication
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
