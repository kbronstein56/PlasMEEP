#!/usr/bin/env python3
"""
High-res solver campaign: uniform FDFD with shifted-Laplacian + geometric MG,
reusing the validated scalar-Hz operator from fdfd_feasibility.

Isolated campaign code — does not modify production PlasMEEP.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, gmres, splu

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "outputs" / "validation" / "fdfd_feasibility"))
sys.path.insert(0, str(ROOT / "scripts" / "validation"))
sys.path.insert(0, str(ROOT / "scripts"))

import fdfd_hz_prototype as proto  # noqa: E402
import sixport_common as sc  # noqa: E402

OUT = Path(__file__).resolve().parents[3] / "outputs" / "validation" / "highres_solver_campaign"


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)


def log(msg: str) -> None:
    print(msg, flush=True)


def true_residual(A: sparse.csr_matrix, x: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(A @ x - b) / (np.linalg.norm(b) + 1e-30))


def build_shifted_operator(
    grid: proto.Grid,
    rho_xx, rho_xy, rho_yx, rho_yy, pec,
    k0: float,
    beta: float,
) -> sparse.csr_matrix:
    """M = same stencil as A but k0^2 -> (1 + 1j*beta)*k0^2 (complex-shifted Helmholtz)."""
    k0_eff = k0 * np.sqrt(1.0 + 1j * beta)
    # build_operator uses k0**2 in mass term; pass complex via abs? No — modify mass.
    # Easiest: build A then replace diagonal mass contribution.
    # Instead rebuild with complex k0^2 by monkey-patching: call with k0 replaced.
    # Our operator uses (k0**2)*Sx*Sy on diagonal — use complex k0_sq.
    return _build_operator_k0sq(grid, rho_xx, rho_xy, rho_yx, rho_yy, pec, (1.0 + 1j * beta) * (k0**2))


def _build_operator_k0sq(grid, rho_xx, rho_xy, rho_yx, rho_yy, pec, k0_sq: complex) -> sparse.csr_matrix:
    """Copy of proto.build_operator with arbitrary complex k0^2 mass coefficient."""
    Nx, Ny, dx = grid.Nx, grid.Ny, grid.dx
    omega = 2 * np.pi * sc.fs_a
    dpml_cells = int(round(sc.dpml_ports * grid.res))
    sigx = proto.pml_sigma_1d(Nx, dpml_cells)
    sigy = proto.pml_sigma_1d(Ny, dpml_cells)
    sx = 1.0 + 1j * sigx / max(omega, 1e-30)
    sy = 1.0 + 1j * sigy / max(omega, 1e-30)
    Sx = np.broadcast_to(sx[None, :], (Ny, Nx)).copy()
    Sy = np.broadcast_to(sy[:, None], (Ny, Nx)).copy()

    rho_xx = rho_xx.copy(); rho_xy = rho_xy.copy(); rho_yx = rho_yx.copy(); rho_yy = rho_yy.copy()
    pec = pec.astype(bool)
    for arr in (rho_xx, rho_xy, rho_yx, rho_yy):
        arr[pec] = 0.0

    inv_dx2 = 1.0 / (dx * dx)
    inv_dx2_4 = 1.0 / (4.0 * dx * dx)

    def harm(a, b):
        return 2 * a * b / (a + b + 1e-30)

    rxx_xp = np.zeros_like(rho_xx); rxx_xm = np.zeros_like(rho_xx)
    ryy_yp = np.zeros_like(rho_yy); ryy_ym = np.zeros_like(rho_yy)
    rxx_xp[:, :-1] = harm(rho_xx[:, :-1], rho_xx[:, 1:])
    rxx_xm[:, 1:] = harm(rho_xx[:, 1:], rho_xx[:, :-1])
    ryy_yp[:-1, :] = harm(rho_yy[:-1, :], rho_yy[1:, :])
    ryy_ym[1:, :] = harm(rho_yy[1:, :], rho_yy[:-1, :])

    rxy_xp = np.zeros_like(rho_xy); rxy_xm = np.zeros_like(rho_xy)
    ryx_yp = np.zeros_like(rho_yx); ryx_ym = np.zeros_like(rho_yx)
    rxy_xp[:, :-1] = 0.5 * (rho_xy[:, :-1] + rho_xy[:, 1:])
    rxy_xm[:, 1:] = 0.5 * (rho_xy[:, 1:] + rho_xy[:, :-1])
    ryx_yp[:-1, :] = 0.5 * (rho_yx[:-1, :] + rho_yx[1:, :])
    ryx_ym[1:, :] = 0.5 * (rho_yx[1:, :] + rho_yx[:-1, :])

    cx = inv_dx2 / Sx
    cy = inv_dx2 / Sy
    diag = -cx * (rxx_xp + rxx_xm) - cy * (ryy_yp + ryy_ym) + k0_sq * Sx * Sy

    data, rows, cols = [], [], []

    def add_diag_array(coeff: np.ndarray, diy: int, dix: int):
        iy0 = max(0, -diy); iy1 = Ny - max(0, diy)
        ix0 = max(0, -dix); ix1 = Nx - max(0, dix)
        c = coeff[iy0:iy1, ix0:ix1].ravel()
        if c.size == 0:
            return
        IY, IX = np.mgrid[iy0:iy1, ix0:ix1]
        r = (IY * Nx + IX).ravel()
        cc = ((IY + diy) * Nx + (IX + dix)).ravel()
        mask = (np.abs(c) > 0) & (~pec[iy0:iy1, ix0:ix1].ravel())
        if not np.any(mask):
            return
        rows.append(r[mask]); cols.append(cc[mask]); data.append(c[mask])

    add_diag_array(diag, 0, 0)
    add_diag_array(cx * rxx_xp, 0, +1)
    add_diag_array(cx * rxx_xm, 0, -1)
    add_diag_array(cy * ryy_yp, +1, 0)
    add_diag_array(cy * ryy_ym, -1, 0)
    if np.max(np.abs(rho_xy)) + np.max(np.abs(rho_yx)) > 0:
        cxy = inv_dx2_4 / Sx
        cyx = inv_dx2_4 / Sy
        add_diag_array(+cxy * rxy_xp, +1, +1)
        add_diag_array(-cxy * rxy_xp, -1, +1)
        add_diag_array(-cxy * rxy_xm, +1, -1)
        add_diag_array(+cxy * rxy_xm, -1, -1)
        add_diag_array(+cyx * ryx_yp, +1, +1)
        add_diag_array(-cyx * ryx_yp, +1, -1)
        add_diag_array(-cyx * ryx_ym, -1, +1)
        add_diag_array(+cyx * ryx_ym, -1, -1)

    pec_idx = np.flatnonzero(pec.ravel())
    if pec_idx.size:
        rows.append(pec_idx); cols.append(pec_idx)
        data.append(np.ones(pec_idx.size, dtype=np.complex128))

    r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(data)
    return sparse.coo_matrix((v, (r, c)), shape=(Nx * Ny, Nx * Ny)).tocsr()


def restrict_inject(fine: np.ndarray, Nx: int, Ny: int) -> Tuple[np.ndarray, int, int]:
    """Full-weighting restriction to even grid (Ny must be even-ish)."""
    F = fine.reshape(Ny, Nx)
    # trim to even
    Ny2, Nx2 = (Ny // 2) * 2, (Nx // 2) * 2
    F = F[:Ny2, :Nx2]
    # 2x2 average
    C = 0.25 * (F[0::2, 0::2] + F[0::2, 1::2] + F[1::2, 0::2] + F[1::2, 1::2])
    return C.ravel(), Nx2 // 2, Ny2 // 2


def prolong_bilinear(coarse: np.ndarray, nxc: int, nyc: int, Nxf: int, Nyf: int) -> np.ndarray:
    C = coarse.reshape(nyc, nxc)
    F = np.zeros((Nyf, Nxf), dtype=np.complex128)
    # nearest + average for odd rows/cols
    F[0:2 * nyc:2, 0:2 * nxc:2] = C
    F[0:2 * nyc:2, 1:2 * nxc:2] = C
    F[1:2 * nyc:2, 0:2 * nxc:2] = C
    F[1:2 * nyc:2, 1:2 * nxc:2] = C
    if Nyf > 2 * nyc:
        F[2 * nyc :, :] = F[2 * nyc - 1 : 2 * nyc, :]
    if Nxf > 2 * nxc:
        F[:, 2 * nxc :] = F[:, 2 * nxc - 1 : 2 * nxc]
    return F[:Nyf, :Nxf].ravel()


class GeometricMGPC:
    """
    Simple geometric V-cycle PC for complex-shifted Helmholtz on uniform grids.
    Coarsest level: SuperLU. Smoother: damped Jacobi.
    """

    def __init__(self, M: sparse.csr_matrix, Nx: int, Ny: int, n_levels: int = 4, damp: float = 0.7):
        self.levels: List[Dict[str, Any]] = []
        A = M.tocsr()
        nx, ny = Nx, Ny
        for lev in range(n_levels):
            diag = A.diagonal()
            invd = 1.0 / np.where(np.abs(diag) > 1e-30, diag, 1.0)
            entry = {"A": A, "Nx": nx, "Ny": ny, "invd": invd, "damp": damp}
            if lev == n_levels - 1 or nx < 40 or ny < 40:
                entry["lu"] = splu(A.tocsc())
                self.levels.append(entry)
                break
            self.levels.append(entry)
            # Galerkin coarse: build by injecting stencil approx via restriction of residual operator
            # Practical shortcut: form coarse matrix by restricting identity columns in batches — too slow.
            # Instead: rebuild is not available; use rediscretization approximation via 2x subsample of
            # diagonal + nearest — crude. Better: use scipy to form P^T A P with sparse prolong.
            P = self._prolong_matrix(nx // 2, ny // 2, nx, ny)
            R = P.T.tocsr()  # injection-like (transpose of bilinear)
            A = (R @ A @ P).tocsr()
            nx, ny = nx // 2, ny // 2
        log(f"  MG levels={len(self.levels)} coarsest={self.levels[-1]['Nx']}x{self.levels[-1]['Ny']}")

    @staticmethod
    def _prolong_matrix(nxc, nyc, nxf, nyf) -> sparse.csr_matrix:
        # Each coarse node maps to up to 4 fine nodes with weight 1 (piecewise constant prolong)
        rows, cols, data = [], [], []
        for jy in range(nyc):
            for jx in range(nxc):
                c = jy * nxc + jx
                for dy in (0, 1):
                    for dx in (0, 1):
                        iy, ix = 2 * jy + dy, 2 * jx + dx
                        if iy < nyf and ix < nxf:
                            rows.append(iy * nxf + ix)
                            cols.append(c)
                            data.append(1.0)
        return sparse.coo_matrix((data, (rows, cols)), shape=(nxf * nyf, nxc * nyc)).tocsr()

    def _smooth(self, lev: int, x: np.ndarray, b: np.ndarray, sweeps: int = 2) -> np.ndarray:
        A = self.levels[lev]["A"]
        invd = self.levels[lev]["invd"]
        damp = self.levels[lev]["damp"]
        for _ in range(sweeps):
            r = b - A @ x
            x = x + damp * (invd * r)
        return x

    def _vcycle(self, lev: int, b: np.ndarray, x: Optional[np.ndarray] = None) -> np.ndarray:
        L = self.levels[lev]
        n = L["A"].shape[0]
        if x is None:
            x = np.zeros(n, dtype=np.complex128)
        if "lu" in L:
            return L["lu"].solve(b)
        x = self._smooth(lev, x, b, sweeps=2)
        r = b - L["A"] @ x
        # restrict
        nxc, nyc = self.levels[lev + 1]["Nx"], self.levels[lev + 1]["Ny"]
        P = self._prolong_matrix(nxc, nyc, L["Nx"], L["Ny"])
        R = 0.25 * P.T.tocsr()
        rc = R @ r
        ec = self._vcycle(lev + 1, rc)
        x = x + P @ ec
        x = self._smooth(lev, x, b, sweeps=2)
        return x

    def apply(self, b: np.ndarray) -> np.ndarray:
        return self._vcycle(0, b)


def try_pyamg_pc(M: sparse.csr_matrix):
    """Attempt PyAMG on shifted operator; return LinearOperator or None."""
    try:
        import pyamg
    except ImportError:
        return None, "pyamg missing"
    try:
        # PyAMG classical SA often assumes real SPD; try anyway on complex
        ml = pyamg.smoothed_aggregation_solver(M.tocsr(), max_levels=6, max_coarse=500)
        def mv(x):
            return ml.solve(x, tol=1e-4, maxiter=5)
        return LinearOperator(M.shape, matvec=mv, dtype=np.complex128), "pyamg_sa"
    except Exception as e:
        return None, f"pyamg_failed:{e}"


def assemble_case(level: str, ppc: float, bias_a: float = 0.0):
    from physical_units import meep_resolution_from_points_per_cm
    res = meep_resolution_from_points_per_cm(ppc, a_m=sc.a)
    grid = proto.make_grid(res)
    k0 = 2 * np.pi * sc.fs_a
    t0 = time.perf_counter()
    maps = proto.material_maps(level, grid, bias_a=bias_a)
    t_mat = time.perf_counter() - t0
    t1 = time.perf_counter()
    A = proto.build_operator(grid, *maps, k0)
    t_op = time.perf_counter() - t1
    b = proto.inject_mode_source(grid, 0, level=level)
    return {
        "grid": grid, "k0": k0, "maps": maps, "A": A, "b": b,
        "timings": {"material": t_mat, "operator": t_op},
        "res": int(res), "ppc": ppc, "level": level, "bias_a": bias_a,
    }


def solve_with_method(
    case: Dict[str, Any],
    method: str,
    beta: float = 0.5,
    tol: float = 1e-4,
    maxiter: int = 200,
    timeout_s: float = 600.0,
) -> Dict[str, Any]:
    A: sparse.csr_matrix = case["A"]
    b: np.ndarray = case["b"]
    grid = case["grid"]
    maps = case["maps"]
    k0 = case["k0"]
    t_setup0 = time.perf_counter()
    hist: List[float] = []
    niter = [0]
    reusable: Dict[str, Any] = {}

    if method == "direct":
        lu = splu(A.tocsc())
        t_setup = time.perf_counter() - t_setup0
        t1 = time.perf_counter()
        x = lu.solve(b)
        t_solve = time.perf_counter() - t1
        reusable = {"lu": lu}
        info = 0
    elif method.startswith("shifted_mg"):
        M = build_shifted_operator(grid, *[a.copy() for a in maps[:4]], maps[4].copy(), k0, beta)
        mg = GeometricMGPC(M, grid.Nx, grid.Ny, n_levels=5)
        t_setup = time.perf_counter() - t_setup0
        def pc(x):
            return mg.apply(x)
        PCop = LinearOperator(A.shape, matvec=pc, dtype=np.complex128)
        reusable = {"mg": mg, "beta": beta}

        def cb(pr):
            niter[0] += 1
            hist.append(float(pr))
            if time.perf_counter() - t_setup0 > timeout_s:
                raise TimeoutError("solve timeout")

        t1 = time.perf_counter()
        try:
            x, info = gmres(A, b, M=PCop, rtol=tol, atol=0, maxiter=maxiter, restart=40,
                            callback=cb, callback_type="pr_norm")
        except TimeoutError:
            x = np.zeros(A.shape[0], dtype=np.complex128)
            info = -99
        t_solve = time.perf_counter() - t1
    elif method.startswith("shifted_ilu"):
        from scipy.sparse.linalg import spilu
        M = build_shifted_operator(grid, *[a.copy() for a in maps[:4]], maps[4].copy(), k0, beta)
        ilu = spilu(M.tocsc(), drop_tol=1e-3, fill_factor=5.0)
        t_setup = time.perf_counter() - t_setup0
        PCop = LinearOperator(A.shape, matvec=ilu.solve, dtype=np.complex128)
        reusable = {"ilu": ilu, "beta": beta}

        def cb(pr):
            niter[0] += 1
            hist.append(float(pr))

        t1 = time.perf_counter()
        x, info = gmres(A, b, M=PCop, rtol=tol, atol=0, maxiter=maxiter, restart=40,
                        callback=cb, callback_type="pr_norm")
        t_solve = time.perf_counter() - t1
    elif method.startswith("shifted_pyamg"):
        M = build_shifted_operator(grid, *[a.copy() for a in maps[:4]], maps[4].copy(), k0, beta)
        PCop, tag = try_pyamg_pc(M)
        t_setup = time.perf_counter() - t_setup0
        if PCop is None:
            return {"method": method, "error": tag, "converged": False}
        reusable = {"tag": tag, "beta": beta}

        def cb(pr):
            niter[0] += 1
            hist.append(float(pr))

        t1 = time.perf_counter()
        x, info = gmres(A, b, M=PCop, rtol=tol, atol=0, maxiter=maxiter, restart=40,
                        callback=cb, callback_type="pr_norm")
        t_solve = time.perf_counter() - t1
    else:
        raise ValueError(method)

    resid = true_residual(A, x, b)
    # port proxies
    pec = maps[4]
    proxies = {}
    n_ports = 1 if case["level"] == "level1" else 6
    for p in range(n_ports):
        proxies[str(p)] = float(np.real(proto.flux_proxy_through_port(grid, x, p, pec=pec)))
    p0 = proxies.get("0", 1.0) + 1e-30
    norms = {k: v / p0 for k, v in proxies.items()}

    return {
        "method": method,
        "beta": beta,
        "info": int(info) if info is not None else 0,
        "iterations": int(niter[0]),
        "true_residual": resid,
        "pr_norm_history": hist[:: max(1, len(hist) // 50)] if hist else [],
        "converged": bool(resid < max(tol * 10, 1e-3) and info == 0),
        "timings_s": {
            "setup": t_setup,
            "solve": t_solve,
            "total": t_setup + t_solve + case["timings"]["material"] + case["timings"]["operator"],
            "material": case["timings"]["material"],
            "operator": case["timings"]["operator"],
        },
        "Nx": grid.Nx, "Ny": grid.Ny, "unknowns": grid.Nx * grid.Ny, "nnz": int(A.nnz),
        "rss_gib": rss_gib(),
        "normalized_power_proxy": norms,
        "reusable_keys": list(reusable.keys()),
        "_reusable": reusable,
        "_x": x,
        "_A": A,
        "_case": case,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", default="level1")
    ap.add_argument("--ppc", type=float, default=10.0)
    ap.add_argument("--method", default="shifted_mg")
    ap.add_argument("--beta", type=float, default=0.5)
    ap.add_argument("--bias-a", type=float, default=0.0)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--maxiter", type=int, default=200)
    ap.add_argument("--timeout", type=float, default=600.0)
    ap.add_argument("--json-out", required=True)
    ap.add_argument("--second-rhs", action="store_true")
    args = ap.parse_args()

    log(f"=== assemble {args.level} ppc={args.ppc} method={args.method} beta={args.beta} ===")
    case = assemble_case(args.level, args.ppc, bias_a=args.bias_a)
    log(f"  N={case['grid'].Nx}x{case['grid'].Ny} nnz={case['A'].nnz} mat={case['timings']['material']:.2f}s op={case['timings']['operator']:.2f}s")
    result = solve_with_method(case, args.method, beta=args.beta, tol=args.tol,
                               maxiter=args.maxiter, timeout_s=args.timeout)
    log(f"  true_resid={result.get('true_residual')} iters={result.get('iterations')} "
        f"setup={result.get('timings_s',{}).get('setup')} solve={result.get('timings_s',{}).get('solve')} "
        f"conv={result.get('converged')}")

    if args.second_rhs and result.get("converged") and "_reusable" in result:
        ru = result["_reusable"]
        A = result["_A"]
        b2 = proto.inject_mode_source(case["grid"], 1, level=args.level)
        t0 = time.perf_counter()
        if "lu" in ru:
            x2 = ru["lu"].solve(b2)
        elif "mg" in ru:
            # reuse MG PC
            def pc(x):
                return ru["mg"].apply(x)
            PCop = LinearOperator(A.shape, matvec=pc, dtype=np.complex128)
            x2, info2 = gmres(A, b2, M=PCop, rtol=args.tol, atol=0, maxiter=args.maxiter, restart=40)
        else:
            x2, info2 = None, -1
        t2 = time.perf_counter() - t0
        result["second_rhs"] = {
            "wall_s": t2,
            "true_residual": true_residual(A, x2, b2) if x2 is not None else None,
        }
        log(f"  second RHS {t2:.3f}s resid={result['second_rhs']['true_residual']}")

    # strip non-json
    out = {k: v for k, v in result.items() if not k.startswith("_")}
    Path(args.json_out).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
    return 0 if result.get("converged") else 1


if __name__ == "__main__":
    raise SystemExit(main())
