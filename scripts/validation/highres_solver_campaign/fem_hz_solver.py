#!/usr/bin/env python3
"""
CG1 FEM scalar-Hz solver on locally refined Delaunay meshes.

Weak form for B=0 / isotropic cells:
  ∫ (1/ε) ∇Hz · ∇v dx  -  k0² ∫ Hz v dx  =  ⟨f, v⟩
with complex-stretched mass/stiffness in PML strips.
PEC walls: natural Neumann for TE-to-z (no essential BC).
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import gmres, splu, LinearOperator

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / ".conda_envs" / "fdfd_solver" / "lib" / "python3.13" / "site-packages"))

import sixport_common as sc  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402
from PMMCirculatorInverse import PMMI  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "highres_solver_campaign" / "phase5_fem"


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)


def log(msg: str) -> None:
    print(msg, flush=True)


def bulb_centers() -> np.ndarray:
    pmm_tmp = PMMI(a=sc.a, res=32, nx=14, ny=12, dpml=sc.dpml, B=np.zeros(3))
    pmm_tmp.Rod_Array_Hexagon_train(
        xy_cen=np.array([7.0, 6.0]), side_dim=6, r=sc.r_plasma, d=sc.d_exp, bulbs=False, uniform=True
    )
    locs = np.asarray(pmm_tmp.train_elem_locs, dtype=float)
    return locs - np.array([7.0, 6.0]) + np.array([sc.nx_ports / 2, sc.ny_ports / 2])


def load_mesh(grade: str, device_mode: str = "full") -> Tuple[np.ndarray, np.ndarray]:
    path = OUT / f"mesh_{grade}_{device_mode}.npz"
    data = np.load(path)
    return data["points"], data["triangles"]


def element_materials(points: np.ndarray, tris: np.ndarray, bias_a: float = 0.0) -> np.ndarray:
    """Return complex ε per triangle (isotropic; gyrotropic handled via ε_perp for B!=0 scalar approx)."""
    cents = points[tris].mean(axis=1)
    eps = np.ones(len(tris), dtype=np.complex128)
    eps_p, _eta = gyrotropic_drude_eps_eta(sc.fs_a, sc.fp_a, sc.gamma_a, bias_a)
    # For full anisotropic FEM we'd need tensor assembly; at B=0 eps_p is enough.
    # For B!=0 smoke, use ε_perp as isotropic proxy ONLY if eta~0; else use det-based scalar
    # effective: for Faraday, scalar reduction uses full ρ — here for B!=0 use ε_eff = eps_perp
    # (documented limitation) OR ρ_xx = eps_perp/det for isotropic-like.
    if abs(bias_a) > 0:
        # Use det-consistent isotropic proxy for magnitude: ε_eff = det/ε_perp = ε_perp - η²/ε_perp
        # Better: store ρ = 1/ε_perp for diagonal-dominant approximation — NOT exact.
        # Exact anisotropic FEM is more work; for B!=0 we still assemble with ε_perp and note limitation.
        pass
    eps_q = 3.8 + 0j
    bulbs = bulb_centers()
    r_q = sc.r_bulb_outer
    r_vac = sc.r_bulb_inner
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    for cx, cy in bulbs:
        d2 = (cents[:, 0] - cx) ** 2 + (cents[:, 1] - cy) ** 2
        eps[d2 <= r_q**2] = eps_q
        eps[(d2 <= r_q**2) & (d2 > r_vac**2)] = eps_q
        # vacuum liner
        eps[(d2 <= r_vac**2) & (d2 > r_p**2)] = 1.0
        eps[d2 <= r_p**2] = eps_p
    return eps


def pml_factor_at(xy: np.ndarray) -> np.ndarray:
    """Complex stretch sx*sy approx as scalar factor on mass/stiffness in PML strips."""
    nx, ny = sc.nx_ports, sc.ny_ports
    dp = sc.dpml_ports
    omega = 2 * np.pi * sc.fs_a
    x, y = xy[:, 0], xy[:, 1]
    sigx = np.zeros_like(x)
    sigy = np.zeros_like(y)
    # left/right
    m = x < dp
    sigx[m] = 2.0 * ((dp - x[m]) / dp) ** 2
    m = x > nx - dp
    sigx[m] = 2.0 * ((x[m] - (nx - dp)) / dp) ** 2
    m = y < dp
    sigy[m] = 2.0 * ((dp - y[m]) / dp) ** 2
    m = y > ny - dp
    sigy[m] = 2.0 * ((y[m] - (ny - dp)) / dp) ** 2
    sx = 1.0 + 1j * sigx / omega
    sy = 1.0 + 1j * sigy / omega
    return sx, sy


def pec_node_mask(points: np.ndarray) -> np.ndarray:
    """Nodes within ~half wall thickness of horn polylines → PEC (Dirichlet Hz unused / identity)."""
    from sixport_common import full_horns

    pec = np.zeros(len(points), dtype=bool)
    nx, ny = sc.nx_ports, sc.ny_ports
    rad = 0.55 * sc.wall_thickness
    rad2 = rad * rad
    for horn in full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], dtype=float) + np.array([nx / 2, ny / 2])
            for k in range(len(poly) - 1):
                p0, p1 = poly[k], poly[k + 1]
                L = float(np.linalg.norm(p1 - p0))
                n = max(2, int(L / (0.5 * rad)))
                for t in np.linspace(0, 1, n):
                    c = p0 * (1 - t) + p1 * t
                    d2 = (points[:, 0] - c[0]) ** 2 + (points[:, 1] - c[1]) ** 2
                    pec |= d2 <= rad2
    return pec


def assemble_hz_fem(
    points: np.ndarray,
    tris: np.ndarray,
    eps_tri: np.ndarray,
    k0: float,
    pec_nodes: Optional[np.ndarray] = None,
) -> sparse.csr_matrix:
    """Vectorized P1 assembly of stretched Helmholtz with isotropic 1/ε."""
    n = len(points)
    cents = points[tris].mean(axis=1)
    sx, sy = pml_factor_at(cents)
    inv_eps = 1.0 / eps_tri

    x = points[tris, 0]
    y = points[tris, 1]
    twice_area = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    area = 0.5 * np.abs(twice_area)
    good = area > 1e-18
    b = np.stack(
        [y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1
    ) / twice_area[:, None]
    c = np.stack(
        [x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1
    ) / twice_area[:, None]

    ax = (sy / sx) * inv_eps
    ay = (sx / sy) * inv_eps
    mass = (k0**2) * sx * sy

    rows = []
    cols = []
    data = []
    for i in range(3):
        for j in range(3):
            Ke = area * (ax * b[:, i] * b[:, j] + ay * c[:, i] * c[:, j])
            Me = mass * area * (1.0 / 6.0 if i == j else 1.0 / 12.0)
            val = np.where(good, Ke - Me, 0.0)
            rows.append(tris[:, i])
            cols.append(tris[:, j])
            data.append(val)

    r = np.concatenate(rows)
    cidx = np.concatenate(cols)
    v = np.concatenate(data)
    A = sparse.coo_matrix((v, (r, cidx)), shape=(n, n)).tocsr()
    if pec_nodes is not None and np.any(pec_nodes):
        # Identity rows/cols for PEC nodes (TE Neumann approx via removing DOF)
        pec_idx = np.flatnonzero(pec_nodes)
        A = A.tolil()
        for i in pec_idx:
            A.rows[i] = [i]
            A.data[i] = [1.0 + 0j]
        A = A.tocsr()
    return A


def inject_source(points: np.ndarray, port: int = 0, level: str = "level3") -> np.ndarray:
    b = np.zeros(len(points), dtype=np.complex128)
    if level == "level1":
        cx, cy = sc.nx_ports / 2, sc.ny_ports / 2
        span = 0.96 * sc.clear_width
        s_line = np.linspace(-span / 2, span / 2, 41)
        for s in s_line:
            xy = np.array([cx, cy + s])
            d = np.linalg.norm(points - xy, axis=1)
            i = int(np.argmin(d))
            amp = max(np.cos(np.pi * s / span), 0.0)
            b[i] += amp
    else:
        from sixport_common import set_geometry_context, horn_for_port, effective_port_dir

        set_geometry_context(res=50, horn_walls="prism")
        horn = horn_for_port(port, 50)
        center = np.asarray(horn["source_center"], dtype=float) + np.array(
            [sc.nx_ports / 2, sc.ny_ports / 2]
        )
        u = np.asarray(effective_port_dir(port), dtype=float)
        u /= np.linalg.norm(u)
        tang = np.array([-u[1], u[0]])
        span = 0.96 * sc.clear_width
        for s in np.linspace(-span / 2, span / 2, 41):
            xy = center + s * tang
            d = np.linalg.norm(points - xy, axis=1)
            i = int(np.argmin(d))
            amp = max(np.cos(np.pi * s / span), 0.0)
            b[i] += amp
    nrm = np.linalg.norm(b)
    if nrm > 0:
        b /= nrm
    else:
        raise RuntimeError("source injection produced zero RHS")
    return b


def port_proxy(points: np.ndarray, u: np.ndarray, port: int) -> float:
    from sixport_common import set_geometry_context, monitor_center_for_port, effective_port_dir

    set_geometry_context(res=50, horn_walls="prism")
    center = np.asarray(monitor_center_for_port(port, 50), dtype=float) + np.array(
        [sc.nx_ports / 2, sc.ny_ports / 2]
    )
    uhat = np.asarray(effective_port_dir(port), dtype=float)
    uhat /= np.linalg.norm(uhat)
    tang = np.array([-uhat[1], uhat[0]])
    span = 0.96 * sc.clear_width
    acc = 0.0
    wsum = 0.0
    for s in np.linspace(-span / 2, span / 2, 41):
        xy = center + s * tang
        d = np.linalg.norm(points - xy, axis=1)
        i = int(np.argmin(d))
        w = max(np.cos(np.pi * s / span), 0.0)
        acc += w * float(np.abs(u[i]) ** 2)
        wsum += w
    return float(acc / (wsum + 1e-30))


def solve_fem(
    grade: str,
    level: str = "level3",
    bias_a: float = 0.0,
    method: str = "direct",
    second_rhs: bool = False,
    tol: float = 1e-4,
    maxiter: int = 300,
) -> Dict[str, Any]:
    device_mode = "horns_only" if level == "level2" else ("full" if level == "level3" else "full")
    if level == "level1":
        # reuse full-domain mesh but vacuum materials except guide — for smoke use level3 mesh
        device_mode = "full"

    log(f"=== FEM {grade} level={level} method={method} bias={bias_a} ===")
    t0 = time.perf_counter()
    points, tris = load_mesh(grade, device_mode if level != "level1" else "full")
    if level == "level1":
        # vacuum everywhere
        eps = np.ones(len(tris), dtype=np.complex128)
    elif level == "level2":
        # try horns_only mesh if present else full vacuum+horns geometry only (eps=1)
        try:
            points, tris = load_mesh(grade, "horns_only")
        except Exception:
            pass
        eps = np.ones(len(tris), dtype=np.complex128)
    else:
        eps = element_materials(points, tris, bias_a=bias_a)
    t_mat = time.perf_counter() - t0

    k0 = 2 * np.pi * sc.fs_a
    pec = pec_node_mask(points) if level != "level1" else np.zeros(len(points), dtype=bool)
    log(f"  pec_frac={pec.mean():.4f}")
    t1 = time.perf_counter()
    A = assemble_hz_fem(points, tris, eps, k0, pec_nodes=pec)
    t_op = time.perf_counter() - t1
    log(f"  nodes={len(points)} tris={len(tris)} nnz={A.nnz} mat={t_mat:.2f}s asm={t_op:.2f}s rss={rss_gib():.2f}")

    b0 = inject_source(points, 0, level=level)
    if level != "level1":
        b0[pec] = 0.0
    t_setup = 0.0
    t_solve = 0.0
    reusable = {}
    niter = -1
    info = 0

    if method == "direct":
        t2 = time.perf_counter()
        lu = splu(A.tocsc())
        t_setup = time.perf_counter() - t2
        t3 = time.perf_counter()
        x = lu.solve(b0)
        t_solve = time.perf_counter() - t3
        reusable = {"lu": lu}
        log(f"  factor={t_setup:.2f}s sub={t_solve:.2f}s")
    else:
        # Jacobi GMRES fallback
        diag = A.diagonal()
        invd = 1.0 / np.where(np.abs(diag) > 1e-30, diag, 1.0)
        M = LinearOperator(A.shape, matvec=lambda v: invd * v, dtype=np.complex128)
        hist = []

        def cb(pr):
            hist.append(float(pr))

        t2 = time.perf_counter()
        x, info = gmres(A, b0, M=M, rtol=tol, atol=0, maxiter=maxiter, restart=50,
                        callback=cb, callback_type="pr_norm")
        t_solve = time.perf_counter() - t2
        niter = len(hist)
        log(f"  gmres iters={niter} info={info} time={t_solve:.2f}s")

    resid = float(np.linalg.norm(A @ x - b0) / (np.linalg.norm(b0) + 1e-30))
    n_ports = 1 if level == "level1" else 6
    proxies = {str(p): port_proxy(points, x, p) for p in range(n_ports)}
    p0 = proxies.get("0", 1.0) + 1e-30
    norms = {k: v / p0 for k, v in proxies.items()}

    out: Dict[str, Any] = {
        "grade": grade,
        "level": level,
        "method": method,
        "bias_a": bias_a,
        "n_nodes": len(points),
        "n_triangles": len(tris),
        "nnz": int(A.nnz),
        "true_residual": resid,
        "iterations": niter,
        "info": int(info),
        "timings_s": {
            "material": t_mat,
            "assemble": t_op,
            "setup": t_setup,
            "solve": t_solve,
            "wall_first_rhs": t_mat + t_op + t_setup + t_solve,
        },
        "normalized_power_proxy": norms,
        "rss_gib": rss_gib(),
        "converged": resid < 1e-3,
    }

    if second_rhs and "lu" in reusable:
        b1 = inject_source(points, 1, level=level)
        t4 = time.perf_counter()
        x1 = reusable["lu"].solve(b1)
        t_sub = time.perf_counter() - t4
        out["second_rhs"] = {
            "wall_s": t_sub,
            "true_residual": float(np.linalg.norm(A @ x1 - b1) / (np.linalg.norm(b1) + 1e-30)),
        }
        out["timings_s"]["estimated_six_source_s"] = (
            out["timings_s"]["wall_first_rhs"] + 5.0 * t_sub
        )
        log(f"  second RHS sub={t_sub:.3f}s resid={out['second_rhs']['true_residual']:.2e}")

    if level == "level1":
        # energy left/right of center
        left = float(np.sum(np.abs(x[points[:, 0] < sc.nx_ports / 2]) ** 2))
        right = float(np.sum(np.abs(x[points[:, 0] >= sc.nx_ports / 2]) ** 2))
        out["level1_R_over_L"] = right / (left + 1e-30)
        log(f"  level1 R/L={out['level1_R_over_L']:.3f}")

    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="FEM-L")
    ap.add_argument("--level", default="level3", choices=["level1", "level2", "level3"])
    ap.add_argument("--method", default="direct", choices=["direct", "gmres"])
    ap.add_argument("--bias-a", type=float, default=0.0)
    ap.add_argument("--second-rhs", action="store_true")
    ap.add_argument("--json-out", required=True)
    args = ap.parse_args()
    result = solve_fem(
        args.grade, level=args.level, bias_a=args.bias_a, method=args.method, second_rhs=args.second_rhs
    )
    Path(args.json_out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if result.get("converged") else 1


if __name__ == "__main__":
    raise SystemExit(main())
