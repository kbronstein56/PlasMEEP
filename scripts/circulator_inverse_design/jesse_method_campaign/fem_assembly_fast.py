#!/usr/bin/env python3
"""Numerically equivalent FEM assembly accelerations for the A2 production path.

Optimizations (mathematics unchanged):
- Precompute frequency-independent mesh geometry / PML stretch / COO index pattern
- Cache frequency-independent device statics (mesh, masks, ports, exterior)
- Assemble dA only on orbit-support elements, or contract dA without building CSR
- Pin empty rows via diagonal COO add (no LIL Python loop)

Reference path remains available via a2_efficient_core.use_fast_kernel=False.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
from scipy import sparse

import sixport_common as sc
from boundary_ports import build_port_from_boundary, modal_c_vector
from diff_fem_core import per_rod_plasma_masks, rho_from_scales
from dtn_port_geometry import port_outer_center
from fem_validated_solver import pml_sx_sy
from full91_dtn_solve import exterior_beyond_ports_mask
from full91_geometry import bulb_centers_device
from full91_solve import load_mesh, material_masks
from modes import mode_constants
from multimode_opt_core import orbits_and_masks
from plasma_sensitivity import drho_ds


def assemble_port_operator_fast(A, points, ports, modes, n_prop: int = 2, n_evan: int = 2):
    """Exact same Robin updates as direct_dtn_s.assemble_port_operator.

    The reference path converts CSR→LIL→CSR once per (port, mode). Here we
    convert to LIL once, apply all updates, then convert back once.
    """
    A = A.tolil()
    for port in ports:
        if port["n_edges"] < 1:
            continue
        for md in modes[: n_prop + n_evan]:
            idx, c = modal_c_vector(points, port, md["m"])
            if len(idx) == 0:
                continue
            beta = complex(md["beta"])
            if abs(beta.imag) > 1e-12 and beta.imag > 0:
                beta_eff = 1j * beta.imag
            else:
                beta_eff = -1j * float(beta.real)
            block = beta_eff * np.outer(c, c)
            for i, p in enumerate(idx):
                for j, q in enumerate(idx):
                    A[p, q] = A[p, q] + block[i, j]
    return A.tocsr()


@dataclass
class MeshGeomCache:
    points: np.ndarray
    tris: np.ndarray
    area: np.ndarray
    good: np.ndarray
    bx: np.ndarray
    by: np.ndarray
    ax: np.ndarray
    ay: np.ndarray
    sxsy: np.ndarray
    rows: np.ndarray  # (9, n_tris) int
    cols: np.ndarray  # (9, n_tris) int
    n: int


@dataclass
class DeviceStaticCache:
    grade: str
    points: np.ndarray
    tris: np.ndarray
    quartz: np.ndarray
    walls: np.ndarray
    rod_masks: list
    exterior: np.ndarray
    ports: list  # geometry-only; measure set by caller
    geom: MeshGeomCache
    orbits: list
    o_masks: list
    fp_ref_ord: float


_GEOM: dict[str, MeshGeomCache] = {}
_DEVICE: dict[tuple, DeviceStaticCache] = {}


def mesh_geom_cache(points: np.ndarray, tris: np.ndarray) -> MeshGeomCache:
    """Frequency-independent geometric quantities for assemble_anisotropic."""
    points = np.asarray(points, float)
    tris = np.asarray(tris, int)
    n = len(points)
    cents = points[tris].mean(axis=1)
    sx, sy = pml_sx_sy(cents)
    x = points[tris, 0]
    y = points[tris, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    area = 0.5 * np.abs(twice)
    good = area > 1e-18
    bx = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1) / twice[:, None]
    by = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1) / twice[:, None]
    ax = sy / sx
    ay = sx / sy
    sxsy = sx * sy
    rows = np.empty((9, len(tris)), dtype=np.int64)
    cols = np.empty((9, len(tris)), dtype=np.int64)
    k = 0
    for i in range(3):
        for j in range(3):
            rows[k] = tris[:, i]
            cols[k] = tris[:, j]
            k += 1
    return MeshGeomCache(
        points=points,
        tris=tris,
        area=area,
        good=good,
        bx=bx,
        by=by,
        ax=ax,
        ay=ay,
        sxsy=sxsy,
        rows=rows,
        cols=cols,
        n=n,
    )


def get_mesh_geom(grade: str) -> MeshGeomCache:
    if grade not in _GEOM:
        points, tris = load_mesh(grade)
        _GEOM[grade] = mesh_geom_cache(points, tris)
    return _GEOM[grade]


def _pin_empty_diag(A: sparse.csr_matrix, empty: np.ndarray) -> sparse.csr_matrix:
    """Exactly replace empty rows by identity when those rows have no entries."""
    if not np.any(empty):
        return A
    idx = np.flatnonzero(empty)
    D = sparse.coo_matrix(
        (np.ones(len(idx), dtype=np.complex128), (idx, idx)),
        shape=A.shape,
    )
    return (A + D).tocsr()


def assemble_anisotropic_cached(
    geo: MeshGeomCache,
    rho_xx,
    rho_xy,
    rho_yx,
    rho_yy,
    k0: float,
    pec: Optional[np.ndarray] = None,
    mass_scale: float = 1.0,
    pin_empty: bool = True,
    elem_idx: Optional[np.ndarray] = None,
) -> sparse.csr_matrix:
    """Same weak form as fem_validated_solver.assemble_anisotropic."""
    if elem_idx is None:
        sel = slice(None)
        n_sel = len(geo.tris)
        rows = geo.rows
        cols = geo.cols
        area = geo.area
        good = geo.good
        bx = geo.bx
        by = geo.by
        ax = geo.ax
        ay = geo.ay
        sxsy = geo.sxsy
        rxx = np.asarray(rho_xx)
        rxy = np.asarray(rho_xy)
        ryx = np.asarray(rho_yx)
        ryy = np.asarray(rho_yy)
    else:
        sel = np.asarray(elem_idx, dtype=np.int64)
        n_sel = len(sel)
        if n_sel == 0:
            return sparse.csr_matrix((geo.n, geo.n), dtype=np.complex128)
        rows = geo.rows[:, sel]
        cols = geo.cols[:, sel]
        area = geo.area[sel]
        good = geo.good[sel]
        bx = geo.bx[sel]
        by = geo.by[sel]
        ax = geo.ax[sel]
        ay = geo.ay[sel]
        sxsy = geo.sxsy[sel]
        rxx = np.asarray(rho_xx)[sel]
        rxy = np.asarray(rho_xy)[sel]
        ryx = np.asarray(rho_yx)[sel]
        ryy = np.asarray(rho_yy)[sel]

    mass = mass_scale * (k0**2) * sxsy
    pec_el = (
        (np.abs(rxx) < 1e-10)
        & (np.abs(ryy) < 1e-10)
        & (np.abs(rxy) < 1e-10)
        & (np.abs(ryx) < 1e-10)
    )
    physical = good & ~pec_el
    data = np.empty((9, n_sel), dtype=np.complex128)
    k = 0
    for i in range(3):
        for j in range(3):
            term = (
                rxx * ax * bx[:, j] * bx[:, i]
                + rxy * by[:, j] * bx[:, i]
                + ryx * bx[:, j] * by[:, i]
                + ryy * ay * by[:, j] * by[:, i]
            )
            Ke = area * term
            Me = mass * area * (1.0 / 6.0 if i == j else 1.0 / 12.0)
            data[k] = np.where(physical, Ke - Me, 0.0)
            k += 1
    A = sparse.coo_matrix(
        (data.ravel(), (rows.ravel(), cols.ravel())),
        shape=(geo.n, geo.n),
    ).tocsr()
    if pin_empty:
        row_sum = np.abs(A).sum(axis=1).A1
        empty = row_sum < 1e-14
        A = _pin_empty_diag(A, empty)
    if pec is not None and np.any(pec):
        # Rare on production A2 path (pec usually all-False); keep exact semantics.
        A = A.tolil()
        for i in np.flatnonzero(pec):
            A.rows[i] = [i]
            A.data[i] = [1.0 + 0j]
        A = A.tocsr()
    return A


def assemble_dA_rod_cached(
    geo: MeshGeomCache,
    rod_mask: np.ndarray,
    f_ord,
    fp_ref_ord,
    gamma,
    fc,
    k0,
    s: float = 1.0,
) -> sparse.csr_matrix:
    """dA/ds for one orbit/rod mask — assemble only support elements."""
    dxx, dxy, dyx, dyy = drho_ds(f_ord, fp_ref_ord, gamma, fc, s)
    mask = np.asarray(rod_mask, dtype=bool) & geo.good
    elem_idx = np.flatnonzero(mask)
    if elem_idx.size == 0:
        return sparse.csr_matrix((geo.n, geo.n), dtype=np.complex128)
    # Constant dρ on the support; broadcast via arrays of length n_sel.
    n_sel = elem_idx.size
    rxx = np.full(n_sel, dxx, dtype=np.complex128)
    rxy = np.full(n_sel, dxy, dtype=np.complex128)
    ryx = np.full(n_sel, dyx, dtype=np.complex128)
    ryy = np.full(n_sel, dyy, dtype=np.complex128)
    # Pass already-sliced rho by temporarily building full? Better special-case.
    # Reuse assemble with elem_idx and rho arrays sized to full mesh:
    z = np.zeros(len(geo.tris), dtype=np.complex128)
    rho = [z.copy(), z.copy(), z.copy(), z.copy()]
    rho[0][elem_idx] = dxx
    rho[1][elem_idx] = dxy
    rho[2][elem_idx] = dyx
    rho[3][elem_idx] = dyy
    return assemble_anisotropic_cached(
        geo, *rho, k0, pec=None, mass_scale=0.0, pin_empty=False, elem_idx=elem_idx
    )


def contract_dA_orbit(
    geo: MeshGeomCache,
    rod_mask: np.ndarray,
    f_ord,
    fp_ref_ord,
    gamma,
    fc,
    k0,
    s: float,
    U: np.ndarray,
    Lam: np.ndarray,
) -> float:
    """2/6 Re Σ_j Λ_j · (dA U_j) without building the sparse dA matrix.

    Mathematically identical to assembling dA then matvec/dot, for mass_scale=0
    and pin_empty=False (the production assemble_dA_rod path).
    """
    dxx, dxy, dyx, dyy = drho_ds(f_ord, fp_ref_ord, gamma, fc, s)
    sel = np.flatnonzero(np.asarray(rod_mask, dtype=bool) & geo.good)
    if sel.size == 0:
        return 0.0
    area = geo.area[sel]
    bx = geo.bx[sel]
    by = geo.by[sel]
    ax = geo.ax[sel]
    ay = geo.ay[sel]
    tris = geo.tris[sel]
    # pec_el never triggers on nonzero dρ support
    acc = np.zeros(U.shape[1], dtype=np.float64)
    for i in range(3):
        for j in range(3):
            term = (
                dxx * ax * bx[:, j] * bx[:, i]
                + dxy * by[:, j] * bx[:, i]
                + dyx * bx[:, j] * by[:, i]
                + dyy * ay * by[:, j] * by[:, i]
            )
            Ke = area * term  # (n_sel,)
            li = Lam[tris[:, i], :]
            uj = U[tris[:, j], :]
            acc += 2.0 * np.real(np.sum(li * (Ke[:, None] * uj), axis=0))
    return float(np.mean(acc))


def get_device_static(grade: str, measure: str = "projected") -> DeviceStaticCache:
    key = (grade, measure)
    if key in _DEVICE:
        return _DEVICE[key]
    geom = get_mesh_geom(grade)
    points, tris = geom.points, geom.tris
    masks = material_masks(points, tris, sc)
    _, quartz, walls, _ = masks
    centers = bulb_centers_device(sc)
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    # Physical mask / ports depend on which elements are excised (walls+exterior),
    # not on the particular nonzero plasma ρ. Build a reference ρ at a fixed freq.
    s_ref = np.ones(len(rod_masks), float)
    rho, _, _, _ = rho_from_scales(
        len(tris),
        rod_masks,
        s_ref,
        quartz,
        walls,
        5.0e9,
        0.05,
        sc.a,
        sc.fp_Hz,
        sc.gamma_Hz,
    )
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    rough = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        rough.append({"center": c, "outward": n_hat})
    exterior = exterior_beyond_ports_mask(points, tris, rough)
    for arr in rho:
        arr[exterior] = 0.0
    ports = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        ps = build_port_from_boundary(
            points,
            tris,
            rho,
            c,
            n_hat,
            tang,
            W=float(sc.clear_width),
            plane_tol=0.12,
            port=p,
        )
        ps["measure"] = measure
        ports.append(ps)
    orbits, _, o_masks = orbits_and_masks(rod_masks=rod_masks)
    dev = DeviceStaticCache(
        grade=grade,
        points=points,
        tris=tris,
        quartz=quartz,
        walls=walls,
        rod_masks=rod_masks,
        exterior=exterior,
        ports=ports,
        geom=geom,
        orbits=orbits,
        o_masks=o_masks,
        fp_ref_ord=sc.fp_Hz * sc.a / 299792458.0,
    )
    _DEVICE[key] = dev
    return dev


def build_system_cached(grade, s_vec, f_hz, b_tesla, measure="projected"):
    """Same return dict as adjoint_direct_dtn.build_system, faster assembly."""
    dev = get_device_static(grade, measure=measure)
    rho, f_ord, gamma, fc = rho_from_scales(
        len(dev.tris),
        dev.rod_masks,
        s_vec,
        dev.quartz,
        dev.walls,
        f_hz,
        b_tesla,
        sc.a,
        sc.fp_Hz,
        sc.gamma_Hz,
    )
    for arr in rho:
        arr[dev.exterior] = 0.0
    k0 = 2.0 * np.pi * (f_hz * sc.a / 299792458.0)
    A0 = assemble_anisotropic_cached(dev.geom, *rho, k0, pec=None, mass_scale=1.0, pin_empty=True)
    modes = mode_constants(k0, float(sc.clear_width), n_modes=4)
    A = assemble_port_operator_fast(A0, dev.points, dev.ports, modes, n_prop=2, n_evan=2)
    return {
        "points": dev.points,
        "tris": dev.tris,
        "rho": rho,
        "A": A,
        "k0": k0,
        "ports": dev.ports,
        "modes": modes,
        "rod_masks": dev.rod_masks,
        "f_ord": f_ord,
        "gamma": gamma,
        "fc": fc,
        "fp_ref_ord": dev.fp_ref_ord,
        "_geom": dev.geom,
        "_o_masks": dev.o_masks,
    }


def clear_caches() -> None:
    _GEOM.clear()
    _DEVICE.clear()
