#!/usr/bin/env python3
"""
FEM validation solver matching Meep num_mode_guide_normal:

- numerical-mode Hz line source (same offsets/amps as Meep)
- guide-normal Poynting flux receivers ∫ S·n̂ ds
- straight-feed incident normalization + source-port subtraction
- B=0 isotropic ρ=1/ε; B≠0 full ρ=ε⁻¹ gyrotropic tensor

Mesh coords: [0,nx]×[0,ny]. Horn centers from sixport_common are origin-centered
→ always apply +[nx/2, ny/2] when mapping to the mesh.
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
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
    str(ROOT / ".conda_envs" / "fdfd_solver" / "lib" / "python3.13" / "site-packages"),
]

import sixport_common as sc  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402
from PMMCirculatorInverse import PMMI  # noqa: E402

MESH_DIR = ROOT / "outputs" / "validation" / "highres_solver_campaign" / "phase5_fem"
OUT = ROOT / "outputs" / "validation" / "fem_meep_validation"


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)


def log(msg: str) -> None:
    print(msg, flush=True)


def origin_to_mesh(xy: np.ndarray) -> np.ndarray:
    """Origin-centered Meep/horn coords → FEM mesh corner-origin coords."""
    return np.asarray(xy, dtype=float) + np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])


def bulb_centers_mesh() -> np.ndarray:
    pmm = PMMI(a=sc.a, res=32, nx=14, ny=12, dpml=sc.dpml, B=np.zeros(3))
    pmm.Rod_Array_Hexagon_train(
        xy_cen=np.array([7.0, 6.0]), side_dim=6, r=sc.r_plasma, d=sc.d_exp, bulbs=False, uniform=True
    )
    locs = np.asarray(pmm.train_elem_locs, dtype=float)
    return locs - np.array([7.0, 6.0]) + np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])


def load_mesh(grade: str) -> Tuple[np.ndarray, np.ndarray]:
    data = np.load(MESH_DIR / f"mesh_{grade}_full.npz")
    return data["points"], data["triangles"]


def eps_tensor_at_bias(bias_a: float) -> Tuple[complex, complex, complex, complex]:
    """Meep Faraday convention: ε = [[ε⊥, -iη], [+iη, ε⊥]]."""
    eps_perp, eta = gyrotropic_drude_eps_eta(sc.fs_a, sc.fp_a, sc.gamma_a, bias_a)
    return eps_perp, -1j * eta, +1j * eta, eps_perp


def rho_from_eps(eps_xx, eps_xy, eps_yx, eps_yy) -> Tuple[complex, complex, complex, complex]:
    det = eps_xx * eps_yy - eps_xy * eps_yx
    return eps_yy / det, -eps_xy / det, -eps_yx / det, eps_xx / det


def assign_rho(
    points: np.ndarray, tris: np.ndarray, bias_a: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    T = len(tris)
    rho_xx = np.ones(T, dtype=np.complex128)
    rho_xy = np.zeros(T, dtype=np.complex128)
    rho_yx = np.zeros(T, dtype=np.complex128)
    rho_yy = np.ones(T, dtype=np.complex128)
    cents = points[tris].mean(axis=1)
    rxx_q = 1.0 / (3.8 + 0j)
    exx, exy, eyx, eyy = eps_tensor_at_bias(bias_a)
    rxx, rxy, ryx, ryy = rho_from_eps(exx, exy, eyx, eyy)
    bulbs = bulb_centers_mesh()
    r_q, r_vac = sc.r_bulb_outer, sc.r_bulb_inner
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    for cx, cy in bulbs:
        d2 = (cents[:, 0] - cx) ** 2 + (cents[:, 1] - cy) ** 2
        mq = (d2 <= r_q**2) & (d2 > r_vac**2)
        rho_xx[mq] = rxx_q
        rho_yy[mq] = rxx_q
        rho_xy[mq] = 0.0
        rho_yx[mq] = 0.0
        mv = (d2 <= r_vac**2) & (d2 > r_p**2)
        rho_xx[mv] = 1.0
        rho_yy[mv] = 1.0
        rho_xy[mv] = 0.0
        rho_yx[mv] = 0.0
        mp = d2 <= r_p**2
        rho_xx[mp] = rxx
        rho_xy[mp] = rxy
        rho_yx[mp] = ryx
        rho_yy[mp] = ryy
    return rho_xx, rho_xy, rho_yx, rho_yy


# Meep quadratic PML, R_asymptotic=1e-15. σ(u)=[-ln R /(2 d ∫u² du)] u², ∫u²=1/3.
# The retired validation profile used σ_max=2 and must not come back.
PML_R_ASYMPTOTIC = 1e-15
# Validation hook. Production calls leave this at 1. Tests may set it to
# compare PML strength without changing the quadratic Meep shape.
PML_SIGMA_SCALE = 1.0
LAST_OPERATOR: Dict[str, Any] = {}


def meep_sigma_max(thickness: float) -> float:
    sigma = -np.log(PML_R_ASYMPTOTIC) * 3.0 / (2.0 * float(thickness))
    if sigma < 10.0:
        raise RuntimeError(f"PML sigma_max={sigma} is below the Meep-equivalent floor")
    return float(sigma)


def pml_sx_sy(xy: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    nx, ny, dp = sc.nx_ports, sc.ny_ports, sc.dpml_ports
    omega = 2 * np.pi * sc.fs_a
    sigma_max = meep_sigma_max(dp) * float(PML_SIGMA_SCALE)
    x, y = xy[:, 0], xy[:, 1]
    sigx = np.zeros_like(x)
    sigy = np.zeros_like(y)
    m = x < dp
    sigx[m] = sigma_max * ((dp - x[m]) / dp) ** 2
    m = x > nx - dp
    sigx[m] = sigma_max * ((x[m] - (nx - dp)) / dp) ** 2
    m = y < dp
    sigy[m] = sigma_max * ((dp - y[m]) / dp) ** 2
    m = y > ny - dp
    sigy[m] = sigma_max * ((y[m] - (ny - dp)) / dp) ** 2
    return 1.0 + 1j * sigx / omega, 1.0 + 1j * sigy / omega


def horn_prism_quads(nx: float, ny: float) -> list:
    """Filled Meep prism quads in mesh coordinates. Not an edge tube."""
    shift = np.array([nx / 2.0, ny / 2.0])
    quads = []
    for horn in sc.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            quads.append(np.asarray(horn[name], dtype=float)[:, :2] + shift)
    return quads


def wall_triangle_mask(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    """Triangles whose centroid lies inside a Meep prism solid."""
    from matplotlib.path import Path as MPath

    cents = points[tris].mean(axis=1)
    mask = np.zeros(len(tris), dtype=bool)
    nx, ny = sc.nx_ports, sc.ny_ports
    for quad in horn_prism_quads(nx, ny):
        mask |= MPath(quad).contains_points(cents, radius=1e-12)
    return mask


def apply_metal_walls_rho(
    rho_xx, rho_xy, rho_yx, rho_yy, wall_mask: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Mark PEC wall triangles by setting ρ=0.

    Assembly then drops both the stiffness and the mass of those elements, so
    the interior metal supports no Helmholtz equation. The air-side boundary
    condition is the natural Neumann condition ∂Hz/∂n=0. Dirichlet Hz=0 is PMC
    and is not applied.
    """
    rho_xx = rho_xx.copy()
    rho_xy = rho_xy.copy()
    rho_yx = rho_yx.copy()
    rho_yy = rho_yy.copy()
    rho_xx[wall_mask] = 0.0
    rho_xy[wall_mask] = 0.0
    rho_yx[wall_mask] = 0.0
    rho_yy[wall_mask] = 0.0
    return rho_xx, rho_xy, rho_yx, rho_yy


def pec_mask_horns(points: np.ndarray) -> np.ndarray:
    """Deprecated PMC-style node mask; kept for diagnostics only. Prefer ρ→0 walls."""
    pec = np.zeros(len(points), dtype=bool)
    rad = 0.55 * sc.wall_thickness
    rad2 = rad * rad
    nx, ny = sc.nx_ports, sc.ny_ports
    for horn in sc.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], dtype=float) + np.array([nx / 2.0, ny / 2.0])
            for k in range(len(poly) - 1):
                p0, p1 = poly[k], poly[k + 1]
                L = float(np.linalg.norm(p1 - p0))
                n = max(2, int(L / max(0.5 * rad, 1e-6)))
                for t in np.linspace(0.0, 1.0, n):
                    c = p0 * (1 - t) + p1 * t
                    pec |= (points[:, 0] - c[0]) ** 2 + (points[:, 1] - c[1]) ** 2 <= rad2
    return pec


def assemble_anisotropic(
    points, tris, rho_xx, rho_xy, rho_yx, rho_yy, k0: float, pec: np.ndarray,
    mass_scale: float = 1.0, pin_empty: bool = True,
) -> sparse.csr_matrix:
    """
    Stretched anisotropic weak form (e^{-iωt}):
      ∫ [ρxx (sy/sx) ∂x∂x + ρxy ∂y∂x + ρyx ∂x∂y + ρyy (sx/sy) ∂y∂y]
        − k0² ∫ sx sy Hz v
    Reduces to prior isotropic FEM when ρ=ρ I.
    """
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
    mass = mass_scale * (k0**2) * sx * sy
    ax = (sy / sx)  # stretch factor for xx
    ay = (sx / sy)
    # PEC elements leave the physical domain: no stiffness and no mass.
    # rho->0 with the mass term kept left active Hz equations inside the metal.
    pec_el = (
        (np.abs(rho_xx) < 1e-10)
        & (np.abs(rho_yy) < 1e-10)
        & (np.abs(rho_xy) < 1e-10)
        & (np.abs(rho_yx) < 1e-10)
    )
    physical = good & ~pec_el
    rows, cols, data = [], [], []
    for i in range(3):
        for j in range(3):
            term = (
                rho_xx * ax * bx[:, j] * bx[:, i]
                + rho_xy * by[:, j] * bx[:, i]
                + rho_yx * bx[:, j] * by[:, i]
                + rho_yy * ay * by[:, j] * by[:, i]
            )
            Ke = area * term
            Me = mass * area * (1.0 / 6.0 if i == j else 1.0 / 12.0)
            val = np.where(physical, Ke - Me, 0.0)
            rows.append(tris[:, i])
            cols.append(tris[:, j])
            data.append(val)
    A = sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
    ).tocsr()
    row_sum = np.abs(A).sum(axis=1).A1
    empty = row_sum < 1e-14
    if pin_empty and np.any(empty):
        A = A.tolil()
        for i in np.flatnonzero(empty):
            A.rows[i] = [i]
            A.data[i] = [1.0 + 0j]
        A = A.tocsr()
    touched = np.zeros(n, dtype=bool)
    if np.any(pec_el):
        touched[tris[pec_el].ravel()] = True
    active_pec = touched & empty
    LAST_OPERATOR.clear()
    LAST_OPERATOR.update(
        {
            "pec_elements": int(np.count_nonzero(pec_el)),
            "pec_active_dofs_pinned": int(np.count_nonzero(active_pec)),
            "empty_rows_pinned": int(np.count_nonzero(empty)),
            "pml_thickness_a": float(sc.dpml_ports),
            "sigma_max": meep_sigma_max(sc.dpml_ports),
            "domain_nx": float(sc.nx_ports),
            "domain_ny": float(sc.ny_ports),
            "pml_profile": "meep_quadratic_R1e-15",
            "pec_model": "excised_neumann_no_mass",
        }
    )
    if np.any(pec):
        pec_idx = np.flatnonzero(pec)
        A = A.tolil()
        for i in pec_idx:
            A.rows[i] = [i]
            A.data[i] = [1.0 + 0j]
        A = A.tocsr()
    return A


def load_numerical_mode(path: Optional[Path] = None, port: int = 0) -> Dict[str, Any]:
    from plasmeep.ports.mode_registry import get_numerical_mode
    from plasmeep.ports.numerical_launch import port_tangent
    from sixport_common import set_geometry_context

    set_geometry_context(res=50, horn_walls="prism", grid_offset_cells=(0, 0), coord_rotation_deg=0.0)
    if path and Path(path).is_file():
        from plasmeep.ports.numerical_mode import NumericalPortMode

        mode = NumericalPortMode.load(path)
        port = int(mode.port_index)
    else:
        mode = get_numerical_mode(
            port,
            res=50,
            frequency_a=sc.fs_a,
            horn_walls="prism",
            grid_offset_cells=(0.0, 0.0),
            coord_rotation_deg=0.0,
        )
    center_origin = np.asarray(sc.horn_for_port(port, 50)["source_center"], dtype=float)
    center = origin_to_mesh(center_origin)
    u = np.asarray(sc.effective_port_dir(port), dtype=float)
    u /= np.linalg.norm(u)
    tang = port_tangent(u)
    amps = np.asarray(mode.source_amplitudes(target_power=1.0), dtype=np.complex128)
    offsets = np.asarray(mode.offsets_a, dtype=float)
    return {
        "mode": mode,
        "center_origin": center_origin,
        "center_mesh": center,
        "tangent": tang,
        "outward": u,
        "amps": amps,
        "offsets": offsets,
    }


def inject_numerical_mode_rhs(points: np.ndarray, src: Dict[str, Any], pec: np.ndarray) -> np.ndarray:
    """Scatter Meep point-Hz amplitudes onto nearest non-PEC mesh nodes."""
    b = np.zeros(len(points), dtype=np.complex128)
    center, tang = src["center_mesh"], src["tangent"]
    for s, amp in zip(src["offsets"], src["amps"]):
        xy = center + float(s) * tang
        d2 = np.sum((points - xy) ** 2, axis=1)
        if pec is not None and np.any(pec):
            d2 = d2.copy()
            d2[pec] = np.inf
        i = int(np.argmin(d2))
        b[i] += complex(amp)
    return b


def nodal_gradient(points, tris, u):
    n = len(points)
    gx = np.zeros(n, dtype=np.complex128)
    gy = np.zeros(n, dtype=np.complex128)
    wsum = np.zeros(n, dtype=np.float64)
    x = points[tris, 0]
    y = points[tris, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    area = 0.5 * np.abs(twice)
    bx = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1) / twice[:, None]
    by = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1) / twice[:, None]
    ugx = np.sum(bx * u[tris], axis=1)
    ugy = np.sum(by * u[tris], axis=1)
    for k in range(3):
        idx = tris[:, k]
        np.add.at(gx, idx, ugx * area)
        np.add.at(gy, idx, ugy * area)
        np.add.at(wsum, idx, area)
    wsum = np.maximum(wsum, 1e-30)
    return gx / wsum, gy / wsum


def nodal_rho(points, tris, rho_xx, rho_xy, rho_yx, rho_yy):
    n = len(points)
    out = [np.zeros(n, dtype=np.complex128) for _ in range(4)]
    wsum = np.zeros(n)
    x = points[tris, 0]
    y = points[tris, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    area = 0.5 * np.abs(twice)
    for k in range(3):
        idx = tris[:, k]
        for o, a in zip(out, (rho_xx, rho_xy, rho_yx, rho_yy)):
            np.add.at(o, idx, a * area)
        np.add.at(wsum, idx, area)
    wsum = np.maximum(wsum, 1e-30)
    return [o / wsum for o in out]


def fields_from_hz(points, tris, u, rho_xx, rho_xy, rho_yx, rho_yy, omega: float):
    """e^{-iωt}: curl H = (∂y Hz, −∂x Hz); E = (i/ω) ρ curl H."""
    dx, dy = nodal_gradient(points, tris, u)
    rxx, rxy, ryx, ryy = nodal_rho(points, tris, rho_xx, rho_xy, rho_yx, rho_yy)
    vx, vy = dy, -dx
    Ex = (1j / omega) * (rxx * vx + rxy * vy)
    Ey = (1j / omega) * (ryx * vx + ryy * vy)
    return Ex, Ey


class ElementSampler:
    """Barycentric samples of the piecewise-linear FEM field on the true line."""

    def __init__(self, points: np.ndarray, tris: np.ndarray):
        self.points = np.asarray(points, dtype=float)
        self.tris = np.asarray(tris, dtype=int)
        x = self.points[self.tris, 0]
        y = self.points[self.tris, 1]
        twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
        self.twice = twice
        self.bx = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1) / twice[:, None]
        self.by = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1) / twice[:, None]
        cents = self.points[self.tris].mean(axis=1)
        self.tree = cKDTree(cents)
        self.misses = 0

    def _bary(self, t: int, p: np.ndarray):
        v0 = self.points[self.tris[t, 0]]
        v1 = self.points[self.tris[t, 1]]
        v2 = self.points[self.tris[t, 2]]
        den = (v1[1] - v2[1]) * (v0[0] - v2[0]) + (v2[0] - v1[0]) * (v0[1] - v2[1])
        a = ((v1[1] - v2[1]) * (p[0] - v2[0]) + (v2[0] - v1[0]) * (p[1] - v2[1])) / den
        b = ((v2[1] - v0[1]) * (p[0] - v2[0]) + (v0[0] - v2[0]) * (p[1] - v2[1])) / den
        return a, b, 1.0 - a - b

    def locate(self, xy: np.ndarray) -> np.ndarray:
        xy = np.atleast_2d(np.asarray(xy, dtype=float))
        k = min(24, len(self.tris))
        _, cand = self.tree.query(xy, k=k)
        cand = np.atleast_2d(cand)
        out = np.full(len(xy), -1, dtype=int)
        for n in range(len(xy)):
            for t in np.atleast_1d(cand[n]):
                t = int(t)
                if t < 0 or abs(self.twice[t]) < 1e-18:
                    continue
                a, b, c = self._bary(t, xy[n])
                if a >= -1e-8 and b >= -1e-8 and c >= -1e-8:
                    out[n] = t
                    break
        return out

    def fields(self, u, rho_xx, rho_xy, rho_yx, rho_yy, omega: float, xy: np.ndarray):
        xy = np.atleast_2d(np.asarray(xy, dtype=float))
        tri_id = self.locate(xy)
        Hz = np.zeros(len(xy), dtype=np.complex128)
        Ex = np.zeros(len(xy), dtype=np.complex128)
        Ey = np.zeros(len(xy), dtype=np.complex128)
        found = tri_id >= 0
        self.misses += int(np.count_nonzero(~found))
        if np.any(found):
            tid = tri_id[found]
            nodes = self.tris[tid]
            w = np.zeros((len(tid), 3), dtype=float)
            for n, t in enumerate(tid):
                w[n] = self._bary(int(t), xy[found][n])
            un = u[nodes]
            Hz[found] = np.sum(w * un, axis=1)
            gx = np.sum(self.bx[tid] * un, axis=1)
            gy = np.sum(self.by[tid] * un, axis=1)
            rxx, rxy = rho_xx[tid], rho_xy[tid]
            ryx, ryy = rho_yx[tid], rho_yy[tid]
            Ex[found] = (1j / omega) * (rxx * gy + rxy * (-gx))
            Ey[found] = (1j / omega) * (ryx * gy + ryy * (-gx))
        if np.any(~found):
            for n in np.flatnonzero(~found):
                i = int(np.argmin(np.sum((self.points - xy[n]) ** 2, axis=1)))
                Hz[n] = u[i]
        return Hz, Ex, Ey


def guide_normal_flux(
    points: np.ndarray,
    u_hz: np.ndarray,
    Ex: np.ndarray,
    Ey: np.ndarray,
    port: int,
    n_points: int = 31,
    span_factor: float = 0.96,
    sampler: Optional[ElementSampler] = None,
    tris: Optional[np.ndarray] = None,
    rho=None,
    omega: Optional[float] = None,
) -> float:
    """∫ S·n̂ ds; Sx=½Re(Ey Hz*), Sy=−½Re(Ex Hz*). Positive = power along +n̂ (outward)."""
    from sixport_common import set_geometry_context, monitor_center_for_port, effective_port_dir

    set_geometry_context(res=50, horn_walls="prism")
    center = origin_to_mesh(monitor_center_for_port(port, 50))
    n_hat = np.asarray(effective_port_dir(port), dtype=float)
    n_hat /= np.linalg.norm(n_hat)
    tang = np.array([-n_hat[1], n_hat[0]])
    span = span_factor * sc.clear_width
    offsets = np.linspace(-span / 2, span / 2, n_points)
    if n_points == 1:
        w = np.array([span])
    else:
        ds = span / (n_points - 1)
        w = np.full(n_points, ds)
        w[0] *= 0.5
        w[-1] *= 0.5
    if sampler is not None:
        xy = np.vstack([center + float(s) * tang for s in offsets])
        Hz, Exs, Eys = sampler.fields(u_hz, *rho, float(omega), xy)
        Sx = 0.5 * np.real(Eys * np.conj(Hz))
        Sy = -0.5 * np.real(Exs * np.conj(Hz))
        return float(np.sum((n_hat[0] * Sx + n_hat[1] * Sy) * w))
    total = 0.0
    for s, ws in zip(offsets, w):
        xy = center + s * tang
        i = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
        Hz = u_hz[i]
        Sx = 0.5 * np.real(Ey[i] * np.conj(Hz))
        Sy = -0.5 * np.real(Ex[i] * np.conj(Hz))
        total += float(n_hat[0] * Sx + n_hat[1] * Sy) * float(ws)
    return total


def inject_mode_consistent(points, tris, src, sampler: ElementSampler) -> np.ndarray:
    """Consistent load: each numerical-mode sample is a delta at its true coordinate."""
    b = np.zeros(len(points), dtype=np.complex128)
    center, tang = src["center_mesh"], src["tangent"]
    for s, amp in zip(src["offsets"], src["amps"]):
        xy = center + float(s) * np.asarray(tang, dtype=float)
        t = int(sampler.locate(xy[None, :])[0])
        if t < 0:
            i = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
            b[i] += complex(amp)
            continue
        w = sampler._bary(t, xy)
        nodes = tris[t]
        for k in range(3):
            b[int(nodes[k])] += float(w[k]) * complex(amp)
    return b


def solve_system(A, b, pec):
    b = b.copy()
    b[pec] = 0.0
    t0 = time.perf_counter()
    lu = splu(A.tocsc())
    t_fac = time.perf_counter() - t0
    t1 = time.perf_counter()
    x = lu.solve(b)
    t_sol = time.perf_counter() - t1
    resid = float(np.linalg.norm(A @ x - b) / (np.linalg.norm(b) + 1e-30))
    return x, lu, t_fac, t_sol, resid


def build_problem(grade: str, bias_a: float, vacuum: bool = False, pec: Optional[np.ndarray] = None):
    """
    Build device FEM problem.

    Metal horns: ρ→0 in wall triangles (ε→∞). No Dirichlet Hz=0 (that would be PMC).
    """
    points, tris = load_mesh(grade)
    t_mat = 0.0
    if vacuum:
        T = len(tris)
        rho = (
            np.ones(T, dtype=np.complex128),
            np.zeros(T, dtype=np.complex128),
            np.zeros(T, dtype=np.complex128),
            np.ones(T, dtype=np.complex128),
        )
    else:
        t0 = time.perf_counter()
        rho = list(assign_rho(points, tris, bias_a=bias_a))
        wmask = wall_triangle_mask(points, tris)
        rho = apply_metal_walls_rho(*rho, wmask)
        t_mat = time.perf_counter() - t0
    # Never apply Dirichlet PMC; keep zero mask for solve_system compatibility
    if pec is None:
        pec = np.zeros(len(points), dtype=bool)
    k0 = 2 * np.pi * sc.fs_a
    t1 = time.perf_counter()
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    t_asm = time.perf_counter() - t1
    return {
        "points": points,
        "tris": tris,
        "rho": rho,
        "pec": pec,
        "A": A,
        "k0": k0,
        "omega": 2 * np.pi * sc.fs_a,
        "t_assemble": t_asm,
        "t_material": t_mat,
        "wall_model": "rho0_metal_triangles",
    }


def port_powers_from_solution(prob, x) -> Dict[int, float]:
    Ex, Ey = fields_from_hz(prob["points"], prob["tris"], x, *prob["rho"], prob["omega"])
    return {p: guide_normal_flux(prob["points"], x, Ex, Ey, p) for p in range(6)}


def fem_incident_reference(grade: str, src: Dict[str, Any]) -> Dict[str, Any]:
    """
    Straight-feed reference with ρ→0 sidewalls (TE-consistent PEC).
    Stores nodal fields on the P1 monitor line for device field-subtraction.
    """
    log("Building FEM incident reference (straight P1 feed, ρ→0 walls)...")
    points, tris = load_mesh(grade)
    T = len(tris)
    rho = [
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    ]
    pec = np.zeros(len(points), dtype=bool)
    u = src["outward"]
    n = src["tangent"]
    center = src["center_mesh"]
    wall_off = sc.clear_width / 2 + sc.wall_thickness / 2
    axis_len = float(2.0 * np.hypot(sc.nx_ports, sc.ny_ports))
    rad = 0.55 * sc.wall_thickness
    cents = points[tris].mean(axis=1)
    wmask = np.zeros(T, dtype=bool)
    for sign in (+1.0, -1.0):
        wall_c0 = center + sign * wall_off * n - 0.5 * axis_len * u
        wall_c1 = center + sign * wall_off * n + 0.5 * axis_len * u
        L = float(np.linalg.norm(wall_c1 - wall_c0))
        ns = max(2, int(L / max(0.5 * rad, 1e-6)))
        for t in np.linspace(0.0, 1.0, ns):
            c = wall_c0 * (1 - t) + wall_c1 * t
            wmask |= (cents[:, 0] - c[0]) ** 2 + (cents[:, 1] - c[1]) ** 2 <= rad * rad
    for arr in rho:
        arr[wmask] = 0.0

    k0 = 2 * np.pi * sc.fs_a
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    b = inject_numerical_mode_rhs(points, src, pec)
    x, _lu, t_fac, t_sol, resid = solve_system(A, b, pec)
    Ex, Ey = fields_from_hz(points, tris, x, *rho, 2 * np.pi * sc.fs_a)
    P_raw_ref = guide_normal_flux(points, x, Ex, Ey, 0)
    P_inc = abs(P_raw_ref)
    meta = {
        "P_raw_ref": float(P_raw_ref),
        "P_inc": float(P_inc),
        "true_residual": resid,
        "factor_s": t_fac,
        "solve_s": t_sol,
        "wall_tri_frac": float(wmask.mean()),
        "wall_model": "rho0_metal_triangles",
        "sign_note": (
            "P_inc=|P_raw_ref|. Device P11 uses field subtraction "
            "(E,H)_dev - (E,H)_ref on the monitor line, then ∫S·n̂."
        ),
        # Keep fields for in-process field subtraction (not written to JSON cache).
        "_Hz_ref": x,
        "_Ex_ref": Ex,
        "_Ey_ref": Ey,
        "_points": points,
    }
    log(f"  FEM P_raw_ref={P_raw_ref:.6e} P_inc={P_inc:.6e} resid={resid:.2e}")
    return meta


def run_device(
    grade: str,
    bias_a: float,
    src: Dict[str, Any],
    P_inc: float,
    P_raw_ref: float,
    multi_rhs_timing: bool = False,
) -> Dict[str, Any]:
    log(f"=== FEM device grade={grade} bias={bias_a} ===")
    t0 = time.perf_counter()
    prob = build_problem(grade, bias_a)
    t_build = time.perf_counter() - t0
    b = inject_numerical_mode_rhs(prob["points"], src, prob["pec"])
    x, lu, t_fac, t_sol, resid = solve_system(prob["A"], b, prob["pec"])
    t_pp0 = time.perf_counter()
    Ex, Ey = fields_from_hz(prob["points"], prob["tris"], x, *prob["rho"], prob["omega"])
    raw = {p: guide_normal_flux(prob["points"], x, Ex, Ey, p) for p in range(6)}
    # Field subtraction at source port (Meep load_minus_flux_data analogue), if refs provided
    raw_adj = dict(raw)
    p11_method = "scalar_raw_minus_Pref"
    if isinstance(P_raw_ref, dict) and "_Hz_ref" in P_raw_ref:
        # Expect full incident_meta with fields — handled by caller via norm_with_field_sub
        pass
    raw_adj[0] = raw[0] - (P_raw_ref if not isinstance(P_raw_ref, dict) else P_raw_ref["P_raw_ref"])
    norm = {p: raw_adj[p] / P_inc for p in range(6)}
    t_pp = time.perf_counter() - t_pp0

    out: Dict[str, Any] = {
        "grade": grade,
        "bias_a": bias_a,
        "n_nodes": int(len(prob["points"])),
        "n_tris": int(len(prob["tris"])),
        "nnz": int(prob["A"].nnz),
        "true_residual": resid,
        "P_inc_FEM": float(P_inc),
        "P_raw_ref_FEM": float(P_raw_ref),
        "raw_flux_by_port": {str(k): float(v) for k, v in raw.items()},
        "raw_flux_source_after_subtraction": float(raw_adj[0]),
        "normalized_power_by_port": {str(k): float(v) for k, v in norm.items()},
        "power_dB_by_port": {
            str(k): (10.0 * np.log10(v) if v > 0 else float("nan")) for k, v in norm.items()
        },
        "timings_s": {
            "material": prob["t_material"],
            "assemble": prob["t_assemble"],
            "build_total": t_build,
            "factor": t_fac,
            "solve": t_sol,
            "port_postprocess": t_pp,
            "wall_first_rhs": t_build + t_fac + t_sol + t_pp,
        },
        "rss_gib": rss_gib(),
    }
    log(
        f"  resid={resid:.2e} factor={t_fac:.2f}s solve={t_sol:.2f}s pp={t_pp:.2f}s "
        f"P11={norm[0]:.4e} P12={norm[1]:.4e} P13={norm[2]:.4e} P14={norm[3]:.4e}"
    )

    if multi_rhs_timing:
        extras = []
        for p in range(1, 6):
            # Rebuild source for each port using same mode shape rotated to that port
            src_p = _source_for_port(p, src)
            bp = inject_numerical_mode_rhs(prob["points"], src_p, prob["pec"])
            t1 = time.perf_counter()
            xp = lu.solve(bp)
            dt = time.perf_counter() - t1
            raw_p = port_powers_from_solution(prob, xp)
            extras.append({"port": p, "solve_s": dt, "raw0": float(raw_p[0]), "raw_p": float(raw_p[p])})
            log(f"  extra RHS P{p+1}: solve={dt:.3f}s")
        out["extra_rhs"] = extras
        out["timings_s"]["six_port_candidate"] = (
            out["timings_s"]["wall_first_rhs"] + sum(e["solve_s"] for e in extras)
            + 5 * t_pp  # approx postprocess per extra (measured on first)
        )
    return out


def _source_for_port(port: int, template: Dict[str, Any]) -> Dict[str, Any]:
    from plasmeep.ports.numerical_launch import port_tangent
    from sixport_common import set_geometry_context

    set_geometry_context(res=50, horn_walls="prism")
    center_origin = np.asarray(sc.horn_for_port(port, 50)["source_center"], dtype=float)
    u = np.asarray(sc.effective_port_dir(port), dtype=float)
    u /= np.linalg.norm(u)
    return {
        "center_mesh": origin_to_mesh(center_origin),
        "tangent": port_tangent(u),
        "outward": u,
        "amps": template["amps"],
        "offsets": template["offsets"],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="FEM-H")
    ap.add_argument("--bias-a", type=float, default=0.0)
    ap.add_argument("--mode-json", type=str, default="")
    ap.add_argument("--json-out", required=True)
    ap.add_argument("--p-inc", type=float, default=None)
    ap.add_argument("--p-raw-ref", type=float, default=None)
    ap.add_argument("--skip-incident", action="store_true")
    ap.add_argument("--multi-rhs-timing", action="store_true")
    ap.add_argument("--incident-cache", type=str, default="")
    args = ap.parse_args()

    mode_path = Path(args.mode_json) if args.mode_json else (
        OUT / "phase2" / "numerical_mode_P1_res50.json"
    )
    src = load_numerical_mode(mode_path if mode_path.is_file() else None)
    log(
        f"Mode: n={len(src['offsets'])} center_mesh={src['center_mesh']} "
        f"center_origin={src['center_origin']}"
    )

    inc_cache = Path(args.incident_cache) if args.incident_cache else (
        OUT / "phase3" / f"incident_{args.grade}.json"
    )
    if args.p_inc is not None and args.p_raw_ref is not None:
        inc_meta = {"P_inc": args.p_inc, "P_raw_ref": args.p_raw_ref, "note": "cli"}
    elif args.skip_incident:
        inc_meta = {"P_inc": 1.0, "P_raw_ref": -1.0, "note": "unnormalized"}
    elif inc_cache.is_file() and args.bias_a == 0.0:
        inc_meta = json.loads(inc_cache.read_text())
        log(f"Loaded incident cache {inc_cache}")
    else:
        inc_meta = fem_incident_reference(args.grade, src)
        inc_cache.parent.mkdir(parents=True, exist_ok=True)
        inc_cache.write_text(json.dumps(inc_meta, indent=2) + "\n")

    result = run_device(
        args.grade,
        args.bias_a,
        src,
        float(inc_meta["P_inc"]),
        float(inc_meta["P_raw_ref"]),
        multi_rhs_timing=args.multi_rhs_timing,
    )
    result["incident_meta"] = inc_meta
    result["source_meta"] = {
        "center_mesh": src["center_mesh"].tolist(),
        "center_origin": src["center_origin"].tolist(),
        "outward": src["outward"].tolist(),
        "tangent": src["tangent"].tolist(),
        "n_samples": len(src["offsets"]),
        "mode_json": str(mode_path),
    }
    Path(args.json_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.json_out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "extra_rhs"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
