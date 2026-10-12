"""Validated Hz FEM bulk operator on local meshes (unity stretch, no domain PML).

Same weak form as fem_validated_solver.assemble_anisotropic with sx=sy=1.
Used for straight-guide DtN validation without Meep/sixport imports.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree


def assemble_anisotropic(
    points,
    tris,
    rho_xx,
    rho_xy,
    rho_yx,
    rho_yy,
    k0: float,
    pec: np.ndarray,
    mass_scale: float = 1.0,
    pin_empty: bool = True,
) -> sparse.csr_matrix:
    n = len(points)
    sx = sy = 1.0
    x = points[tris, 0]
    y = points[tris, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    area = 0.5 * np.abs(twice)
    good = area > 1e-18
    bx = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1) / twice[:, None]
    by = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1) / twice[:, None]
    mass = mass_scale * (k0**2) * sx * sy
    ax = sy / sx
    ay = sx / sy
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
    if np.any(pec):
        pec_idx = np.flatnonzero(pec)
        A = A.tolil()
        for i in pec_idx:
            A.rows[i] = [i]
            A.data[i] = [1.0 + 0j]
        A = A.tocsr()
    return A


class ElementSampler:
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
