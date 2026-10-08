#!/usr/bin/env python3
"""Shared material-only differentiability helpers for the validated Hz FEM.

Fixed geometry. Design variables are per-rod density scales s_i with
fp_i^2 = s_i * fp_ref^2. Forward sensitivity: A dx/ds = -(dA/ds) x.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.sparse.linalg import SuperLU, splu

from fem_validated_solver import ElementSampler, assemble_anisotropic
from gyrotropic_tensor import ordinary_from_si, rho_xy, tensor_ordinary
from plasma_sensitivity import element_drho


@dataclass
class ForwardState:
    points: np.ndarray
    tris: np.ndarray
    rho: list
    A: object
    lu: SuperLU
    x: np.ndarray
    k0: float
    f_ord: float
    gamma: float
    fc: float
    fp_ref_ord: float
    b_tesla: float
    s_vec: np.ndarray
    rod_masks: list
    quartz: np.ndarray
    walls: np.ndarray


def per_rod_plasma_masks(points, tris, centers, r_plasma: float) -> list[np.ndarray]:
    cents = points[tris].mean(1)
    out = []
    for c in centers:
        d2 = (cents[:, 0] - c[0]) ** 2 + (cents[:, 1] - c[1]) ** 2
        out.append(d2 <= r_plasma**2)
    return out


def rho_from_scales(
    n_tris: int,
    rod_masks: Sequence[np.ndarray],
    s_vec: Sequence[float],
    quartz: np.ndarray,
    walls: np.ndarray,
    f_hz: float,
    b_tesla: float,
    a_m: float,
    fp_hz: float,
    gamma_hz: float,
):
    rho = [
        np.ones(n_tris, np.complex128),
        np.zeros(n_tris, np.complex128),
        np.zeros(n_tris, np.complex128),
        np.ones(n_tris, np.complex128),
    ]
    rho[0][quartz] = 1.0 / 3.8
    rho[3][quartz] = 1.0 / 3.8
    f = gamma = fc = None
    for mask, s in zip(rod_masks, s_vec):
        f, fp, gamma, fc = ordinary_from_si(f_hz, fp_hz * np.sqrt(float(s)), gamma_hz, b_tesla, a_m)
        exx, exy, eyx, eyy, _, _ = tensor_ordinary(f, fp, gamma, fc)
        rxx, rxy, ryx, ryy, _ = rho_xy(exx, exy, eyx, eyy)
        rho[0][mask] = rxx
        rho[1][mask] = rxy
        rho[2][mask] = ryx
        rho[3][mask] = ryy
    for arr in rho:
        arr[walls] = 0.0
    return rho, f, gamma, fc


def assemble_dA_rod(points, tris, rod_mask, f_ord, fp_ref_ord, gamma, fc, k0, s: float = 1.0):
    dr = element_drho(rod_mask, f_ord, fp_ref_ord, gamma, fc, s)
    return assemble_anisotropic(
        points, tris, *dr, k0, np.zeros(len(points), dtype=bool), mass_scale=0.0, pin_empty=False
    )


def solve_forward(
    points,
    tris,
    rod_masks,
    quartz,
    walls,
    s_vec,
    bvec,
    f_hz: float,
    b_tesla: float,
    a_m: float,
    fp_hz: float,
    gamma_hz: float,
) -> ForwardState:
    k0 = 2.0 * np.pi * (f_hz * a_m / 299792458.0)
    rho, f, gamma, fc = rho_from_scales(
        len(tris), rod_masks, s_vec, quartz, walls, f_hz, b_tesla, a_m, fp_hz, gamma_hz
    )
    A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    lu = splu(A.tocsc())
    x = lu.solve(bvec)
    fp_ref_ord = fp_hz * a_m / 299792458.0
    return ForwardState(
        points=points,
        tris=tris,
        rho=rho,
        A=A,
        lu=lu,
        x=x,
        k0=k0,
        f_ord=f,
        gamma=gamma,
        fc=fc,
        fp_ref_ord=fp_ref_ord,
        b_tesla=b_tesla,
        s_vec=np.asarray(s_vec, float),
        rod_masks=list(rod_masks),
        quartz=quartz,
        walls=walls,
    )


def sensitivity_column(state: ForwardState, rod_index: int) -> np.ndarray:
    """dx/ds_k using the existing LU of A. One RHS, no refactor."""
    dA = assemble_dA_rod(
        state.points,
        state.tris,
        state.rod_masks[rod_index],
        state.f_ord,
        state.fp_ref_ord,
        state.gamma,
        state.fc,
        state.k0,
        float(state.s_vec[rod_index]),
    )
    return state.lu.solve(-(dA @ state.x))


def monitor_amp(sampler: ElementSampler, uh, rho, k0, xy, weights) -> complex:
    hz, _, _ = sampler.fields(uh, *rho, k0, xy)
    return complex(np.dot(weights, hz))
