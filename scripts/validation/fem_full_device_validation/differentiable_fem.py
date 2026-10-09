#!/usr/bin/env python3
"""Clean material-only differentiable FEM API over the validated assembler.

    solve_forward(parameters, frequency, B, source_b)
    compute_outputs(state, monitors)
    compute_jacobian(state, monitors, rod_indices=None)
    compute_objective(state, objective)
    compute_gradient(state, objective)

Geometry/mesh fixed. Parameters are per-rod density scales s_i.
No giant dx arrays are written to disk by default.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from diff_fem_core import (
    ForwardState,
    assemble_dA_rod,
    monitor_amp,
    sensitivity_column,
    solve_forward as _solve_forward_core,
)
from fem_validated_solver import ElementSampler


@dataclass
class Monitor:
    name: str
    xy: np.ndarray
    weights: np.ndarray


@dataclass
class Objective:
    """Real scalar J = |L x|^2 for a linear complex monitor L."""

    name: str
    L: sparse.csr_matrix  # shape (1, n)


def solve_forward(
    points,
    tris,
    rod_masks,
    quartz,
    walls,
    parameters: Sequence[float],
    frequency_hz: float,
    B_tesla: float,
    source_b: np.ndarray,
    a_m: float,
    fp_hz: float,
    gamma_hz: float,
) -> ForwardState:
    return _solve_forward_core(
        points, tris, rod_masks, quartz, walls, parameters, source_b,
        frequency_hz, B_tesla, a_m, fp_hz, gamma_hz,
    )


def compute_outputs(state: ForwardState, monitors: Sequence[Monitor]) -> dict:
    sampler = ElementSampler(state.points, state.tris)
    out = {}
    for m in monitors:
        out[m.name] = monitor_amp(sampler, state.x, state.rho, state.k0, m.xy, m.weights)
    return out


def compute_jacobian(
    state: ForwardState,
    monitors: Sequence[Monitor],
    rod_indices: Optional[Sequence[int]] = None,
) -> dict:
    """Return complex Jacobian dict without storing full dx fields.

    J[name][k] = d(monitor)/ds_k for requested rods (default: all).
    Reuses state.lu.
    """
    idx = list(range(len(state.rod_masks))) if rod_indices is None else list(rod_indices)
    sampler = ElementSampler(state.points, state.tris)
    J = {m.name: np.zeros(len(idx), np.complex128) for m in monitors}
    for j, k in enumerate(idx):
        dx = sensitivity_column(state, k)
        for m in monitors:
            J[m.name][j] = monitor_amp(sampler, dx, state.rho, state.k0, m.xy, m.weights)
    return {"rod_indices": idx, "J": J}


def compute_objective(state: ForwardState, objective: Objective) -> float:
    y = complex((objective.L @ state.x)[0])
    return float(abs(y) ** 2)


def compute_gradient(state: ForwardState, objective: Objective) -> np.ndarray:
    """Adjoint gradient of J=|L x|^2 w.r.t. all s_k. One A^H solve."""
    y = complex((objective.L @ state.x)[0])
    rhs = np.asarray(objective.L.conj().T @ np.array([y])).ravel()
    try:
        lam = state.lu.solve(rhs, trans="H")
    except TypeError:
        lam = splu(state.A.conj().T.tocsc()).solve(rhs)
    g = np.zeros(len(state.rod_masks), float)
    for k in range(len(state.rod_masks)):
        dA = assemble_dA_rod(
            state.points, state.tris, state.rod_masks[k],
            state.f_ord, state.fp_ref_ord, state.gamma, state.fc, state.k0,
            float(state.s_vec[k]),
        )
        g[k] = -2.0 * np.real(np.vdot(lam, dA @ state.x))
    return g


def parameter_bounds(n: int, s_max: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    lo = np.zeros(n, float)
    hi = np.full(n, s_max, float)
    return lo, hi
