#!/usr/bin/env python3
"""Direct DtN scattering: S = I + 2i √(β_m β_n) c_m^T A^{-1} c_n.

Frozen FEM pair (−iβ Robin, +2iβ incident RHS) yields 2iβ (c^T A^{-1} c) = −1
on a matched port, so S = I + 2i Q gives R=0, |T|=1 on the empty guide.
(The −I + 2i Q form belongs to the opposite RHS sign.)
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from boundary_ports import (  # noqa: E402
    add_incident_rhs_from_c,
    add_modal_robin_from_c,
    build_port_from_boundary,
    modal_c_vector,
)
from modes import mode_constants, power_scale  # noqa: E402


def channel_list(n_ports: int, n_prop: int = 2) -> list[tuple[int, int]]:
    """(port, mode) for propagating channels."""
    return [(p, m) for p in range(n_ports) for m in range(n_prop)]


def assemble_port_operator(
    A: sparse.spmatrix,
    points: np.ndarray,
    ports: Sequence[dict],
    modes: list[dict],
    n_prop: int = 2,
    n_evan: int = 2,
) -> sparse.spmatrix:
    for port in ports:
        if port["n_edges"] < 1:
            continue
        for md in modes[: n_prop + n_evan]:
            idx, c = modal_c_vector(points, port, md["m"])
            if len(idx) == 0:
                continue
            A = add_modal_robin_from_c(A, idx, c, md["beta"])
    return A


def build_Q_columns(
    points: np.ndarray,
    ports: Sequence[dict],
    modes: list[dict],
    n_prop: int = 2,
) -> tuple[list[np.ndarray], list[np.ndarray], list[complex]]:
    """For each propagating channel: (idx, c, beta)."""
    cols = []
    for p, port in enumerate(ports):
        for m in range(n_prop):
            idx, c = modal_c_vector(points, port, m)
            cols.append((idx, c, complex(modes[m]["beta"])))
    return cols


def direct_S_from_factor(
    lu,
    n_dof: int,
    channel_cols: list[tuple[np.ndarray, np.ndarray, complex]],
) -> np.ndarray:
    """S_mn = δ_mn + 2i √(β_m β_n) c_m^T A^{-1} c_n."""
    n = len(channel_cols)
    C = np.zeros((n_dof, n), dtype=np.complex128)
    betas = np.zeros(n, dtype=np.complex128)
    for j, (idx, c, beta) in enumerate(channel_cols):
        betas[j] = beta
        if len(idx):
            C[idx, j] = c
    X = np.zeros_like(C)
    for j in range(n):
        X[:, j] = lu.solve(C[:, j])
    S = np.zeros((n, n), dtype=np.complex128)
    for m in range(n):
        idx_m, c_m, _ = channel_cols[m]
        for n_ in range(n):
            if len(idx_m) == 0:
                cmT_xn = 0.0 + 0.0j
            else:
                cmT_xn = complex(np.dot(np.asarray(c_m, np.complex128), X[idx_m, n_]))
            sbeta = np.sqrt(betas[m] * betas[n_])
            S[m, n_] = 2j * sbeta * cmT_xn
            if m == n_:
                S[m, n_] += 1.0
    return S


def onsager_metrics(Sp: np.ndarray, Sm: np.ndarray) -> dict:
    D = Sp - Sm.T
    fro = float(np.linalg.norm(D) / max(np.linalg.norm(Sp), 1e-30))
    worst = float(np.max(np.abs(D)))
    return {"r_F": fro, "worst_abs_entry": worst, "PASS": fro <= 5e-3}
