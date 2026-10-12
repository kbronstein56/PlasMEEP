#!/usr/bin/env python3
"""Modal Robin (first-order DtN) contributions for parallel-plate ports.

Weak form (e^{-iωt}), after integration by parts on the scalar Hz equation:

    ∫ ∇Hz·∇v − k0² Hz v  +  ∫_Γ (∂n Hz) v  = 0

First-order absorbing boundary for mode m with outward normal n and
propagation constant β_m (>0 for propagating):

    ∂n Hz − i β_m Hz = 0     (pure outgoing at that port)

moved to the left-hand side gives a boundary contribution

    + ∫_Γ (−i β_m) Hz v ds   (frozen sign; see MULTIMODE_DTN_WEAK_FORM.md)

Modal Galerkin on a port segment Γ (piecewise-linear Hz):

    Hz|_Γ ≈ Σ_n u_n N_n ,   test v = N_p
    add  (i β_m) ∫_Γ φ_m N_p ds ∫_Γ φ_m N_q ds  to A[p,q]

Prescribed unit-power incident mode m at a port with outward n and
wave traveling into the domain along +t (guide axis):

    At x = x_p with axis +x and outward n = −x̂ (left end):
        incident travels +x :  H_inc = φ_m e^{+iβ x}
        ∂n H = −∂x H = −iβ H_inc  at the boundary
    First-order:  ∂n H + iβ H = 2 iβ a_m φ_m  (with a_m the incident modal amp)

contributes to RHS (weak form):  b_p += ∫ 2 iβ a_m φ_m N_p ds  (sign folded into API)

Sign convention frozen in MULTIMODE_DTN_WEAK_FORM.md — see `side` argument.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
from scipy import sparse


def port_line_integrals(
    node_indices: Sequence[int],
    y_coords: np.ndarray,
    weights: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    """c_n = ∫ φ(y) N_n ds with hat functions on a 1D port mesh (2-node elements)."""
    node_indices = list(node_indices)
    y_coords = np.asarray(y_coords, float)
    phi = np.asarray(phi, float)
    weights = np.asarray(weights, float)
    c = np.zeros(len(node_indices), dtype=np.complex128)
    # Assume ordered nodes along port; each weight segment between consecutive samples
    # belongs to both adjacent hats with trapezoid split.
    if len(node_indices) != len(y_coords):
        raise ValueError("node_indices and y_coords length mismatch")
    n = len(node_indices)
    if n == 1:
        c[0] = weights[0] * phi[0]
        return c
    for s in range(n - 1):
        ds = float(weights[s]) if s < len(weights) else float(weights[-1])
        ph = 0.5 * (phi[s] + phi[s + 1])
        c[s] += 0.5 * ds * ph
        c[s + 1] += 0.5 * ds * ph
    return c


def add_modal_robin(
    A: sparse.spmatrix,
    node_indices: Sequence[int],
    y_coords: np.ndarray,
    weights: np.ndarray,
    phi: np.ndarray,
    beta: complex,
) -> sparse.spmatrix:
    """Add i β (c c^T) to A for one propagating mode."""
    c = port_line_integrals(node_indices, y_coords, weights, phi)
    beta = complex(beta)
    if abs(beta.imag) > 1e-12 and beta.imag > 0:
        # Evanescent: use iα with decay — first-order outgoing: ∂n Hz + α Hz = 0
        beta_eff = 1j * beta.imag if beta.imag > 0 else beta
    else:
        beta_eff = -1j * float(beta.real)
    block = beta_eff * np.outer(c, c)
    idx = np.asarray(node_indices, dtype=int)
    A = A.tolil()
    for i, p in enumerate(idx):
        for j, q in enumerate(idx):
            A[p, q] = A[p, q] + block[i, j]
    return A.tocsr()


def add_modal_incident_rhs(
    b: np.ndarray,
    node_indices: Sequence[int],
    y_coords: np.ndarray,
    weights: np.ndarray,
    phi: np.ndarray,
    beta: complex,
    amplitude: complex,
    *,
    side: str,
) -> np.ndarray:
    """Add incident-mode forcing for unit incoming modal amplitude.

    side='left'  : port at min-x, outward n = −x̂, incident travels +x into domain
    side='right' : port at max-x, outward n = +x̂, incident travels −x into domain
    """
    c = port_line_integrals(node_indices, y_coords, weights, phi)
    beta = complex(beta)
    if side == "left":
        # Prescribed +x-traveling incident at x=0 (outward n = −x̂); see MULTIMODE_DTN_WEAK_FORM.md
        coeff = 2j * beta * complex(amplitude)
    elif side == "right":
        # Prescribed −x-traveling incident at x=L (outward n = +x̂)
        coeff = 2j * beta * complex(amplitude)
    else:
        raise ValueError(f"unknown side {side}")
    idx = np.asarray(node_indices, dtype=int)
    for i, p in enumerate(idx):
        b[p] += coeff * c[i]
    return b
