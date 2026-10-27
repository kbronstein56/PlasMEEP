#!/usr/bin/env python3
"""Engineering-aligned Objective A2 (proposed stage-2).

Mode-preserving TE0 circulator FoM that cannot hide T_des→0 by also
shrinking reverse power.

    J_A2 = log(T_des + ε) + log((T_des + ε)/(T_rev + ε))
         - W_THR T_thr - W_LEAK T_other - W_REFL R - W_TE1 TE1_total

Equivalent:
    J_A2 = 2 log(T_des + ε) - log(T_rev + ε) - linear leakage terms.

ε = 1e-6 (power units). No hard IL/isolation dB thresholds.

Weights on thr/other/R/TE1 match production Objective A so those channels
remain discouraged; the log terms replace the linear ±T_des/T_rev pair.

Maximize J_A2 (same convention as A).
"""
from __future__ import annotations

import numpy as np

# Match production A for non-desired/reverse absolute channels
W_THR, W_LEAK, W_REFL, W_TE1 = 1.5, 1.0, 0.25, 1.0
EPS = 1e-6

A2_FREEZE = {
    "name": "A2_hybrid",
    "eps": EPS,
    "formula": "2*log(Tdes+eps)-log(Trev+eps)-W_thr*Tthr-W_leak*Tother-W_refl*R-W_te1*TE1",
    "W_THR": W_THR,
    "W_LEAK": W_LEAK,
    "W_REFL": W_REFL,
    "W_TE1": W_TE1,
    "note": (
        "No IL/isolation hard thresholds. log(Tdes) forces useful transmission; "
        "log(Tdes/Trev) enforces isolation. Pathological Tdes≈Trev≈0 → J→−∞."
    ),
}


def powers_from_S_column(S: np.ndarray, j_src: int):
    col = 2 * j_src
    des, rev, thr = (j_src + 1) % 6, (j_src - 1) % 6, (j_src + 3) % 6
    Tdes = abs(S[2 * des, col]) ** 2
    Trev = abs(S[2 * rev, col]) ** 2
    Tthr = abs(S[2 * thr, col]) ** 2
    R = abs(S[col, col]) ** 2
    other = 0.0
    for p in range(6):
        if p in (j_src, des, rev, thr):
            continue
        other += abs(S[2 * p, col]) ** 2
    te1 = float(sum(abs(S[2 * i + 1, col]) ** 2 for i in range(6)))
    return {
        "Tdes": Tdes,
        "Trev": Trev,
        "Tthr": Tthr,
        "other": other,
        "R": R,
        "te1": te1,
        "des": des,
        "rev": rev,
        "thr": thr,
    }


def J_j_A2(S: np.ndarray, j_src: int, eps: float = EPS) -> float:
    p = powers_from_S_column(S, j_src)
    Td, Tr = p["Tdes"] + eps, p["Trev"] + eps
    return float(
        np.log(Td)
        + np.log(Td / Tr)
        - W_THR * p["Tthr"]
        - W_LEAK * p["other"]
        - W_REFL * p["R"]
        - W_TE1 * p["te1"]
    )


def objective_A2_from_S(S: np.ndarray, eps: float = EPS) -> float:
    return float(np.mean([J_j_A2(S, j, eps=eps) for j in range(6)]))


def channel_weights_A2(j_src: int, S: np.ndarray, eps: float = EPS) -> np.ndarray:
    """∂J_j/∂|S_m,col|² for the TE0 drive column — state-dependent.

    ∂/∂Tdes [2 log(Tdes+ε)] = 2/(Tdes+ε)
    ∂/∂Trev [-log(Trev+ε)] = -1/(Trev+ε)
    """
    p = powers_from_S_column(S, j_src)
    w = np.zeros(12, float)
    des, rev, thr = p["des"], p["rev"], p["thr"]
    w[2 * des] = 2.0 / (p["Tdes"] + eps)
    w[2 * rev] = -1.0 / (p["Trev"] + eps)
    w[2 * thr] = -W_THR
    w[2 * j_src] = -W_REFL
    for port in range(6):
        if port not in (j_src, des, rev, thr):
            w[2 * port] = -W_LEAK
        w[2 * port + 1] = -W_TE1
    return w


def scalar_adjoint_rhs_A2(C, betas, S, j_src: int, eps: float = EPS):
    """RHS for Aᵀ λ with A2 channel weights (same bilinear convention as production)."""
    n = 2 * j_src
    w = channel_weights_A2(j_src, S, eps=eps)
    r = np.zeros(C.shape[0], dtype=np.complex128)
    for m in range(12):
        if abs(w[m]) < 1e-30:
            continue
        sbeta = np.sqrt(betas[m] * betas[n])
        alpha = w[m] * np.conj(S[m, n]) * (-2j * sbeta)
        r += alpha * C[:, m]
    return r
