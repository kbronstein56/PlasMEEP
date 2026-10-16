#!/usr/bin/env python3
"""Efficient frequency-domain adjoint for the frozen TE0-drive multimode objective.

NEW campaign code. Does not modify frozen multimode_opt_core / port package.

Target per frequency:
  1 LU factorization
  6 TE0 forward RHS (multi-RHS when possible)
  6 scalar-objective transpose-adjoint RHS
  S TE0 columns built from U (no second pass through direct_S_from_factor)

Gradient mathematics equivalent to the validated channel-wise adjoint.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from adjoint_direct_dtn import (  # noqa: E402
    build_system,
    channel_weights,
    objective_from_S,
)
from diff_fem_core import assemble_dA_rod  # noqa: E402
from direct_dtn_s import build_Q_columns  # noqa: E402
from multimode_opt_core import (  # noqa: E402
    B_PLUS,
    F_HZ,
    expand_tied,
    orbits_and_masks,
    te0_metrics_from_S,
)


def _pack_C(points, cols):
    n_dof = len(points)
    n_ch = len(cols)
    C = np.zeros((n_dof, n_ch), dtype=np.complex128)
    betas = np.zeros(n_ch, dtype=np.complex128)
    for j, (idx, c, beta) in enumerate(cols):
        betas[j] = beta
        if len(idx):
            C[idx, j] = c
    return C, betas


def _solve_multi(lu, B: np.ndarray, trans: str = "N") -> np.ndarray:
    """Solve A X = B or A.T X = B. Prefer multi-RHS; fall back to columns."""
    try:
        X = lu.solve(B, trans=trans)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        return np.asarray(X, dtype=np.complex128)
    except Exception:
        X = np.zeros_like(B, dtype=np.complex128)
        for j in range(B.shape[1]):
            try:
                X[:, j] = lu.solve(B[:, j], trans=trans)
            except TypeError:
                if trans != "N":
                    raise
                X[:, j] = lu.solve(B[:, j])
        return X


def S_from_U(C, U_te0, betas, te0_cols):
    """Build full 12×12 S with only TE0 incident columns filled (others 0).

    S_mn = δ_mn + 2i √(β_m β_n) C_m · U_n  for n in te0_cols.
    Uses the same bilinear (non-sesquilinear) product as direct_S_from_factor.
    """
    n_ch = C.shape[1]
    S = np.zeros((n_ch, n_ch), dtype=np.complex128)
    for j_loc, n in enumerate(te0_cols):
        for m in range(n_ch):
            sbeta = np.sqrt(betas[m] * betas[n])
            q = complex(np.dot(C[:, m], U_te0[:, j_loc]))
            S[m, n] = (1.0 if m == n else 0.0) + 2j * sbeta * q
    return S


def scalar_adjoint_rhs(C, betas, S, j_src: int) -> np.ndarray:
    """RHS for Aᵀ λ = r with r = Σ_m α_m C[:,m], α chosen so

        ∂J_j/∂s contribution = 2 Re(λ · dA U_n) / 6
    matches the validated channel-wise formula.
    """
    n = 2 * j_src
    w = channel_weights(j_src)
    r = np.zeros(C.shape[0], dtype=np.complex128)
    for m in range(12):
        if abs(w[m]) < 1e-30:
            continue
        # α_m such that 2 Re(α_m V_m · dA U) = 2 Re(w_m conj(S) dS)
        # dS = -2i √(βm βn) V_m · dA U
        # so α_m = w_m conj(S_mn) (-2i √βmβn)
        sbeta = np.sqrt(betas[m] * betas[n])
        alpha = w[m] * np.conj(S[m, n]) * (-2j * sbeta)
        r += alpha * C[:, m]
    return r


def eval_and_grad_efficient(
    grade: str,
    q,
    orbits,
    f_hz: float = F_HZ,
    b_tesla: float = B_PLUS,
    measure: str = "projected",
    profile: bool = False,
):
    """Efficient TE0-drive objective + C6 gradient."""
    timings = {}
    s = expand_tied(np.asarray(q, float), orbits)

    t0 = time.perf_counter()
    st = build_system(grade, s, f_hz, b_tesla, measure)
    timings["assemble_s"] = time.perf_counter() - t0

    _, _, o_masks = orbits_and_masks(rod_masks=st["rod_masks"])
    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    C_te0 = C[:, te0_cols]

    t0 = time.perf_counter()
    lu = splu(A.tocsc())
    timings["factor_lu"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    U_te0 = _solve_multi(lu, C_te0, trans="N")
    timings["forward_te0_rhs"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    S = S_from_U(C, U_te0, betas, te0_cols)
    J = objective_from_S(S)
    met = te0_metrics_from_S(S)
    timings["S_and_objective"] = time.perf_counter() - t0

    # scalar adjoint RHS for each TE0 source → 6 transpose solves
    t0 = time.perf_counter()
    R = np.column_stack([scalar_adjoint_rhs(C, betas, S, j) for j in range(6)])
    Lam = _solve_multi(lu, R, trans="T")
    timings["adjoint_te0_rhs"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    g = np.zeros(len(orbits), float)
    for m, o in enumerate(orbits):
        s_m = float(s[o["rods"][0]])
        dA = assemble_dA_rod(
            points,
            st["tris"],
            o_masks[m],
            st["f_ord"],
            st["fp_ref_ord"],
            st["gamma"],
            st["fc"],
            st["k0"],
            s=s_m,
        )
        acc = 0.0
        for j_loc in range(6):
            dAu = dA @ U_te0[:, j_loc]
            # Bilinear form matching validated np.dot(V, dAu), not vdot.
            acc += 2.0 * np.real(np.dot(Lam[:, j_loc], dAu))
        g[m] = float(acc / 6.0)
    timings["dA_contraction"] = time.perf_counter() - t0
    timings["total"] = sum(timings.values())

    bundle = {
        "st": st,
        "lu": lu,
        "C": C,
        "U_te0": U_te0,
        "Lam": Lam,
        "betas": betas,
        "S": S,
        "metrics": met,
        "J": J,
        "timings": timings if profile else None,
        "n_forward_rhs": 6,
        "n_adjoint_rhs": 6,
        "n_lu": 1,
    }
    return bundle, g, s


def baseline_eval_and_grad(grade, q, orbits, f_hz=F_HZ, b_tesla=B_PLUS):
    """Frozen path for comparison (imports multimode_opt_core)."""
    from multimode_opt_core import eval_and_grad

    return eval_and_grad(grade, q, orbits, f_hz, b_tesla, with_absorption=False)
