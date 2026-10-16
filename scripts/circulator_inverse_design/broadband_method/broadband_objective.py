#!/usr/bin/env python3
"""Frozen differentiable broadband aggregation of per-frequency J(f).

See BROADBAND_DIFFERENTIABLE_OBJECTIVE.md for the freeze.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from efficient_adjoint import eval_and_grad_efficient  # noqa: E402
from multimode_opt_core import B_PLUS, expand_tied  # noqa: E402

# Frozen hyperparameters
TAU = 0.05
ALPHA = 0.25
F_OBJ_GHZ = [4.0 + 0.2 * i for i in range(11)]  # 4.0 .. 6.0
F_AUDIT_GHZ = [4.0 + 0.05 * i for i in range(41)]  # 4.0 .. 6.0


def softmin_bb(Jf: np.ndarray, tau: float = TAU, alpha: float = ALPHA):
    """Return J_BB and weights w_i = ∂J_BB/∂J_i.

    J_softmin = -τ log(mean_i exp(-J_i/τ)), with max-shift for stability.
    Softmax p_i = ∂J_softmin/∂J_i.
    """
    Jf = np.asarray(Jf, float)
    n = len(Jf)
    m = np.max(-Jf / tau)
    e = np.exp(-Jf / tau - m)
    J_soft = float(-tau * (m + np.log(np.mean(e))))
    p = e / np.sum(e)
    J_mean = float(np.mean(Jf))
    J_bb = (1.0 - alpha) * J_soft + alpha * J_mean
    w = (1.0 - alpha) * p + alpha / n
    return float(J_bb), w, J_soft, J_mean, p


def eval_bb_and_grad(
    grade: str,
    q,
    orbits,
    freqs_hz,
    b_tesla: float = B_PLUS,
    tau: float = TAU,
    alpha: float = ALPHA,
):
    """Broadband J_BB and ∂J_BB/∂q via chain rule on efficient per-freq grads."""
    Js = []
    gs = []
    mets = []
    for f in freqs_hz:
        bun, g, s = eval_and_grad_efficient(grade, q, orbits, f_hz=float(f), b_tesla=b_tesla)
        Js.append(bun["J"])
        gs.append(g)
        mets.append(bun["metrics"])
    Js = np.asarray(Js, float)
    J_bb, w, J_soft, J_mean, p = softmin_bb(Js, tau=tau, alpha=alpha)
    g_bb = sum(wi * gi for wi, gi in zip(w, gs))
    return {
        "J_BB": float(J_bb),
        "J_softmin": J_soft,
        "J_mean": J_mean,
        "J_per_freq": Js.tolist(),
        "weights": w.tolist(),
        "softmax_p": p.tolist(),
        "freqs_Hz": [float(f) for f in freqs_hz],
        "metrics_per_freq": mets,
        "grad_q": np.asarray(g_bb, float),
        "tau": tau,
        "alpha": alpha,
        "worst_J": float(np.min(Js)),
        "worst_f_Hz": float(freqs_hz[int(np.argmin(Js))]),
    }


def dense_audit(grade, q, orbits, freqs_hz=None, b_tesla=B_PLUS):
    """Forward-only J(f) on a dense grid (no gradients)."""
    from efficient_adjoint import eval_and_grad_efficient as eg

    if freqs_hz is None:
        freqs_hz = [f * 1e9 for f in F_AUDIT_GHZ]
    # cheaper: only need J — still uses efficient path but we could add eval-only
    rows = []
    for f in freqs_hz:
        bun, g, s = eg(grade, q, orbits, f_hz=float(f), b_tesla=b_tesla)
        m = bun["metrics"]
        rows.append(
            {
                "f_Hz": float(f),
                "f_GHz": float(f) / 1e9,
                "J": bun["J"],
                "T_des": m["T_des_avg"],
                "IL_dB": m["IL_avg_dB"],
                "isolation_dB": m["isolation_avg_dB"],
                "TE1": m["TE1_avg"],
                "RL_dB": m["RL_avg_dB"],
            }
        )
    Js = np.array([r["J"] for r in rows])
    return {
        "rows": rows,
        "worst_J": float(np.min(Js)),
        "worst_f_GHz": rows[int(np.argmin(Js))]["f_GHz"],
        "mean_J": float(np.mean(Js)),
    }
