#!/usr/bin/env python3
"""Forward-only S / metrics eval (no adjoint, no dA) for spectral audits."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
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

from adjoint_direct_dtn import build_system  # noqa: E402
from direct_dtn_s import build_Q_columns  # noqa: E402
from efficient_adjoint import S_from_U, _pack_C, _solve_multi  # noqa: E402
from mode_objectives import mode_resolved_metrics  # noqa: E402
from multimode_opt_core import B_PLUS, expand_tied, te0_metrics_from_S  # noqa: E402


def eval_forward(
    grade: str,
    s: np.ndarray,
    f_hz: float,
    b_tesla: float = B_PLUS,
    measure: str = "projected",
):
    t0 = time.perf_counter()
    st = build_system(grade, s, f_hz, b_tesla, measure)
    t_asm = time.perf_counter() - t0
    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    t0 = time.perf_counter()
    lu = splu(A.tocsc())
    t_lu = time.perf_counter() - t0
    t0 = time.perf_counter()
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    t_rhs = time.perf_counter() - t0
    S = S_from_U(C, U, betas, te0_cols)
    met = te0_metrics_from_S(S)
    mres = mode_resolved_metrics(S)
    return {
        "S": S,
        "metrics": met,
        "mode_resolved": mres,
        "J_A": mres["J_A"],
        "J_B": mres["J_B"],
        "timings": {
            "assemble_s": t_asm,
            "factor_lu": t_lu,
            "forward_te0_rhs": t_rhs,
            "total": t_asm + t_lu + t_rhs,
        },
        "st": st,
        "lu": lu,
        "C": C,
        "U_te0": U,
        "betas": betas,
    }


def eval_q_forward(grade, q, orbits, f_hz, b_tesla=B_PLUS):
    s = expand_tied(np.asarray(q, float), orbits)
    return eval_forward(grade, s, f_hz, b_tesla), s
