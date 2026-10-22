#!/usr/bin/env python3
"""Efficient adjoint with pluggable channel_weights / objective_from_S."""
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
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from adjoint_direct_dtn import build_system  # noqa: E402
from diff_fem_core import assemble_dA_rod  # noqa: E402
from direct_dtn_s import build_Q_columns  # noqa: E402
from efficient_adjoint import S_from_U, _pack_C, _solve_multi  # noqa: E402
from multimode_opt_core import B_PLUS, F_HZ, orbits_and_masks  # noqa: E402


def scalar_adjoint_rhs_generic(C, betas, S, j_src, channel_weights_fn):
    n = 2 * j_src
    w = channel_weights_fn(j_src)
    r = np.zeros(C.shape[0], dtype=np.complex128)
    for m in range(12):
        if abs(w[m]) < 1e-30:
            continue
        sbeta = np.sqrt(betas[m] * betas[n])
        alpha = w[m] * np.conj(S[m, n]) * (-2j * sbeta)
        r += alpha * C[:, m]
    return r


def eval_and_grad_generic(
    grade,
    q,
    orbits,
    objective_from_S_fn,
    channel_weights_fn,
    expand_tied_fn,
    f_hz=F_HZ,
    b_tesla=B_PLUS,
    measure="projected",
):
    s = expand_tied_fn(np.asarray(q, float), orbits)
    st = build_system(grade, s, f_hz, b_tesla, measure)
    _, _, o_masks = orbits_and_masks(rod_masks=st["rod_masks"])
    # rebuild orbit masks for possibly different orbit sets
    from multimode_opt_core import orbits_and_masks as oam

    # o_masks from production assumes C6 orbits; rebuild manually
    o_masks = []
    for o in orbits:
        mask = np.zeros(len(st["tris"]), dtype=bool)
        for i in o["rods"]:
            mask |= st["rod_masks"][i]
        o_masks.append(mask)

    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    lu = splu(A.tocsc())
    U_te0 = _solve_multi(lu, C[:, te0_cols], trans="N")
    S = S_from_U(C, U_te0, betas, te0_cols)
    J = float(objective_from_S_fn(S))
    R = np.column_stack(
        [scalar_adjoint_rhs_generic(C, betas, S, j, channel_weights_fn) for j in range(6)]
    )
    Lam = _solve_multi(lu, R, trans="T")
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
            acc += 2.0 * np.real(np.dot(Lam[:, j_loc], dAu))
        g[m] = float(acc / 6.0)
    return {"J": J, "S": S, "st": st, "lu": lu}, g, s
