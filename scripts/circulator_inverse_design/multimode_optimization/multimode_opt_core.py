#!/usr/bin/env python3
"""Frozen multimode direct-DtN objective, metrics, and C6 adjoint gradient.

Do not change weights, S definition, or projected measure without a new
validation campaign. See MULTIMODE_SINGLE_FREQ_OBJECTIVE_FROZEN.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
DES = ROOT / "outputs" / "circulator_inverse_design"
OUT = DES / "multimode_optimization"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from adjoint_direct_dtn import (  # noqa: E402
    S_HI,
    S_LO,
    W_DES,
    W_LEAK,
    W_REFL,
    W_REV,
    W_TE1,
    W_THR,
    build_system,
    channel_weights,
    objective_from_S,
)
from boundary_ports import add_incident_rhs_from_c, modal_c_vector  # noqa: E402
from diff_fem_core import assemble_dA_rod  # noqa: E402
from direct_dtn_s import build_Q_columns, direct_S_from_factor, onsager_metrics  # noqa: E402
from full91_direct_dtn import load_s  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402
from modes import power_scale  # noqa: E402
from sixfold_orbits import build_orbits, expand_tied, reduce_gradient  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402

F_HZ = 3.85e9
B_PLUS = 0.05
ONSAGER_BAR = 5e-3
BALANCE_BAR = 0.05
WEIGHTS = {
    "des": W_DES,
    "rev": W_REV,
    "thr": W_THR,
    "leak": W_LEAK,
    "refl": W_REFL,
    "te1": W_TE1,
}


def orbits_and_masks(points=None, tris=None, rod_masks=None):
    centers = bulb_centers_device(sc)
    orbits, origin = build_orbits(centers)
    o_masks = None
    if rod_masks is not None:
        o_masks = []
        for o in orbits:
            m = np.zeros(len(rod_masks[0]), dtype=bool)
            for i in o["rods"]:
                m |= rod_masks[i]
            o_masks.append(m)
    return orbits, origin, o_masks


def q_from_s(s: np.ndarray, orbits) -> np.ndarray:
    q = np.zeros(len(orbits), float)
    for m, o in enumerate(orbits):
        vals = [float(s[i]) for i in o["rods"]]
        if max(vals) - min(vals) > 1e-8:
            # free-91 or broken tying — take mean of orbit
            q[m] = float(np.mean(vals))
        else:
            q[m] = vals[0]
    return q


def load_start(tag: str, orbits):
    if tag == "uniform":
        q = np.ones(len(orbits), float)
        s = expand_tied(q, orbits)
        return q, s, tag
    s = load_s(tag)
    q = q_from_s(s, orbits)
    # re-expand so tying is exact
    s = expand_tied(q, orbits)
    return q, s, tag


def te0_metrics_from_S(S: np.ndarray) -> dict:
    """Per-source and aggregate TE0 engineering metrics from 12×12 S."""
    per = []
    for j in range(6):
        col = 2 * j
        des, rev, thr = (j + 1) % 6, (j - 1) % 6, (j + 3) % 6
        Tdes = float(abs(S[2 * des, col]) ** 2)
        Trev = float(abs(S[2 * rev, col]) ** 2)
        Tthr = float(abs(S[2 * thr, col]) ** 2)
        R = float(abs(S[col, col]) ** 2)
        other = 0.0
        for p in range(6):
            if p in (j, des, rev, thr):
                continue
            other += float(abs(S[2 * p, col]) ** 2)
        te1 = float(sum(abs(S[2 * i + 1, col]) ** 2 for i in range(6)))
        P_out = float(np.sum(np.abs(S[:, col]) ** 2))
        IL = float(-10.0 * np.log10(max(Tdes, 1e-30)))
        iso = float(10.0 * np.log10(max(Tdes, 1e-30) / max(Trev, 1e-30)))
        RL = float(-10.0 * np.log10(max(R, 1e-30)))
        Jj = (
            W_DES * Tdes
            - W_REV * Trev
            - W_THR * Tthr
            - W_LEAK * other
            - W_REFL * R
            - W_TE1 * te1
        )
        per.append(
            {
                "source": j,
                "T_des": Tdes,
                "T_rev": Trev,
                "T_thr": Tthr,
                "T_other": other,
                "R": R,
                "TE1_total": te1,
                "P_out": P_out,
                "IL_dB": IL,
                "isolation_dB": iso,
                "RL_dB": RL,
                "J_j": float(Jj),
            }
        )

    def avg(key):
        return float(np.mean([r[key] for r in per]))

    def mn(key):
        return float(np.min([r[key] for r in per]))

    def mx(key):
        return float(np.max([r[key] for r in per]))

    return {
        "per_source": per,
        "T_des_avg": avg("T_des"),
        "T_des_min": mn("T_des"),
        "T_rev_avg": avg("T_rev"),
        "T_thr_avg": avg("T_thr"),
        "T_other_avg": avg("T_other"),
        "R_avg": avg("R"),
        "TE1_avg": avg("TE1_total"),
        "P_out_avg": avg("P_out"),
        "IL_avg_dB": avg("IL_dB"),
        "IL_worst_dB": mx("IL_dB"),
        "isolation_avg_dB": avg("isolation_dB"),
        "isolation_min_dB": mn("isolation_dB"),
        "RL_avg_dB": avg("RL_dB"),
        "RL_min_dB": mn("RL_dB"),
        "J": float(np.mean([r["J_j"] for r in per])),
        "J_std": float(np.std([r["J_j"] for r in per])),
    }


def _channel_basis(st, need_adjoint: bool = True):
    points, A, ports, modes = st["points"], st["A"], st["ports"], st["modes"]
    lu = splu(A.tocsc())
    cols = build_Q_columns(points, ports, modes, n_prop=2)
    n_dof = len(points)
    n_ch = 12
    C = np.zeros((n_dof, n_ch), dtype=np.complex128)
    betas = np.zeros(n_ch, dtype=np.complex128)
    for j, (idx, c, beta) in enumerate(cols):
        betas[j] = beta
        if len(idx):
            C[idx, j] = c
    U = np.zeros_like(C)
    for j in range(n_ch):
        U[:, j] = lu.solve(C[:, j])
    V = None
    if need_adjoint:
        # Reuse the same SuperLU factor for Aᵀ solves (avoid a second factorization).
        V = np.zeros_like(C)
        for m in range(n_ch):
            try:
                V[:, m] = lu.solve(C[:, m], trans="T")
            except TypeError:
                if "_luT_cache" not in st:
                    st["_luT_cache"] = splu(A.T.tocsc())
                V[:, m] = st["_luT_cache"].solve(C[:, m])
    S = direct_S_from_factor(lu, n_dof, cols)
    return lu, cols, C, U, V, betas, S


def te0_absorption_mean(st, lu) -> float:
    """Mean material absorption over 6 unit-power TE0 drives."""
    points, tris, rho, ports, modes, k0 = (
        st["points"],
        st["tris"],
        st["rho"],
        st["ports"],
        st["modes"],
        st["k0"],
    )
    vals = []
    for j in range(6):
        idx, c = modal_c_vector(points, ports[j], 0)
        scale = power_scale(modes[0]["beta"], k0)
        a_field = 1.0 / max(abs(scale), 1e-30)
        b = np.zeros(len(points), dtype=np.complex128)
        b = add_incident_rhs_from_c(b, idx, c, modes[0]["beta"], a_field)
        x = lu.solve(b)
        vals.append(float(material_absorption(points, tris, x, rho, k0)))
    return float(np.mean(vals))


def eval_design(grade: str, s_vec: np.ndarray, f_hz: float = F_HZ, b_tesla: float = B_PLUS,
                measure: str = "projected", with_absorption: bool = True):
    st = build_system(grade, s_vec, f_hz, b_tesla, measure)
    lu, cols, C, U, V, betas, S = _channel_basis(st, need_adjoint=False)
    met = te0_metrics_from_S(S)
    abs_mean = te0_absorption_mean(st, lu) if with_absorption else None
    if abs_mean is not None:
        met["P_material_avg"] = abs_mean
        met["rem_avg"] = float(1.0 - met["P_out_avg"] - abs_mean)
    return {
        "st": st,
        "lu": lu,
        "cols": cols,
        "C": C,
        "U": U,
        "V": V,
        "betas": betas,
        "S": S,
        "metrics": met,
        "J": met["J"],
    }


def gradient_q_adjoint(bundle, orbits, o_masks, s_vec) -> np.ndarray:
    """∂J/∂q_m via one dA per orbit (validated adjoint formula)."""
    st = bundle["st"]
    S, U, V, betas = bundle["S"], bundle["U"], bundle["V"], bundle["betas"]
    points = st["points"]
    g = np.zeros(len(orbits), float)
    for m, o in enumerate(orbits):
        s_m = float(s_vec[o["rods"][0]])
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
        dAu = {}
        for j_src in range(6):
            n = 2 * j_src
            dAu[n] = dA @ U[:, n]

        def dS_adj(m_ch, n, _dAu=dAu, _V=V, _b=betas):
            sbeta = np.sqrt(_b[m_ch] * _b[n])
            return -2j * sbeta * complex(np.dot(_V[:, m_ch], _dAu[n]))

        acc = 0.0
        for j_src in range(6):
            n = 2 * j_src
            w = channel_weights(j_src)
            for mch in range(12):
                if abs(w[mch]) < 1e-30:
                    continue
                dS = dS_adj(mch, n)
                acc += w[mch] * 2.0 * np.real(np.conj(S[mch, n]) * dS)
        g[m] = float(acc / 6.0)
    return g


def eval_and_grad(grade, q, orbits, f_hz=F_HZ, b_tesla=B_PLUS, measure="projected",
                  with_absorption=True):
    s = expand_tied(np.asarray(q, float), orbits)
    # build system once; need o_masks from rod_masks
    st0 = build_system(grade, s, f_hz, b_tesla, measure)
    _, _, o_masks = orbits_and_masks(rod_masks=st0["rod_masks"])
    lu, cols, C, U, V, betas, S = _channel_basis(st0)
    met = te0_metrics_from_S(S)
    if with_absorption:
        met["P_material_avg"] = te0_absorption_mean(st0, lu)
        met["rem_avg"] = float(1.0 - met["P_out_avg"] - met["P_material_avg"])
    bundle = {
        "st": st0,
        "lu": lu,
        "cols": cols,
        "C": C,
        "U": U,
        "V": V,
        "betas": betas,
        "S": S,
        "metrics": met,
        "J": met["J"],
    }
    gq = gradient_q_adjoint(bundle, orbits, o_masks, s)
    return bundle, gq, s


def dump_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


__all__ = [
    "F_HZ",
    "B_PLUS",
    "S_LO",
    "S_HI",
    "ONSAGER_BAR",
    "BALANCE_BAR",
    "WEIGHTS",
    "OUT",
    "DES",
    "ROOT",
    "orbits_and_masks",
    "load_start",
    "expand_tied",
    "q_from_s",
    "te0_metrics_from_S",
    "eval_design",
    "eval_and_grad",
    "objective_from_S",
    "dump_json",
    "load_s",
    "onsager_metrics",
]
