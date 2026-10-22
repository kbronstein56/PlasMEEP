#!/usr/bin/env python3
"""Frozen mode-preserving (A) vs port-power (B) objectives.

Weights frozen BEFORE any optimization results. Do not retune after seeing J.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# FROZEN WEIGHTS (derived from production multimode objective magnitudes)
# ---------------------------------------------------------------------------
# Shared port-role weights (same numerical values as production Objective A):
W_DES = 1.0
W_REV = 2.0
W_THR = 1.5
W_LEAK = 1.0
W_REFL = 0.25

# Objective A only: penalize ALL outgoing TE1 (mode preservation)
W_TE1_A = 1.0

# Objective B: no separate TE1 penalty; TE1 at desired port is rewarded via total
W_TE1_B = 0.0

OBJECTIVE_FREEZE = {
    "W_DES": W_DES,
    "W_REV": W_REV,
    "W_THR": W_THR,
    "W_LEAK": W_LEAK,
    "W_REFL": W_REFL,
    "W_TE1_A": W_TE1_A,
    "W_TE1_B": W_TE1_B,
    "note": (
        "A = production TE0-preserving FoM. "
        "B applies the same port-role weights to TOTAL propagating power "
        "(TE0+TE1) per port role; TE1 at desired is not penalized."
    ),
}


def _port_powers(S: np.ndarray, j: int):
    """Return TE0/TE1 powers for des/rev/thr/other/refl for TE0 source j."""
    col = 2 * j
    des, rev, thr = (j + 1) % 6, (j - 1) % 6, (j + 3) % 6

    def te01(port):
        return float(abs(S[2 * port, col]) ** 2), float(abs(S[2 * port + 1, col]) ** 2)

    d0, d1 = te01(des)
    r0, r1 = te01(rev)
    t0, t1 = te01(thr)
    R0, R1 = te01(j)
    o0 = o1 = 0.0
    for p in range(6):
        if p in (j, des, rev, thr):
            continue
        a, b = te01(p)
        o0 += a
        o1 += b
    return {
        "des0": d0,
        "des1": d1,
        "rev0": r0,
        "rev1": r1,
        "thr0": t0,
        "thr1": t1,
        "other0": o0,
        "other1": o1,
        "R0": R0,
        "R1": R1,
    }


def objective_A_from_S(S: np.ndarray) -> float:
    """Mode-preserving: reward desired TE0; penalize all TE1 and bad ports."""
    js = []
    for j in range(6):
        p = _port_powers(S, j)
        te1_all = p["des1"] + p["rev1"] + p["thr1"] + p["other1"] + p["R1"]
        js.append(
            W_DES * p["des0"]
            - W_REV * p["rev0"]
            - W_THR * p["thr0"]
            - W_LEAK * p["other0"]
            - W_REFL * p["R0"]
            - W_TE1_A * te1_all
        )
    return float(np.mean(js))


def objective_B_from_S(S: np.ndarray) -> float:
    """Port-power: reward total desired power; no TE1-as-TE1 penalty."""
    js = []
    for j in range(6):
        p = _port_powers(S, j)
        js.append(
            W_DES * (p["des0"] + p["des1"])
            - W_REV * (p["rev0"] + p["rev1"])
            - W_THR * (p["thr0"] + p["thr1"])
            - W_LEAK * (p["other0"] + p["other1"])
            - W_REFL * (p["R0"] + p["R1"])
        )
    return float(np.mean(js))


def channel_weights_A(j_src: int) -> np.ndarray:
    w = np.zeros(12, float)
    des, rev, thr = (j_src + 1) % 6, (j_src - 1) % 6, (j_src + 3) % 6
    w[2 * des] = W_DES
    w[2 * rev] = -W_REV
    w[2 * thr] = -W_THR
    w[2 * j_src] = -W_REFL
    for p in range(6):
        if p not in (j_src, des, rev, thr):
            w[2 * p] = -W_LEAK
        w[2 * p + 1] = -W_TE1_A
    return w


def channel_weights_B(j_src: int) -> np.ndarray:
    """Same port roles on TE0 and TE1 channels; desired TE1 gets +W_DES."""
    w = np.zeros(12, float)
    des, rev, thr = (j_src + 1) % 6, (j_src - 1) % 6, (j_src + 3) % 6
    w[2 * des] = W_DES
    w[2 * des + 1] = W_DES
    w[2 * rev] = -W_REV
    w[2 * rev + 1] = -W_REV
    w[2 * thr] = -W_THR
    w[2 * thr + 1] = -W_THR
    w[2 * j_src] = -W_REFL
    w[2 * j_src + 1] = -W_REFL
    for p in range(6):
        if p not in (j_src, des, rev, thr):
            w[2 * p] = -W_LEAK
            w[2 * p + 1] = -W_LEAK
    return w


def mode_resolved_metrics(S: np.ndarray) -> dict:
    """Spectra-friendly aggregates with TE0/TE1 split at desired and reverse."""
    rows = []
    for j in range(6):
        p = _port_powers(S, j)
        rows.append(
            {
                "source": j,
                "des0": p["des0"],
                "des1": p["des1"],
                "des_tot": p["des0"] + p["des1"],
                "rev0": p["rev0"],
                "rev1": p["rev1"],
                "rev_tot": p["rev0"] + p["rev1"],
                "thr_tot": p["thr0"] + p["thr1"],
                "R_tot": p["R0"] + p["R1"],
                "iso_tot_dB": float(
                    10.0
                    * np.log10(
                        max(p["des0"] + p["des1"], 1e-30)
                        / max(p["rev0"] + p["rev1"], 1e-30)
                    )
                ),
                "iso_te0_dB": float(
                    10.0 * np.log10(max(p["des0"], 1e-30) / max(p["rev0"], 1e-30))
                ),
            }
        )

    def avg(k):
        return float(np.mean([r[k] for r in rows]))

    return {
        "per_source": rows,
        "des0_avg": avg("des0"),
        "des1_avg": avg("des1"),
        "des_tot_avg": avg("des_tot"),
        "rev0_avg": avg("rev0"),
        "rev1_avg": avg("rev1"),
        "rev_tot_avg": avg("rev_tot"),
        "iso_tot_avg_dB": avg("iso_tot_dB"),
        "iso_te0_avg_dB": avg("iso_te0_dB"),
        "J_A": objective_A_from_S(S),
        "J_B": objective_B_from_S(S),
    }


def write_freeze(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "MODE_OBJECTIVE_WEIGHTS_FROZEN.json").write_text(
        json.dumps(OBJECTIVE_FREEZE, indent=2) + "\n"
    )
