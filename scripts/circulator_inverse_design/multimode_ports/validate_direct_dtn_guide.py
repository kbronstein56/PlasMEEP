#!/usr/bin/env python3
"""Straight-guide validation of true-boundary + direct DtN S."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_ports"
sys.path[:0] = [str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports")]

from boundary_ports import (  # noqa: E402
    add_modal_robin_from_c,
    build_port_from_boundary,
    modal_c_vector,
)
from direct_dtn_s import assemble_port_operator, build_Q_columns, direct_S_from_factor  # noqa: E402
from feed_constants import clear_width, fs_a, wall_thickness  # noqa: E402
from fem_local import assemble_anisotropic  # noqa: E402
from modes import mode_constants  # noqa: E402
from straight_guide_dtn import build_mesh, pec_wall_mask  # noqa: E402

THRESH = {"T_min": 0.995, "R_max": 0.015, "X_max": 0.02}


def solve_S(L=24.0, h=0.025, n_evan=2):
    W = float(clear_width)
    wall = float(wall_thickness)
    k0 = 2 * np.pi * fs_a
    points, tris = build_mesh(L, W, wall, h)
    T = len(tris)
    rho = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    wmask = pec_wall_mask(points, tris, W, wall)
    for arr in rho:
        arr[wmask] = 0.0
    A0 = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    modes = mode_constants(k0, W, n_modes=2 + n_evan)

    left = build_port_from_boundary(
        points, tris, rho, np.array([0.0, 0.0]), np.array([-1.0, 0.0]), np.array([0.0, 1.0]),
        W=W, plane_tol=max(1.5 * h, 0.02), port=0,
    )
    right = build_port_from_boundary(
        points, tris, rho, np.array([L, 0.0]), np.array([1.0, 0.0]), np.array([0.0, 1.0]),
        W=W, plane_tol=max(1.5 * h, 0.02), port=1,
    )
    ports = [left, right]
    A = assemble_port_operator(A0.copy(), points, ports, modes, n_prop=2, n_evan=n_evan)
    lu = splu(A.tocsc())
    cols = build_Q_columns(points, ports, modes, n_prop=2)
    S = direct_S_from_factor(lu, len(points), cols)
    geo = {
        "left": {"n_edges": left["n_edges"], "n_nodes": left["n_nodes"], "length": left["total_length"]},
        "right": {"n_edges": right["n_edges"], "n_nodes": right["n_nodes"], "length": right["total_length"]},
        "W": W,
    }
    return S, geo


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    S, geo = solve_S()
    # channels: (L,TE0),(L,TE1),(R,TE0),(R,TE1)
    tests = {
        "TE0_T": abs(S[2, 0]),
        "TE0_R": abs(S[0, 0]),
        "TE0_X1": max(abs(S[1, 0]), abs(S[3, 0])),
        "TE1_T": abs(S[3, 1]),
        "TE1_R": abs(S[1, 1]),
        "TE1_X0": max(abs(S[0, 1]), abs(S[2, 1])),
        "recip_B0": float(np.linalg.norm(S - S.T) / max(np.linalg.norm(S), 1e-30)),
    }
    ok = (
        tests["TE0_T"] >= THRESH["T_min"]
        and tests["TE0_R"] <= THRESH["R_max"]
        and tests["TE0_X1"] <= THRESH["X_max"]
        and tests["TE1_T"] >= THRESH["T_min"]
        and tests["TE1_R"] <= THRESH["R_max"]
        and tests["TE1_X0"] <= THRESH["X_max"]
        and tests["recip_B0"] <= 5e-3
        and abs(geo["left"]["length"] - geo["W"]) < 0.05 * geo["W"]
        and abs(geo["right"]["length"] - geo["W"]) < 0.05 * geo["W"]
    )
    report = {
        "thresholds": THRESH,
        "geometry": geo,
        "tests": tests,
        "S_abs": np.abs(S).tolist(),
        "PASS": bool(ok),
    }
    (OUT / "direct_dtn_guide_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
