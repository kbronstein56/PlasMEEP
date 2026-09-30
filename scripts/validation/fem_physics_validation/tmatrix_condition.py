#!/usr/bin/env python3
"""Raw and row-column-scaled condition numbers of the T-matrix reference.

Scaling is a reference-solver check. It is not fitted to FEM.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, solve_clusters  # noqa: E402
from analytic_sweeps import R_CORE, R_GAP, R_SHELL, SRC, plasma  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def seven():
    return np.array(
        [[0.0, 0.0], [-0.86602540378, 0.5], [-0.86602540378, -0.5], [0.86602540378, 0.5], [0.86602540378, -0.5], [0.0, 1.0], [0.0, -1.0]]
    )


def three():
    return np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])


def equilibrate(M, rhs):
    row = 1.0 / np.maximum(np.max(np.abs(M), axis=1), 1e-30)
    col = 1.0 / np.maximum(np.max(np.abs(M), axis=0), 1e-30)
    Ms = (row[:, None] * M) * col[None, :]
    ys = np.linalg.solve(Ms, row * rhs)
    xs = col * ys
    return Ms, xs


def one_case(name, centers, coated, nmaxes):
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    coat = ([R_CORE, R_GAP, R_SHELL], [eps, 1.0 + 0j, 3.8 + 0j]) if coated else None
    outer = R_SHELL if coated else R_CORE
    probes = {"forward": np.array([[2.8, 0.0]]), "near": np.array([[outer + 0.06, 0.0]])}
    rows = []
    prev = None
    for nmax in nmaxes:
        b = solve_clusters(centers, k0, R_CORE, eps, SRC, nmax, coated=coat)
        M = solve_clusters.last_M
        rhs = solve_clusters.last_rhs
        Ms, xs = equilibrate(M, rhs)
        b_s = xs.reshape(b.shape)
        obs = {}
        obs_s = {}
        for key, xy in probes.items():
            obs[key] = complex(cluster_field(xy, centers, k0, outer, eps, SRC, b, nmax)[0])
            obs_s[key] = complex(cluster_field(xy, centers, k0, outer, eps, SRC, b_s, nmax)[0])
        rec = {
            "case": name,
            "m_max": nmax,
            "dim": int(M.shape[0]),
            "cond_raw": float(np.linalg.cond(M)),
            "cond_scaled": float(np.linalg.cond(Ms)),
            "residual_raw": solve_clusters.last_residual,
            "residual_scaled": float(np.linalg.norm(M @ xs - rhs) / (np.linalg.norm(rhs) + 1e-30)),
            "forward": [obs["forward"].real, obs["forward"].imag],
            "near": [obs["near"].real, obs["near"].imag],
            "scale_changes_forward": float(abs(obs["forward"] - obs_s["forward"])),
            "scale_changes_near": float(abs(obs["near"] - obs_s["near"])),
        }
        if prev is not None:
            rec["d_forward"] = float(abs(obs["forward"] - prev[0]))
            rec["d_near"] = float(abs(obs["near"] - prev[1]))
        prev = (obs["forward"], obs["near"])
        rows.append(rec)
        print(name, nmax, "cond", f"{rec['cond_raw']:.3e}", "scaled", f"{rec['cond_scaled']:.3e}", "dF", rec.get("d_forward"), "resid", rec["residual_raw"], flush=True)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nmaxes = [2, 3, 4, 5, 6, 8, 10, 12, 14]
    rows = []
    rows += one_case("bare_1", np.zeros((1, 2)), False, nmaxes)
    rows += one_case("bare_3", three(), False, nmaxes)
    rows += one_case("bare_7", seven(), False, nmaxes)
    rows += one_case("coated_1", np.zeros((1, 2)), True, nmaxes)
    rows += one_case("coated_7", seven(), True, nmaxes)
    (OUT / "tmatrix_condition.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("COND_DONE", flush=True)


if __name__ == "__main__":
    main()
