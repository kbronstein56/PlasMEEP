#!/usr/bin/env python3
"""Multipole convergence and resonance sensitivity of the coated cylinder at 5.30 GHz."""
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
from analytic_maxwell import cluster_field, coated_T, solve_clusters  # noqa: E402
from analytic_sweeps import R_CORE, R_GAP, R_SHELL, plasma, probes  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
SRC = np.array([-4.5, 0.0])


def forward_at(ghz: float, nmax: int):
    k0 = 2 * np.pi * float(sc.fs_a) * ghz / 3.85
    eps = plasma(float(sc.fs_a) * ghz / 3.85)
    layers = [eps, 1.0 + 0j, 3.8 + 0j]
    radii = [R_CORE, R_GAP, R_SHELL]
    b = solve_clusters(np.zeros((1, 2)), k0, R_CORE, eps, SRC, nmax, coated=(radii, layers))
    xy = np.array([[2.8, 0.0]])
    from scipy.special import hankel1

    field = cluster_field(xy, np.zeros((1, 2)), k0, R_SHELL, eps, SRC, b, nmax)[0]
    vac = hankel1(0, k0 * np.linalg.norm(xy[0] - SRC))
    ratio = field / vac
    T = coated_T(k0, radii, layers, nmax)
    return {
        "ghz": ghz,
        "nmax": nmax,
        "eps": [eps.real, eps.imag],
        "k0": k0,
        "forward": [ratio.real, ratio.imag],
        "forward_abs": float(abs(ratio)),
        "T_max_abs": float(np.max(np.abs(T))),
        "T_last_abs": float(abs(T[-1])),
        "residual": float(solve_clusters.last_residual),
        "cond": float(solve_clusters.last_cond),
    }


def main():
    rows = [forward_at(5.30, n) for n in (4, 6, 8, 10, 12, 14, 16)]
    base = complex(*rows[-1]["forward"])
    for rec, prev in zip(rows, [None] + rows[:-1]):
        cur = complex(*rec["forward"])
        rec["d_forward_from_nmax16"] = float(abs(cur - base))
        rec["d_forward_from_previous"] = None if prev is None else float(abs(cur - complex(*prev["forward"])))
    # Local frequency slope of the forward phase. A large slope means the phase
    # tolerance is sensitive, which is reported and not used to change the threshold.
    neigh = [forward_at(f, 12) for f in (5.20, 5.25, 5.30, 5.35, 5.40)]
    phases = [np.angle(complex(*r["forward"])) * 180 / np.pi for r in neigh]
    mags_db = [20 * np.log10(r["forward_abs"]) for r in neigh]
    dph = (phases[3] - phases[1]) / 0.10
    ddb = (mags_db[3] - mags_db[1]) / 0.10
    out = {
        "multipole": rows,
        "frequency": neigh,
        "phase_deg_per_GHz": float(dph),
        "db_per_GHz": float(ddb),
        "note": "One-cylinder T is diagonal. cond is 1. Residual is the dense solve residual.",
    }
    (OUT / "coated_530_analytic.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"phase_deg_per_GHz": dph, "db_per_GHz": ddb, "d16": rows[-1]["d_forward_from_nmax16"], "Tlast16": rows[-1]["T_last_abs"], "cond8": rows[2]["cond"]}, indent=2))
    print("COATED_530_ANALYTIC_DONE", flush=True)


if __name__ == "__main__":
    main()
