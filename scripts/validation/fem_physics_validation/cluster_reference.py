#!/usr/bin/env python3
"""Compare the independent cylinder T-matrix with saved body-fitted FEM fields.

The FEM files are the coated-bulb cluster solves from the earlier campaign.
This script does not refit anything. Truncation is reported by increasing m_max.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import hankel1

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, coated_T, mie_T, solve_clusters  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
RING = ROOT / "outputs" / "validation" / "fem_meep_validation" / "plasma_ring_campaign"
DISK = ROOT / "outputs" / "validation" / "fem_meep_validation" / "plasma_interface"
SRC = np.array([-6.0, 0.0])
R_CORE = 0.230
R_GAP = float(sc.r_bulb_inner)
R_SHELL = float(sc.r_bulb_outer)


def plasma(f_a: float) -> complex:
    eps, _ = gyrotropic_drude_eps_eta(float(f_a), float(sc.fp_a), float(sc.gamma_a), 0.0)
    return complex(eps)


def ratio_at(xy, centers, k0, radius, eps, nmax, coated=None):
    b = solve_clusters(centers, k0, radius, eps, SRC, nmax, coated=coated)
    field = cluster_field(xy, centers, k0, radius if coated is None else coated[0][-1], eps, SRC, b, nmax)
    vac = hankel1(0, k0 * np.linalg.norm(xy - SRC, axis=1))
    return field / vac, b


def metrics(fem, analytic):
    d = fem - analytic
    return {
        "abs_err": float(abs(d)),
        "db": float(20 * np.log10(abs(fem) / abs(analytic))) if abs(analytic) > 1e-8 else None,
        "phase_deg": float(np.angle(fem / analytic) * 180 / np.pi) if abs(analytic) > 1e-8 else None,
        "fem_abs": float(abs(fem)),
        "analytic_abs": float(abs(analytic)),
    }


def compare_cut(npz_path: Path, centers, coated, nmax: int) -> dict:
    z = np.load(npz_path)
    x = np.asarray(z["x"], float)
    fem = np.asarray(z["Hz_disk"], np.complex128) / np.asarray(z["Hz_vac"], np.complex128)
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    xy = np.column_stack([x, np.zeros_like(x)])
    analytic, _ = ratio_at(xy, centers, k0, R_CORE, eps, nmax, coated=coated)
    outside = np.isfinite(analytic)
    err = np.abs(fem[outside] - analytic[outside])
    scale = np.mean(np.abs(analytic[outside]))
    def at(xv):
        i = int(np.argmin(np.abs(x - xv)))
        return metrics(fem[i], analytic[i]) | {"x": float(x[i])}
    return {
        "file": npz_path.name,
        "nmax": nmax,
        "outside_rel_L2": float(np.sqrt(np.mean(err**2)) / scale),
        "outside_max_abs": float(np.max(err)),
        "forward": at(3.0),
        "backward": at(-3.0),
    }


def truncation(centers, coated, nmaxes):
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    xy = np.array([[3.0, 0.0], [-3.0, 0.0]])
    prev = None
    rows = []
    for nmax in nmaxes:
        analytic, _ = ratio_at(xy, centers, k0, R_CORE, eps, nmax, coated=coated)
        rec = {"nmax": nmax, "T": [analytic[0].real, analytic[0].imag], "R": [analytic[1].real, analytic[1].imag]}
        if prev is not None:
            rec["dT"] = float(abs(analytic[0] - prev[0]))
            rec["dR"] = float(abs(analytic[1] - prev[1]))
        prev = analytic
        rows.append(rec)
        print("trunc", rec, flush=True)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    centers7 = np.array(json.loads((RING / "fem_cluster7.json").read_text())[-1]["centers_a"])
    centers3 = np.array(json.loads((RING / "fem_cluster3.json").read_text())[-1]["centers_a"])
    coated = ([R_CORE, R_GAP, R_SHELL], [eps, 1.0 + 0j, 3.8 + 0j])
    # Confirm radii.
    print("radii", R_CORE, R_GAP, R_SHELL, "eps", eps, flush=True)
    bare_one = np.array([[0.0, 0.0]])
    # Frequency sweep of one bare cylinder, analytic only.
    sweep = []
    for f_GHz in np.linspace(3.4, 4.4, 21):
        f_a = f_GHz * 1e9 * 0.02 / 2.99792458e8
        kk = 2 * np.pi * f_a
        ee = plasma(f_a)
        T = mie_T(kk, ee, R_CORE, 10)
        sweep.append({"f_GHz": float(f_GHz), "max_abs_T": float(np.max(np.abs(T))), "eps_re": ee.real})
    out = {
        "eps_fs": [eps.real, eps.imag],
        "bare_truncation_one": truncation(bare_one, None, [2, 4, 6, 8, 10]),
        "coated_truncation_one": truncation(bare_one, coated, [2, 4, 6, 8, 10, 12]),
        "coated_truncation_seven": truncation(centers7, coated, [2, 4, 6, 8, 10, 12]),
        "bare_disk_fem": compare_cut(DISK / "fem_disk_FEM-F.npz", bare_one, None, 10),
        "coated_cluster3_fem": compare_cut(RING / "fem_cluster3_FEM-F.npz", centers3, coated, 10),
        "coated_cluster7_fem_C": compare_cut(RING / "fem_cluster7_FEM-C.npz", centers7, coated, 10),
        "coated_cluster7_fem_F": compare_cut(RING / "fem_cluster7_FEM-F.npz", centers7, coated, 10),
        "mie_sweep_bare": sweep,
    }
    # Two-cylinder orientation. Physical pair, two angles. Analytic only.
    pair0 = np.array([[0.0, 0.0], [1.0, 0.0]])
    ang = np.deg2rad(30.0)
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    pair30 = pair0 @ rot.T
    xy = np.array([[3.0, 0.0]])
    a0, _ = ratio_at(xy, pair0, k0, R_CORE, eps, 10, coated=None)
    a30, _ = ratio_at(xy, pair30, k0, R_CORE, eps, 10, coated=None)
    out["two_bare_orientation"] = {
        "axis_T": [a0[0].real, a0[0].imag],
        "rotated_30deg_T": [a30[0].real, a30[0].imag],
        "note": "These differ because the pair orientation relative to the source changed. A grid artifact would be an extra difference at fixed physical orientation.",
    }
    (OUT / "cluster_reference.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in out if "fem" in k or k.startswith("two")}, indent=2))
    print("CLUSTER_REF_DONE", flush=True)


if __name__ == "__main__":
    main()
