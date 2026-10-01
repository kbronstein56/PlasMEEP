#!/usr/bin/env python3
"""Two- and three-cylinder gyrotropic scattering against an independent T-matrix."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, solve_clusters  # noqa: E402
from fem_scatterers import mesh_circles, to_mesh  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from gyrotropic_fem import exterior_field, gyrotropic_T, rho_tuple  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_gyrotropic_validation"


def solve_gyro(centers, k0, radius, rho_xx, rho_xy, source, nmax):
    centers = np.asarray(centers, float)
    n_obj = len(centers)
    ns = np.arange(-nmax, nmax + 1)
    n_mode = len(ns)
    T = gyrotropic_T(k0, radius, rho_xx, rho_xy, nmax)
    dim = n_obj * n_mode
    M = np.eye(dim, dtype=np.complex128)
    rhs = np.zeros(dim, dtype=np.complex128)
    for l in range(n_obj):
        d = source - centers[l]
        R = np.linalg.norm(d)
        phi = np.arctan2(d[1], d[0])
        for im, m in enumerate(ns):
            ext = hankel_coeff(int(m), k0, R, phi)
            rhs[l * n_mode + im] = T[im] * ext
    for l in range(n_obj):
        for j in range(n_obj):
            if j == l:
                continue
            d = centers[j] - centers[l]
            R = np.linalg.norm(d)
            phi = np.arctan2(d[1], d[0])
            for im, m in enumerate(ns):
                row = l * n_mode + im
                for inn, n in enumerate(ns):
                    G = hankel_coeff(int(m - n), k0, R, phi)
                    M[row, j * n_mode + inn] -= T[im] * G
    b = np.linalg.solve(M, rhs).reshape(n_obj, n_mode)
    resid = float(np.linalg.norm(M @ b.ravel() - rhs) / (np.linalg.norm(rhs) + 1e-30))
    return b, resid


def hankel_coeff(n, k0, R, phi):
    from scipy.special import hankel1

    return hankel1(n, k0 * R) * np.exp(-1j * n * phi)


def fem_ratio(pts, tris, k0, centers, radius, rho, src, xy):
    origin = pts[tris].mean(1) - np.array([7.0, 5.0])
    inside = np.zeros(len(tris), dtype=bool)
    for c in centers:
        inside |= np.linalg.norm(origin - c, axis=1) <= radius
    T = len(tris)
    rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    vac = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    for k in range(4):
        rho_e[k][inside] = rho[k]
    sampler = ElementSampler(pts, tris)
    src_m = to_mesh(src)
    bvec = np.zeros(len(pts), np.complex128)
    t_src = int(sampler.locate(src_m[None, :])[0])
    w = sampler._bary(t_src, src_m)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    out = []
    for rr in (vac, rho_e):
        A = assemble_anisotropic(pts, tris, *rr, k0, np.zeros(len(pts), dtype=bool))
        uh = splu(A.tocsc()).solve(bvec)
        hz, _, _ = sampler.fields(uh, *rr, k0, to_mesh(xy))
        out.append(hz)
    return out[1] / out[0]


def main():
    from scipy.special import hankel1

    OUT.mkdir(parents=True, exist_ok=True)
    f, fp, gamma = 0.45, 0.15, 1.0e-4
    radius = 0.230
    k0 = 2 * np.pi * f
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 14.0, 10.0, 1.2
    sc.fs_a = f
    src = np.array([-4.5, 0.0])
    probes = {
        "forward": np.array([2.8, 0.0]),
        "backward": np.array([-2.8, 0.0]),
        "side_p60": np.array([2.8 * np.cos(np.pi / 3), 2.8 * np.sin(np.pi / 3)]),
        "side_m60": np.array([2.8 * np.cos(-np.pi / 3), 2.8 * np.sin(-np.pi / 3)]),
    }
    names = list(probes)
    xy = np.array(list(probes.values()), float)
    geometries = {
        "two": np.array([[0.0, 0.0], [1.15, 0.35]]),
        "three": np.array([[0.0, 0.0], [1.15, 0.35], [-0.70, 0.95]]),
        "two_flipped": np.array([[0.0, 0.0], [1.15, -0.35]]),
    }
    # Analytic Onsager and B=0 reduction, no FEM.
    checks = []
    for name, centers in geometries.items():
        rho0, _, _ = rho_tuple(f, fp, gamma, 0.0)
        eps0 = 1.0 / rho0[0]
        b_new, resid_new = solve_gyro(centers, k0, radius, rho0[0], rho0[1], src, 8)
        b_old = solve_clusters(centers, k0, radius, complex(eps0), src, 8)
        field_new = cluster_field(xy, centers, k0, radius, eps0, src, b_new, 8)
        field_old = cluster_field(xy, centers, k0, radius, eps0, src, b_old, 8)
        rho_p, _, _ = rho_tuple(f, fp, gamma, 0.08)
        rho_m, _, _ = rho_tuple(f, fp, gamma, -0.08)
        b_p, _ = solve_gyro(centers, k0, radius, rho_p[0], rho_p[1], src, 10)
        # Green: field at probe from src at +B versus field at src from probe at -B.
        ons = []
        for p in xy:
            fp_b, _ = solve_gyro(centers, k0, radius, rho_p[0], rho_p[1], src, 10)
            fm_b, _ = solve_gyro(centers, k0, radius, rho_m[0], rho_m[1], p, 10)
            g_ps = cluster_field(p[None, :], centers, k0, radius, 1.0, src, fp_b, 10)[0]
            g_sp = cluster_field(src[None, :], centers, k0, radius, 1.0, p, fm_b, 10)[0]
            ons.append(abs(g_ps - g_sp) / abs(g_ps))
        field_p = cluster_field(xy, centers, k0, radius, 1.0, src, b_p, 10)
        checks.append({
            "geometry": name,
            "b0_cluster_resid": float(np.max(np.abs(field_new - field_old))),
            "tmatrix_residual": resid_new,
            "onsager_rel_max": float(max(ons)),
            "side_abs_p60": abs(field_p[2]),
            "side_abs_m60": abs(field_p[3]),
        })
        print("analytic", checks[-1], flush=True)
    # Single-cylinder side swap, analytic.
    rho_p, _, _ = rho_tuple(f, fp, gamma, 0.05)
    rho_m, _, _ = rho_tuple(f, fp, gamma, -0.05)
    hp = exterior_field(xy, np.zeros(2), k0, radius, rho_p[0], rho_p[1], src, 12)
    hm = exterior_field(xy, np.zeros(2), k0, radius, rho_m[0], rho_m[1], src, 12)
    h0 = exterior_field(xy, np.zeros(2), k0, radius, rho_tuple(f, fp, gamma, 0.0)[0][0], 0.0, src, 12)
    side = {
        "p60_plus_vs_m60_minus": abs(hp[2] - hm[3]) / abs(hp[2]),
        "m60_plus_vs_p60_minus": abs(hp[3] - hm[2]) / abs(hp[3]),
        "b0_side_symmetry": abs(h0[2] - h0[3]) / abs(h0[2]),
        "plus_side_asymmetry": abs(hp[2] - hp[3]) / abs(hp[2]),
    }
    print("side", side, flush=True)
    rows = []
    jobs = [
        ("two", 0.0, 0.04, 0.012),
        ("two", 0.08, 0.04, 0.012),
        ("two", -0.08, 0.04, 0.012),
        ("two", 0.08, 0.02, 0.008),
        ("three", 0.0, 0.04, 0.012),
        ("three", 0.08, 0.04, 0.012),
        ("three", -0.08, 0.04, 0.012),
        ("two_flipped", 0.08, 0.04, 0.012),
    ]
    for name, fc, h, he in jobs:
        centers = geometries[name]
        rho, _, _ = rho_tuple(f, fp, gamma, fc)
        print(f"mesh {name} fc={fc} h={h}", flush=True)
        pts, tris = mesh_circles(centers, [(radius,)] * len(centers), h, he)
        print(" nodes", len(pts), flush=True)
        ratio = fem_ratio(pts, tris, k0, centers, radius, rho, src, xy)
        b, _ = solve_gyro(centers, k0, radius, rho[0], rho[1], src, 10)
        analytic = cluster_field(xy, centers, k0, radius, 1.0, src, b, 10)
        inc = hankel1(0, k0 * np.linalg.norm(xy - src, axis=1))
        ar = analytic / inc
        per = {}
        for i, probe in enumerate(names):
            per[probe] = {
                "db": float(20 * np.log10(abs(ratio[i]) / abs(ar[i]))),
                "phase_deg": float(np.angle(ratio[i] / ar[i]) * 180 / np.pi),
            }
        rows.append({"geometry": name, "fc": fc, "h": h, "dofs": int(len(pts)), "probes": per})
        print("fem", name, fc, h, per["forward"], per["side_p60"], flush=True)
    payload = {"analytic": checks, "single_side": side, "rows": rows}
    (OUT / "gyrotropic_cluster.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
