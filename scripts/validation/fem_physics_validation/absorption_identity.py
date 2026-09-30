#!/usr/bin/env python3
"""Compare volume absorption to the closed-contour flux for a point source.

The contour does not enclose the source. For the e^{-iωt} convention the
identity on a source-free contour is

    ∮ S·n_out dl + (ω/2) ∫ Im(ε) |E|^2 dA = 0.

E = (i/(ωε)) (∂y Hz, -∂x Hz) and ω = k0. Quartz and vacuum contribute no
absorption. Thresholds are not changed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu
from scipy.special import hankel1, jv, jvp

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
import analytic_sweeps as ans  # noqa: E402
import fem_scatterers as fsc  # noqa: E402
from analytic_maxwell import _coated_one, cluster_field, mie_T, solve_clusters  # noqa: E402
from analytic_sweeps import R_CORE, R_GAP, R_SHELL, plasma  # noqa: E402
from fem_scatterers import assign, mesh_circles, to_mesh  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def incident_coeff(k0, src, center, nmax):
    d = np.asarray(src, float) - np.asarray(center, float)
    radius = np.linalg.norm(d)
    phi = np.arctan2(d[1], d[0])
    ns = np.arange(-nmax, nmax + 1)
    a = np.array([hankel1(int(n), k0 * radius) * np.exp(-1j * int(n) * phi) for n in ns])
    return ns, a


def grad_hz_polar(ns, coeff, k_in, rho, theta):
    """Cartesian gradient of sum coeff_n J_n(k rho) exp(i n theta)."""
    acc_r = np.zeros(rho.shape, np.complex128)
    acc_t = np.zeros(rho.shape, np.complex128)
    for c, n in zip(coeff, ns):
        if abs(c) < 1e-18:
            continue
        n = int(n)
        z = k_in * rho
        phase = np.exp(1j * n * theta)
        acc_r += c * k_in * jvp(n, z, 1) * phase
        # (1/rho) d/dtheta. At rho=0 this term is sampled only for rho>0.
        acc_t += c * (1j * n) * jv(n, z) * phase / np.maximum(rho, 1e-14)
    gr = acc_r * np.cos(theta) - acc_t * np.sin(theta)
    gy = acc_r * np.sin(theta) + acc_t * np.cos(theta)
    return gr, gy


def volume_from_coeff(ns, coeff, k0, eps, radius, n_rho=480, n_th=720):
    rho = np.linspace(radius / n_rho, radius, n_rho)
    theta = np.linspace(0.0, 2 * np.pi, n_th, endpoint=False)
    rr, tt = np.meshgrid(rho, theta, indexing="ij")
    k_in = k0 * np.sqrt(complex(eps))
    gx, gy = grad_hz_polar(ns, coeff, k_in, rr, tt)
    e2 = (np.abs(gx) ** 2 + np.abs(gy) ** 2) / abs(k0 * eps) ** 2
    dr = rho[1] - rho[0]
    dth = theta[1] - theta[0]
    # (ω/2) Im(ε) ∫ |E|^2 dA
    return float(0.5 * k0 * np.imag(eps) * np.sum(e2 * rr * dr * dth))


def contour_flux(field_fn, radius, n=1440):
    ang = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    xy = np.column_stack([radius * np.cos(ang), radius * np.sin(ang)])
    step = 1e-5
    hz = field_fn(xy)
    dx = (field_fn(xy + np.array([step, 0.0])) - field_fn(xy - np.array([step, 0.0]))) / (2 * step)
    dy = (field_fn(xy + np.array([0.0, step])) - field_fn(xy - np.array([0.0, step]))) / (2 * step)
    ex = (1j / 1.0) * dy  # placeholder, k0 applied by caller through field scaling
    return hz, dx, dy, xy, ang


def air_flux(field_fn, k0, radius, n=1440):
    ang = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    xy = np.column_stack([radius * np.cos(ang), radius * np.sin(ang)])
    step = 1e-5
    hz = field_fn(xy)
    dx = (field_fn(xy + np.array([step, 0.0])) - field_fn(xy - np.array([step, 0.0]))) / (2 * step)
    dy = (field_fn(xy + np.array([0.0, step])) - field_fn(xy - np.array([0.0, step]))) / (2 * step)
    ex = (1j / k0) * dy
    ey = (1j / k0) * (-dx)
    sx = 0.5 * np.real(ey * np.conj(hz))
    sy = -0.5 * np.real(ex * np.conj(hz))
    nrm = xy / radius
    dl = radius * (ang[1] - ang[0])
    return float(np.sum((sx * nrm[:, 0] + sy * nrm[:, 1]) * dl))


def bare_coeff(k0, eps, radius, src, nmax):
    ns, a = incident_coeff(k0, src, np.zeros(2), nmax)
    T = mie_T(k0, eps, radius, nmax)
    b = T * a
    kp = k0 * np.sqrt(complex(eps))
    c = np.zeros_like(a)
    for i, n in enumerate(ns):
        c[i] = (a[i] * jv(int(n), k0 * radius) + b[i] * hankel1(int(n), k0 * radius)) / jv(int(n), kp * radius)
    return ns, c, b


def coated_core_coeff(k0, radii, layers, src, nmax):
    ns, a = incident_coeff(k0, src, np.zeros(2), nmax)
    c = np.zeros_like(a)
    for i, n in enumerate(ns):
        _t, sol = _coated_one(k0, radii, layers, int(n))
        c[i] = sol[0] * a[i]
    return ns, c


def fem_volume(k0, eps, h, h_edge, lx, ly, dpml, src):
    fsc.LX, fsc.LY, fsc.DPML = lx, ly, dpml
    ans.SRC[:] = src
    fsc.SRC[:] = src
    centers = np.zeros((1, 2))
    pts, tris = mesh_circles(centers, [(R_CORE,)], h, h_edge)
    areas = {}
    rho = assign(pts, tris, centers, R_CORE, eps, False, areas)
    sampler = ElementSampler(pts, tris)
    bvec = np.zeros(len(pts), np.complex128)
    src_m = to_mesh(src)
    t_src = int(sampler.locate(src_m[None, :])[0])
    w = sampler._bary(t_src, src_m)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    uh = splu(A.tocsc()).solve(bvec)
    twice = (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1]) - (
        pts[tris[:, 2], 0] - pts[tris[:, 0], 0]
    ) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
    area = 0.5 * np.abs(twice)
    bx = np.stack(
        [
            pts[tris[:, 1], 1] - pts[tris[:, 2], 1],
            pts[tris[:, 2], 1] - pts[tris[:, 0], 1],
            pts[tris[:, 0], 1] - pts[tris[:, 1], 1],
        ],
        axis=1,
    ) / twice[:, None]
    by = np.stack(
        [
            pts[tris[:, 2], 0] - pts[tris[:, 1], 0],
            pts[tris[:, 0], 0] - pts[tris[:, 2], 0],
            pts[tris[:, 1], 0] - pts[tris[:, 0], 0],
        ],
        axis=1,
    ) / twice[:, None]
    dux = np.sum(bx * uh[tris], axis=1)
    duy = np.sum(by * uh[tris], axis=1)
    origin = pts[tris].mean(1) - np.array([lx / 2.0, ly / 2.0])
    inside = np.linalg.norm(origin, axis=1) <= R_CORE
    rho_e = 1.0 / eps
    ex = (1j / k0) * rho_e * duy
    ey = (1j / k0) * rho_e * (-dux)
    e2 = np.abs(ex) ** 2 + np.abs(ey) ** 2
    p_abs = float(0.5 * k0 * np.imag(eps) * np.sum(e2[inside] * area[inside]))
    return {
        "h": h,
        "dofs": int(len(pts)),
        "volume_absorption": p_abs,
        "plasma_area": areas["plasma"],
        "plasma_exact": areas["plasma_exact"],
        "area_rel": float(areas["plasma"] / areas["plasma_exact"] - 1.0),
    }


def main():
    k0 = 2 * np.pi * float(sc.fs_a)
    eps = plasma(float(sc.fs_a))
    src = np.array([-4.5, 0.0])
    nmax = 16
    ns, coeff, _b = bare_coeff(k0, eps, R_CORE, src, nmax)
    vol = volume_from_coeff(ns, coeff, k0, eps, R_CORE)

    def bare_field(xy):
        b = solve_clusters(np.zeros((1, 2)), k0, R_CORE, eps, src, nmax)
        return cluster_field(xy, np.zeros((1, 2)), k0, R_CORE, eps, src, b, nmax)

    flux = air_flux(bare_field, k0, 1.2)
    # Lossless dielectric must give a null contour.
    def die_field(xy):
        b = solve_clusters(np.zeros((1, 2)), k0, R_CORE, 3.8 + 0j, src, nmax)
        return cluster_field(xy, np.zeros((1, 2)), k0, R_CORE, 3.8 + 0j, src, b, nmax)

    flux_die = air_flux(die_field, k0, 1.2)
    radii = [R_CORE, R_GAP, R_SHELL]
    layers = [eps, 1.0 + 0j, 3.8 + 0j]
    ns_c, coeff_c = coated_core_coeff(k0, radii, layers, src, nmax)
    vol_c = volume_from_coeff(ns_c, coeff_c, k0, eps, R_CORE)

    def coated_field(xy):
        b = solve_clusters(np.zeros((1, 2)), k0, R_CORE, eps, src, nmax, coated=(radii, layers))
        return cluster_field(xy, np.zeros((1, 2)), k0, R_SHELL, eps, src, b, nmax)

    flux_c = air_flux(coated_field, k0, 1.2)
    # Two coated cylinders. Absorption is minus the analytic exterior flux.
    centers = np.array([[-0.5, 0.0], [0.5, 0.0]])

    def pair_field(xy):
        b = solve_clusters(centers, k0, R_CORE, eps, src, 12, coated=(radii, layers))
        return cluster_field(xy, centers, k0, R_SHELL, eps, src, b, 12)

    flux_pair = air_flux(pair_field, k0, 1.8)
    rec = {
        "k0": k0,
        "eps": [eps.real, eps.imag],
        "bare_volume": vol,
        "bare_outward_flux": flux,
        "bare_balance": flux + vol,
        "bare_rel_balance": (flux + vol) / vol,
        "dielectric_outward_flux": flux_die,
        "coated_volume": vol_c,
        "coated_outward_flux": flux_c,
        "coated_balance": flux_c + vol_c,
        "coated_rel_balance": (flux_c + vol_c) / vol_c,
        "coated_pair_outward_flux": flux_pair,
        "coated_pair_absorption_from_flux": -flux_pair,
        "fem": [],
    }
    print("analytic", {k: rec[k] for k in rec if k != "fem"}, flush=True)
    saved = (fsc.LX, fsc.LY, fsc.DPML, ans.SRC.copy())
    try:
        for h, he in ((0.02, 0.006), (0.01, 0.004)):
            row = fem_volume(k0, eps, h, he, 8.0, 6.0, 0.8, np.array([-2.6, 0.0]))
            # Analytic absorption for this closer source, same cylinder.
            ns2, c2, _ = bare_coeff(k0, eps, R_CORE, np.array([-2.6, 0.0]), nmax)
            row["analytic_volume"] = volume_from_coeff(ns2, c2, k0, eps, R_CORE)
            row["rel_vs_analytic"] = row["volume_absorption"] / row["analytic_volume"] - 1.0
            rec["fem"].append(row)
            print("fem", row, flush=True)
    finally:
        fsc.LX, fsc.LY, fsc.DPML = saved[0], saved[1], saved[2]
        ans.SRC[:] = saved[3]
        fsc.SRC[:] = saved[3]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "absorption_identity.json").write_text(json.dumps(rec, indent=2) + "\n")
    print("WROTE absorption_identity.json", flush=True)


if __name__ == "__main__":
    main()
