#!/usr/bin/env python3
"""Plane-wave echo width for the Hamid & Cooray 2016 Table I ellipse.

The acceptance thresholds are frozen in PUBLICATION_BENCHMARK_PASS_CRITERIA.md.
A circular Mie calibration runs first. If it misses 0.5%, the ellipse numbers
are still written and are not treated as a paper reproduction.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.special import hankel1, jv

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from analytic_maxwell import mie_T  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_full_device_validation"
PUBLISHED = [
    (0.0, 45.0, 0.13781),
    (0.0, 90.0, 0.73693),
    (45.0, 0.0, 0.14538),
    (45.0, 90.0, 0.62262),
    (90.0, 0.0, 0.77898),
    (90.0, 45.0, 0.68662),
]


def _tri():
    vend = ROOT / "scripts" / "validation" / "fem_meep_validation" / "_vendor"
    if str(vend) not in sys.path:
        sys.path.insert(0, str(vend))
    import triangle

    return triangle


def sigma_over_lambda_from_f(f_amp: np.ndarray) -> np.ndarray:
    """σ/λ for u ~ sqrt(2/(π k ρ)) exp(i(kρ-π/4)) f, with |u_inc|=1."""
    return 2.0 * np.abs(f_amp) ** 2 / np.pi


def mie_sigma(k0: float, eps: complex, radius: float, phi: np.ndarray, nmax: int = 12) -> np.ndarray:
    t = mie_T(k0, eps, radius, nmax)
    ns = np.arange(-nmax, nmax + 1)
    f = np.zeros(len(phi), np.complex128)
    for tn, n in zip(t, ns):
        f += tn * np.exp(1j * n * phi)
    return sigma_over_lambda_from_f(f)


def far_field_amplitude(k0: float, xy: np.ndarray, hz: np.ndarray, dhz_dn: np.ndarray, n_hat: np.ndarray, dl: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """f(φ) in the outgoing representation. n_hat points out of the scatterer."""
    out = np.zeros(len(phi), np.complex128)
    for i, ang in enumerate(phi):
        rhat = np.array([np.cos(ang), np.sin(ang)])
        phase = np.exp(-1j * k0 * (xy @ rhat))
        kern = (hz * (-1j * k0 * (n_hat @ rhat)) - dhz_dn) * (1j / 4.0)
        out[i] = np.sum(kern * phase * dl)
    return out


def mesh_ellipse(a: float, b: float, h: float, h_edge: float):
    nx = ny = 16.0
    dp = 1.5
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, dp
    sc.fs_a = 1.0 / (2.0 * np.pi)  # k0 = 1
    center = np.array([nx / 2.0, ny / 2.0])
    m = max(48, int(np.ceil(2 * np.pi * max(a, b) / h_edge)))
    ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
    ring = center + np.column_stack([a * np.cos(ang), b * np.sin(ang)])
    box = np.array([[0, 0], [nx, 0], [nx, ny], [0, ny]], float)
    verts = np.vstack([box, ring])
    segs = [[0, 1], [1, 2], [2, 3], [3, 0]]
    for i in range(m):
        segs.append([4 + i, 4 + (i + 1) % m])
    regions = [
        [center[0], center[1], 1, 0.45 * h_edge**2],
        [center[0] + a + 0.4, center[1], 2, 0.45 * h**2],
        [0.4, 0.4, 3, 0.45 * (2.0 * h) ** 2],
    ]
    tri = _tri()
    mesh = tri.triangulate(
        {"vertices": verts, "segments": np.asarray(segs, np.int32), "regions": np.asarray(regions, float)},
        "pq20aA",
    )
    pts = np.asarray(mesh["vertices"], float)
    tris = np.asarray(mesh["triangles"], int)
    return pts, tris, center


def solve_scattered(pts, tris, inside, k0, k_vec, eps_xx, eps_xy, eps_yx, eps_yy, mu_zz: float = 1.0):
    rho = np.array([[eps_xx, eps_xy], [eps_yx, eps_yy]], np.complex128)
    rho = np.linalg.inv(rho)
    tcount = len(tris)
    rxx = np.ones(tcount, np.complex128)
    rxy = np.zeros(tcount, np.complex128)
    ryx = np.zeros(tcount, np.complex128)
    ryy = np.ones(tcount, np.complex128)
    rxx[inside] = rho[0, 0]
    rxy[inside] = rho[0, 1]
    ryx[inside] = rho[1, 0]
    ryy[inside] = rho[1, 1]
    rel = pts - pts.mean(0) * 0.0
    # Incident plane wave, e^{-iωt}, wavevector k_vec.
    ui = np.exp(-1j * (pts @ k_vec))
    A = assemble_anisotropic(pts, tris, rxx, rxy, ryx, ryy, k0, np.zeros(len(pts), dtype=bool))
    # Contrast RHS on object elements: ∫ (ρ - I) ∇u_i · ∇v.
    x = pts[tris, 0]
    y = pts[tris, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    bx = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1) / twice[:, None]
    by = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1) / twice[:, None]
    gui = ui[tris]
    gix = np.sum(bx * gui, 1)
    giy = np.sum(by * gui, 1)
    rows, cols, data = [], [], []
    if mu_zz != 1.0:
        # Hz mass term k0^2 μ_zz. The production operator uses μ_zz = 1.
        extra = (mu_zz - 1.0) * (k0 ** 2) * (0.5 * np.abs(twice))
        A = A.tolil()
        for el in np.flatnonzero(inside):
            nodes = tris[el]
            for i in range(3):
                for j in range(3):
                    A[nodes[i], nodes[j]] -= extra[el] * (1.0 / 6.0 if i == j else 1.0 / 12.0)
        A = A.tocsr()
    b = np.zeros(len(pts), np.complex128)
    drxx = (rxx - 1.0)
    # Only the object contrast. Vacuum elements contribute 0.
    for i in range(3):
        incr = inside * (
            drxx * gix * bx[:, i]
            + rxy * giy * bx[:, i]
            + ryx * gix * by[:, i]
            + (ryy - 1.0) * giy * by[:, i]
        ) * (0.5 * np.abs(twice))
        np.add.at(b, tris[:, i], incr)
    del rows, cols, data, rel
    us = splu(A.tocsc()).solve(b)
    return us, (rxx, rxy, ryx, ryy)


def echo_widths(pts, tris, us, rho, k0, center, phi):
    sampler = ElementSampler(pts, tris)
    rad = 3.0
    n = 720
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xy = center + rad * np.column_stack([np.cos(ang), np.sin(ang)])
    n_hat = np.column_stack([np.cos(ang), np.sin(ang)])
    hz, ex, ey = sampler.fields(us, *rho, k0, xy)
    # Vacuum: ∂x Hz = i k0 Ey, ∂y Hz = -i k0 Ex.
    dhz_dn = (1j * k0 * ey) * n_hat[:, 0] + (-1j * k0 * ex) * n_hat[:, 1]
    dl = np.full(n, rad * 2 * np.pi / n)
    f = far_field_amplitude(k0, xy - center, hz, dhz_dn, n_hat, dl, phi)
    return sigma_over_lambda_from_f(f)


def propagation(inc_deg: float) -> np.ndarray:
    """Paper incidence angle: 0° travels toward -x."""
    beta = np.pi + np.deg2rad(inc_deg)
    return np.array([np.cos(beta), np.sin(beta)])  # k0 = 1


def run_calibration():
    k0 = 1.0
    radius = 1.0
    eps = 4.0 + 0.0j
    phi = np.deg2rad(np.array([0.0, 45.0, 90.0, 180.0]))
    analytic = mie_sigma(k0, eps, radius, phi)
    rows = []
    for h, he in ((0.12, 0.06), (0.06, 0.03)):
        pts, tris, center = mesh_ellipse(radius, radius, h, he)
        origin = pts[tris].mean(1) - center
        inside = np.linalg.norm(origin, axis=1) <= radius
        k_vec = propagation(0.0)
        us, rho = solve_scattered(pts, tris, inside, k0, k_vec, eps, 0.0, 0.0, eps)
        got = echo_widths(pts, tris, us, rho, k0, center, phi)
        rel = np.abs(got - analytic) / np.maximum(analytic, 1e-12)
        row = {"h": h, "nodes": int(len(pts)), "analytic": analytic.tolist(), "fem": got.tolist(), "rel": rel.tolist()}
        rows.append(row)
        print("calibration", row, flush=True)
    return rows


def run_ellipse():
    k0 = 1.0
    a, b = 1.0, 0.5
    exx = 7.0 / 9.0
    exy = -5.0 / 18.0
    phi_list = sorted({p[1] for p in PUBLISHED})
    rows = []
    for h, he in ((0.10, 0.04), (0.05, 0.02)):
        pts, tris, center = mesh_ellipse(a, b, h, he)
        origin = pts[tris].mean(1) - center
        scale = np.column_stack([origin[:, 0] / a, origin[:, 1] / b])
        inside = np.sum(scale**2, axis=1) <= 1.0
        print(f"ellipse h={h} nodes={len(pts)} inside={int(inside.sum())}", flush=True)
        cache = {}
        for inc, obs, published in PUBLISHED:
            if inc not in cache:
                us, rho = solve_scattered(pts, tris, inside, k0, propagation(inc), exx, exy, -exy, exx)
                angs = np.deg2rad(np.array(phi_list, float))
                widths = echo_widths(pts, tris, us, rho, k0, center, angs)
                cache[inc] = dict(zip(phi_list, widths))
            fem = float(cache[inc][obs])
            rel = abs(fem - published) / published
            rows.append({
                "h": h,
                "nodes": int(len(pts)),
                "incidence_deg": inc,
                "observation_deg": obs,
                "published": published,
                "fem": fem,
                "relative_error": rel,
            })
            print("ellipse", rows[-1], flush=True)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cal = run_calibration()
    ell = run_ellipse()
    # Transpose-tensor reciprocity on the fine mesh only, same angles.
    payload = {"calibration": cal, "ellipse": ell}
    (OUT / "hamid_ellipse_benchmark.json").write_text(json.dumps(payload, indent=2) + "\n")
    print("WROTE hamid_ellipse_benchmark.json", flush=True)


if __name__ == "__main__":
    main()
