#!/usr/bin/env python3
"""Oblique Voigt interface, passivity, and a two-port Onsager scattering test."""
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
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from gyrotropic_fem import rho_tuple  # noqa: E402
from gyrotropic_tensor import dissipation_matrix, matrix_of, tensor_ordinary  # noqa: E402
from planar_fem import dirichlet, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_gyrotropic_validation"


def _root_outgoing(q2: complex, flux_of) -> complex:
    q = np.sqrt(q2 + 0j)
    cands = [q, -q]
    return max(cands, key=lambda z: (flux_of(z), float(np.imag(z))))


def oblique_coeffs(k0, ky, rho):
    rxx, rxy, ryx, ryy = rho
    eps_eff = 1.0 / rxx
    kx2 = k0 * k0 - ky * ky
    kx = np.sqrt(kx2 + 0j)
    if np.real(kx) < 0:
        kx = -kx
    q2 = (k0 * k0) * eps_eff - ky * ky

    def flux(q):
        # Ey/Hz = -(1/k0) (ρ_yx ky - ρ_yy q), Sx sign is Re(Ey/Hz).
        ey = -(1.0 / k0) * (ryx * ky - ryy * q)
        return float(np.real(ey))

    q = _root_outgoing(q2, flux)
    Z = rxx * q + rxy * ky
    R = (kx - Z) / (kx + Z)
    Tcoef = 1.0 + R
    return {"R": complex(R), "T": complex(Tcoef), "kx": complex(kx), "q": complex(q), "Z": complex(Z)}


def field_oblique(x, y, k0, ky, coef):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    phase = np.exp(1j * ky * y)
    left = np.exp(1j * coef["kx"] * x) + coef["R"] * np.exp(-1j * coef["kx"] * x)
    right = coef["T"] * np.exp(1j * coef["q"] * x)
    return np.where(x <= 0.0, left, right) * phase


def run_oblique():
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 10.0, 8.0, 1.0
    f, fp, gamma = 0.45, 0.15, 1.0e-4
    k0 = 2 * np.pi * f
    sc.fs_a = f
    rows = []
    for ky_frac in (0.0, 0.25, -0.25, 0.45):
        ky = ky_frac * k0
        for fc in (0.0, 0.05, -0.05):
            rho, _, _ = rho_tuple(f, fp, gamma, fc)
            coef = oblique_coeffs(k0, ky, rho)
            x0, y0 = 3.0, 2.5
            h = 0.02
            length, width = 2.4, 1.0
            pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
            local_x = pts[:, 0] - (x0 + length / 2)
            local_y = pts[:, 1] - y0
            cents = pts[tris].mean(1)
            cx = cents[:, 0] - (x0 + length / 2)
            inside = cx >= 0.0
            T = len(tris)
            rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
            for k in range(4):
                rho_e[k][inside] = rho[k]
            exact = field_oblique(local_x, local_y, k0, ky, coef)
            bound = (
                (np.abs(pts[:, 0] - x0) < 1e-10)
                | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
                | (np.abs(pts[:, 1] - y0) < 1e-10)
                | (np.abs(pts[:, 1] - (y0 + width)) < 1e-10)
            )
            A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
            b = np.zeros(len(pts), np.complex128)
            A, b = dirichlet(A, b, bound, exact)
            uh = splu(A.tocsc()).solve(b)
            rel = float(np.linalg.norm(uh - exact) / np.linalg.norm(exact))
            rows.append({
                "ky_over_k0": ky_frac,
                "fc": fc,
                "h": h,
                "dofs": int(len(pts)),
                "field_rel": rel,
                "R_abs": abs(coef["R"]),
                "R_phase_deg": float(np.angle(coef["R"]) * 180 / np.pi),
            })
            print("oblique", rows[-1], flush=True)
    # Sign identities on the analytic R.
    ident = []
    for ky_frac in (0.25, 0.45):
        ky = ky_frac * k0
        rp = oblique_coeffs(k0, ky, rho_tuple(f, fp, gamma, 0.05)[0])["R"]
        rm = oblique_coeffs(k0, -ky, rho_tuple(f, fp, gamma, 0.05)[0])["R"]
        r_neg_ky_neg_b = oblique_coeffs(k0, -ky, rho_tuple(f, fp, gamma, -0.05)[0])["R"]
        r0 = oblique_coeffs(k0, ky, rho_tuple(f, fp, gamma, 0.0)[0])["R"]
        r0m = oblique_coeffs(k0, -ky, rho_tuple(f, fp, gamma, 0.0)[0])["R"]
        ident.append({
            "ky_over_k0": ky_frac,
            "R_ky_plusB_vs_minusky_minusB": abs(rp - r_neg_ky_neg_b),
            "R_ky_plusB_vs_minusky_plusB": abs(rp - rm),
            "B0_even_in_ky": abs(r0 - r0m),
        })
    payload = {"rows": rows, "analytic_identities": ident}
    (OUT / "gyrotropic_oblique.json").write_text(json.dumps(payload, indent=2) + "\n")


def run_passivity():
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    f, fp = 0.45, 0.15
    k0 = 2 * np.pi * f
    sc.fs_a = f
    rows = []
    for gamma in (1.0e-4, 0.01, 1.0e-8):
        for fc in (0.05, -0.05, 0.0):
            rho, _, eperp = rho_tuple(f, fp, gamma, fc)
            exx, exy, eyx, eyy, _, _ = tensor_ordinary(f, fp, gamma, fc)
            diss = dissipation_matrix(matrix_of(exx, exy, eyx, eyy))
            eig = np.linalg.eigvalsh(0.5 * (diss + diss.conj().T))
            h = 0.02
            x0, y0 = 2.2, 2.2
            pts, tris, _, _ = rect_mesh(x0, x0 + 2.4, y0, y0 + 2.0, h)
            cents = pts[tris].mean(1)
            center = np.array([x0 + 1.2, y0 + 1.0])
            disk = np.sum((cents - center) ** 2, axis=1) <= 0.23 ** 2
            T = len(tris)
            rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
            for k in range(4):
                rho_e[k][disk] = rho[k]
            # Source in vacuum, left of the disk.
            src = np.array([x0 + 0.35, y0 + 1.0])
            sampler = ElementSampler(pts, tris)
            t_src = int(sampler.locate(src[None, :])[0])
            w = sampler._bary(t_src, src)
            b = np.zeros(len(pts), np.complex128)
            for k in range(3):
                b[int(tris[t_src, k])] += complex(w[k])
            A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
            uh = splu(A.tocsc()).solve(b)
            # Volume absorption on disk elements.
            areas = 0.5 * np.abs(
                (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1])
                - (pts[tris[:, 2], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
            )
            twice = (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1]) - (
                pts[tris[:, 2], 0] - pts[tris[:, 0], 0]
            ) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
            bx = np.stack(
                [pts[tris[:, 1], 1] - pts[tris[:, 2], 1], pts[tris[:, 2], 1] - pts[tris[:, 0], 1], pts[tris[:, 0], 1] - pts[tris[:, 1], 1]],
                axis=1,
            ) / twice[:, None]
            by = np.stack(
                [pts[tris[:, 2], 0] - pts[tris[:, 1], 0], pts[tris[:, 0], 0] - pts[tris[:, 2], 0], pts[tris[:, 1], 0] - pts[tris[:, 0], 0]],
                axis=1,
            ) / twice[:, None]
            gx = np.sum(bx * uh[tris], axis=1)
            gy = np.sum(by * uh[tris], axis=1)
            rxx, rxy, ryx, ryy = rho
            ex = (1j / k0) * (rxx * gy + rxy * (-gx))
            ey = (1j / k0) * (ryx * gy + ryy * (-gx))
            # Only inside the disk, where this rho applies.
            evec = np.stack([ex, ey], axis=1)
            absorb = np.real(np.einsum("ti,ij,tj->t", evec.conj(), diss, evec))
            p_vol = float(np.sum(areas[disk] * (k0 / 2.0) * absorb[disk]))
            # Contour r=0.45 around the disk, in vacuum.
            nphi = 720
            phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
            rad = 0.45
            xy = center + rad * np.column_stack([np.cos(phi), np.sin(phi)])
            hz, ex_c, ey_c = sampler.fields(uh, rho_e[0], rho_e[1], rho_e[2], rho_e[3], k0, xy)
            sx = 0.5 * np.real(ey_c * np.conj(hz))
            sy = -0.5 * np.real(ex_c * np.conj(hz))
            nr = np.column_stack([np.cos(phi), np.sin(phi)])
            dphi = 2 * np.pi / nphi
            p_out = float(np.sum((sx * nr[:, 0] + sy * nr[:, 1]) * rad * dphi))
            rows.append({
                "gamma": gamma,
                "fc": fc,
                "diss_min_eig": float(np.min(np.real(eig))),
                "volume_absorption": p_vol,
                "contour_outward_flux": p_out,
                "flux_balance": p_vol + p_out,
                "balance_rel": abs(p_vol + p_out) / (abs(p_vol) + 1e-30),
            })
            print("passivity", rows[-1], flush=True)
    (OUT / "gyrotropic_passivity.json").write_text(json.dumps(rows, indent=2) + "\n")


def run_guide():
    """Offset magnetized rod in a parallel-plate channel. Ends sit in PML."""
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 14.0, 10.0, 1.2
    f, fp, gamma = 0.45, 0.15, 1.0e-3
    k0 = 2 * np.pi * f
    sc.fs_a = f
    y0, y1 = 4.0, 6.0
    rows = []
    for h in (0.03, 0.015):
        pts, tris, xs, ys = rect_mesh(0.0, 14.0, y0, y1, h)
        cents = pts[tris].mean(1)
        center = np.array([7.0, 4.55])
        disk = np.sum((cents - center) ** 2, axis=1) <= 0.23 ** 2
        for fc in (0.0, 0.08, -0.08):
            rho, _, _ = rho_tuple(f, fp, gamma, fc)
            T = len(tris)
            rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
            for k in range(4):
                rho_e[k][disk] = rho[k]
            A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
            lu = splu(A.tocsc())
            sampler = ElementSampler(pts, tris)

            def sheet(x_src):
                b = np.zeros(len(pts), np.complex128)
                # A vertical line of nodal loads, normalized by the transverse spacing.
                col = np.where(np.abs(pts[:, 0] - x_src) < 0.51 * h)[0]
                # Keep one x-column: the nodes whose x is closest to x_src.
                xcol = xs[int(np.argmin(np.abs(xs - x_src)))]
                col = np.where(np.abs(pts[:, 0] - xcol) < 1e-10)[0]
                b[col] = 1.0
                return lu.solve(b)

            def flux_at(uh, x_mon):
                yy = np.linspace(y0 + 0.05, y1 - 0.05, 81)
                xy = np.column_stack([np.full(len(yy), x_mon), yy])
                hz, _, ey = sampler.fields(uh, rho_e[0], rho_e[1], rho_e[2], rho_e[3], k0, xy)
                sx = 0.5 * np.real(ey * np.conj(hz))
                return float(np.trapezoid(sx, yy)), complex(np.trapezoid(hz, yy))

            left = sheet(2.4)
            right = sheet(11.6)
            # Monitors on the far side of each source, in the vacuum channel.
            p_right, h_right = flux_at(left, 10.2)
            p_left, h_left = flux_at(right, 3.8)
            # Complex modal amplitudes.
            rows.append({
                "h": h,
                "fc": fc,
                "dofs": int(len(pts)),
                "S21_hz": [h_right.real, h_right.imag],
                "S12_hz": [h_left.real, h_left.imag],
                "P21": p_right,
                "P12": p_left,
            })
            print("guide", rows[-1], flush=True)
    # Pair +B and -B.
    def find(h, fc):
        return next(r for r in rows if r["h"] == h and r["fc"] == fc)

    summary = []
    for h in (0.03, 0.015):
        zp = np.array(find(h, 0.08)["S21_hz"], float).view(np.complex128)[0] if False else complex(*find(h, 0.08)["S21_hz"])
        zm = complex(*find(h, -0.08)["S12_hz"])
        zp12 = complex(*find(h, 0.08)["S12_hz"])
        z0_21 = complex(*find(h, 0.0)["S21_hz"])
        z0_12 = complex(*find(h, 0.0)["S12_hz"])
        summary.append({
            "h": h,
            "onsager_rel": abs(zp - zm) / abs(zp),
            "same_B_nonrecip_rel": abs(zp - zp12) / abs(zp),
            "B0_recip_rel": abs(z0_21 - z0_12) / abs(z0_21),
        })
        print("guide_summary", summary[-1], flush=True)
    (OUT / "gyrotropic_guide.json").write_text(json.dumps({"rows": rows, "summary": summary}, indent=2) + "\n")


def run_passivity_fine():
    """One magnetized case on a finer mesh. The coarse-contour residual is absolute."""
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    f, fp, gamma, fc = 0.45, 0.15, 1.0e-4, 0.05
    k0 = 2 * np.pi * f
    sc.fs_a = f
    rho, _, _ = rho_tuple(f, fp, gamma, fc)
    exx, exy, eyx, eyy, _, _ = tensor_ordinary(f, fp, gamma, fc)
    diss = dissipation_matrix(matrix_of(exx, exy, eyx, eyy))
    rows = []
    for h in (0.02, 0.01):
        x0, y0 = 2.2, 2.2
        pts, tris, _, _ = rect_mesh(x0, x0 + 2.4, y0, y0 + 2.0, h)
        cents = pts[tris].mean(1)
        center = np.array([x0 + 1.2, y0 + 1.0])
        disk = np.sum((cents - center) ** 2, axis=1) <= 0.23 ** 2
        T = len(tris)
        rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
        for k in range(4):
            rho_e[k][disk] = rho[k]
        src = np.array([x0 + 0.35, y0 + 1.0])
        sampler = ElementSampler(pts, tris)
        t_src = int(sampler.locate(src[None, :])[0])
        w = sampler._bary(t_src, src)
        b = np.zeros(len(pts), np.complex128)
        for k in range(3):
            b[int(tris[t_src, k])] += complex(w[k])
        A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
        uh = splu(A.tocsc()).solve(b)
        areas = 0.5 * np.abs(
            (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1])
            - (pts[tris[:, 2], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
        )
        twice = (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1]) - (
            pts[tris[:, 2], 0] - pts[tris[:, 0], 0]
        ) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
        bx = np.stack(
            [pts[tris[:, 1], 1] - pts[tris[:, 2], 1], pts[tris[:, 2], 1] - pts[tris[:, 0], 1], pts[tris[:, 0], 1] - pts[tris[:, 1], 1]],
            axis=1,
        ) / twice[:, None]
        by = np.stack(
            [pts[tris[:, 2], 0] - pts[tris[:, 1], 0], pts[tris[:, 0], 0] - pts[tris[:, 2], 0], pts[tris[:, 1], 0] - pts[tris[:, 0], 0]],
            axis=1,
        ) / twice[:, None]
        gx = np.sum(bx * uh[tris], axis=1)
        gy = np.sum(by * uh[tris], axis=1)
        ex = (1j / k0) * (rho[0] * gy + rho[1] * (-gx))
        ey = (1j / k0) * (rho[2] * gy + rho[3] * (-gx))
        evec = np.stack([ex, ey], axis=1)
        absorb = np.real(np.einsum("ti,ij,tj->t", evec.conj(), diss, evec))
        p_vol = float(np.sum(areas[disk] * (k0 / 2.0) * absorb[disk]))
        nphi = 1440
        phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
        rad = 0.45
        xy = center + rad * np.column_stack([np.cos(phi), np.sin(phi)])
        hz, ex_c, ey_c = sampler.fields(uh, rho_e[0], rho_e[1], rho_e[2], rho_e[3], k0, xy)
        sx = 0.5 * np.real(ey_c * np.conj(hz))
        sy = -0.5 * np.real(ex_c * np.conj(hz))
        dphi = 2 * np.pi / nphi
        p_out = float(np.sum((sx * np.cos(phi) + sy * np.sin(phi)) * rad * dphi))
        rows.append({
            "h": h,
            "dofs": int(len(pts)),
            "volume_absorption": p_vol,
            "contour_outward_flux": p_out,
            "balance": p_vol + p_out,
        })
        print("fine", rows[-1], flush=True)
    (OUT / "gyrotropic_passivity_fine.json").write_text(json.dumps(rows, indent=2) + "\n")


def run_ports():
    """Nodal ports around one offset magnetized disk. Unit nodal loads, not a modal projection."""
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 6.0, 6.0, 1.0
    f, fp, gamma = 0.45, 0.15, 1.0e-3
    k0 = 2 * np.pi * f
    sc.fs_a = f
    rows = []
    for h in (0.025, 0.0125):
        x0, y0 = 2.0, 2.0
        pts, tris, _, _ = rect_mesh(x0, x0 + 1.5, y0, y0 + 1.2, h)
        cents = pts[tris].mean(1)
        center = np.array([x0 + 0.72, y0 + 0.48])
        disk = np.sum((cents - center) ** 2, axis=1) <= 0.20 ** 2
        bound = (
            (np.abs(pts[:, 0] - x0) < 1e-10)
            | (np.abs(pts[:, 0] - (x0 + 1.5)) < 1e-10)
            | (np.abs(pts[:, 1] - y0) < 1e-10)
            | (np.abs(pts[:, 1] - (y0 + 1.2)) < 1e-10)
        )
        targets = {
            "a": np.array([x0 + 0.22, y0 + 0.85]),
            "b": np.array([x0 + 1.22, y0 + 0.90]),
            "c": np.array([x0 + 1.15, y0 + 0.22]),
        }
        nodes = {name: int(np.argmin(np.sum((pts - p) ** 2, axis=1))) for name, p in targets.items()}
        assert len(set(nodes.values())) == 3
        assert not bound[list(nodes.values())].any()

        def assemble(fc):
            rho, _, _ = rho_tuple(f, fp, gamma, fc)
            T = len(tris)
            rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
            for k in range(4):
                rho_e[k][disk] = rho[k]
            A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
            A, _ = dirichlet(A, np.zeros(len(pts), np.complex128), bound, np.zeros(len(pts)))
            return splu(A.tocsc())

        for fc in (0.0, 0.08, -0.08):
            lu = assemble(fc)
            resp = {}
            for src, i in nodes.items():
                b = np.zeros(len(pts), np.complex128)
                b[i] = 1.0
                uh = lu.solve(b)
                uh2 = lu.solve(b)
                resp[src] = {dst: complex(uh[j]) for dst, j in nodes.items()}
                resp[src]["repeat_rel"] = abs(uh[nodes["b"]] - uh2[nodes["b"]]) / (abs(uh[nodes["b"]]) + 1e-30)
            rows.append({"h": h, "fc": fc, "dofs": int(len(pts)), "G": {
                f"{s}{d}": [resp[s][d].real, resp[s][d].imag] for s in nodes for d in nodes if s != d
            }, "repeat_rel": resp["a"]["repeat_rel"]})
            print("ports", h, fc, rows[-1]["repeat_rel"], flush=True)
    def g(h, fc, key):
        rec = next(r for r in rows if r["h"] == h and r["fc"] == fc)
        return complex(*rec["G"][key])

    summary = []
    for h in (0.025, 0.0125):
        pairs = [("ab", "ba"), ("ac", "ca"), ("bc", "cb")]
        item = {"h": h}
        for p, q in pairs:
            gp = g(h, 0.08, p)
            gq_m = g(h, -0.08, q)
            gq = g(h, 0.08, q)
            g0p = g(h, 0.0, p)
            g0q = g(h, 0.0, q)
            item[f"onsager_{p}"] = abs(gp - gq_m) / abs(gp)
            item[f"sameB_{p}"] = abs(gp - gq) / abs(gp)
            item[f"B0_{p}"] = abs(g0p - g0q) / abs(g0p)
        summary.append(item)
        print("port_summary", item, flush=True)
    (OUT / "gyrotropic_ports.json").write_text(json.dumps({"rows": rows, "summary": summary}, indent=2) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1]
    {"oblique": run_oblique, "passivity": run_passivity, "guide": run_guide, "fine": run_passivity_fine, "ports": run_ports}[cmd]()


if __name__ == "__main__":
    main()
