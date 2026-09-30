#!/usr/bin/env python3
"""Convergence of Hz, its gradient, and reconstructed E/H.

P1 elements give an O(h^2) field and an O(h) gradient in L2. A one-sided
nodal difference is the same order as that gradient. Element-constant
gradients and the area-weighted nodal recovery are reported separately.
Thresholds are not changed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_maxwell import kz_of  # noqa: E402
from fem_canonical import plasma  # noqa: E402
from fem_validated_solver import assemble_anisotropic, nodal_gradient  # noqa: E402
from planar_fem import assign_rho, dirichlet, l2_error, rect_mesh, solve_dirichlet  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def element_grad(pts, tris, uh):
    xy = pts[tris]
    twice = (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1]) - (xy[:, 2, 0] - xy[:, 0, 0]) * (xy[:, 1, 1] - xy[:, 0, 1])
    bx = np.stack([xy[:, 1, 1] - xy[:, 2, 1], xy[:, 2, 1] - xy[:, 0, 1], xy[:, 0, 1] - xy[:, 1, 1]], axis=1) / twice[:, None]
    by = np.stack([xy[:, 2, 0] - xy[:, 1, 0], xy[:, 0, 0] - xy[:, 2, 0], xy[:, 1, 0] - xy[:, 0, 0]], axis=1) / twice[:, None]
    dux = np.sum(bx * uh[tris], axis=1)
    duy = np.sum(by * uh[tris], axis=1)
    area = 0.5 * np.abs(twice)
    cents = xy.mean(axis=1)
    return dux, duy, area, cents


def rel_l2(num, exact, weight):
    num = np.asarray(num)
    exact = np.asarray(exact)
    if num.ndim == 1:
        err = np.sum(np.abs(num - exact) ** 2 * weight)
        den = np.sum(np.abs(exact) ** 2 * weight)
    else:
        err = np.sum(np.sum(np.abs(num - exact) ** 2, axis=1) * weight)
        den = np.sum(np.sum(np.abs(exact) ** 2, axis=1) * weight)
    return float(np.sqrt(err / (den + 1e-30)))


def orders(rows, key):
    out = []
    for a, b in zip(rows, rows[1:]):
        if a[key] is None or b[key] is None or a[key] <= 0 or b[key] <= 0:
            out.append(None)
            continue
        out.append(float(np.log(a[key] / b[key]) / np.log(a["h"] / b["h"])))
    return out


def homogeneous():
    saved = float(sc.fs_a)
    k0 = 2 * np.pi * saved
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    rows = []
    cases = {"eps1": 1.0 + 0j, "plasma_fs": plasma(saved)}
    for name, eps in cases.items():
        kx = kz_of(k0, eps)
        series = []
        for h in (0.05, 0.025, 0.0125, 0.00625, 0.003125):
            pts, tris, _, _ = rect_mesh(1.5, 2.5, 1.5, 2.2, h)

            def ufn(x, y, kx=kx):
                return np.exp(1j * kx * (x - 1.5))

            bound = (
                (np.abs(pts[:, 0] - 1.5) < 1e-10) | (np.abs(pts[:, 0] - 2.5) < 1e-10)
                | (np.abs(pts[:, 1] - 1.5) < 1e-10) | (np.abs(pts[:, 1] - 2.2) < 1e-10)
            )
            rho = assign_rho(pts[tris].mean(1), lambda x, y, eps=eps: eps)
            uh = solve_dirichlet(pts, tris, rho, k0, ufn(pts[:, 0], pts[:, 1]), bound)
            _, rel = l2_error(pts, tris, uh, ufn)
            dux, duy, area, cents = element_grad(pts, tris, uh)
            u_c = ufn(cents[:, 0], cents[:, 1])
            g_exact = np.column_stack([1j * kx * u_c, np.zeros_like(u_c)])
            g_num = np.column_stack([dux, duy])
            grad_l2 = rel_l2(g_num, g_exact, area)
            ey = (1j / k0) * (1.0 / eps) * (-dux)
            exact_ratio = kx / (k0 * eps)
            ratio = ey / u_c
            # Interior elements, away from the Dirichlet boundary.
            interior = (cents[:, 0] > 1.5 + 2 * h) & (cents[:, 0] < 2.5 - 2 * h) & (cents[:, 1] > 1.5 + 2 * h) & (cents[:, 1] < 2.2 - 2 * h)
            point_rel = np.abs(ratio[interior] - exact_ratio) / np.abs(exact_ratio)
            gx, gy = nodal_gradient(pts, tris, uh)
            u_n = ufn(pts[:, 0], pts[:, 1])
            nodal_rel = rel_l2(np.column_stack([gx, gy]), np.column_stack([1j * kx * u_n, np.zeros_like(u_n)]), np.ones(len(pts)))
            # One-sided difference at an interior grid node, the previous sample.
            xq, yq = 2.0, 1.85
            i0 = np.where((np.abs(pts[:, 0] - xq) < 1e-8) & (np.abs(pts[:, 1] - yq) < 0.51 * h))[0]
            i1 = np.where((np.abs(pts[:, 0] - (xq + h)) < 1e-8) & (np.abs(pts[:, 1] - yq) < 0.51 * h))[0]
            onesided = None
            if len(i0) and len(i1) and abs(uh[i0[0]]) > 1e-8:
                dhz = (uh[i1[0]] - uh[i0[0]]) / h
                ey1 = (1j / k0) * (1.0 / eps) * (-dhz)
                onesided = float(abs(ey1 / uh[i0[0]] - exact_ratio) / abs(exact_ratio))
            # Element power on the line of centroids nearest x=2.
            on_cut = np.abs(cents[:, 0] - 2.0) < 0.51 * h
            sx = 0.5 * np.real(ey * np.conj(u_c))
            # Each cut element represents a vertical slice of height ~ h * (area / (0.5*h*h)) roughly.
            # Use the analytic strip integral of the element Sx values whose centroids fall in the cut.
            sx_fem = float(np.sum(sx[on_cut]) * h)
            # Analytic Sx = 0.5 Re(Ey conj Hz) = 0.5 Re(exact_ratio) |Hz|^2, |Hz|=1.
            height = 0.7
            sx_exact = float(0.5 * np.real(exact_ratio) * height)
            rec = {
                "case": name,
                "h": h,
                "dofs": int(len(pts)),
                "hz_l2": rel,
                "grad_l2": grad_l2,
                "nodal_grad_l2": nodal_rel,
                "element_EH_median": float(np.median(point_rel)),
                "element_EH_max": float(np.max(point_rel)),
                "onesided_EH": onesided,
                "power_rel": abs(sx_fem - sx_exact) / (abs(sx_exact) + 1e-30),
            }
            series.append(rec)
            print("eh homog", rec, flush=True)
        for key in ("hz_l2", "grad_l2", "nodal_grad_l2", "element_EH_median", "onesided_EH"):
            od = orders(series, key)
            for rec, order in zip(series[1:], od):
                rec[key + "_order"] = order
        rows.extend(series)
    return rows


def guide():
    rows = []
    width, k0 = 1.0, 5.0
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    for h in (0.02, 0.01, 0.005, 0.0025):
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
        x0, y0, length = 2.2, 2.2, 1.6

        def ufn(x, y, ky=ky, beta=beta, x0=x0, y0=y0):
            return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

        pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
        bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
        rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
        values = ufn(pts[:, 0], pts[:, 1])
        A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
        b = np.zeros(len(pts), dtype=np.complex128)
        A, b = dirichlet(A, b, bound, values)
        uh = splu(A.tocsc()).solve(b)
        _, rel = l2_error(pts, tris, uh, ufn)
        dux, duy, area, cents = element_grad(pts, tris, uh)
        u_c = ufn(cents[:, 0], cents[:, 1])
        g_exact_x = 1j * beta * u_c
        g_exact_y = -ky * np.sin(ky * (cents[:, 1] - y0)) * np.exp(1j * beta * (cents[:, 0] - x0))
        grad_l2 = rel_l2(np.column_stack([dux, duy]), np.column_stack([g_exact_x, g_exact_y]), area)
        ey = (1j / k0) * (-dux)
        exact_ratio = beta / k0
        # Avoid the cosine node.
        band = (np.abs(cents[:, 1] - (y0 + 0.25 * width)) < 1.5 * h) & (cents[:, 0] > x0 + 3 * h) & (cents[:, 0] < x0 + length - 3 * h)
        ratio = ey[band] / u_c[band]
        point_rel = np.abs(ratio - exact_ratio) / np.abs(exact_ratio)
        yq = y0 + 0.25 * width
        xq = x0 + 0.8
        i0 = np.where((np.abs(pts[:, 0] - xq) < 1e-8) & (np.abs(pts[:, 1] - yq) < 0.51 * h))[0]
        i1 = np.where((np.abs(pts[:, 0] - (xq + h)) < 1e-8) & (np.abs(pts[:, 1] - yq) < 0.51 * h))[0]
        onesided = None
        if len(i0) and len(i1) and abs(uh[i0[0]]) > 1e-8:
            dhz = (uh[i1[0]] - uh[i0[0]]) / h
            ey1 = (1j / k0) * (-dhz)
            onesided = float(abs(ey1 / uh[i0[0]] - exact_ratio) / abs(exact_ratio))
        gx, gy = nodal_gradient(pts, tris, uh)
        u_n = ufn(pts[:, 0], pts[:, 1])
        gnx = 1j * beta * u_n
        gny = -ky * np.sin(ky * (pts[:, 1] - y0)) * np.exp(1j * beta * (pts[:, 0] - x0))
        nodal_rel = rel_l2(np.column_stack([gx, gy]), np.column_stack([gnx, gny]), np.ones(len(pts)))
        on_cut = np.abs(cents[:, 0] - (x0 + 0.8)) < 0.51 * h
        sx = 0.5 * np.real(ey * np.conj(u_c))
        # Analytic guide power: 0.5 Re(beta) * (width/2) / k0, for |cos| integrated.
        p_exact = 0.5 * np.real(beta) * (width / 2.0) / k0
        p_fem = float(np.sum(sx[on_cut]) * h)
        rec = {
            "case": "guide_w1_k0_5",
            "h": h,
            "dofs": int(len(pts)),
            "hz_l2": rel,
            "grad_l2": grad_l2,
            "nodal_grad_l2": nodal_rel,
            "element_EH_median": float(np.median(point_rel)) if len(point_rel) else None,
            "element_EH_max": float(np.max(point_rel)) if len(point_rel) else None,
            "onesided_EH": onesided,
            "power_rel": abs(p_fem - p_exact) / abs(p_exact),
            "beta": [float(beta.real), float(beta.imag)],
        }
        rows.append(rec)
        print("eh guide", rec, flush=True)
    for key in ("hz_l2", "grad_l2", "nodal_grad_l2", "element_EH_median", "onesided_EH"):
        od = orders(rows, key)
        for rec, order in zip(rows[1:], od):
            rec[key + "_order"] = order
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = {"homogeneous": homogeneous(), "guide": guide()}
    (OUT / "eh_orders.json").write_text(json.dumps(out, indent=2) + "\n")
    print("EH_ORDERS_DONE", flush=True)


if __name__ == "__main__":
    main()
