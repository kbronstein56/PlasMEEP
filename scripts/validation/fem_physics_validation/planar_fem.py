#!/usr/bin/env python3
"""Homogeneous, PEC-guide, Fresnel, and slab tests of the 2D P1 Hz solver.

Analytic fields come from analytic_maxwell.py. The FEM side calls the
validation assembler. Pass criteria are in PASS_CRITERIA.md and are not
adjusted here.
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
from analytic_maxwell import fresnel_ht, kz_of, layer_field, stack_response  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402
from fem_validated_solver import assemble_anisotropic, pml_sx_sy  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def rect_mesh(x0, x1, y0, y1, h):
    nx = int(round((x1 - x0) / h))
    ny = int(round((y1 - y0) / h))
    xs = x0 + np.linspace(0.0, x1 - x0, nx + 1)
    ys = y0 + np.linspace(0.0, y1 - y0, ny + 1)
    pts = np.array([(x, y) for y in ys for x in xs], float)
    tris = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            b = a + 1
            c = a + (nx + 1)
            d = c + 1
            tris.append((a, b, d))
            tris.append((a, d, c))
    return pts, np.asarray(tris, int), xs, ys


def dirichlet(A, b, bound, values):
    A = A.tolil()
    b = b.copy()
    for i in np.flatnonzero(bound):
        A.rows[i] = [int(i)]
        A.data[i] = [1.0 + 0j]
        b[i] = values[i]
    return A.tocsr(), b


def l2_error(pts, tris, uh, u_exact_fn):
    acc = 0.0
    norm = 0.0
    for tri in tris:
        xy = pts[tri]
        area = 0.5 * abs(
            (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
        )
        mid = xy.mean(0)
        u = u_exact_fn(mid[0], mid[1])
        uh_m = np.mean(uh[tri])
        acc += area * abs(uh_m - u) ** 2
        norm += area * abs(u) ** 2
    return float(np.sqrt(acc)), float(np.sqrt(acc / norm))


def solve_dirichlet(pts, tris, rho_e, k0, values, bound):
    A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
    b = np.zeros(len(pts), dtype=np.complex128)
    A, b = dirichlet(A, b, bound, values)
    uh = splu(A.tocsc()).solve(b)
    return uh


def assign_rho(tris_xy, eps_of_xy):
    T = len(tris_xy)
    rho = [np.ones(T, dtype=np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    for i, xy in enumerate(tris_xy):
        eps = complex(eps_of_xy(xy[0], xy[1]))
        rho[0][i] = 1.0 / eps
        rho[3][i] = 1.0 / eps
    return rho


def plasma_eps(f_a: float) -> complex:
    eps, _ = gyrotropic_drude_eps_eta(float(f_a), float(sc.fp_a), float(sc.gamma_a), 0.0)
    return complex(eps)


def check_fabry_perot():
    k0 = 1.7
    d = 0.35
    eps = 3.8 + 0j
    resp = stack_response(k0, [(d, eps)])
    k2 = kz_of(k0, eps)
    fr = fresnel_ht(k0, 1.0, eps)
    r01 = fr["r"]
    # From the dielectric side back to air. Swap.
    r10 = fresnel_ht(k0, eps, 1.0)["r"]
    t01 = fr["t"]
    t10 = fresnel_ht(k0, eps, 1.0)["t"]
    phase = np.exp(2j * k2 * d)
    r = r01 + t01 * t10 * r10 * phase / (1.0 - r10 * r10 * phase)
    # Downstream vacuum-normalized t: field at x>d relative to e^{i k0 x}.
    # At the back face the transmitted Hz is t01 * e^{i k2 d} / (1 - r10^2 e^{2 i k2 d}) times extra transmissions.
    denom = 1.0 - r10 * r10 * phase
    hz_back = t01 * np.exp(1j * k2 * d) * t10 / denom
    t_over_vac = hz_back * np.exp(-1j * k0 * d)
    return {
        "r_err": abs(resp["r"] - r),
        "t_err": abs(resp["t_over_vacuum"] - t_over_vac),
    }


def homogeneous(levels):
    """Dirichlet recovery of a uniform-medium plane wave, plus the E/H ratio."""
    saved = float(sc.fs_a)
    rows = []
    cases = {
        "eps1": 1.0 + 0j,
        "eps3.8": 3.8 + 0j,
        "eps_complex": 2.0 + 0.3j,
        "plasma_fs": plasma_eps(saved),
    }
    k0 = 2 * np.pi * saved
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    for name, eps in cases.items():
        kx = kz_of(k0, eps)

        def ufn(x, y, kx=kx):
            return np.exp(1j * kx * (x - 1.5))

        for h in levels:
            pts, tris, _, _ = rect_mesh(1.5, 2.5, 1.5, 2.2, h)
            sx, sy = pml_sx_sy(pts)
            if np.max(np.abs(sx - 1)) > 1e-8 or np.max(np.abs(sy - 1)) > 1e-8:
                raise RuntimeError("homogeneous box entered the PML")
            local_x = pts[:, 0]
            bound = (
                (np.abs(pts[:, 0] - 1.5) < 1e-10)
                | (np.abs(pts[:, 0] - 2.5) < 1e-10)
                | (np.abs(pts[:, 1] - 1.5) < 1e-10)
                | (np.abs(pts[:, 1] - 2.2) < 1e-10)
            )
            cents = pts[tris].mean(1)
            rho = assign_rho(cents, lambda x, y, eps=eps: eps)
            uh = solve_dirichlet(pts, tris, rho, k0, ufn(pts[:, 0], pts[:, 1]), bound)
            l2, rel = l2_error(pts, tris, uh, ufn)
            # Ey/Hz at an interior element. Ey = kx/(k0 ε) Hz for this wave, since ω=k0 in these units.
            tri = tris[len(tris) // 2]
            xy = pts[tri]
            twice = (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
            by = np.array([xy[2, 0] - xy[1, 0], xy[0, 0] - xy[2, 0], xy[1, 0] - xy[0, 0]]) / twice
            bx = np.array([xy[1, 1] - xy[2, 1], xy[2, 1] - xy[0, 1], xy[0, 1] - xy[1, 1]]) / twice
            dHz_dx = bx @ uh[tri]
            hz_m = np.mean(uh[tri])
            ey = (1j / k0) * (1.0 / eps) * (-dHz_dx)
            ratio = ey / hz_m
            exact = kx / (k0 * eps)
            rows.append({
                "case": name,
                "h": h,
                "L2": l2,
                "rel_L2": rel,
                "EH_rel": float(abs(ratio - exact) / abs(exact)),
                "phase_deg_center": float(np.angle(uh[len(pts) // 2] / ufn(pts[len(pts) // 2, 0], pts[len(pts) // 2, 1])) * 180 / np.pi),
            })
            print("homogeneous", rows[-1], flush=True)
    sc.fs_a = saved
    return rows


def pec_guide(levels):
    """m=1 parallel-plate mode. Plates y=const are natural Neumann, which is PEC for Hz."""
    k0 = 2.4
    width = 1.0
    length = 1.2
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 6.0, 6.0, 1.0
    x0, y0 = 2.0, 2.0

    def ufn(x, y):
        return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

    rows = []
    for h in levels:
        pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
        # Dirichlet only on the two ends. Top and bottom stay natural Neumann.
        bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
        cents = pts[tris].mean(1)
        rho = assign_rho(cents, lambda x, y: 1.0)
        uh = solve_dirichlet(pts, tris, rho, k0, ufn(pts[:, 0], pts[:, 1]), bound)
        l2, rel = l2_error(pts, tris, uh, ufn)
        rows.append({"h": h, "beta": [beta.real, beta.imag], "L2": l2, "rel_L2": rel, "max_nodal": float(np.max(np.abs(uh - ufn(pts[:, 0], pts[:, 1]))))})
        print("pec", rows[-1], flush=True)
    return rows


def fresnel_box(levels):
    k0 = 2 * np.pi * float(sc.fs_a)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    rows = []
    specs = [
        ("air_to_quartz", 1.0 + 0j, 3.8 + 0j, [0.0, 20.0, 40.0]),
        ("quartz_to_air", 3.8 + 0j, 1.0 + 0j, [0.0, 15.0, 40.0]),
        ("air_to_complex", 1.0 + 0j, 2.0 + 0.3j, [0.0, 25.0, 50.0]),
        ("air_to_plasma", 1.0 + 0j, plasma_eps(float(sc.fs_a)), [0.0, 20.0, 40.0]),
    ]
    for name, e1, e2, angles in specs:
        for ang in angles:
            ky = kz_of(k0, e1) * np.sin(np.deg2rad(ang))
            fr = fresnel_ht(k0, e1, e2, ky)
            kx1, kx2 = fr["kx1"], fr["kx2"]
            r, t = fr["r"], fr["t"]

            def ufn(x, y, ky=ky, kx1=kx1, kx2=kx2, r=r, t=t):
                # Interface at local x=0. Incident Hz = exp(i kx1 x + i ky y).
                yy = y - 1.5
                x_local = (x - 1.5) - 0.5
                phase = np.exp(1j * ky * yy)
                if x_local <= 0.0:
                    return phase * (np.exp(1j * kx1 * x_local) + r * np.exp(-1j * kx1 * x_local))
                return phase * t * np.exp(1j * kx2 * x_local)

            # The expression above has a cancelled dummy term. ufn is the interface-referenced field.
            for h in levels:
                pts, tris, _, _ = rect_mesh(1.5, 2.5, 1.5, 2.3, h)
                bound = (
                    (np.abs(pts[:, 0] - 1.5) < 1e-10)
                    | (np.abs(pts[:, 0] - 2.5) < 1e-10)
                    | (np.abs(pts[:, 1] - 1.5) < 1e-10)
                    | (np.abs(pts[:, 1] - 2.3) < 1e-10)
                )
                cents = pts[tris].mean(1)

                def eps_of(x, y, e1=e1, e2=e2):
                    return e1 if x - 1.5 <= 0.5 else e2

                rho = assign_rho(cents, eps_of)
                values = np.array([ufn(x, y) for x, y in pts])
                uh = solve_dirichlet(pts, tris, rho, k0, values, bound)
                l2, rel = l2_error(pts, tris, uh, ufn)
                rows.append({
                    "case": name,
                    "angle_deg": ang,
                    "h": h,
                    "rel_L2": rel,
                    "L2": l2,
                    "analytic_r_abs": abs(r),
                    "analytic_T_power": fr["T_power"],
                    "analytic_power_sum": fr["power_sum"],
                })
                print("fresnel", rows[-1], flush=True)
    return rows


def driven_slabs(levels):
    """Line-source strip with PML. Compare Hz/Hz_vacuum to the transfer-matrix field."""
    saved_fs = float(sc.fs_a)
    freqs = [saved_fs, saved_fs * 0.92, saved_fs * 1.08]
    stacks = {
        "quartz_0.20": [(0.20, 3.8 + 0j)],
        "quartz_0.40": [(0.40, 3.8 + 0j)],
        "plasma_0.20": None,  # filled per frequency
        "plasma_0.40": None,
        "quartz_plasma_quartz": None,
    }
    rows = []
    for f_a in freqs:
        sc.fs_a = f_a
        k0 = 2 * np.pi * f_a
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = 10.0, 4.0, 1.0
        eps_p = plasma_eps(f_a)
        built = {
            "quartz_0.20": [(0.20, 3.8 + 0j)],
            "quartz_0.40": [(0.40, 3.8 + 0j)],
            "plasma_0.20": [(0.20, eps_p)],
            "plasma_0.40": [(0.40, eps_p)],
            "q_p_q": [(0.10, 3.8 + 0j), (0.40, eps_p), (0.10, 3.8 + 0j)],
        }
        for h in levels:
            pts, tris, xs, ys = rect_mesh(0.0, 10.0, 1.6, 2.4, h)
            sx, sy = pml_sx_sy(pts[tris].mean(1))
            if np.max(np.abs(sy - 1)) > 1e-8:
                raise RuntimeError("strip entered the y PML")
            cents = pts[tris].mean(1)
            src_col = np.where(np.abs(xs - 2.2) < 1e-10)[0]
            if len(src_col) != 1:
                raise RuntimeError("source plane is not a grid line")
            b = np.zeros(len(pts), dtype=np.complex128)
            # Uniform sheet at x=2.2.
            on_src = np.abs(pts[:, 0] - 2.2) < 1e-10
            b[on_src] = 1.0
            pec = np.zeros(len(pts), dtype=bool)
            rho_vac = assign_rho(cents, lambda x, y: 1.0)
            A = assemble_anisotropic(pts, tris, *rho_vac, k0, pec)
            vac = splu(A.tocsc()).solve(b)
            for name, stack in built.items():
                x_front = 4.0

                def eps_of(x, y, stack=stack, x_front=x_front):
                    local = x - x_front
                    cursor = 0.0
                    for thickness, eps in stack:
                        if cursor <= local < cursor + thickness:
                            return eps
                        cursor += thickness
                    return 1.0 + 0j

                rho = assign_rho(cents, eps_of)
                A = assemble_anisotropic(pts, tris, *rho, k0, pec)
                hz = splu(A.tocsc()).solve(b)
                # Probes: average over the y column.
                def column(x_probe, field):
                    m = np.abs(pts[:, 0] - x_probe) < 1e-10
                    return np.mean(field[m])

                x_r = 3.2
                L = sum(t for t, _ in stack)
                x_t = 4.0 + L + 0.6
                # Snap probes onto the grid.
                x_t = xs[np.argmin(np.abs(xs - x_t))]
                ratio_t = column(x_t, hz) / column(x_t, vac)
                ratio_r = column(x_r, hz) / column(x_r, vac)
                analytic = layer_field(np.array([x_r - x_front, x_t - x_front]), k0, stack, 2.2 - x_front)
                vac_a = np.exp(1j * k0 * (np.array([x_r, x_t]) - 2.2))
                a_r = analytic[0] / vac_a[0]
                a_t = analytic[1] / vac_a[1]
                def db(z, a):
                    return float(20 * np.log10(abs(z) / abs(a)))
                def ph(z, a):
                    return float(np.angle(z / a) * 180 / np.pi)
                rec = {
                    "stack": name,
                    "f_GHz": float(f_a * 2.99792458e8 / 0.02 / 1e9),
                    "h": h,
                    "T_db": db(ratio_t, a_t),
                    "T_phase_deg": ph(ratio_t, a_t),
                    "R_db": db(ratio_r, a_r),
                    "R_phase_deg": ph(ratio_r, a_r),
                    "T_abs_err": float(abs(ratio_t - a_t)),
                    "R_abs_err": float(abs(ratio_r - a_r)),
                    "analytic_T_abs": float(abs(a_t)),
                    "analytic_R_abs": float(abs(a_r)),
                    "eps_plasma": [eps_p.real, eps_p.imag],
                }
                rows.append(rec)
                print("slab", rec, flush=True)
    sc.fs_a = saved_fs
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fabry = check_fabry_perot()
    print("fabry", fabry, flush=True)
    if fabry["r_err"] > 1e-10 or fabry["t_err"] > 1e-10:
        raise RuntimeError(f"transfer matrix disagrees with Fabry-Perot: {fabry}")
    fine = [0.05, 0.025]
    slab_h = [0.04, 0.02]
    out = {
        "fabry_perot_error": fabry,
        "homogeneous": homogeneous(fine),
        "pec_guide": pec_guide(fine),
        "fresnel_dirichlet": fresnel_box(fine),
        "driven_slabs": driven_slabs(slab_h),
    }
    (OUT / "planar_fem.json").write_text(json.dumps(out, indent=2) + "\n")
    print("PLANAR_DONE", flush=True)


if __name__ == "__main__":
    main()
