#!/usr/bin/env python3
"""Propagating PEC mode, impedance refinement, PML padding, reciprocity, absorption sign."""
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
from analytic_maxwell import fresnel_ht, kz_of, stack_response  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402
from fem_validated_solver import assemble_anisotropic, pml_sx_sy  # noqa: E402
from planar_fem import assign_rho, dirichlet, l2_error, plasma_eps, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def propagating_guide():
    k0 = 5.0
    width = 1.0
    length = 1.2
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 6.0, 6.0, 1.0
    x0, y0 = 2.0, 2.0

    def ufn(x, y):
        return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

    rows = []
    for h in (0.05, 0.025, 0.0125):
        pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
        bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
        rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0)
        from planar_fem import solve_dirichlet
        uh = solve_dirichlet(pts, tris, rho, k0, ufn(pts[:, 0], pts[:, 1]), bound)
        l2, rel = l2_error(pts, tris, uh, ufn)
        rows.append({"h": h, "beta": [beta.real, beta.imag], "rel_L2": rel, "max_nodal": float(np.max(np.abs(uh - ufn(pts[:, 0], pts[:, 1]))))})
        print("guide", rows[-1], flush=True)
    return rows


def impedance_refine():
    saved = float(sc.fs_a)
    k0 = 2 * np.pi * saved
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    rows = []
    for name, eps in (("eps1", 1.0 + 0j), ("plasma_fs", plasma_eps(saved))):
        kx = kz_of(k0, eps)

        def ufn(x, y, kx=kx):
            return np.exp(1j * kx * (x - 1.5))

        h = 0.0125
        pts, tris, _, _ = rect_mesh(1.5, 2.5, 1.5, 2.2, h)
        bound = (
            (np.abs(pts[:, 0] - 1.5) < 1e-10) | (np.abs(pts[:, 0] - 2.5) < 1e-10)
            | (np.abs(pts[:, 1] - 1.5) < 1e-10) | (np.abs(pts[:, 1] - 2.2) < 1e-10)
        )
        rho = assign_rho(pts[tris].mean(1), lambda x, y, eps=eps: eps)
        from planar_fem import solve_dirichlet
        uh = solve_dirichlet(pts, tris, rho, k0, ufn(pts[:, 0], pts[:, 1]), bound)
        l2, rel = l2_error(pts, tris, uh, ufn)
        tri = tris[len(tris) // 3]
        xy = pts[tri]
        twice = (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[2, 0] - xy[0, 0]) * (xy[1, 1] - xy[0, 1])
        bx = np.array([xy[1, 1] - xy[2, 1], xy[2, 1] - xy[0, 1], xy[0, 1] - xy[1, 1]]) / twice
        ey = (1j / k0) * (1.0 / eps) * (-(bx @ uh[tri]))
        exact = kx / (k0 * eps)
        rows.append({"case": name, "h": h, "rel_L2": rel, "EH_rel": float(abs(ey / np.mean(uh[tri]) - exact) / abs(exact))})
        print("EH", rows[-1], flush=True)
    sc.fs_a = saved
    return rows


def pml_padding():
    """Quartz slab transmission versus PML thickness. Physical probe stays fixed."""
    f_a = float(sc.fs_a)
    k0 = 2 * np.pi * f_a
    rows = []
    for dp in (0.5, 1.0, 1.5):
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = 10.0, 4.0, dp
        h = 0.02
        pts, tris, xs, ys = rect_mesh(0.0, 10.0, 1.6, 2.4, h)
        if np.max(np.abs(pml_sx_sy(pts[tris].mean(1))[1] - 1)) > 1e-8:
            raise RuntimeError("y PML entered the strip")
        cents = pts[tris].mean(1)
        b = np.zeros(len(pts), dtype=np.complex128)
        b[np.abs(pts[:, 0] - 2.2) < 1e-10] = 1.0
        pec = np.zeros(len(pts), dtype=bool)

        def eps_of(x, y, dp=dp):
            return 3.8 + 0j if 4.0 <= x < 4.4 else 1.0 + 0j

        results = {}
        for tag, fn in (("vac", lambda x, y: 1.0), ("slab", eps_of)):
            rho = assign_rho(cents, fn)
            A = assemble_anisotropic(pts, tris, *rho, k0, pec)
            results[tag] = splu(A.tocsc()).solve(b)
        x_t = 5.2
        m = np.abs(pts[:, 0] - x_t) < 1e-10
        ratio = np.mean(results["slab"][m]) / np.mean(results["vac"][m])
        rows.append({"dpml": dp, "T_re": ratio.real, "T_im": ratio.imag, "T_abs": abs(ratio)})
        print("pml", rows[-1], flush=True)
    ref = rows[1]["T_re"] + 1j * rows[1]["T_im"]
    for rec, raw in zip(rows, rows):
        z = rec["T_re"] + 1j * rec["T_im"]
        rec["db_vs_dpml_1"] = float(20 * np.log10(abs(z) / abs(ref)))
        rec["phase_vs_dpml_1_deg"] = float(np.angle(z / ref) * 180 / np.pi)
    return rows


def reciprocity():
    """Green's function symmetry for a lossless quartz disk in a Dirichlet box."""
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    k0 = 2 * np.pi * float(sc.fs_a)
    h = 0.02
    pts, tris, _, _ = rect_mesh(1.4, 2.6, 1.4, 2.6, h)
    bound = (
        (np.abs(pts[:, 0] - 1.4) < 1e-10) | (np.abs(pts[:, 0] - 2.6) < 1e-10)
        | (np.abs(pts[:, 1] - 1.4) < 1e-10) | (np.abs(pts[:, 1] - 2.6) < 1e-10)
    )
    cents = pts[tris].mean(1)
    rho = assign_rho(cents, lambda x, y: 3.8 if (x - 2.0) ** 2 + (y - 2.0) ** 2 <= 0.15**2 else 1.0)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    A, _ = dirichlet(A, np.zeros(len(pts), dtype=np.complex128), bound, np.zeros(len(pts), dtype=np.complex128))
    lu = splu(A.tocsc())
    # Interior nodes away from the disk.
    a = int(np.argmin((pts[:, 0] - 1.7) ** 2 + (pts[:, 1] - 2.0) ** 2))
    b = int(np.argmin((pts[:, 0] - 2.3) ** 2 + (pts[:, 1] - 2.15) ** 2))
    ea = np.zeros(len(pts), dtype=np.complex128); ea[a] = 1
    eb = np.zeros(len(pts), dtype=np.complex128); eb[b] = 1
    ga = lu.solve(ea)
    gb = lu.solve(eb)
    return {
        "G_ab": [ga[b].real, ga[b].imag],
        "G_ba": [gb[a].real, gb[a].imag],
        "abs_diff": float(abs(ga[b] - gb[a])),
        "rel": float(abs(ga[b] - gb[a]) / abs(ga[b])),
    }


def absorption_sign():
    """For the lossy plasma slab, reflected plus transmitted power is below 1."""
    f_a = float(sc.fs_a)
    k0 = 2 * np.pi * f_a
    eps, _ = gyrotropic_drude_eps_eta(f_a, float(sc.fp_a), float(sc.gamma_a), 0.0)
    resp = stack_response(k0, [(0.40, complex(eps))])
    # Power reflection |r|^2. Downstream vacuum-normalized t has the same | | as the
    # transmitted/incident field ratio, and both sides are air, so T_power = |t|^2.
    R = abs(resp["r"]) ** 2
    T = abs(resp["t_over_vacuum"]) ** 2
    return {"R": R, "T": T, "R_plus_T": R + T, "absorption_proxy": 1 - (R + T), "eps": [eps.real, eps.imag]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = {
        "propagating_pec": propagating_guide(),
        "impedance_h_0125": impedance_refine(),
        "pml_padding": pml_padding(),
        "reciprocity": reciprocity(),
        "plasma_slab_power": absorption_sign(),
        "lossless_fresnel_power": fresnel_ht(2 * np.pi * float(sc.fs_a), 1.0, 3.8)["power_sum"],
    }
    print(json.dumps(out, indent=2), flush=True)
    (OUT / "conservation.json").write_text(json.dumps(out, indent=2) + "\n")
    print("EXTRAS_DONE", flush=True)


if __name__ == "__main__":
    main()
