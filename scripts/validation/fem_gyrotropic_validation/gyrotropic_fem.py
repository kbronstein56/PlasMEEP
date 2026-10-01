#!/usr/bin/env python3
"""FEM checks of the gyrotropic Hz operator. Subcommands: mms, homo, slab, onsager, cylinder."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu
from scipy.special import hankel1, h1vp, jv, jvp

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
from gyrotropic_tensor import rho_xy, tensor_ordinary, voigt_eps_eff  # noqa: E402
from mms import run_case  # noqa: E402
from planar_fem import assign_rho, dirichlet, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_gyrotropic_validation"


def rho_tuple(f, fp, gamma, fc):
    exx, exy, eyx, eyy, _, eta = tensor_ordinary(f, fp, gamma, fc)
    rxx, rxy, ryx, ryy, _ = rho_xy(exx, exy, eyx, eyy)
    return (rxx, rxy, ryx, ryy), eta, exx


def run_mms():
    # Propagating-like and production-frequency tensors. The load is manufactured.
    cases = {
        "prod_B0": (0.2568443533025771, 0.5337025523170433, 6.671281903963041e-05, 0.0),
        "prod_Bp": (0.2568443533025771, 0.5337025523170433, 6.671281903963041e-05, 0.08),
        "prod_Bm": (0.2568443533025771, 0.5337025523170433, 6.671281903963041e-05, -0.08),
        "prop_Bp": (0.45, 0.15, 1.0e-4, 0.05),
    }
    out = {}
    for name, spec in cases.items():
        rho, eta, _ = rho_tuple(*spec)
        out[name] = run_case(name, rho, [16, 32, 64])
        out[name]["eta_abs"] = abs(eta)
    (OUT / "gyrotropic_mms.json").write_text(json.dumps(out, indent=2) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.6))
    for name, rec in out.items():
        h = [r["h"] for r in rec["levels"]]
        ax[0].loglog(h, [r["L2"] for r in rec["levels"]], "o-", label=name)
        ax[1].loglog(h, [r["H1"] for r in rec["levels"]], "o-", label=name)
    href = np.array([1 / 16, 1 / 64])
    ax[0].loglog(href, 0.2 * href**2, "k--", lw=0.8)
    ax[1].loglog(href, 0.8 * href, "k--", lw=0.8)
    ax[0].invert_xaxis()
    ax[1].invert_xaxis()
    ax[0].set_xlabel("h")
    ax[1].set_xlabel("h")
    ax[0].set_ylabel("L2")
    ax[1].set_ylabel("H1")
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "gyrotropic_mms.png", dpi=140)
    print("MMS_DONE", flush=True)


def _beta(k0, eps_eff):
    k = k0 * np.sqrt(eps_eff)
    if np.real(k) < 0:
        k = -k
    return complex(k)


def run_homo():
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    f, fp, gamma = 0.45, 0.15, 1.0e-4
    rows = []
    for fc in (0.0, 0.05, -0.05):
        rho, eta, eperp = rho_tuple(f, fp, gamma, fc)
        k0 = 2 * np.pi * f
        sc.fs_a = f
        beta = _beta(k0, voigt_eps_eff(eperp, eta))
        for h in (0.02, 0.01):
            x0, y0, length, width = 2.2, 2.2, 1.2, 0.8
            pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
            bound = (
                (np.abs(pts[:, 0] - x0) < 1e-10)
                | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
                | (np.abs(pts[:, 1] - y0) < 1e-10)
                | (np.abs(pts[:, 1] - (y0 + width)) < 1e-10)
            )

            def ufn(x, y, beta=beta, x0=x0):
                return np.exp(1j * beta * (x - x0))

            T = len(tris)
            rho_e = [np.full(T, z) for z in rho]
            A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
            b = np.zeros(len(pts), np.complex128)
            A, b = dirichlet(A, b, bound, ufn(pts[:, 0], pts[:, 1]))
            uh = splu(A.tocsc()).solve(b)
            exact = ufn(pts[:, 0], pts[:, 1])
            rel = float(np.linalg.norm(uh - exact) / np.linalg.norm(exact))
            cents = pts[tris].mean(1)
            twice = (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1]) - (
                pts[tris[:, 2], 0] - pts[tris[:, 0], 0]
            ) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
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
            gx = np.sum(bx * uh[tris], axis=1)
            gy = np.sum(by * uh[tris], axis=1)
            rxx, rxy, ryx, ryy = rho
            ex = (1j / k0) * (rxx * gy + rxy * (-gx))
            ey = (1j / k0) * (ryx * gy + ryy * (-gx))
            u_c = ufn(cents[:, 0], cents[:, 1])
            ex_a = (beta / k0) * rxy * u_c
            ey_a = (beta / k0) * ryy * u_c
            e_rel = float(np.linalg.norm(np.column_stack([ex - ex_a, ey - ey_a])) / (np.linalg.norm(np.column_stack([ex_a, ey_a])) + 1e-30))
            rows.append({
                "fc": fc,
                "h": h,
                "dofs": int(len(pts)),
                "hz_rel": rel,
                "e_rel": e_rel,
                "beta": [beta.real, beta.imag],
                "ex_over_hz_fem": [ (ex/u_c).mean().real, (ex/u_c).mean().imag ],
                "ex_over_hz_exact": [(ex_a/u_c).mean().real, (ex_a/u_c).mean().imag],
            })
            print("homo", rows[-1], flush=True)
    (OUT / "gyrotropic_homo.json").write_text(json.dumps(rows, indent=2) + "\n")


def slab_response(k0, eps_eff, thickness):
    """Normal-incidence Hz slab. eps_eff is the Voigt permittivity 1/ρ_xx."""
    k2 = _beta(k0, eps_eff)
    phi = k2 * thickness
    c, s = np.cos(phi), np.sin(phi)

    def back(hz0, flux0):
        hz = c * hz0 + (eps_eff / k2) * s * flux0
        flux = c * flux0 - (k2 / eps_eff) * s * hz0
        return flux - 1j * k0 * hz, hz

    b0, _ = back(1.0, 1j * k0)
    br, _ = back(1.0, -1j * k0)
    r = -b0 / br
    hz0 = 1.0 + r
    flux0 = 1j * k0 * (1.0 - r)
    hz_d = c * hz0 + (eps_eff / k2) * s * flux0
    return r, hz_d


def run_slab():
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 12.0, 6.0, 1.0
    f, fp, gamma = 0.45, 0.15, 1.0e-4
    k0 = 2 * np.pi * f
    sc.fs_a = f
    thickness = 0.40
    rows = []
    for fc in (0.0, 0.05, -0.05):
        rho, eta, eperp = rho_tuple(f, fp, gamma, fc)
        eff = voigt_eps_eff(eperp, eta)
        r_a, t_a = slab_response(k0, eff, thickness)
        # Box 0..6 in local coordinates, placed at x=2.2 so PML is idle.
        x0, y0 = 2.2, 2.2
        h = 0.02
        length, width = 6.0, 0.8
        pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
        local_x = pts[:, 0] - x0
        # Slab occupies local x in [2.0, 2.4].
        cents = pts[tris].mean(1)
        lx = cents[:, 0] - x0
        inside = (lx >= 2.0) & (lx <= 2.4)
        T = len(tris)
        rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
        rho_e[0][inside] = rho[0]
        rho_e[1][inside] = rho[1]
        rho_e[2][inside] = rho[2]
        rho_e[3][inside] = rho[3]
        # Analytic field on the boundary: incident+reflected for x<2, interior not needed on boundary if boundary is outside.
        def field_at(x_local):
            if x_local <= 2.0:
                return np.exp(1j * k0 * (x_local - 2.0)) + r_a * np.exp(-1j * k0 * (x_local - 2.0))
            if x_local >= 2.4:
                return t_a * np.exp(1j * k0 * (x_local - 2.4))
            # interior reconstruction from the left state
            k2 = _beta(k0, eff)
            # state at slab entrance
            hz_l = 1.0 + r_a
            flux_l = 1j * k0 * (1.0 - r_a)
            d = x_local - 2.0
            kd = k2 * d
            return np.cos(kd) * hz_l + (eff / k2) * np.sin(kd) * flux_l

        bound = (
            (np.abs(local_x) < 1e-10)
            | (np.abs(local_x - length) < 1e-10)
            | (np.abs(pts[:, 1] - y0) < 1e-10)
            | (np.abs(pts[:, 1] - (y0 + width)) < 1e-10)
        )
        values = np.array([field_at(x) for x in local_x], dtype=np.complex128)
        A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
        b = np.zeros(len(pts), np.complex128)
        A, b = dirichlet(A, b, bound, values)
        uh = splu(A.tocsc()).solve(b)
        rel = float(np.linalg.norm(uh - values) / np.linalg.norm(values))
        # Transmission sample on the downstream Dirichlet edge, compared with analytic t.
        down = np.abs(local_x - length) < 1e-10
        # The boundary is exact by construction. Sample an interior station x_local=4.0.
        station = np.argmin((local_x - 4.0) ** 2 + (pts[:, 1] - (y0 + 0.4)) ** 2)
        t_fem = uh[station] / np.exp(1j * k0 * (local_x[station] - 2.4))
        rows.append({
            "fc": fc,
            "h": h,
            "dofs": int(len(pts)),
            "field_rel": rel,
            "T_db": float(20 * np.log10(abs(t_fem) / abs(t_a))),
            "T_phase_deg": float(np.angle(t_fem / t_a) * 180 / np.pi),
            "analytic_R_abs": abs(r_a),
            "analytic_T_abs": abs(t_a),
        })
        print("slab", rows[-1], flush=True)
    (OUT / "gyrotropic_slab.json").write_text(json.dumps(rows, indent=2) + "\n")


def run_onsager():
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 6.0, 6.0, 1.0
    f, fp, gamma = 0.45, 0.15, 1.0e-3
    k0 = 2 * np.pi * f
    sc.fs_a = f
    h = 0.025
    x0, y0 = 2.0, 2.0
    pts, tris, _, _ = rect_mesh(x0, x0 + 1.5, y0, y0 + 1.2, h)
    cents = pts[tris].mean(1)
    disk = (cents[:, 0] - (x0 + 0.7)) ** 2 + (cents[:, 1] - (y0 + 0.55)) ** 2 <= 0.18 ** 2
    bound = (
        (np.abs(pts[:, 0] - x0) < 1e-10)
        | (np.abs(pts[:, 0] - (x0 + 1.5)) < 1e-10)
        | (np.abs(pts[:, 1] - y0) < 1e-10)
        | (np.abs(pts[:, 1] - (y0 + 1.2)) < 1e-10)
    )

    def assemble(fc):
        rho, _, _ = rho_tuple(f, fp, gamma, fc)
        T = len(tris)
        rho_e = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
        for k in range(4):
            rho_e[k][disk] = rho[k]
        A = assemble_anisotropic(pts, tris, *rho_e, k0, np.zeros(len(pts), dtype=bool))
        A, _ = dirichlet(A, np.zeros(len(pts), np.complex128), bound, np.zeros(len(pts)))
        return A.tocsr()

    free = np.flatnonzero(~bound)

    def interior(A):
        return A[free][:, free]

    def resid(fc):
        # Dirichlet row replacement makes the full matrix unsymmetric.
        # The interior block is the operator under test.
        num = interior(assemble(fc)) - interior(assemble(-fc)).T
        den = interior(assemble(fc))
        return float(np.linalg.norm(num.data) / (np.linalg.norm(den.data) + 1e-30))

    A0 = interior(assemble(0.0))
    sym0 = float(np.linalg.norm((A0 - A0.T).data) / (np.linalg.norm(A0.data) + 1e-30))
    # Green exchange.
    Ap = assemble(0.05)
    Am = assemble(-0.05)
    lup = splu(Ap.tocsc())
    lum = splu(Am.tocsc())
    a = int(np.argmin((pts[:, 0] - (x0 + 0.25)) ** 2 + (pts[:, 1] - (y0 + 0.4)) ** 2))
    b = int(np.argmin((pts[:, 0] - (x0 + 1.15)) ** 2 + (pts[:, 1] - (y0 + 0.85)) ** 2))
    ea = np.zeros(len(pts), np.complex128); ea[a] = 1
    eb = np.zeros(len(pts), np.complex128); eb[b] = 1
    gap = lup.solve(ea)
    gbm = lum.solve(eb)
    # G(b,a; +B) versus G(a,b; -B)
    onsager_g = abs(gap[b] - gbm[a]) / abs(gap[b])
    # At +B, G(a,b) need not equal G(b,a).
    gbp = lup.solve(eb)
    nonrecip = abs(gap[b] - gbp[a]) / abs(gap[b])
    rec = {
        "dofs": int(len(pts)),
        "matrix_onsager_fc_0.05": resid(0.05),
        "matrix_symmetric_B0": sym0,
        "green_onsager_rel": float(onsager_g),
        "green_same_B_nonrecip_rel": float(nonrecip),
    }
    print("onsager", rec, flush=True)
    (OUT / "gyrotropic_onsager.json").write_text(json.dumps(rec, indent=2) + "\n")


def gyrotropic_T(k0, radius, rho_xx, rho_xy, nmax):
    kp = k0 / np.sqrt(rho_xx)
    ns = np.arange(-nmax, nmax + 1)
    T = np.zeros(len(ns), dtype=np.complex128)
    for i, n in enumerate(ns):
        n = int(n)
        j = jv(n, k0 * radius)
        jp = jvp(n, k0 * radius, 1)
        h = hankel1(n, k0 * radius)
        hp = h1vp(n, k0 * radius, 1)
        jin = jv(n, kp * radius)
        jnp = jvp(n, kp * radius, 1)
        alpha = rho_xx * kp * jnp + rho_xy * (1j * n / radius) * jin
        right = alpha / jin
        T[i] = (right * j - k0 * jp) / (k0 * hp - right * h)
    return T


def exterior_field(xy, center, k0, radius, rho_xx, rho_xy, src, nmax):
    T = gyrotropic_T(k0, radius, rho_xx, rho_xy, nmax)
    xy = np.atleast_2d(np.asarray(xy, float))
    rel = xy - center
    rho = np.linalg.norm(rel, axis=1)
    phi = np.arctan2(rel[:, 1], rel[:, 0])
    srel = np.asarray(src, float) - center
    rho_s = np.linalg.norm(srel)
    phi_s = np.arctan2(srel[1], srel[0])
    inc = hankel1(0, k0 * np.linalg.norm(xy - src, axis=1))
    scat = np.zeros(len(xy), np.complex128)
    ns = np.arange(-nmax, nmax + 1)
    for i, n in enumerate(ns):
        alpha = hankel1(int(n), k0 * rho_s) * np.exp(-1j * int(n) * phi_s)
        scat += T[i] * alpha * hankel1(int(n), k0 * rho) * np.exp(1j * int(n) * phi)
    return inc + scat


def run_cylinder():
    from analytic_maxwell import mie_T
    from fem_scatterers import mesh_circles, to_mesh

    f, fp, gamma = 0.45, 0.15, 1.0e-4
    radius = 0.230
    src = np.array([-4.5, 0.0])
    probes = {
        "forward": np.array([2.8, 0.0]),
        "backward": np.array([-2.8, 0.0]),
        "side_p60": np.array([2.8 * np.cos(np.pi / 3), 2.8 * np.sin(np.pi / 3)]),
        "side_m60": np.array([2.8 * np.cos(-np.pi / 3), 2.8 * np.sin(-np.pi / 3)]),
    }
    names = list(probes)
    xy = np.array(list(probes.values()), float)
    rows = []
    # B=0 reduction of the new T against isotropic Mie.
    k0 = 2 * np.pi * f
    rho0, _, e0 = rho_tuple(f, fp, gamma, 0.0)
    eps0 = 1.0 / rho0[0]
    t_new = gyrotropic_T(k0, radius, rho0[0], rho0[1], 8)
    t_old = mie_T(k0, complex(eps0), radius, 8)
    print("mie_b0_diff", float(np.max(np.abs(t_new - t_old))), flush=True)
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 14.0, 10.0, 1.2
    sc.fs_a = f
    for fc in (0.0, 0.05, -0.05):
        rho, eta, _ = rho_tuple(f, fp, gamma, fc)
        for h, he in ((0.04, 0.012), (0.02, 0.008)):
            if abs(fc) > 0 and h < 0.03 and fc < 0:
                continue
            pts, tris = mesh_circles(np.zeros((1, 2)), [(radius,)], h, he)
            print(f"cylinder fc={fc} nodes {len(pts)}", flush=True)
            origin = pts[tris].mean(1) - np.array([7.0, 5.0])
            inside = np.linalg.norm(origin, axis=1) <= radius
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
            fields = {}
            for tag, rr in (("vac", vac), ("obj", rho_e)):
                A = assemble_anisotropic(pts, tris, *rr, k0, np.zeros(len(pts), dtype=bool))
                uh = splu(A.tocsc()).solve(bvec)
                fields[tag] = uh
            ratio = []
            xy_m = to_mesh(xy)
            for tag in ("vac", "obj"):
                hz, _, _ = sampler.fields(fields[tag], * (vac if tag == "vac" else rho_e), k0, xy_m)
                ratio.append(hz)
            fem_ratio = ratio[1] / ratio[0]
            analytic = exterior_field(xy, np.zeros(2), k0, radius, rho[0], rho[1], src, 10)
            inc = hankel1(0, k0 * np.linalg.norm(xy - src, axis=1))
            ar = analytic / inc
            per = {}
            for i, name in enumerate(names):
                per[name] = {
                    "db": float(20 * np.log10(abs(fem_ratio[i]) / abs(ar[i]))),
                    "phase_deg": float(np.angle(fem_ratio[i] / ar[i]) * 180 / np.pi),
                    "abs_err": float(abs(fem_ratio[i] - ar[i])),
                }
            rows.append({"fc": fc, "h": h, "dofs": int(len(pts)), "probes": per, "mie_b0_max_abs": float(np.max(np.abs(t_new - t_old)))})
            print("cylinder", fc, h, per["forward"], per["side_p60"], per["side_m60"], flush=True)
    (OUT / "gyrotropic_cylinder.json").write_text(json.dumps(rows, indent=2) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1]
    {"mms": run_mms, "homo": run_homo, "slab": run_slab, "onsager": run_onsager, "cylinder": run_cylinder}[cmd]()


if __name__ == "__main__":
    main()
