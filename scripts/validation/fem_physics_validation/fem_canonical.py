#!/usr/bin/env python3
"""Guide, port, Poynting, PML-strength, and reciprocity sweeps.

Thresholds are the frozen ones in PASS_CRITERIA.md.
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
import fem_validated_solver as fv  # noqa: E402
from analytic_maxwell import fresnel_ht, kz_of  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402
from planar_fem import assign_rho, dirichlet, l2_error, rect_mesh, solve_dirichlet  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def plasma(f_a: float) -> complex:
    eps, _ = gyrotropic_drude_eps_eta(float(f_a), float(sc.fp_a), float(sc.gamma_a), 0.0)
    return complex(eps)


def residual(A, x, b):
    num = np.linalg.norm(A @ x - b)
    den = np.linalg.norm(b) + 1e-30
    return float(num / den)


def poynting_cut(pts, uh, k0, eps, x_cut):
    """Integrated Sx on a vertical grid line. Ey = (i/k0/eps) (-dHz/dx) from the local triangle."""
    m = np.abs(pts[:, 0] - x_cut) < 1e-10
    # Use neighboring nodes one step to the right for a one-sided difference along the row.
    ys = np.unique(np.round(pts[m, 1], 10))
    if len(ys) < 2:
        return None
    dy = float(ys[1] - ys[0])
    acc = 0.0
    for y in ys:
        i0 = np.where(m & (np.abs(pts[:, 1] - y) < 1e-8))[0]
        i1 = np.where((np.abs(pts[:, 0] - (x_cut + dy)) < 1e-8) & (np.abs(pts[:, 1] - y) < 1e-8))[0]
        if len(i0) != 1 or len(i1) != 1:
            # h may differ in x and y. Find the next x node.
            return None
        dhz = (uh[i1[0]] - uh[i0[0]]) / dy
        ey = (1j / k0) * (1.0 / eps) * (-dhz)
        sx = 0.5 * np.real(ey * np.conj(uh[i0[0]]))
        acc += sx * dy
    return float(acc)


def guide_matrix():
    """Several widths and frequencies. Plates are natural Neumann."""
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    rows = []
    specs = []
    for width in (0.8, 1.0, 1.4):
        for k0, tag in ((3.2, "low"), (5.0, "mid"), (7.5, "high")):
            specs.append((width, k0, tag))
    for width, k0, tag in specs:
        ky = np.pi / width
        beta = kz_of(k0, 1.0, ky)
        x0, y0 = 2.2, 2.2
        length = 1.6

        def ufn(x, y, ky=ky, beta=beta, x0=x0, y0=y0):
            return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

        for h in (0.04, 0.02, 0.01):
            if abs(width / h - round(width / h)) > 1e-8 or abs(length / h - round(length / h)) > 1e-8:
                continue
            pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
            bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
            rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0)
            values = ufn(pts[:, 0], pts[:, 1])
            A = fv.assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
            b = np.zeros(len(pts), dtype=np.complex128)
            A, b = dirichlet(A, b, bound, values)
            uh = splu(A.tocsc()).solve(b)
            rel = l2_error(pts, tris, uh, ufn)[1]
            # Analytic guide-normal power. |cos|^2 integrates to width/2. eps=1, omega=k0.
            p_exact = 0.5 * np.real(beta) * (width / 2.0) / k0
            # Sample power at three interior stations by interpolating the exact nodal difference
            # only when h equals the x spacing, which it does.
            powers = {}
            for station, xv in (("near", x0 + 0.4), ("mid", x0 + 0.8), ("far", x0 + 1.2)):
                if abs((xv - x0) / h - round((xv - x0) / h)) > 1e-8:
                    continue
                # Rebuild a local difference using node spacing h.
                col = np.abs(pts[:, 0] - xv) < 1e-10
                ys = pts[col, 1]
                order = np.argsort(ys)
                idx = np.flatnonzero(col)[order]
                idx_r = []
                ok = True
                for y, i in zip(ys[order], idx):
                    j = np.where((np.abs(pts[:, 0] - (xv + h)) < 1e-10) & (np.abs(pts[:, 1] - y) < 1e-8))[0]
                    if len(j) != 1:
                        ok = False
                        break
                    idx_r.append(int(j[0]))
                if not ok:
                    continue
                idx_r = np.array(idx_r)
                dhz = (uh[idx_r] - uh[idx]) / h
                ey = (1j / k0) * (-dhz)
                sx = 0.5 * np.real(ey * np.conj(uh[idx]))
                dy = h
                powers[station] = float(np.sum(sx) * dy)
            # E/H at the guide center: Ey/Hz should be beta/k0 for eps=1.
            center = np.array([x0 + 0.8, y0 + 0.5 * width])
            ic = int(np.argmin((pts[:, 0] - center[0]) ** 2 + (pts[:, 1] - center[1]) ** 2))
            # Use the mid-station derivative if present.
            eh_rel = None
            if "mid" in powers:
                exact_ratio = beta / k0
                # Reconstruct ratio from the same mid column's mean |ey/hz| is noisy at a node of cos.
                # Use a node near y = y0 + width/4 where cos is not zero.
                yq = y0 + 0.25 * width
                i0 = np.where((np.abs(pts[:, 0] - (x0 + 0.8)) < 1e-10) & (np.abs(pts[:, 1] - yq) < 0.51 * h))[0]
                i1 = np.where((np.abs(pts[:, 0] - (x0 + 0.8 + h)) < 1e-10) & (np.abs(pts[:, 1] - yq) < 0.51 * h))[0]
                if len(i0) and len(i1):
                    dhz = (uh[i1[0]] - uh[i0[0]]) / h
                    ey = (1j / k0) * (-dhz)
                    if abs(uh[i0[0]]) > 1e-8:
                        eh_rel = float(abs(ey / uh[i0[0]] - exact_ratio) / abs(exact_ratio))
            rec = {
                "width": width,
                "k0": k0,
                "tag": tag,
                "h": h,
                "beta": [beta.real, beta.imag],
                "cutoff": bool(abs(beta.real) < 1e-8),
                "rel_L2": rel,
                "power_exact": p_exact if abs(beta.real) > 1e-8 else 0.0,
                "powers": powers,
                "EH_rel": eh_rel,
                "residual": residual(A, uh, b),
                "dofs": int(len(pts)),
            }
            if powers and abs(p_exact) > 1e-8:
                rec["power_rel_mid"] = abs(powers.get("mid", np.nan) - p_exact) / abs(p_exact)
            rows.append(rec)
            print("guide", width, k0, h, "L2", f"{rel:.3e}", "Pres", rec.get("power_rel_mid"), "EH", eh_rel, flush=True)
    return rows


def homogeneous_power():
    saved = float(sc.fs_a)
    k0 = 2 * np.pi * saved
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    rows = []
    cases = {
        "eps1": 1.0 + 0j,
        "eps3.8": 3.8 + 0j,
        "eps_complex": 2.0 + 0.3j,
        "plasma_fs": plasma(saved),
    }
    for name, eps in cases.items():
        kx = kz_of(k0, eps)

        def ufn(x, y, kx=kx):
            return np.exp(1j * kx * (x - 1.5))

        for h in (0.05, 0.025, 0.0125):
            pts, tris, _, _ = rect_mesh(1.5, 2.5, 1.5, 2.2, h)
            bound = (
                (np.abs(pts[:, 0] - 1.5) < 1e-10) | (np.abs(pts[:, 0] - 2.5) < 1e-10)
                | (np.abs(pts[:, 1] - 1.5) < 1e-10) | (np.abs(pts[:, 1] - 2.2) < 1e-10)
            )
            rho = assign_rho(pts[tris].mean(1), lambda x, y, eps=eps: eps)
            values = ufn(pts[:, 0], pts[:, 1])
            A = fv.assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
            b = np.zeros(len(pts), dtype=np.complex128)
            A, b = dirichlet(A, b, bound, values)
            uh = splu(A.tocsc()).solve(b)
            rel = l2_error(pts, tris, uh, ufn)[1]
            # Phase per meter from two interior nodes.
            x_a, x_b = 1.8, 2.2
            ia = int(np.argmin((pts[:, 0] - x_a) ** 2 + (pts[:, 1] - 1.85) ** 2))
            ib = int(np.argmin((pts[:, 0] - x_b) ** 2 + (pts[:, 1] - 1.85) ** 2))
            dphi = np.angle(uh[ib] / uh[ia])
            dx = pts[ib, 0] - pts[ia, 0]
            beta_fem = dphi / dx
            # Power on the line x=2.0, height 0.7. Analytic Sx = 0.5 Re(kx/(k0 eps)) |Hz|^2.
            col = np.abs(pts[:, 0] - 2.0) < 1e-10
            # one-sided dx = h
            acc = 0.0
            nsum = 0
            for i in np.flatnonzero(col):
                j = np.where((np.abs(pts[:, 0] - (2.0 + h)) < 1e-10) & (np.abs(pts[:, 1] - pts[i, 1]) < 1e-8))[0]
                if len(j) != 1:
                    continue
                dhz = (uh[j[0]] - uh[i]) / h
                ey = (1j / k0) * (1.0 / eps) * (-dhz)
                sx = 0.5 * np.real(ey * np.conj(uh[i]))
                acc += sx
                nsum += 1
            sx_fem = acc / max(nsum, 1)
            hz_mean = np.mean(uh[col])
            sx_exact = 0.5 * np.real(kx / (k0 * eps)) * abs(hz_mean) ** 2
            rows.append({
                "case": name,
                "h": h,
                "rel_L2": rel,
                "beta_exact": [kx.real, kx.imag],
                "beta_fem": float(beta_fem),
                "beta_rel": float(abs(beta_fem - kx) / (abs(kx) + 1e-30)),
                "Sx_fem": float(sx_fem),
                "Sx_exact_using_fem_amplitude": float(sx_exact),
                "Sx_rel": float(abs(sx_fem - sx_exact) / (abs(sx_exact) + 1e-30)),
                "residual": residual(A, uh, b),
            })
            print("homog", name, h, "L2", f"{rel:.3e}", "beta", rows[-1]["beta_rel"], "Sx", rows[-1]["Sx_rel"], flush=True)
    sc.fs_a = saved
    return rows


def pml_strength_and_padding():
    """Quartz slab transmission versus PML strength scale and air gap before the PML."""
    f_a = float(sc.fs_a)
    k0 = 2 * np.pi * f_a
    rows = []
    # Reference geometry: source 2.2, slab [4.0, 4.4], probe 5.2, domain [0,10]x, dpml=1.
    def one(dp, scale, x_shift):
        fv.PML_SIGMA_SCALE = scale
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = 10.0, 4.0, dp
        h = 0.02
        pts, tris, xs, _ = rect_mesh(0.0, 10.0, 1.6, 2.4, h)
        sy = fv.pml_sx_sy(pts[tris].mean(1))[1]
        if np.max(np.abs(sy - 1)) > 1e-6:
            raise RuntimeError("y PML in the strip")
        b = np.zeros(len(pts), dtype=np.complex128)
        xsnd = 2.2 + x_shift
        # Keep the source and the slab on grid lines. x_shift is 0 or -0.4.
        b[np.abs(pts[:, 0] - xsnd) < 1e-10] = 1.0
        pec = np.zeros(len(pts), dtype=bool)
        front = 4.0 + x_shift

        def eps_of(x, y, front=front):
            return 3.8 + 0j if front <= x < front + 0.4 else 1.0

        fields = {}
        for tag, fn in (("vac", lambda x, y: 1.0), ("slab", eps_of)):
            rho = assign_rho(pts[tris].mean(1), fn)
            A = fv.assemble_anisotropic(pts, tris, *rho, k0, pec)
            fields[tag] = splu(A.tocsc()).solve(b)
        xt = front + 0.4 + 0.8
        xt = xs[np.argmin(np.abs(xs - xt))]
        m = np.abs(pts[:, 0] - xt) < 1e-10
        ratio = np.mean(fields["slab"][m]) / np.mean(fields["vac"][m])
        fv.PML_SIGMA_SCALE = 1.0
        return ratio

    ref = one(1.0, 1.0, 0.0)
    for dp in (0.6, 1.0, 1.4):
        for scale in (0.5, 1.0, 2.0):
            if dp != 1.0 and scale != 1.0:
                continue
            z = one(dp, scale, 0.0)
            rows.append({
                "dpml": dp,
                "sigma_scale": scale,
                "db_vs_ref": float(20 * np.log10(abs(z) / abs(ref))),
                "phase_deg_vs_ref": float(np.angle(z / ref) * 180 / np.pi),
                "T_abs": float(abs(z)),
            })
            print("pml", rows[-1], flush=True)
    # Move the whole experiment 0.4 toward the PML (less air padding on the left).
    z = one(1.0, 1.0, -0.4)
    rows.append({
        "dpml": 1.0,
        "sigma_scale": 1.0,
        "air_pad_shift": -0.4,
        "db_vs_ref": float(20 * np.log10(abs(z) / abs(ref))),
        "phase_deg_vs_ref": float(np.angle(z / ref) * 180 / np.pi),
    })
    print("pml pad", rows[-1], flush=True)
    return rows


def reciprocity_pairs():
    """Several source/receiver pairs on a quartz inclusion. Direct Gab versus Gba."""
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 4.0, 4.0, 1.0
    k0 = 2 * np.pi * float(sc.fs_a)
    h = 0.025
    pts, tris, _, _ = rect_mesh(1.3, 2.7, 1.3, 2.7, h)
    bound = (
        (np.abs(pts[:, 0] - 1.3) < 1e-10) | (np.abs(pts[:, 0] - 2.7) < 1e-10)
        | (np.abs(pts[:, 1] - 1.3) < 1e-10) | (np.abs(pts[:, 1] - 2.7) < 1e-10)
    )
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 3.8 if (x - 2.0) ** 2 + (y - 2.0) ** 2 <= 0.18 ** 2 else 1.0)
    A = fv.assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    A, _ = dirichlet(A, np.zeros(len(pts), np.complex128), bound, np.zeros(len(pts), np.complex128))
    lu = splu(A.tocsc())
    targets = [(1.55, 1.7), (2.45, 1.7), (1.55, 2.35), (2.45, 2.35), (2.0, 1.5), (2.0, 2.5)]
    idx = [int(np.argmin((pts[:, 0] - x) ** 2 + (pts[:, 1] - y) ** 2)) for x, y in targets]
    rows = []
    for i, a in enumerate(idx):
        ea = np.zeros(len(pts), np.complex128)
        ea[a] = 1.0
        ga = lu.solve(ea)
        for b in idx[i + 1:]:
            eb = np.zeros(len(pts), np.complex128)
            eb[b] = 1.0
            gb = lu.solve(eb)
            diff = abs(ga[b] - gb[a])
            rows.append({"a": targets[i], "b": targets[idx.index(b)], "abs_diff": float(diff), "rel": float(diff / (abs(ga[b]) + 1e-30))})
    print("reciprocity pairs", len(rows), "max rel", max(r["rel"] for r in rows), flush=True)
    return rows


def mirror_symmetry():
    """A symmetric two-circle layout must satisfy Hz(x,y)=Hz(x,-y) on one solve."""
    # Deferred to the scatterer driver, which has the body-fitted mesh.
    return {"note": "measured in fem_scatterers mirror case"}


def fresnel_power_table():
    k0 = 2 * np.pi * float(sc.fs_a)
    eps_p = plasma(float(sc.fs_a))
    rows = []
    for name, e1, e2, angles in (
        ("air_quartz", 1.0, 3.8 + 0j, (0.0, 20.0, 40.0)),
        ("quartz_air", 3.8 + 0j, 1.0, (0.0, 15.0, 40.0)),
        ("air_complex", 1.0, 2.0 + 0.3j, (0.0, 25.0, 50.0)),
        ("air_plasma", 1.0, eps_p, (0.0, 20.0, 40.0)),
    ):
        for ang in angles:
            ky = kz_of(k0, e1) * np.sin(np.deg2rad(ang))
            fr = fresnel_ht(k0, e1, e2, ky)
            rows.append({
                "case": name,
                "angle_deg": ang,
                "r": [fr["r"].real, fr["r"].imag],
                "t": [fr["t"].real, fr["t"].imag],
                "R_power": fr["R_power"],
                "T_power": fr["T_power"],
                "power_sum": fr["power_sum"],
            })
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fv.PML_SIGMA_SCALE = 1.0
    out = {
        "guide": guide_matrix(),
        "homogeneous_power": homogeneous_power(),
        "pml": pml_strength_and_padding(),
        "reciprocity_pairs": reciprocity_pairs(),
        "fresnel_coefficients": fresnel_power_table(),
        "mirror": mirror_symmetry(),
    }
    (OUT / "canonical_suite.json").write_text(json.dumps(out, indent=2) + "\n")
    print("CANONICAL_DONE", flush=True)


if __name__ == "__main__":
    main()
