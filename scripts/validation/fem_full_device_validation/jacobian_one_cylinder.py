#!/usr/bin/env python3
"""One-cylinder forward Jacobian vs centered FD, including ε/ρ FD and B=±0.05 T.

Frozen thresholds: JACOBIAN_PASS_CRITERIA.md (rel ≤ 1e-4 with truncation visible).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from gyrotropic_tensor import ordinary_from_si, rho_xy, tensor_ordinary  # noqa: E402
from plasma_sensitivity import deps_ds, drho_ds, element_drho  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_full_device_validation"
HS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
REL_TOL = 1e-4
C = 299792458.0


def build_mesh():
    vend = ROOT / "scripts" / "validation" / "fem_meep_validation" / "_vendor"
    sys.path.insert(0, str(vend))
    import triangle

    nx = ny = 8.0
    dp = 1.0
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, dp
    center = np.array([4.0, 4.0])
    radius = 0.40
    m = 80
    ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
    ring = center + radius * np.column_stack([np.cos(ang), np.sin(ang)])
    box = np.array([[0, 0], [nx, 0], [nx, ny], [0, ny]], float)
    verts = np.vstack([box, ring])
    segs = [[0, 1], [1, 2], [2, 3], [3, 0]] + [[4 + i, 4 + (i + 1) % m] for i in range(m)]
    regions = [
        [center[0], center[1], 1, 0.45 * 0.04**2],
        [center[0] + 1.2, center[1], 2, 0.45 * 0.08**2],
        [0.3, 0.3, 3, 0.45 * 0.16**2],
    ]
    mesh = triangle.triangulate(
        {"vertices": verts, "segments": np.asarray(segs, np.int32), "regions": np.asarray(regions, float)},
        "pq20aA",
    )
    pts = np.asarray(mesh["vertices"], float)
    tris = np.asarray(mesh["triangles"], int)
    origin = pts[tris].mean(1) - center
    inside = np.linalg.norm(origin, axis=1) <= radius
    return pts, tris, inside, center


def constitutive_fd_check(b_tesla: float, s0: float = 1.0):
    """Independent FD of ε(s) and ρ(s) before FEM assembly blame."""
    f, fp_ref_ord, gamma, fc = ordinary_from_si(sc.fs_Hz, sc.fp_Hz, sc.gamma_Hz, b_tesla, sc.a)
    # analytic at s0
    de = deps_ds(f, fp_ref_ord, gamma, fc, s0)
    dr = drho_ds(f, fp_ref_ord, gamma, fc, s0)
    rows = []
    for h in HS:
        def eps_at(s):
            fp = fp_ref_ord * np.sqrt(s)
            return tensor_ordinary(f, fp, gamma, fc)[:4]

        def rho_at(s):
            e = eps_at(s)
            return rho_xy(*e)[:4]

        ep = eps_at(s0 + h)
        em = eps_at(s0 - h)
        rp = rho_at(s0 + h)
        rm = rho_at(s0 - h)
        de_fd = [(ep[i] - em[i]) / (2 * h) for i in range(4)]
        dr_fd = [(rp[i] - rm[i]) / (2 * h) for i in range(4)]
        names = ["exx", "exy", "eyx", "eyy"]
        eps_rel = []
        rho_rel = []
        for i in range(4):
            a = de[i] if i < 2 or True else de[i]
            # de tuple is (dexx, dexy, deyx, deyy, fp) — first 4
            ana_e = de[i]
            ana_r = dr[i]
            eps_rel.append(abs(de_fd[i] - ana_e) / max(abs(ana_e), 1e-30))
            rho_rel.append(abs(dr_fd[i] - ana_r) / max(abs(ana_r), 1e-30))
        rows.append({
            "h": h,
            "eps_rel_max": float(max(eps_rel)),
            "rho_rel_max": float(max(rho_rel)),
            "eps_rel": {names[i]: float(eps_rel[i]) for i in range(4)},
            "rho_rel": {names[i]: float(rho_rel[i]) for i in range(4)},
        })
    emin = min(r["eps_rel_max"] for r in rows)
    rmin = min(r["rho_rel_max"] for r in rows)
    return {
        "B_T": b_tesla,
        "eps_min_rel": emin,
        "rho_min_rel": rmin,
        "eps_pass": emin <= 1e-6,
        "rho_pass": rmin <= 1e-6,
        "sweep": rows,
    }


def rho_elements(inside, s, b_tesla):
    f, fp, gamma, fc = ordinary_from_si(sc.fs_Hz, sc.fp_Hz * np.sqrt(s), sc.gamma_Hz, b_tesla, sc.a)
    exx, exy, eyx, eyy, _, _ = tensor_ordinary(f, fp, gamma, fc)
    rxx, rxy, ryx, ryy, _ = rho_xy(exx, exy, eyx, eyy)
    n = len(inside)
    rho = [np.ones(n, np.complex128), np.zeros(n, np.complex128), np.zeros(n, np.complex128), np.ones(n, np.complex128)]
    rho[0][inside] = rxx
    rho[1][inside] = rxy
    rho[2][inside] = ryx
    rho[3][inside] = ryy
    return rho, f, fp, gamma, fc


def solve_state(pts, tris, inside, s, b_tesla, k0, bvec):
    rho, f, fp, gamma, fc = rho_elements(inside, s, b_tesla)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    x = splu(A.tocsc()).solve(bvec)
    return x, A, (f, fp, gamma, fc), rho


def gate_sweep(errs):
    imin = int(np.argmin(errs))
    minimum = errs[imin]
    larger = max(errs[: max(imin, 1)]) if imin > 0 else (errs[1] if len(errs) > 1 else errs[0])
    # Prefer a larger-h truncation marker among h > h_min when possible
    if imin > 0:
        larger = max(errs[:imin])
    else:
        larger = max(errs[1:]) if len(errs) > 1 else errs[0]
    return imin, minimum, larger, bool(minimum <= REL_TOL and larger >= 10.0 * minimum)


def run_case(pts, tris, inside, center, b_tesla, k0, bvec, sampler, probe):
    t0 = time.time()
    x, A, ords, rho0 = solve_state(pts, tris, inside, 1.0, b_tesla, k0, bvec)
    f, fp, gamma, fc = ords
    fp_ref_ord = sc.fp_Hz * sc.a / C
    dr = element_drho(inside, f, fp_ref_ord, gamma, fc, 1.0)
    dA = assemble_anisotropic(
        pts, tris, *dr, k0, np.zeros(len(pts), dtype=bool), mass_scale=0.0, pin_empty=False
    )
    lu = splu(A.tocsc())
    dx = lu.solve(-(dA @ x))
    hz, ex, ey = sampler.fields(x, *rho0, k0, probe[None, :])
    dhz, _, _ = sampler.fields(dx, *([np.zeros(len(tris), np.complex128)] * 4), k0, probe[None, :])
    direct_sample = complex(dhz[0])
    # product-rule dE would need dρ; here check Hz and |Hz|^2 only as frozen criteria
    j = abs(hz[0]) ** 2
    dj = 2.0 * np.real(np.conj(hz[0]) * direct_sample)
    # complex scattering proxy: transmitted probe
    sweep = []
    for h in HS:
        xp, _, _, _ = solve_state(pts, tris, inside, 1.0 + h, b_tesla, k0, bvec)
        xm, _, _, _ = solve_state(pts, tris, inside, 1.0 - h, b_tesla, k0, bvec)
        dx_fd = (xp - xm) / (2.0 * h)
        rel = float(np.linalg.norm(dx_fd - dx) / np.linalg.norm(dx))
        hp, _, _ = sampler.fields(xp, *([np.zeros(len(tris))] * 4), k0, probe[None, :])
        hm, _, _ = sampler.fields(xm, *([np.zeros(len(tris))] * 4), k0, probe[None, :])
        d_sample = (hp[0] - hm[0]) / (2.0 * h)
        sample_rel = abs(d_sample - direct_sample) / max(abs(direct_sample), 1e-30)
        # real/imag separately
        sample_re_rel = abs(d_sample.real - direct_sample.real) / max(abs(direct_sample.real), 1e-30)
        sample_im_rel = abs(d_sample.imag - direct_sample.imag) / max(abs(direct_sample.imag), 1e-30)
        mag = abs(hz[0])
        dmag = np.real(np.conj(hz[0]) * direct_sample) / max(mag, 1e-30)
        mag_p, mag_m = abs(hp[0]), abs(hm[0])
        dmag_fd = (mag_p - mag_m) / (2.0 * h)
        mag_rel = abs(dmag_fd - dmag) / max(abs(dmag), 1e-30)
        phase = np.angle(hz[0])
        # d(arg)/ds = Im(conj(h) dh) / |h|^2
        dphase = np.imag(np.conj(hz[0]) * direct_sample) / max(j, 1e-30)
        dphase_fd = (np.angle(hp[0]) - np.angle(hm[0])) / (2.0 * h)
        # unwrap small jumps
        if abs(dphase_fd) > np.pi / h:
            dphase_fd = ((np.angle(hp[0]) - np.angle(hm[0]) + np.pi) % (2 * np.pi) - np.pi) / (2.0 * h)
        phase_rel = abs(dphase_fd - dphase) / max(abs(dphase), 1e-30)
        jp, jm = abs(hp[0]) ** 2, abs(hm[0]) ** 2
        dj_fd = (jp - jm) / (2.0 * h)
        j_rel = abs(dj_fd - dj) / max(abs(dj), 1e-30)
        sweep.append({
            "h": h,
            "field_rel": rel,
            "sample_rel": float(sample_rel),
            "sample_re_rel": float(sample_re_rel),
            "sample_im_rel": float(sample_im_rel),
            "mag_rel": float(mag_rel),
            "phase_rel": float(phase_rel),
            "power_rel": float(j_rel),
        })
        print(
            f"B={b_tesla:+.2f} h={h:.1e} field {rel:.3e} sample {sample_rel:.3e} "
            f"power {j_rel:.3e} mag {mag_rel:.3e} phase {phase_rel:.3e}",
            flush=True,
        )
    field_errs = [r["field_rel"] for r in sweep]
    imin, minimum, larger, passed = gate_sweep(field_errs)
    sample_errs = [r["sample_rel"] for r in sweep]
    _, smin, slarger, spass = gate_sweep(sample_errs)
    power_errs = [r["power_rel"] for r in sweep]
    _, pmin, plarger, ppass = gate_sweep(power_errs)
    return {
        "B_T": b_tesla,
        "nodes": int(len(pts)),
        "plasma_tris": int(inside.sum()),
        "minimum_field_rel": minimum,
        "h_at_minimum": sweep[imin]["h"],
        "truncation_visible": bool(larger >= 10.0 * minimum),
        "field_pass": bool(passed),
        "sample_pass": bool(spass),
        "power_pass": bool(ppass),
        "pass": bool(passed and spass and ppass),
        "sample_min_rel": smin,
        "power_min_rel": pmin,
        "sweep": sweep,
        "sample_abs": abs(complex(hz[0])),
        "dj": float(dj),
        "wall_s": time.time() - t0,
    }


def main():
    pts, tris, inside, center = build_mesh()
    sc.fs_a = sc.fs_Hz * sc.a / C
    k0 = 2 * np.pi * sc.fs_a
    sampler = ElementSampler(pts, tris)
    src = center + np.array([-1.5, 0.2])
    bvec = np.zeros(len(pts), np.complex128)
    t_src = int(sampler.locate(src[None, :])[0])
    w = sampler._bary(t_src, src)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    probe = center + np.array([1.3, 0.4])
    print(f"nodes={len(pts)} inside={int(inside.sum())}", flush=True)

    const_rows = []
    for b in (0.0, 0.05, -0.05):
        c = constitutive_fd_check(b)
        const_rows.append(c)
        print(f"CONST B={b} eps_min={c['eps_min_rel']:.3e} rho_min={c['rho_min_rel']:.3e} "
              f"pass={c['eps_pass'] and c['rho_pass']}", flush=True)

    fem_rows = []
    for b in (0.0, 0.05, -0.05):
        row = run_case(pts, tris, inside, center, b, k0, bvec, sampler, probe)
        fem_rows.append(row)
        print("CASE", row["B_T"], "PASS" if row["pass"] else "FAIL",
              "min", row["minimum_field_rel"], flush=True)

    out = {
        "criteria": {"rel_tol": REL_TOL, "hs": list(HS)},
        "constitutive_fd": const_rows,
        "fem_jacobian": fem_rows,
        "overall_pass": all(c["eps_pass"] and c["rho_pass"] for c in const_rows)
        and all(r["pass"] for r in fem_rows),
    }
    (OUT / "jacobian_one_cylinder.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE jacobian_one_cylinder.json overall", out["overall_pass"], flush=True)


if __name__ == "__main__":
    main()
