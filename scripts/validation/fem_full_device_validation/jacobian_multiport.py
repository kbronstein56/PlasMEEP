#!/usr/bin/env python3
"""Reduced multiport: derivatives of complex amplitudes, phase, and Poynting power.

Two open sides (source + receive), one plasma disk parameter s.
"""
from __future__ import annotations

import json
import sys
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
from plasma_sensitivity import element_drho  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_full_device_validation"
HS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
REL_TOL = 1e-4
C = 299792458.0


def build_mesh():
    vend = ROOT / "scripts" / "validation" / "fem_meep_validation" / "_vendor"
    sys.path.insert(0, str(vend))
    import triangle

    nx, ny, dp = 12.0, 6.0, 1.0
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, dp
    center = np.array([6.0, 3.0])
    radius = 0.45
    m = 72
    ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
    ring = center + radius * np.column_stack([np.cos(ang), np.sin(ang)])
    box = np.array([[0, 0], [nx, 0], [nx, ny], [0, ny]], float)
    verts = np.vstack([box, ring])
    segs = [[0, 1], [1, 2], [2, 3], [3, 0]] + [[4 + i, 4 + (i + 1) % m] for i in range(m)]
    regions = [
        [center[0], center[1], 1, 0.45 * 0.05**2],
        [1.5, 3.0, 2, 0.45 * 0.10**2],
        [0.3, 0.3, 3, 0.45 * 0.18**2],
    ]
    mesh = triangle.triangulate(
        {"vertices": verts, "segments": np.asarray(segs, np.int32), "regions": np.asarray(regions, float)},
        "pq20aA",
    )
    pts = np.asarray(mesh["vertices"], float)
    tris = np.asarray(mesh["triangles"], int)
    inside = np.linalg.norm(pts[tris].mean(1) - center, axis=1) <= radius
    return pts, tris, inside, center, (nx, ny, dp)


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


def line_weights(n=31, span=1.8):
    off = np.linspace(-span / 2, span / 2, n)
    w = np.full(n, span / (n - 1))
    w[0] *= 0.5
    w[-1] *= 0.5
    return off, w


def inject(pts, tris, sampler, xy, weights):
    b = np.zeros(len(pts), np.complex128)
    tid = sampler.locate(xy)
    for n in range(len(xy)):
        t = int(tid[n])
        if t < 0:
            continue
        bw = sampler._bary(t, xy[n])
        for k in range(3):
            b[int(tris[t, k])] += float(bw[k]) * float(weights[n])
    return b


def amp(sampler, uh, rho, k0, xy, weights):
    hz, _, _ = sampler.fields(uh, *rho, k0, xy)
    return complex(np.dot(weights, hz))


def flux_x(sampler, uh, rho, k0, xy, coord):
    """Guide-normal Poynting through a vertical line (n = +x)."""
    hz, ex, ey = sampler.fields(uh, *rho, k0, xy)
    # S_x = 0.5 Re(Ey Hz*)
    sx = 0.5 * np.real(ey * np.conj(hz))
    return float(np.trapezoid(sx, coord))


def solve_state(pts, tris, inside, s, b_tesla, k0, bvec):
    rho, f, fp, gamma, fc = rho_elements(inside, s, b_tesla)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    x = splu(A.tocsc()).solve(bvec)
    return x, A, rho, (f, gamma, fc)


def gate(errs):
    imin = int(np.argmin(errs))
    minimum = errs[imin]
    larger = max(errs[:imin]) if imin > 0 else (max(errs[1:]) if len(errs) > 1 else errs[0])
    return imin, minimum, larger, bool(minimum <= REL_TOL and larger >= 10.0 * minimum)


def main():
    pts, tris, inside, center, (nx, ny, dp) = build_mesh()
    sc.fs_a = sc.fs_Hz * sc.a / C
    k0 = 2 * np.pi * sc.fs_a
    sampler = ElementSampler(pts, tris)
    off, w = line_weights()
    # source line near left interior; receive near right
    src_xy = np.column_stack([np.full(len(off), 2.0), 3.0 + off])
    mon_xy = np.column_stack([np.full(len(off), 10.0), 3.0 + off])
    bvec = inject(pts, tris, sampler, src_xy, w)
    print(f"nodes={len(pts)} plasma={int(inside.sum())}", flush=True)

    cases = []
    for b_tesla in (0.0, 0.05):
        x, A, rho, ords = solve_state(pts, tris, inside, 1.0, b_tesla, k0, bvec)
        f, gamma, fc = ords
        fp_ref = sc.fp_Hz * sc.a / C
        dr = element_drho(inside, f, fp_ref, gamma, fc, 1.0)
        dA = assemble_anisotropic(
            pts, tris, *dr, k0, np.zeros(len(pts), dtype=bool), mass_scale=0.0, pin_empty=False
        )
        lu = splu(A.tocsc())
        dx = lu.solve(-(dA @ x))
        a0 = amp(sampler, x, rho, k0, mon_xy, w)
        # d amp / ds from Hz only (linear in x); rho at monitor is vacuum so no dρ term
        da = amp(sampler, dx, rho, k0, mon_xy, w)
        p0 = flux_x(sampler, x, rho, k0, mon_xy, 3.0 + off)
        # power is real; differentiate via complex calculus on fields at the line
        hz, ex, ey = sampler.fields(x, *rho, k0, mon_xy)
        dhz, dex, dey = sampler.fields(dx, *rho, k0, mon_xy)  # vacuum ρ at monitor
        # d(Sx)/ds = 0.5 Re( dEy Hz* + Ey dHz* ) = 0.5 Re( dey conj(hz) + ey conj(dhz) )
        dsx = 0.5 * np.real(dey * np.conj(hz) + ey * np.conj(dhz))
        dp = float(np.trapezoid(dsx, 3.0 + off))
        # scalar objective J = |a|^2
        j = abs(a0) ** 2
        dj = 2.0 * np.real(np.conj(a0) * da)
        phase = np.angle(a0)
        dphase = np.imag(np.conj(a0) * da) / max(j, 1e-30)

        sweep = []
        for h in HS:
            xp, _, rhop, _ = solve_state(pts, tris, inside, 1.0 + h, b_tesla, k0, bvec)
            xm, _, rhom, _ = solve_state(pts, tris, inside, 1.0 - h, b_tesla, k0, bvec)
            ap = amp(sampler, xp, rhop, k0, mon_xy, w)
            am = amp(sampler, xm, rhom, k0, mon_xy, w)
            da_fd = (ap - am) / (2.0 * h)
            amp_rel = abs(da_fd - da) / max(abs(da), 1e-30)
            amp_re = abs(da_fd.real - da.real) / max(abs(da.real), 1e-30)
            amp_im = abs(da_fd.imag - da.imag) / max(abs(da.imag), 1e-30)
            pp = flux_x(sampler, xp, rhop, k0, mon_xy, 3.0 + off)
            pm = flux_x(sampler, xm, rhom, k0, mon_xy, 3.0 + off)
            dp_fd = (pp - pm) / (2.0 * h)
            power_rel = abs(dp_fd - dp) / max(abs(dp), 1e-30)
            jp, jm = abs(ap) ** 2, abs(am) ** 2
            dj_fd = (jp - jm) / (2.0 * h)
            j_rel = abs(dj_fd - dj) / max(abs(dj), 1e-30)
            dph_fd = (np.angle(ap) - np.angle(am)) / (2.0 * h)
            # unwrap
            raw = np.angle(ap) - np.angle(am)
            raw = (raw + np.pi) % (2 * np.pi) - np.pi
            dph_fd = raw / (2.0 * h)
            phase_rel = abs(dph_fd - dphase) / max(abs(dphase), 1e-30)
            field_rel = float(
                np.linalg.norm(((xp - xm) / (2.0 * h)) - dx) / np.linalg.norm(dx)
            )
            sweep.append({
                "h": h,
                "field_rel": field_rel,
                "amp_rel": float(amp_rel),
                "amp_re_rel": float(amp_re),
                "amp_im_rel": float(amp_im),
                "power_rel": float(power_rel),
                "obj_rel": float(j_rel),
                "phase_rel": float(phase_rel),
            })
            print(
                f"B={b_tesla} h={h:.1e} field {field_rel:.3e} amp {amp_rel:.3e} "
                f"power {power_rel:.3e} obj {j_rel:.3e} phase {phase_rel:.3e}",
                flush=True,
            )

        def pack(key):
            errs = [r[key] for r in sweep]
            imin, minimum, larger, passed = gate(errs)
            return {
                "min": minimum,
                "h": sweep[imin]["h"],
                "truncation": bool(larger >= 10.0 * minimum),
                "pass": bool(passed),
            }

        gates = {k: pack(k) for k in ("field_rel", "amp_rel", "power_rel", "obj_rel", "phase_rel")}
        ok = all(g["pass"] for g in gates.values())
        cases.append({
            "B_T": b_tesla,
            "nodes": int(len(pts)),
            "a0_abs": abs(a0),
            "p0": p0,
            "dj": float(dj),
            "dp": float(dp),
            "gates": gates,
            "pass": bool(ok),
            "sweep": sweep,
        })
        print("CASE", b_tesla, "PASS" if ok else "FAIL", flush=True)

    out = {"criteria_rel": REL_TOL, "cases": cases, "overall_pass": all(c["pass"] for c in cases)}
    (OUT / "jacobian_multiport.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE jacobian_multiport.json", out["overall_pass"], flush=True)


if __name__ == "__main__":
    main()
