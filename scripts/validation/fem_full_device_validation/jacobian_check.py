#!/usr/bin/env python3
"""Direct density Jacobian versus centered differences on one plasma disk.

Thresholds are frozen in JACOBIAN_PASS_CRITERIA.md.
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
HS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
REL_TOL = 1e-4


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
    return x, A, (f, fp, gamma, fc)


def main():
    pts, tris, inside, center = build_mesh()
    sc.fs_a = sc.fs_Hz * sc.a / 299792458.0
    k0 = 2 * np.pi * sc.fs_a
    sampler = ElementSampler(pts, tris)
    src = center + np.array([-1.5, 0.2])
    bvec = np.zeros(len(pts), np.complex128)
    t_src = int(sampler.locate(src[None, :])[0])
    w = sampler._bary(t_src, src)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    probe = center + np.array([1.3, 0.4])
    rows = []
    print(f"nodes={len(pts)} inside={int(inside.sum())}", flush=True)
    for b_tesla in (0.0, 0.05):
        x, A, ords = solve_state(pts, tris, inside, 1.0, b_tesla, k0, bvec)
        f, fp, gamma, fc = ords
        # fp here is fp_ref*sqrt(s)=fp_ref. element_drho wants fp_ref.
        dr = element_drho(inside, f, sc.fp_Hz * sc.a / 299792458.0, gamma, fc, 1.0)
        dA = assemble_anisotropic(
            pts, tris, *dr, k0, np.zeros(len(pts), dtype=bool), mass_scale=0.0, pin_empty=False
        )
        lu = splu(A.tocsc())
        dx = lu.solve(-(dA @ x))
        hz, _, _ = sampler.fields(x, *rho_elements(inside, 1.0, b_tesla)[0], k0, probe[None, :])
        # derivative of the sample from the same basis as a nodal interpolation of dx
        dhz, _, _ = sampler.fields(dx, *([np.zeros(len(tris), np.complex128)] * 4), k0, probe[None, :])
        # fields() reconstructs E from rho; Hz itself is the interpolant and does not use rho.
        direct_sample = complex(dhz[0])
        j = abs(hz[0]) ** 2
        dj = 2.0 * np.real(np.conj(hz[0]) * direct_sample)
        sweep = []
        for h in HS:
            xp, _, _ = solve_state(pts, tris, inside, 1.0 + h, b_tesla, k0, bvec)
            xm, _, _ = solve_state(pts, tris, inside, 1.0 - h, b_tesla, k0, bvec)
            dx_fd = (xp - xm) / (2.0 * h)
            rel = float(np.linalg.norm(dx_fd - dx) / np.linalg.norm(dx))
            hp, _, _ = sampler.fields(xp, *([np.zeros(len(tris))] * 4), k0, probe[None, :])
            hm, _, _ = sampler.fields(xm, *([np.zeros(len(tris))] * 4), k0, probe[None, :])
            d_sample = (hp[0] - hm[0]) / (2.0 * h)
            sample_rel = abs(d_sample - direct_sample) / max(abs(direct_sample), 1e-30)
            jp = abs(hp[0]) ** 2
            jm = abs(hm[0]) ** 2
            dj_fd = (jp - jm) / (2.0 * h)
            j_rel = abs(dj_fd - dj) / max(abs(dj), 1e-30)
            sweep.append({"h": h, "field_rel": rel, "sample_rel": float(sample_rel), "power_rel": float(j_rel)})
            print(f"B={b_tesla} h={h:.1e} field {rel:.3e} sample {sample_rel:.3e} power {j_rel:.3e}", flush=True)
        field_errs = [r["field_rel"] for r in sweep]
        imin = int(np.argmin(field_errs))
        minimum = field_errs[imin]
        larger = max(field_errs[: max(imin, 1)])
        passed = minimum <= REL_TOL and larger >= 10.0 * minimum
        rows.append({
            "B_T": b_tesla,
            "nodes": int(len(pts)),
            "minimum_field_rel": minimum,
            "h_at_minimum": sweep[imin]["h"],
            "truncation_visible": bool(larger >= 10.0 * minimum),
            "pass": bool(passed),
            "sweep": sweep,
            "sample_abs": abs(complex(hz[0])),
            "dj": float(dj),
        })
        print("CASE", rows[-1]["B_T"], "PASS" if passed else "FAIL", "min", minimum, flush=True)
    (OUT / "jacobian_fd_comparison.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("WROTE jacobian_fd_comparison.json", flush=True)


if __name__ == "__main__":
    main()
