#!/usr/bin/env python3
"""Multi-cylinder Jacobian: independent s_i per rod, column indexing check.

J_mk = dy_m / dp_k compared to centered FD. Frozen rel ≤ 1e-4 with truncation.
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

    nx = ny = 10.0
    dp = 1.0
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, dp
    centers = [np.array([4.0, 5.0]), np.array([5.5, 5.0]), np.array([7.0, 5.0])]
    radius = 0.35
    m = 64
    box = np.array([[0, 0], [nx, 0], [nx, ny], [0, ny]], float)
    verts = [box]
    segs = [[0, 1], [1, 2], [2, 3], [3, 0]]
    n0 = 4
    for c in centers:
        ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
        ring = c + radius * np.column_stack([np.cos(ang), np.sin(ang)])
        verts.append(ring)
        segs += [[n0 + i, n0 + (i + 1) % m] for i in range(m)]
        n0 += m
    verts = np.vstack(verts)
    regions = [[c[0], c[1], i + 1, 0.45 * 0.05**2] for i, c in enumerate(centers)]
    regions.append([0.4, 0.4, 10, 0.45 * 0.18**2])
    mesh = triangle.triangulate(
        {"vertices": verts, "segments": np.asarray(segs, np.int32), "regions": np.asarray(regions, float)},
        "pq20aA",
    )
    pts = np.asarray(mesh["vertices"], float)
    tris = np.asarray(mesh["triangles"], int)
    cents = pts[tris].mean(1)
    masks = [np.linalg.norm(cents - c, axis=1) <= radius for c in centers]
    return pts, tris, masks, centers


def rho_from_s(n_tris, masks, s_vec, b_tesla):
    rho = [
        np.ones(n_tris, np.complex128),
        np.zeros(n_tris, np.complex128),
        np.zeros(n_tris, np.complex128),
        np.ones(n_tris, np.complex128),
    ]
    f = gamma = fc = None
    for mask, s in zip(masks, s_vec):
        f, fp, gamma, fc = ordinary_from_si(sc.fs_Hz, sc.fp_Hz * np.sqrt(s), sc.gamma_Hz, b_tesla, sc.a)
        exx, exy, eyx, eyy, _, _ = tensor_ordinary(f, fp, gamma, fc)
        rxx, rxy, ryx, ryy, _ = rho_xy(exx, exy, eyx, eyy)
        rho[0][mask] = rxx
        rho[1][mask] = rxy
        rho[2][mask] = ryx
        rho[3][mask] = ryy
    return rho, f, gamma, fc


def solve(pts, tris, masks, s_vec, b_tesla, k0, bvec):
    rho, f, gamma, fc = rho_from_s(len(tris), masks, s_vec, b_tesla)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    x = splu(A.tocsc()).solve(bvec)
    return x, A, rho, (f, gamma, fc)


def gate(errs):
    imin = int(np.argmin(errs))
    minimum = errs[imin]
    larger = max(errs[:imin]) if imin > 0 else (max(errs[1:]) if len(errs) > 1 else errs[0])
    return imin, minimum, larger, bool(minimum <= REL_TOL and larger >= 10.0 * minimum)


def main():
    pts, tris, masks, centers = build_mesh()
    sc.fs_a = sc.fs_Hz * sc.a / C
    k0 = 2 * np.pi * sc.fs_a
    sampler = ElementSampler(pts, tris)
    src = np.array([1.5, 5.0])
    bvec = np.zeros(len(pts), np.complex128)
    t_src = int(sampler.locate(src[None, :])[0])
    w = sampler._bary(t_src, src)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    # three probes, one near each rod (downstream side)
    probes = [c + np.array([0.55, 0.15]) for c in centers]
    n_p = len(masks)
    s0 = np.ones(n_p)
    print(f"nodes={len(pts)} rods={[int(m.sum()) for m in masks]}", flush=True)

    results = []
    for b_tesla in (0.0, 0.05):
        x, A, rho, ords = solve(pts, tris, masks, s0, b_tesla, k0, bvec)
        f, gamma, fc = ords
        fp_ref = sc.fp_Hz * sc.a / C
        lu = splu(A.tocsc())
        # direct Jacobian columns: dx/ds_k and probe samples
        dx_cols = []
        d_samples = np.zeros((len(probes), n_p), np.complex128)
        for k in range(n_p):
            dr = element_drho(masks[k], f, fp_ref, gamma, fc, 1.0)
            dA = assemble_anisotropic(
                pts, tris, *dr, k0, np.zeros(len(pts), dtype=bool), mass_scale=0.0, pin_empty=False
            )
            dx = lu.solve(-(dA @ x))
            dx_cols.append(dx)
            for m, pr in enumerate(probes):
                dhz, _, _ = sampler.fields(dx, *([np.zeros(len(tris), np.complex128)] * 4), k0, pr[None, :])
                d_samples[m, k] = complex(dhz[0])
        # FD columns
        columns = []
        for k in range(n_p):
            sweep = []
            for h in HS:
                sp = s0.copy(); sp[k] = 1.0 + h
                sm = s0.copy(); sm[k] = 1.0 - h
                xp, _, _, _ = solve(pts, tris, masks, sp, b_tesla, k0, bvec)
                xm, _, _, _ = solve(pts, tris, masks, sm, b_tesla, k0, bvec)
                dx_fd = (xp - xm) / (2.0 * h)
                field_rel = float(np.linalg.norm(dx_fd - dx_cols[k]) / np.linalg.norm(dx_cols[k]))
                sample_rels = []
                for m, pr in enumerate(probes):
                    hp, _, _ = sampler.fields(xp, *([np.zeros(len(tris))] * 4), k0, pr[None, :])
                    hm, _, _ = sampler.fields(xm, *([np.zeros(len(tris))] * 4), k0, pr[None, :])
                    d_fd = (hp[0] - hm[0]) / (2.0 * h)
                    sample_rels.append(
                        abs(d_fd - d_samples[m, k]) / max(abs(d_samples[m, k]), 1e-30)
                    )
                # cross-talk: perturbing rod k should not match column j≠k better than column k
                cross = []
                for j in range(n_p):
                    if j == k:
                        continue
                    cross.append(float(np.linalg.norm(dx_fd - dx_cols[j]) / np.linalg.norm(dx_cols[j])))
                sweep.append({
                    "h": h,
                    "field_rel": field_rel,
                    "sample_rel_max": float(max(sample_rels)),
                    "cross_min_field_rel": float(min(cross)) if cross else None,
                })
                print(
                    f"B={b_tesla} col={k} h={h:.1e} field {field_rel:.3e} "
                    f"sample {max(sample_rels):.3e} cross_min {min(cross):.3e}",
                    flush=True,
                )
            errs = [r["field_rel"] for r in sweep]
            imin, minimum, larger, passed = gate(errs)
            # indexing: at best h, own column must beat every other column by ≥10×
            best = sweep[imin]
            index_ok = best["cross_min_field_rel"] is None or best["cross_min_field_rel"] >= 10.0 * best["field_rel"]
            columns.append({
                "rod": k,
                "minimum_field_rel": minimum,
                "h_at_minimum": sweep[imin]["h"],
                "truncation_visible": bool(larger >= 10.0 * minimum),
                "index_ok": bool(index_ok),
                "pass": bool(passed and index_ok),
                "sweep": sweep,
            })
            print("COL", k, "PASS" if columns[-1]["pass"] else "FAIL", minimum, "index", index_ok, flush=True)
        results.append({
            "B_T": b_tesla,
            "nodes": int(len(pts)),
            "n_params": n_p,
            "columns": columns,
            "pass": all(c["pass"] for c in columns),
        })
    out = {"criteria_rel": REL_TOL, "cases": results, "overall_pass": all(c["pass"] for c in results)}
    (OUT / "jacobian_multicylinder.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE jacobian_multicylinder.json", out["overall_pass"], flush=True)


if __name__ == "__main__":
    main()
