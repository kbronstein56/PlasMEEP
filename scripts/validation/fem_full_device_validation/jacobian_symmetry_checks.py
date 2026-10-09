#!/usr/bin/env python3
"""Physics-based derivative checks on an exactly y-mirrored two-rod mesh.

At B=0 with source and probe on the mirror axis, dJ/ds for mirror-pair rods
must agree (reciprocal isotropic operator + exact geometric symmetry).

At ±B, ε_xy → −ε_xy under B→−B is equivalent to reflecting the gyrotropy.
For this axial setup the pairing is dJ/ds0(+B) = dJ/ds1(−B).
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
C = 299792458.0
REL_TOL = 1e-4  # same order as Jacobian gate once geometry is exact


def build_exact_mirror_mesh():
    """Mesh the upper half, then reflect across y=mid to force exact symmetry."""
    vend = ROOT / "scripts" / "validation" / "fem_meep_validation" / "_vendor"
    sys.path.insert(0, str(vend))
    import triangle

    nx, ny, dp = 10.0, 8.0, 1.0
    mid = ny / 2.0
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, dp
    c_up = np.array([5.0, mid + 1.2])
    radius = 0.35
    m = 40
    # upper-half domain including the midline
    box = np.array([[0, mid], [nx, mid], [nx, ny], [0, ny]], float)
    ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
    ring = c_up + radius * np.column_stack([np.cos(ang), np.sin(ang)])
    # keep only vertices with y >= mid - 1e-12; clip ring crossings
    verts = np.vstack([box, ring])
    segs = [[0, 1], [1, 2], [2, 3], [3, 0]] + [[4 + i, 4 + (i + 1) % m] for i in range(m)]
    regions = [[c_up[0], c_up[1], 1, 0.45 * 0.05**2], [1.0, mid + 0.5, 2, 0.45 * 0.12**2]]
    mesh = triangle.triangulate(
        {"vertices": verts, "segments": np.asarray(segs, np.int32), "regions": np.asarray(regions, float)},
        "pq20aA",
    )
    pts_u = np.asarray(mesh["vertices"], float)
    tris_u = np.asarray(mesh["triangles"], int)
    # drop degenerate / below mid
    cents_u = pts_u[tris_u].mean(1)
    keep = cents_u[:, 1] >= mid - 1e-10
    tris_u = tris_u[keep]

    # reflect: new points, map upper midline points to themselves
    # Build full vertex list: upper verts + reflected strict-upper verts
    on_mid = np.abs(pts_u[:, 1] - mid) < 1e-10
    pts_u[on_mid, 1] = mid
    refl_idx = -np.ones(len(pts_u), int)
    new_pts = [pts_u]
    next_id = len(pts_u)
    for i, p in enumerate(pts_u):
        if on_mid[i]:
            refl_idx[i] = i
        else:
            refl_idx[i] = next_id
            new_pts.append(np.array([[p[0], 2 * mid - p[1]]]))
            next_id += 1
    pts = np.vstack(new_pts)
    # lower triangles with reflected node ids (flip orientation)
    tris_l = []
    for t in tris_u:
        a, b, c = int(refl_idx[t[0]]), int(refl_idx[t[1]]), int(refl_idx[t[2]])
        tris_l.append([a, c, b])
    tris = np.vstack([tris_u, np.asarray(tris_l, int)])
    centers = [c_up, np.array([c_up[0], 2 * mid - c_up[1]])]
    cents = pts[tris].mean(1)
    masks = [np.linalg.norm(cents - c, axis=1) <= radius + 1e-9 for c in centers]
    return pts, tris, masks, centers, mid


def rho_from(masks, s_vec, b_tesla, n):
    rho = [np.ones(n, np.complex128), np.zeros(n, np.complex128), np.zeros(n, np.complex128), np.ones(n, np.complex128)]
    f = gamma = fc = None
    for mask, s in zip(masks, s_vec):
        f, fp, gamma, fc = ordinary_from_si(sc.fs_Hz, sc.fp_Hz * np.sqrt(s), sc.gamma_Hz, b_tesla, sc.a)
        e = tensor_ordinary(f, fp, gamma, fc)[:4]
        r = rho_xy(*e)[:4]
        for i in range(4):
            rho[i][mask] = r[i]
    return rho, f, gamma, fc


def dJ_ds(pts, tris, masks, b_tesla, k0, bvec, probe):
    s0 = np.ones(len(masks))
    rho, f, gamma, fc = rho_from(masks, s0, b_tesla, len(tris))
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    lu = splu(A.tocsc())
    x = lu.solve(bvec)
    sampler = ElementSampler(pts, tris)
    hz, _, _ = sampler.fields(x, *rho, k0, probe[None, :])
    y = complex(hz[0])
    fp_ref = sc.fp_Hz * sc.a / C
    grads = []
    for k in range(len(masks)):
        dr = element_drho(masks[k], f, fp_ref, gamma, fc, 1.0)
        dA = assemble_anisotropic(
            pts, tris, *dr, k0, np.zeros(len(pts), dtype=bool), mass_scale=0.0, pin_empty=False
        )
        dx = lu.solve(-(dA @ x))
        dhz, _, _ = sampler.fields(dx, *([np.zeros(len(tris), np.complex128)] * 4), k0, probe[None, :])
        dJ = 2.0 * np.real(np.conj(y) * complex(dhz[0]))
        grads.append(float(dJ))
    return grads, abs(y)


def main():
    pts, tris, masks, centers, mid = build_exact_mirror_mesh()
    sc.fs_a = sc.fs_Hz * sc.a / C
    k0 = 2 * np.pi * sc.fs_a
    sampler = ElementSampler(pts, tris)
    src = np.array([1.8, mid])
    probe = np.array([8.2, mid])
    bvec = np.zeros(len(pts), np.complex128)
    t = int(sampler.locate(src[None, :])[0])
    if t < 0:
        # nudge off the boundary into the mesh
        src = np.array([1.8, mid + 1e-4])
        t = int(sampler.locate(src[None, :])[0])
    w = sampler._bary(t, src)
    for k in range(3):
        bvec[int(tris[t, k])] += complex(w[k])
    # symmetrize source: also inject mirror of source if source not exactly on mid
    src_m = np.array([src[0], 2 * mid - src[1]])
    tm = int(sampler.locate(src_m[None, :])[0])
    if tm >= 0 and tm != t:
        wm = sampler._bary(tm, src_m)
        for k in range(3):
            bvec[int(tris[tm, k])] += complex(wm[k])
        bvec *= 0.5

    print(
        f"nodes={len(pts)} tris={len(tris)} rod_tris={[int(m.sum()) for m in masks]} mid={mid}",
        flush=True,
    )

    g0, _ = dJ_ds(pts, tris, masks, 0.0, k0, bvec, probe)
    sym_rel = abs(g0[0] - g0[1]) / max(0.5 * (abs(g0[0]) + abs(g0[1])), 1e-30)
    print(f"B=0 mirror grads {g0} rel_diff {sym_rel:.3e}", flush=True)

    gp, _ = dJ_ds(pts, tris, masks, 0.05, k0, bvec, probe)
    gm, _ = dJ_ds(pts, tris, masks, -0.05, k0, bvec, probe)
    pair_rel = abs(gp[0] - gm[1]) / max(0.5 * (abs(gp[0]) + abs(gm[1])), 1e-30)
    pair_rel2 = abs(gp[1] - gm[0]) / max(0.5 * (abs(gp[1]) + abs(gm[0])), 1e-30)
    print(f"+B grads {gp} -B grads {gm}", flush=True)
    print(f"pair dJ0(+B) vs dJ1(-B) {pair_rel:.3e}; dJ1(+B) vs dJ0(-B) {pair_rel2:.3e}", flush=True)

    out = {
        "B0_mirror_grads": g0,
        "B0_mirror_rel": float(sym_rel),
        "B0_pass": bool(sym_rel <= REL_TOL),
        "Bp_grads": gp,
        "Bm_grads": gm,
        "B_reversal_pair_rel": [float(pair_rel), float(pair_rel2)],
        "B_reversal_pass": bool(pair_rel <= REL_TOL and pair_rel2 <= REL_TOL),
        "criteria_rel": REL_TOL,
        "note": "Exact reflected mesh. Identity: axial FoM, mirror rods, B→−B swaps rod indices.",
    }
    out["pass"] = bool(out["B0_pass"] and out["B_reversal_pass"])
    (OUT / "jacobian_symmetry_checks.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE jacobian_symmetry_checks.json", out["pass"], flush=True)


if __name__ == "__main__":
    main()
