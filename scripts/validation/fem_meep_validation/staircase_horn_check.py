#!/usr/bin/env python3
"""
One horns-only solve whose PEC boundary is the Meep res-50 pixel staircase.

Pixel (i, j) is metal when its center lies inside a prism quad. That is the
Yee material test Meep uses for eps=-1e20 with subpixel averaging skipped.
The result is a diagnosis of the polygon-vs-raster gap, not a new canonical model.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from b0_campaign import (  # noqa: E402
    MEEP,
    _load_triangle,
    _seg,
    feed_rho,
    horn_prism_quads,
)
from fem_validated_solver import (  # noqa: E402
    LAST_OPERATOR,
    ElementSampler,
    apply_metal_walls_rho,
    assemble_anisotropic,
    guide_normal_flux,
    inject_mode_consistent,
    load_numerical_mode,
    solve_system,
)

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "b0_campaign"


def pixel_mask(nx: float, ny: float, res: int = 50):
    nxp = int(round(nx * res))
    nyp = int(round(ny * res))
    dx, dy = nx / nxp, ny / nyp
    ii, jj = np.meshgrid(np.arange(nxp), np.arange(nyp), indexing="ij")
    centers = np.column_stack([((ii + 0.5) * dx).ravel(), ((jj + 0.5) * dy).ravel()])
    inside = np.zeros(len(centers), dtype=bool)
    for quad in horn_prism_quads(nx, ny):
        inside |= MPath(quad).contains_points(centers, radius=0.0)
    return inside.reshape(nxp, nyp), dx, dy


def boundary_edges(inside: np.ndarray, dx: float, dy: float):
    nxp, nyp = inside.shape
    segs = []
    # Vertical pixel edges.
    for i in range(nxp + 1):
        left = inside[i - 1] if i > 0 else np.zeros(nyp, dtype=bool)
        right = inside[i] if i < nxp else np.zeros(nyp, dtype=bool)
        diff = left != right
        j = 0
        while j < nyp:
            if not diff[j]:
                j += 1
                continue
            j0 = j
            while j < nyp and diff[j]:
                j += 1
            x = i * dx
            segs.append(([x, j0 * dy], [x, j * dy]))
    # Horizontal pixel edges.
    for j in range(nyp + 1):
        below = inside[:, j - 1] if j > 0 else np.zeros(nxp, dtype=bool)
        above = inside[:, j] if j < nyp else np.zeros(nxp, dtype=bool)
        diff = below != above
        i = 0
        while i < nxp:
            if not diff[i]:
                i += 1
                continue
            i0 = i
            while i < nxp and diff[i]:
                i += 1
            y = j * dy
            segs.append(([i0 * dx, y], [i * dx, y]))
    return segs


def mesh_upper(h_bulk: float):
    tr = _load_triangle()
    nx, ny = float(sc.nx_ports), float(sc.ny_ports)
    y0 = ny / 2.0
    inside, dx, dy = pixel_mask(nx, ny)
    edges = boundary_edges(inside, dx, dy)
    verts = []
    index = {}
    segs = []
    segset = set()

    def add_vertex(p):
        x = float(p[0])
        y = float(p[1])
        if x < -1e-8 or x > nx + 1e-8 or y < y0 - 1e-8 or y > ny + 1e-8:
            return None
        if abs(y - y0) <= 1e-8:
            y = y0
        if abs(x) <= 1e-8:
            x = 0.0
        if abs(x - nx) <= 1e-8:
            x = nx
        if abs(y - ny) <= 1e-8:
            y = ny
        key = (round(x, 8), round(y, 8))
        j = index.get(key)
        if j is None:
            j = len(verts)
            index[key] = j
            verts.append([key[0], key[1]])
        return j

    def add_seg(a, b):
        if a is None or b is None or a == b:
            return
        key = (a, b) if a < b else (b, a)
        if key in segset:
            return
        segset.add(key)
        segs.append([a, b])

    for p0, p1 in edges:
        a = np.asarray(p0, float)
        b = np.asarray(p1, float)
        if a[1] < y0 - 1e-8 and b[1] < y0 - 1e-8:
            continue
        if min(a[1], b[1]) < y0 - 1e-8:
            # Crosses the symmetry axis. Keep the upper piece.
            t = (y0 - a[1]) / (b[1] - a[1])
            hit = a + t * (b - a)
            if a[1] >= y0:
                b = hit
            else:
                a = hit
        ia, ib = add_vertex(a), add_vertex(b)
        on_axis = abs(a[1] - y0) <= 1e-8 and abs(b[1] - y0) <= 1e-8
        on_box = (
            (abs(a[0]) <= 1e-8 and abs(b[0]) <= 1e-8)
            or (abs(a[0] - nx) <= 1e-8 and abs(b[0] - nx) <= 1e-8)
            or (abs(a[1] - ny) <= 1e-8 and abs(b[1] - ny) <= 1e-8)
        )
        if not (on_axis or on_box):
            add_seg(ia, ib)

    xs = np.arange(0.0, nx + 0.5 * h_bulk, h_bulk)
    ys = np.arange(y0, ny + 0.5 * h_bulk, h_bulk)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    for p in np.column_stack([xx.ravel(), yy.ravel()]):
        add_vertex(p)
    for p0, p1 in (
        ([0.0, y0], [nx, y0]),
        ([nx, y0], [nx, ny]),
        ([nx, ny], [0.0, ny]),
        ([0.0, ny], [0.0, y0]),
    ):
        for p in _seg(np.asarray(p0, float), np.asarray(p1, float), h_bulk):
            add_vertex(p)

    def boundary_chain(selector, sort_key):
        ids = [i for i, v in enumerate(verts) if selector(v)]
        ids.sort(key=lambda i: sort_key(verts[i]))
        for i in range(len(ids) - 1):
            add_seg(ids[i], ids[i + 1])

    boundary_chain(lambda v: abs(v[1] - y0) <= 1e-8, lambda v: v[0])
    boundary_chain(lambda v: abs(v[0] - nx) <= 1e-8, lambda v: v[1])
    boundary_chain(lambda v: abs(v[1] - ny) <= 1e-8, lambda v: -v[0])
    boundary_chain(lambda v: abs(v[0]) <= 1e-8, lambda v: -v[1])
    amax = 0.45 * h_bulk * h_bulk
    mesh = tr.triangulate(
        {"vertices": np.asarray(verts, float), "segments": np.asarray(segs, np.int32)},
        f"pq20a{amax:.8e}",
    )
    pts = np.asarray(mesh["vertices"], float)
    upper = np.asarray(mesh["triangles"], int)
    on_axis = pts[:, 1] <= y0 + 1e-8
    below = np.flatnonzero(~on_axis)
    mirrored = np.column_stack([pts[below, 0], 2 * y0 - pts[below, 1]])
    full = np.vstack([pts, mirrored])
    mirror_of = np.arange(len(pts))
    mirror_of[below] = np.arange(len(pts), len(pts) + len(below))
    tris = np.vstack([upper, mirror_of[upper]])
    x = full[tris, 0]
    y = full[tris, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    tris = tris[np.abs(twice) > 1e-14]
    return full, tris, inside, dx, dy


def metal_rho(points, tris, inside, dx, dy):
    T = len(tris)
    rho = [
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    ]
    cents = points[tris].mean(1)
    nxp, nyp = inside.shape
    i = np.clip(np.floor(cents[:, 0] / dx).astype(int), 0, nxp - 1)
    j = np.clip(np.floor(cents[:, 1] / dy).astype(int), 0, nyp - 1)
    mask = inside[i, j]
    return list(apply_metal_walls_rho(*rho, mask)), mask


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    points, tris, inside, dx, dy = mesh_upper(0.08)
    print(f"nodes {len(points)} tris {len(tris)} metal pixels {int(inside.sum())}", flush=True)
    k0 = 2 * np.pi * float(sc.fs_a)
    pec = np.zeros(len(points), dtype=bool)
    rho, mask = metal_rho(points, tris, inside, dx, dy)
    sampler = ElementSampler(points, tris)
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    op = dict(LAST_OPERATOR)
    print("OPERATOR", json.dumps(op), flush=True)
    src = load_numerical_mode()
    b = inject_mode_consistent(points, tris, src, sampler)
    x, _lu, t_fac, _ts, resid = solve_system(A, b, pec)
    raw = np.array(
        [
            guide_normal_flux(points, x, None, None, p, sampler=sampler, rho=rho, omega=k0)
            for p in range(6)
        ]
    )
    rho_f = feed_rho(points, tris)
    Aref = assemble_anisotropic(points, tris, *rho_f, k0, pec)
    xref, *_rest = solve_system(Aref, b, pec)
    pinc = abs(guide_normal_flux(points, xref, None, None, 0, sampler=sampler, rho=rho_f, omega=k0))
    norm = raw / pinc
    delta = [
        None if not (np.isfinite(MEEP[i]) and norm[i] > 0) else float(10 * np.log10(norm[i] / MEEP[i]))
        for i in range(6)
    ]
    rec = {
        "model": "meep_res50_pixel_staircase",
        "n_nodes": int(len(points)),
        "n_tris": int(len(tris)),
        "metal_pixels": int(inside.sum()),
        "pixel_dx": dx,
        "factor_s": t_fac,
        "residual": resid,
        "P_inc": float(pinc),
        "normalized": norm.tolist(),
        "delta_dB": delta,
        "operator": op,
        "sampler_misses": int(sampler.misses),
    }
    (OUT / "staircase_horn.json").write_text(json.dumps(rec, indent=2) + "\n")
    print("DELTA", [None if v is None else round(v, 3) for v in delta], "Pinc", round(pinc, 5), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
