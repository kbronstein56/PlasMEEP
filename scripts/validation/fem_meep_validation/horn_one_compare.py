#!/usr/bin/env python3
"""
Export the exact prism polygons Meep receives (horn_walls='prism') and
compare a one-horn B=0 solve: Meep PEC prisms vs FEM rho=0 on those polygons.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath
from scipy.spatial import Delaunay
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import meep as mp  # noqa: E402
import sixport_common as sc  # noqa: E402
from fem_validated_solver import assemble_anisotropic, fields_from_hz  # noqa: E402
from sixport_common import (  # noqa: E402
    effective_port_dir,
    horn_for_port,
    monitor_center_for_port,
    set_geometry_context,
)

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_localization"
WALLS = ("left_flare", "right_flare", "left_feed", "right_feed")


def export_polygons() -> dict:
    set_geometry_context(res=50, horn_walls="prism", grid_offset_cells=(0.0, 0.0), coord_rotation_deg=0.0)
    a = float(sc.a)
    horns = []
    for p in range(6):
        h = horn_for_port(p, 50)
        n = np.asarray(effective_port_dir(p), dtype=float)
        n = n / np.linalg.norm(n)
        t = np.array([-n[1], n[0]])
        walls = []
        for name in WALLS:
            xy = np.asarray(h[name], dtype=float)[:, :2]
            walls.append(
                {
                    "name": name,
                    "vertices_a": xy.tolist(),
                    "vertices_m": (xy * a).tolist(),
                    "kind": "flare" if "flare" in name else "straight_feed",
                }
            )
        horns.append(
            {
                "port": p + 1,
                "outward": n.tolist(),
                "tangent": t.tolist(),
                "throat_a": np.asarray(h["throat_center"], dtype=float).tolist(),
                "source_center_a": np.asarray(h["source_center"], dtype=float).tolist(),
                "monitor_center_a": np.asarray(monitor_center_for_port(p, 50), dtype=float).tolist(),
                "span_a": float(0.96 * sc.clear_width),
                "walls": walls,
                "note": "Vertices are exactly those passed to Add_Prism when horn_walls='prism'. Meep cell is origin-centered.",
            }
        )
    doc = {
        "a_m": a,
        "horn_walls": "prism",
        "coordinate_frame": "Meep origin-centered a-units (cell center = array center)",
        "pec": "mp.perfect_electric_conductor via Add_Prism(PEC=True)",
        "horns": horns,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "meep_horn_polygons.json").write_text(json.dumps(doc, indent=2) + "\n")
    return doc


def all_quads(doc, port=None):
    quads = []
    for h in doc["horns"]:
        if port is not None and h["port"] != port:
            continue
        for w in h["walls"]:
            quads.append(np.asarray(w["vertices_a"], dtype=float))
    return quads


def raster_masks(quads_a, quads_b, pitch=0.01):
    pts = np.vstack(quads_a + quads_b)
    lo = pts.min(0) - 0.3
    hi = pts.max(0) + 0.3
    xs = np.arange(lo[0], hi[0], pitch)
    ys = np.arange(lo[1], hi[1], pitch)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    grid = np.column_stack([xx.ravel(), yy.ravel()])

    def mask(quads):
        m = np.zeros(len(grid), dtype=bool)
        for q in quads:
            m |= MPath(q).contains_points(grid)
        return m

    A = mask(quads_a)
    B = mask(quads_b)
    inter = np.logical_and(A, B).sum()
    union = np.logical_or(A, B).sum()
    cell = pitch * pitch
    return {
        "pitch_a": pitch,
        "intersection_area_a2": float(inter * cell),
        "union_area_a2": float(union * cell),
        "jaccard": float(inter / union) if union else 1.0,
        "symmetric_difference_area_a2": float((union - inter) * cell),
    }


def old_polyline_mask_quads(quads, rad):
    """Approximate the OLD FEM mask: disks along edges, as polygons for raster via points."""
    # Return a point predicate by sampling the same grid later. Here build many small segments' buffers
    # as a set of sample centers; caller uses contains via distance.
    segs = []
    for q in quads:
        for i in range(len(q)):
            segs.append((q[i], q[(i + 1) % len(q)]))
    return segs, rad


def raster_polyline_vs_fill(quads, pitch=0.02, rad=None):
    """Old FEM mask is a tube of radius ~0.55*thickness around prism edges."""
    if rad is None:
        rad = 0.55 * float(sc.wall_thickness)
    pts = np.vstack(quads)
    lo = pts.min(0) - rad - 0.05
    hi = pts.max(0) + rad + 0.05
    xs = np.arange(lo[0], hi[0], pitch)
    ys = np.arange(lo[1], hi[1], pitch)
    ny, nx = len(ys), len(xs)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    grid = np.column_stack([xx.ravel(), yy.ravel()])
    filled = np.zeros(ny * nx, dtype=bool)
    for q in quads:
        filled |= MPath(q).contains_points(grid)
    filled = filled.reshape(ny, nx)
    tube = np.zeros((ny, nx), dtype=bool)
    rpx = max(1, int(np.ceil(rad / pitch)))
    yy_i, xx_i = np.ogrid[-rpx : rpx + 1, -rpx : rpx + 1]
    disk = (xx_i * pitch) ** 2 + (yy_i * pitch) ** 2 <= rad * rad
    for q in quads:
        for i in range(len(q)):
            p0, p1 = q[i], q[(i + 1) % len(q)]
            L = float(np.linalg.norm(p1 - p0))
            n = max(2, int(L / pitch))
            for t in np.linspace(0, 1, n):
                c = p0 * (1 - t) + p1 * t
                ix = int(round((c[0] - lo[0]) / pitch))
                iy = int(round((c[1] - lo[1]) / pitch))
                x0, x1 = max(0, ix - rpx), min(nx, ix + rpx + 1)
                y0, y1 = max(0, iy - rpx), min(ny, iy + rpx + 1)
                dx0, dx1 = x0 - (ix - rpx), rpx + (x1 - ix)
                dy0, dy1 = y0 - (iy - rpx), rpx + (y1 - iy)
                tube[y0:y1, x0:x1] |= disk[dy0:dy1, dx0:dx1]
    inter = np.logical_and(filled, tube).sum()
    union = np.logical_or(filled, tube).sum()
    cell = pitch * pitch
    return {
        "pitch_a": pitch,
        "jaccard_polyline_tube_vs_filled_prism": float(inter / union) if union else 1.0,
        "filled_area_a2": float(filled.sum() * cell),
        "tube_area_a2": float(tube.sum() * cell),
        "symmetric_difference_a2": float((union - inter) * cell),
        "capture_radius_a": rad,
    }


def line_flux(points, Hz, Ex, Ey, center, tangent, n_hat, span, n_points=31):
    offs = np.linspace(-span / 2, span / 2, n_points)
    ds = span / (n_points - 1)
    w = np.full(n_points, ds)
    w[0] *= 0.5
    w[-1] *= 0.5
    total = 0.0
    samples = []
    for s, wi in zip(offs, w):
        xy = center + s * tangent
        i = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
        hx, ex, ey = Hz[i], Ex[i], Ey[i]
        Sx = 0.5 * np.real(ey * np.conj(hx))
        Sy = -0.5 * np.real(ex * np.conj(hx))
        sn = float(n_hat[0] * Sx + n_hat[1] * Sy)
        total += sn * float(wi)
        samples.append({"s": float(s), "Hz_abs": float(abs(hx)), "Hz_phase": float(np.angle(hx)), "Sdotn": sn})
    return total, samples


def run_fem_one_horn(doc) -> dict:
    h = doc["horns"][0]
    quads = all_quads(doc, port=1)
    span = h["span_a"]
    src = np.asarray(h["source_center_a"], dtype=float)
    mon = np.asarray(h["monitor_center_a"], dtype=float)
    throat = np.asarray(h["throat_a"], dtype=float)
    n_hat = np.asarray(h["outward"], dtype=float)
    tang = np.asarray(h["tangent"], dtype=float)
    # local domain around this horn
    pts_w = np.vstack(quads)
    pad = 3.0
    pml_th = 1.0
    lo0 = pts_w.min(0) - pad
    hi0 = pts_w.max(0) + pad
    # Map the local box onto [0, L] so assemble_anisotropic's PML bands
    # (x<dpml, x>L-dpml) absorb outgoing power. A closed cavity gives Re(S)=0.
    shift = -lo0
    quads = [q + shift for q in quads]
    src = src + shift
    mon = mon + shift
    throat = throat + shift
    lo = lo0 + shift
    hi = hi0 + shift
    saved_box = (sc.nx_ports, sc.ny_ports, sc.dpml_ports)
    sc.nx_ports = float(hi[0] - lo[0])
    sc.ny_ports = float(hi[1] - lo[1])
    sc.dpml_ports = pml_th
    hmesh = 0.04
    xs = np.arange(lo[0], hi[0] + hmesh, hmesh)
    ys = np.arange(lo[1], hi[1] + hmesh, hmesh)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    extra = np.vstack(quads)
    points = np.vstack([np.column_stack([xx.ravel(), yy.ravel()]), extra])
    # Polygon vertices can land on grid nodes. Duplicate coordinates produce
    # zero-area triangles and empty matrix rows (exactly singular).
    _, keep = np.unique(np.round(points, 8), axis=0, return_index=True)
    points = points[np.sort(keep)]
    tris = Delaunay(points).simplices
    T = len(tris)
    rho = [np.ones(T, dtype=np.complex128), np.zeros(T, dtype=np.complex128),
           np.zeros(T, dtype=np.complex128), np.ones(T, dtype=np.complex128)]
    cents = points[tris].mean(1)
    metal = np.zeros(T, dtype=bool)
    for q in quads:
        metal |= MPath(q).contains_points(cents)
    for a in rho:
        a[metal] = 0.0
    # Local PML: imaginary stretch is inside assemble_anisotropic and keys off the
    # big device box. Our local coords are near x~10, y~0, outside that PML. Good.
    pec = np.zeros(len(points), dtype=bool)
    k0 = 2 * np.pi * float(sc.fs_a)
    # But assemble_anisotropic PML uses sc.nx_ports. Points with x<2 or x>28 get PML.
    # Horn port0 lives near x=10..14, y~0 — NOT in PML. Outer boundary of this mesh
    # is a hard truncation (Neumann). Keep monitors interior.
    try:
        A = assemble_anisotropic(points, tris, *rho, k0, pec)
    finally:
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = saved_box
    row_sum = np.abs(A).sum(axis=1).A1
    empty = row_sum < 1e-14
    if np.any(empty):
        A = A.tolil()
        for i in np.flatnonzero(empty):
            A.rows[i] = [i]
            A.data[i] = [1.0 + 0j]
        A = A.tocsr()
    # cosine source across the feed
    s = np.linspace(-span / 2, span / 2, 41)
    amp = np.clip(np.cos(np.pi * s / span), 0, None)
    amp = amp / np.sqrt(np.sum(amp ** 2))
    b = np.zeros(len(points), dtype=np.complex128)
    for si, ai in zip(s, amp):
        xy = src + si * tang
        i = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
        b[i] += ai
    t0 = time.perf_counter()
    lu = splu(A.tocsc())
    tf = time.perf_counter() - t0
    x = lu.solve(b)
    Ex, Ey = fields_from_hz(points, tris, x, *rho, k0)
    p_mon, samp_mon = line_flux(points, x, Ex, Ey, mon, tang, n_hat, span)
    # inward plane: 1.2 a inside the throat (toward array, -outward)
    inner = throat - 1.2 * n_hat
    p_in, samp_in = line_flux(points, x, Ex, Ey, inner, tang, n_hat, span)
    # outward of source
    outer = src + 1.2 * n_hat
    p_out, _ = line_flux(points, x, Ex, Ey, outer, tang, n_hat, span)
    return {
        "n_nodes": int(len(points)),
        "metal_tri_frac": float(metal.mean()),
        "factor_s": tf,
        "P_monitor_outward": float(p_mon),
        "P_inner_outward_normal": float(p_in),
        "P_outer_outward_normal": float(p_out),
        "monitor_samples": samp_mon[::4],
        "inner_center_a": inner.tolist(),
        "h_a": hmesh,
    }


def run_meep_one_horn(doc) -> dict:
    h = doc["horns"][0]
    quads = all_quads(doc, port=1)
    span = h["span_a"]
    src = np.asarray(h["source_center_a"], dtype=float)
    mon = np.asarray(h["monitor_center_a"], dtype=float)
    throat = np.asarray(h["throat_a"], dtype=float)
    n_hat = np.asarray(h["outward"], dtype=float)
    tang = np.asarray(h["tangent"], dtype=float)
    pts = np.vstack(quads)
    pad = 3.0
    lo = pts.min(0) - pad
    hi = pts.max(0) + pad
    center = 0.5 * (lo + hi)
    size = hi - lo
    geom = []
    for q in quads:
        verts = [mp.Vector3(float(v[0] - center[0]), float(v[1] - center[1]), 0) for v in q]
        geom.append(mp.Prism(vertices=verts, height=mp.inf, axis=mp.Vector3(0, 0, 1), material=mp.perfect_electric_conductor))

    def to_cell(xy):
        return mp.Vector3(float(xy[0] - center[0]), float(xy[1] - center[1]), 0)

    s = np.linspace(-span / 2, span / 2, 41)
    amp = np.clip(np.cos(np.pi * s / span), 0, None)
    amp = amp / np.sqrt(np.sum(amp ** 2))
    sources = []
    for si, ai in zip(s, amp):
        sources.append(
            mp.Source(
                mp.GaussianSource(frequency=float(sc.fs_a), fwidth=0.10 * float(sc.fs_a)),
                component=mp.Hz,
                center=to_cell(src + si * tang),
                amplitude=float(ai),
            )
        )
    sim = mp.Simulation(
        cell_size=mp.Vector3(float(size[0]), float(size[1]), 0),
        geometry=geom,
        sources=sources,
        resolution=40,
        boundary_layers=[mp.PML(0.8)],
    )

    def add_flux(xy):
        c = to_cell(xy)
        # X flux weighted by n_x, Y by n_y via two regions if needed.
        # Port 0 outward is +x, so a vertical flux line (size in y) of X-flux is correct.
        regions = []
        if abs(n_hat[0]) > 0.2:
            regions.append(mp.FluxRegion(center=c, size=mp.Vector3(0, span, 0), weight=float(n_hat[0])))
        if abs(n_hat[1]) > 0.2:
            regions.append(mp.FluxRegion(center=c, size=mp.Vector3(span, 0, 0), weight=float(n_hat[1])))
        return sim.add_flux(float(sc.fs_a), 0, 1, *regions)

    f_mon = add_flux(mon)
    inner = throat - 1.2 * n_hat
    f_in = add_flux(inner)
    outer = src + 1.2 * n_hat
    f_out = add_flux(outer)
    t0 = time.perf_counter()
    sim.run(until_after_sources=30)
    wall = time.perf_counter() - t0
    return {
        "P_monitor_outward": float(mp.get_fluxes(f_mon)[0]),
        "P_inner_outward_normal": float(mp.get_fluxes(f_in)[0]),
        "P_outer_outward_normal": float(mp.get_fluxes(f_out)[0]),
        "wall_s": wall,
        "res": 40,
        "inner_center_a": inner.tolist(),
    }


def db_ratio(a, b):
    if a == 0 or b == 0 or a * b < 0:
        return float("nan")
    return float(10 * np.log10(abs(a) / abs(b)))


def main() -> int:
    print("export polygons", flush=True)
    doc = export_polygons()
    quads = all_quads(doc)
    cmp = raster_polyline_vs_fill(quads, pitch=0.015)
    # filled vs filled (same polygons) must be 1
    ident = raster_masks(quads, quads, pitch=0.02)
    print("polyline vs fill", json.dumps(cmp, indent=2), flush=True)
    print("identical jaccard", ident["jaccard"], flush=True)
    print("FEM one horn...", flush=True)
    fem = run_fem_one_horn(doc)
    print({k: fem[k] for k in fem if k != "monitor_samples"}, flush=True)
    print("Meep one horn...", flush=True)
    meep = run_meep_one_horn(doc)
    print(meep, flush=True)

    def rel_db(fa, fb, ma, mb):
        # compare power ratios within each solver
        rf = fa / fb if fb else float("nan")
        rm = ma / mb if mb else float("nan")
        if rf <= 0 or rm <= 0:
            return {"fem_ratio": rf, "meep_ratio": rm, "delta_dB": float("nan")}
        return {"fem_ratio": float(rf), "meep_ratio": float(rm), "delta_dB": float(10 * np.log10(rf / rm))}

    out = {
        "polygons": str(OUT / "meep_horn_polygons.json"),
        "why_jaccard_0.32": cmp,
        "same_polygon_jaccard": ident["jaccard"],
        "fem": {k: v for k, v in fem.items() if k != "monitor_samples"},
        "meep": meep,
        "ratio_inner_over_monitor": rel_db(
            abs(fem["P_inner_outward_normal"]), abs(fem["P_monitor_outward"]),
            abs(meep["P_inner_outward_normal"]), abs(meep["P_monitor_outward"]),
        ),
        "sign_monitor": {"fem": fem["P_monitor_outward"], "meep": meep["P_monitor_outward"]},
    }
    (OUT / "one_horn_compare.json").write_text(json.dumps(out, indent=2) + "\n")
    (OUT / "mask_jaccard.json").write_text(json.dumps({"polyline_vs_fill": cmp, "identical": ident}, indent=2) + "\n")
    print(json.dumps(out["ratio_inner_over_monitor"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
