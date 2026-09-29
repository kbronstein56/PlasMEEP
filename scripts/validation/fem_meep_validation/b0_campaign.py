#!/usr/bin/env python3
"""
B=0 validation campaign: excised PEC, Meep PML, mirror-symmetric meshes.

Horn-only gate is 0.10 dB and symmetry splitting <= 0.10 dB.
Later phases run only if that gate passes.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from fem_validated_solver import (  # noqa: E402
    LAST_OPERATOR,
    ElementSampler,
    apply_metal_walls_rho,
    assemble_anisotropic,
    guide_normal_flux,
    horn_prism_quads,
    inject_mode_consistent,
    load_numerical_mode,
    solve_system,
    wall_triangle_mask,
)
import sixport_common as sc  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "b0_campaign"
MEEP = np.array(
    [
        np.nan,
        0.03763160465218182,
        0.08553915041606872,
        0.696684297884412,
        0.08553915041606859,
        0.037631604652181856,
    ]
)
LEVELS = [
    ("FEM-M", 0.10, 0.040, 0.070, 0.08),
    ("FEM-H", 0.070, 0.025, 0.045, 0.05),
    ("FEM-H+", 0.050, 0.016, 0.030, 0.035),
    ("FEM-VH", 0.040, 0.012, 0.022, 0.028),
]


def _unique(pts: np.ndarray) -> np.ndarray:
    _, keep = np.unique(np.round(pts, 8), axis=0, return_index=True)
    return pts[np.sort(keep)]


def _seg(p0, p1, h):
    L = float(np.linalg.norm(p1 - p0))
    n = max(2, int(np.ceil(L / max(h, 1e-6))))
    t = np.linspace(0.0, 1.0, n)[:, None]
    return p0 * (1 - t) + p1 * t


def _load_triangle():
    vend = ROOT / "scripts" / "validation" / "fem_meep_validation" / "_vendor"
    if str(vend) not in sys.path:
        sys.path.insert(0, str(vend))
    import triangle as tr

    return tr


def clip_above(poly, y0: float) -> np.ndarray:
    """Sutherland–Hodgman clip of a polygon to y >= y0."""
    poly = np.asarray(poly, dtype=float)
    if np.all(poly[:, 1] >= y0 - 1e-9):
        return poly.copy()
    if np.all(poly[:, 1] <= y0 + 1e-9):
        return np.zeros((0, 2))
    out = []
    n = len(poly)
    for i in range(n):
        cur = poly[i]
        prev = poly[(i - 1) % n]
        cur_in = cur[1] >= y0 - 1e-12
        prev_in = prev[1] >= y0 - 1e-12
        if cur_in != prev_in:
            t = (y0 - prev[1]) / (cur[1] - prev[1])
            out.append(prev + t * (cur - prev))
        if cur_in:
            out.append(cur.copy())
    if not out:
        return np.zeros((0, 2))
    return np.asarray(out, dtype=float)


def _on_outer_edge(a, b, nx, ny, y0) -> bool:
    def on_side(p):
        return (
            abs(p[0]) <= 1e-8
            or abs(p[0] - nx) <= 1e-8
            or abs(p[1] - y0) <= 1e-8
            or abs(p[1] - ny) <= 1e-8
        )

    if not (on_side(a) and on_side(b)):
        return False
    return abs(a[0] - b[0]) <= 1e-8 or abs(a[1] - b[1]) <= 1e-8


def build_symmetric_mesh(h_bulk, h_wall, h_center, h_pml, h_quartz=None, h_plasma=None):
    """
    Upper-half constrained triangulation, mirrored across the P1–P4 axis.

    Horn prism edges are constrained segments, so the PEC boundary is the
    polygon itself. The mirror makes P2/P6 and P3/P5 meshes identical.
    """
    tr = _load_triangle()
    nx, ny = float(sc.nx_ports), float(sc.ny_ports)
    y0 = ny / 2.0
    verts: list = []
    index: dict = {}
    segs: list = []
    segset: set = set()

    def add_vertex(p):
        x = float(p[0])
        y = float(p[1])
        if x < -1e-8 or x > nx + 1e-8 or y < y0 - 1e-8 or y > ny + 1e-8:
            return None
        if abs(y - y0) <= 1e-8:
            y = y0
        if abs(y - ny) <= 1e-8:
            y = ny
        if abs(x) <= 1e-8:
            x = 0.0
        if abs(x - nx) <= 1e-8:
            x = nx
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

    def add_interior_chain(samples):
        ids = [add_vertex(p) for p in samples]
        clean = []
        for i in ids:
            if i is None:
                continue
            if not clean or clean[-1] != i:
                clean.append(i)
        for i in range(len(clean) - 1):
            add_seg(clean[i], clean[i + 1])

    for quad in horn_prism_quads(nx, ny):
        clipped = clip_above(quad, y0)
        if len(clipped) < 3:
            continue
        for i in range(len(clipped)):
            a = clipped[i]
            b = clipped[(i + 1) % len(clipped)]
            if _on_outer_edge(a, b, nx, ny, y0):
                add_vertex(a)
                add_vertex(b)
                continue
            add_interior_chain(_seg(np.asarray(a, float), np.asarray(b, float), h_wall))

    quartz_centers = None
    if h_quartz is not None:
        from fem_validated_solver import bulb_centers_mesh

        quartz_centers = bulb_centers_mesh()
        r_in = float(sc.r_bulb_inner)
        r_out = float(sc.r_bulb_outer)
        n_span = max(3, int(np.ceil((r_out - r_in) / float(h_quartz))))
        radii = np.linspace(r_in, r_out, n_span + 1)
        for cx, cy in quartz_centers:
            if cy + r_out < y0 - 1e-8:
                continue
            for radius in radii:
                if cy - radius >= y0 - 1e-9:
                    n = max(16, int(np.ceil(2.0 * np.pi * radius / float(h_quartz))))
                    ang = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
                    pts_c = np.column_stack([cx + radius * np.cos(ang), cy + radius * np.sin(ang)])
                    add_interior_chain(np.vstack([pts_c, pts_c[:1]]))
                else:
                    s = float(np.clip((y0 - cy) / radius, -1.0, 1.0))
                    a1 = float(np.arcsin(s))
                    a2 = float(np.pi - a1)
                    n = max(8, int(np.ceil(radius * abs(a2 - a1) / float(h_quartz))))
                    ang = np.linspace(a1, a2, n)
                    pts_c = np.column_stack([cx + radius * np.cos(ang), cy + radius * np.sin(ang)])
                    add_interior_chain(pts_c)
        if h_plasma is not None:
            r_p = 4.6 * float(sc.r_bulb_inner) / 6.5
            for cx, cy in quartz_centers:
                if cy + r_p < y0 - 1e-8:
                    continue
                if cy - r_p >= y0 - 1e-9:
                    n = max(16, int(np.ceil(2.0 * np.pi * r_p / float(h_plasma))))
                    ang = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
                    pts_c = np.column_stack([cx + r_p * np.cos(ang), cy + r_p * np.sin(ang)])
                    add_interior_chain(np.vstack([pts_c, pts_c[:1]]))
                else:
                    s = float(np.clip((y0 - cy) / r_p, -1.0, 1.0))
                    a1 = float(np.arcsin(s))
                    a2 = float(np.pi - a1)
                    n = max(8, int(np.ceil(r_p * abs(a2 - a1) / float(h_plasma))))
                    ang = np.linspace(a1, a2, n)
                    pts_c = np.column_stack([cx + r_p * np.cos(ang), cy + r_p * np.sin(ang)])
                    add_interior_chain(pts_c)

    from sixport_common import effective_port_dir, horn_for_port, monitor_center_for_port

    doc_shift = np.array([nx / 2.0, ny / 2.0])
    free = []
    xs = np.arange(0.0, nx + 0.5 * h_bulk, h_bulk)
    ys = np.arange(y0, ny + 0.5 * h_bulk, h_bulk)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    free.append(np.column_stack([xx.ravel(), yy.ravel()]))
    xs = np.arange(nx / 2.0 - 8.0, nx / 2.0 + 8.0 + h_center, h_center)
    ys = np.arange(y0, ny / 2.0 + 8.0 + h_center, h_center)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    free.append(np.column_stack([xx.ravel(), yy.ravel()]))
    dp = float(sc.dpml_ports)
    xb = np.arange(0.0, nx + h_pml, h_pml)
    yb = np.arange(y0, ny + h_pml, h_pml)
    free.append(np.column_stack([np.full_like(yb, dp), yb]))
    free.append(np.column_stack([np.full_like(yb, nx - dp), yb]))
    free.append(np.column_stack([xb, np.full_like(xb, ny - dp)]))
    free.append(np.column_stack([xb, np.full_like(xb, y0)]))
    for p in range(6):
        horn = horn_for_port(p, 50)
        u = np.asarray(effective_port_dir(p), dtype=float)
        u = u / np.linalg.norm(u)
        tang = np.array([-u[1], u[0]])
        span = 0.96 * float(sc.clear_width)
        for key in ("source_center", "throat_center"):
            c = np.asarray(horn[key], dtype=float)[:2] + doc_shift
            free.append(_seg(c - 0.5 * span * tang, c + 0.5 * span * tang, h_wall))
        c = np.asarray(monitor_center_for_port(p, 50), dtype=float) + doc_shift
        free.append(_seg(c - 0.5 * span * tang, c + 0.5 * span * tang, h_wall))
    free_pts = np.vstack(free)
    if quartz_centers is not None:
        r_in = float(sc.r_bulb_inner)
        r_out = float(sc.r_bulb_outer)
        band = 0.35 * float(h_quartz)
        keep = np.ones(len(free_pts), dtype=bool)
        for cx, cy in quartz_centers:
            dist = np.hypot(free_pts[:, 0] - cx, free_pts[:, 1] - cy)
            keep &= ~((dist > r_in - band) & (dist < r_out + band))
        free_pts = free_pts[keep]
    for p in free_pts:
        add_vertex(p)
    # Closed outer boundary. arange() does not reliably land on nx/ny.
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

    amax = 0.45 * float(h_bulk) * float(h_bulk)
    mesh = tr.triangulate(
        {"vertices": np.asarray(verts, dtype=float), "segments": np.asarray(segs, dtype=np.int32)},
        f"pq20a{amax:.8e}",
    )
    pts = np.asarray(mesh["vertices"], dtype=float)
    upper = np.asarray(mesh["triangles"], dtype=int)
    xu = pts[upper, 0]
    yu = pts[upper, 1]
    twice_u = (xu[:, 1] - xu[:, 0]) * (yu[:, 2] - yu[:, 0]) - (xu[:, 2] - xu[:, 0]) * (yu[:, 1] - yu[:, 0])
    upper_area = 0.5 * np.abs(twice_u).sum()
    if abs(upper_area - 0.5 * nx * ny) / (0.5 * nx * ny) > 0.02:
        OUT.mkdir(parents=True, exist_ok=True)
        np.savez(
            OUT / "pslg_debug.npz",
            vertices=np.asarray(verts, dtype=float),
            segments=np.asarray(segs, dtype=np.int32),
            tri_vertices=pts,
            triangles=upper,
        )
        raise RuntimeError(
            f"upper mesh area {upper_area:.4f} != half-domain {0.5*nx*ny:.4f} "
            f"(nvert {len(pts)} ntri {len(upper)} nseg {len(segs)})"
        )
    pts = np.asarray(mesh["vertices"], dtype=float)
    upper = np.asarray(mesh["triangles"], dtype=int)
    if np.any(pts[:, 1] < y0 - 1e-6):
        raise RuntimeError("constrained mesh left the upper half-domain")
    on_axis = pts[:, 1] <= y0 + 1e-8
    below_idx = np.flatnonzero(~on_axis)
    mirrored_pts = np.column_stack([pts[below_idx, 0], 2.0 * y0 - pts[below_idx, 1]])
    full = np.vstack([pts, mirrored_pts])
    mirror_of = np.arange(len(pts))
    mirror_of[below_idx] = np.arange(len(pts), len(pts) + len(below_idx))
    tris = np.vstack([upper, mirror_of[upper]])
    x = full[tris, 0]
    y = full[tris, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    tris = tris[np.abs(twice) > 1e-14]
    area = 0.5 * np.abs(twice[np.abs(twice) > 1e-14]).sum()
    expect = nx * ny
    if abs(area - expect) / expect > 0.02:
        raise RuntimeError(f"mirrored mesh area {area:.4f} != domain {expect:.4f}")
    return full, tris


def edge_stats(points, tris):
    lens = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        lens.append(np.linalg.norm(points[tris[:, a]] - points[tris[:, b]], axis=1))
    e = np.concatenate(lens)
    return float(e.min()), float(np.median(e)), float(e.max())


def feed_rho(points, tris):
    T = len(tris)
    rho = [np.ones(T, dtype=np.complex128), np.zeros(T, dtype=np.complex128), np.zeros(T, dtype=np.complex128), np.ones(T, dtype=np.complex128)]
    cents = points[tris].mean(1)
    half = 0.5 * float(sc.clear_width)
    th = float(sc.wall_thickness)
    c = cents - np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    wall = (np.abs(c[:, 1]) >= half) & (np.abs(c[:, 1]) <= half + th)
    for a in rho:
        a[wall] = 0.0
    return rho


def horn_rho(points, tris):
    T = len(tris)
    rho = [np.ones(T, dtype=np.complex128), np.zeros(T, dtype=np.complex128), np.zeros(T, dtype=np.complex128), np.ones(T, dtype=np.complex128)]
    mask = wall_triangle_mask(points, tris)
    return list(apply_metal_walls_rho(*rho, mask)), mask


def run_level(name, h_bulk, h_wall, h_center, h_pml, src):
    print(f"\n=== {name} h_bulk={h_bulk} h_wall={h_wall} ===", flush=True)
    t_mesh0 = time.perf_counter()
    points, tris = build_symmetric_mesh(h_bulk, h_wall, h_center, h_pml)
    t_mesh = time.perf_counter() - t_mesh0
    emin, emed, emax = edge_stats(points, tris)
    k0 = 2 * np.pi * float(sc.fs_a)
    pec = np.zeros(len(points), dtype=bool)
    rho, mask = horn_rho(points, tris)
    sampler = ElementSampler(points, tris)
    t0 = time.perf_counter()
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    t_asm = time.perf_counter() - t0
    op = dict(LAST_OPERATOR)
    print("OPERATOR", json.dumps(op), flush=True)
    if op["pec_elements"] < 100:
        raise RuntimeError(f"PEC excision did not see the prisms: {op}")
    if op["sigma_max"] < 20.0:
        raise RuntimeError(f"PML fell back to a weak profile: {op}")
    b = inject_mode_consistent(points, tris, src, sampler)
    x, _lu, t_fac, t_sol, resid = solve_system(A, b, pec)
    raw = np.array(
        [
            guide_normal_flux(
                points, x, None, None, p, sampler=sampler, rho=rho, omega=k0
            )
            for p in range(6)
        ]
    )
    # Incident straight feed, same mesh and same excised-PML operator.
    rho_f = feed_rho(points, tris)
    Aref = assemble_anisotropic(points, tris, *rho_f, k0, pec)
    xref, _lu2, t_fac_ref, _ts, _r = solve_system(Aref, b, pec)
    pinc = abs(
        guide_normal_flux(points, xref, None, None, 0, sampler=sampler, rho=rho_f, omega=k0)
    )
    if sampler.misses:
        raise RuntimeError(f"line samples missed the mesh: {sampler.misses}")
    norm = raw / pinc
    dB = [float(10 * np.log10(v)) if v > 0 else None for v in norm]
    delta = [None if not (np.isfinite(MEEP[i]) and norm[i] > 0) else float(10 * np.log10(norm[i] / MEEP[i])) for i in range(6)]
    rec = {
        "level": name,
        "h_bulk": h_bulk,
        "h_wall": h_wall,
        "h_center": h_center,
        "n_nodes": int(len(points)),
        "n_tris": int(len(tris)),
        "edge_min_a": emin,
        "edge_median_a": emed,
        "edge_max_a": emax,
        "mesh_s": t_mesh,
        "assemble_s": t_asm,
        "factor_s": t_fac,
        "solve_s": t_sol,
        "incident_factor_s": t_fac_ref,
        "residual": resid,
        "operator": op,
        "metal_tri_frac": float(mask.mean()),
        "sampler_misses": int(sampler.misses),
        "mesh": "constrained_prism_edges_mirrored",
        "flux": "element_barycentric",
        "P_inc": float(pinc),
        "normalized": norm.tolist(),
        "dB": dB,
        "delta_dB": delta,
        "P2_minus_P6_dB": float(10 * np.log10(norm[1] / norm[5])) if norm[1] > 0 and norm[5] > 0 else None,
        "P3_minus_P5_dB": float(10 * np.log10(norm[2] / norm[4])) if norm[2] > 0 and norm[4] > 0 else None,
    }
    print(
        name,
        "DOFs", rec["n_nodes"],
        "sym", rec["P2_minus_P6_dB"], rec["P3_minus_P5_dB"],
        "dB", [None if v is None else round(v, 3) for v in delta],
        flush=True,
    )
    return rec


def plot_convergence(rows):
    dofs = [r["n_nodes"] for r in rows]
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
    for p, label in ((1, "P2"), (2, "P3"), (3, "P4"), (4, "P5"), (5, "P6")):
        ax[0].plot(dofs, [r["normalized"][p] for r in rows], marker="o", label=label)
    ax[0].set_xlabel("DOFs")
    ax[0].set_ylabel("normalized power")
    ax[0].legend(fontsize=7)
    for p, label in ((1, "P2"), (2, "P3"), (3, "P4"), (4, "P5"), (5, "P6")):
        ax[1].plot(dofs, [r["delta_dB"][p] for r in rows], marker="o", label=label)
    ax[1].axhline(0.1, color="k", lw=0.6)
    ax[1].axhline(-0.1, color="k", lw=0.6)
    ax[1].set_xlabel("DOFs")
    ax[1].set_ylabel("FEM − Meep (dB)")
    ax[1].legend(fontsize=7)
    ax[2].plot(dofs, [abs(r["P2_minus_P6_dB"]) for r in rows], marker="o", label="|P2−P6|")
    ax[2].plot(dofs, [abs(r["P3_minus_P5_dB"]) for r in rows], marker="o", label="|P3−P5|")
    ax[2].axhline(0.10, color="k", lw=0.6)
    ax[2].set_xlabel("DOFs")
    ax[2].set_ylabel("symmetry split (dB)")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "horn_convergence.png", dpi=140)
    plt.close(fig)


def gate(rows):
    last = rows[-1]
    errs = [abs(v) for v in last["delta_dB"][1:] if v is not None]
    sym = max(abs(last["P2_minus_P6_dB"]), abs(last["P3_minus_P5_dB"]))
    return max(errs) <= 0.10 and sym <= 0.10, max(errs), sym


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    src = load_numerical_mode()
    rows = []
    stop_at = None
    for spec in LEVELS:
        rec = run_level(*spec, src)
        rows.append(rec)
        (OUT / "horn_convergence.json").write_text(json.dumps(rows, indent=2) + "\n")
        errs = [abs(v) for v in rec["delta_dB"][1:] if v is not None]
        sym = max(abs(rec["P2_minus_P6_dB"]), abs(rec["P3_minus_P5_dB"]))
        if max(errs) > 0.25 and rec["level"] in ("FEM-H", "FEM-H+", "FEM-VH"):
            # Stay on horns, but still finish the ladder unless the error is growing badly.
            print("still above 0.25 dB", rec["level"], max(errs), flush=True)
        if len(rows) >= 2:
            # Stop the ladder only if two successive fine meshes are both outside 0.25
            # and not improving. Otherwise continue.
            prev = max(abs(v) for v in rows[-2]["delta_dB"][1:] if v is not None)
            if rec["level"] == "FEM-VH" or (max(errs) <= 0.10 and sym <= 0.10 and rec["level"] != "FEM-M"):
                if max(errs) <= 0.10 and sym <= 0.10:
                    print("GATE PASS", rec["level"], flush=True)
                    break
    plot_convergence(rows)
    ok, worst, sym = gate(rows)
    summary = {"horn_gate_pass": ok, "worst_dB": worst, "symmetry_dB": sym, "levels": [r["level"] for r in rows]}
    (OUT / "horn_gate.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("GATE", summary, flush=True)
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
