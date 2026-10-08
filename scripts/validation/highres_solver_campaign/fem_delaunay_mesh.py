#!/usr/bin/env python3
"""
Lightweight locally-refined triangular mesh (no gmsh, no Meep device build).

Uses hex lattice of Plasma Mirror Module bulb centers from PlasMEEP geometry
constants + horn polylines from sixport_common.full_horns.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.spatial import Delaunay, cKDTree

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))
sys.path.insert(0, str(ROOT / "scripts"))

import sixport_common as sc  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "highres_solver_campaign" / "phase5_fem"
A_MM = sc.a * 1e3

MESH_GRADES = {
    "FEM-L": dict(h_crit_mm=0.8, h_air_mm=6.0),
    "FEM-M": dict(h_crit_mm=0.4, h_air_mm=3.5),
    "FEM-H": dict(h_crit_mm=0.2, h_air_mm=2.5),
}


def mm_to_a(mm: float) -> float:
    return mm / A_MM


def bulb_centers_cell() -> np.ndarray:
    """91 bulb centers in full-domain cell coords."""
    from PMMCirculatorInverse import PMMI

    # Match sixport_common layout: hex train, then shift into device frame
    pmm_tmp = PMMI(
        a=sc.a,
        res=32,
        nx=14,
        ny=12,
        dpml=sc.dpml,
        B=np.array([0.0, 0.0, 0.0]),
    )
    pmm_tmp.Rod_Array_Hexagon_train(
        xy_cen=np.array([14 / 2, 12 / 2]),
        side_dim=6,
        r=sc.r_plasma,
        d=sc.d_exp,
        bulbs=False,
        uniform=True,
    )
    locs = np.asarray(pmm_tmp.train_elem_locs, dtype=float)
    # locs are in the 14x12 frame; shift so array center -> (nx_ports/2, ny_ports/2)
    locs = locs - np.array([14 / 2, 12 / 2]) + np.array([sc.nx_ports / 2, sc.ny_ports / 2])
    assert len(locs) == 91, len(locs)
    return locs


def horn_polylines_cell() -> List[np.ndarray]:
    nx, ny = sc.nx_ports, sc.ny_ports
    out = []
    for horn in sc.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], dtype=float) + np.array([nx / 2, ny / 2])
            out.append(poly)
    return out


def build_points(grade: str, device_mode: str = "full") -> Tuple[np.ndarray, Dict[str, Any]]:
    cfg = MESH_GRADES[grade]
    h_c = mm_to_a(cfg["h_crit_mm"])
    h_a = mm_to_a(cfg["h_air_mm"])
    nx, ny = float(sc.nx_ports), float(sc.ny_ports)

    pts: List[List[float]] = []

    # Domain boundary
    nb = max(20, int(max(nx, ny) / h_a))
    for x in np.linspace(0, nx, nb):
        pts.append([float(x), 0.0]); pts.append([float(x), ny])
    for y in np.linspace(0, ny, nb):
        pts.append([0.0, float(y)]); pts.append([nx, float(y)])

    # Background lattice at air scale
    xs = np.arange(0.5 * h_a, nx, h_a)
    ys = np.arange(0.5 * h_a, ny, h_a)
    X, Y = np.meshgrid(xs, ys)
    pts.extend(np.column_stack([X.ravel(), Y.ravel()]).tolist())

    crit: List[List[float]] = []
    if device_mode == "full":
        bulbs = bulb_centers_cell()
        r_q = sc.r_bulb_outer
        r_p = 4.6 * sc.r_bulb_inner / 6.5
        for cx, cy in bulbs:
            crit.append([float(cx), float(cy)])
            for r in (r_p * 0.5, r_p, 0.5 * (r_p + r_q), r_q, r_q + 1.5 * h_c):
                n = max(16, int(2 * np.pi * r / h_c))
                th = np.linspace(0, 2 * np.pi, n, endpoint=False)
                ring = np.column_stack([cx + r * np.cos(th), cy + r * np.sin(th)])
                pts.extend(ring.tolist())
                crit.extend(ring.tolist())
            # local fine lattice in bulb neighborhood
            rad = r_q + 4 * h_c
            xs2 = np.arange(cx - rad, cx + rad + 1e-12, h_c)
            ys2 = np.arange(cy - rad, cy + rad + 1e-12, h_c)
            for x in xs2:
                for y in ys2:
                    if (x - cx) ** 2 + (y - cy) ** 2 <= rad ** 2:
                        if 0 <= x <= nx and 0 <= y <= ny:
                            pts.append([float(x), float(y)])

    # Horn walls
    for poly in horn_polylines_cell():
        for k in range(len(poly) - 1):
            p0, p1 = poly[k], poly[k + 1]
            L = float(np.linalg.norm(p1 - p0))
            n = max(2, int(L / h_c))
            for t in np.linspace(0, 1, n):
                p = p0 * (1 - t) + p1 * t
                pts.append([float(p[0]), float(p[1])])
                crit.append([float(p[0]), float(p[1])])
            # band around segment
            if L < 1e-12:
                continue
            tang = (p1 - p0) / L
            nrm = np.array([-tang[1], tang[0]])
            for s in np.linspace(0, L, max(2, int(L / (1.5 * h_c)))):
                base = p0 + s * tang
                for w in (-2 * h_c, -h_c, h_c, 2 * h_c):
                    q = base + w * nrm
                    if 0 <= q[0] <= nx and 0 <= q[1] <= ny:
                        pts.append([float(q[0]), float(q[1])])

    P = np.asarray(pts, dtype=float)
    P = P[(P[:, 0] >= 0) & (P[:, 0] <= nx) & (P[:, 1] >= 0) & (P[:, 1] <= ny)]
    # Round-dedup
    P = np.unique(np.round(P, decimals=6), axis=0)
    meta = {
        "h_crit_mm": cfg["h_crit_mm"],
        "h_air_mm": cfg["h_air_mm"],
        "h_crit_a": h_c,
        "h_air_a": h_a,
        "n_crit_samples": len(crit),
        "n_points_raw": int(len(P)),
    }
    return P, meta


def build(grade: str, device_mode: str = "full") -> Dict[str, Any]:
    t0 = time.perf_counter()
    print(f"  sampling {grade}...", flush=True)
    P, meta = build_points(grade, device_mode)
    print(f"  points={len(P)} sample_s={time.perf_counter()-t0:.2f}", flush=True)
    t1 = time.perf_counter()
    tri = Delaunay(P)
    simplices = tri.simplices
    # Drop triangles with very large circumradius (outside domain artifacts)
    keep = []
    for i, t in enumerate(simplices):
        p = P[t]
        # area
        area = 0.5 * abs(
            p[0, 0] * (p[1, 1] - p[2, 1])
            + p[1, 0] * (p[2, 1] - p[0, 1])
            + p[2, 0] * (p[0, 1] - p[1, 1])
        )
        if area < 1e-14:
            continue
        # edge max
        e = [
            np.linalg.norm(p[0] - p[1]),
            np.linalg.norm(p[1] - p[2]),
            np.linalg.norm(p[2] - p[0]),
        ]
        if max(e) > 3.5 * meta["h_air_a"]:
            continue
        keep.append(i)
    simplices = simplices[keep]
    print(f"  tris={len(simplices)} tri_s={time.perf_counter()-t1:.2f}", flush=True)

    # edge stats subsample
    edges = []
    step = max(1, len(simplices) // 5000)
    for t in simplices[::step]:
        p = P[t]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            edges.append(np.linalg.norm(p[a] - p[b]))
    e = np.asarray(edges)

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"mesh_{grade}_{device_mode}.npz", points=P, triangles=simplices)
    stats = {
        "grade": grade,
        "device_mode": device_mode,
        "generator": "Delaunay_light",
        "n_nodes": int(len(P)),
        "n_triangles": int(len(simplices)),
        "scalar_unknowns_cg1": int(len(P)),
        "min_edge_mm": float(e.min() * A_MM),
        "median_edge_mm": float(np.median(e) * A_MM),
        "max_edge_mm": float(e.max() * A_MM),
        "vs_uniform_50ppc_8_4M": float(len(P)) / 8.4e6,
        "memory_est_GiB_csr": float(len(P)) * 12 * 16 / (1024**3),
        "wall_s": time.perf_counter() - t0,
        **meta,
    }
    return stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="all")
    ap.add_argument("--device-mode", default="full")
    args = ap.parse_args()
    grades = list(MESH_GRADES) if args.grade == "all" else [args.grade]
    all_stats = []
    for g in grades:
        print(f"Building {g}", flush=True)
        st = build(g, args.device_mode)
        print(json.dumps(st, indent=2), flush=True)
        (OUT / f"stats_{g}_{args.device_mode}.json").write_text(json.dumps(st, indent=2) + "\n")
        all_stats.append(st)
    (OUT / f"mesh_stats_{args.device_mode}.json").write_text(json.dumps(all_stats, indent=2) + "\n")
    # key question
    for st in all_stats:
        print(
            f"KEY: {st['grade']} unknowns={st['n_nodes']} "
            f"({100*st['vs_uniform_50ppc_8_4M']:.2f}% of 8.4M) "
            f"min_edge={st['min_edge_mm']:.3f} mm",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
