#!/usr/bin/env python3
"""
Locally refined FEM mesh for PlasMEEP six-port + 91-bulb geometry.

Uses gmsh + meshio. Element size ~0.2 mm near critical interfaces (FEM-H),
coarser elsewhere. Coordinates in Meep a-units (a=20 mm).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))
sys.path.insert(0, str(ROOT / "scripts"))

import sixport_common as sc  # noqa: E402
from sixport_common import (  # noqa: E402
    build_circulator_device,
    default_uniform_rho,
    set_geometry_context,
)

OUT = ROOT / "outputs" / "validation" / "highres_solver_campaign" / "phase5_fem"


# Physical: a = 20 mm = 0.02 m. 0.2 mm = 0.01 a. 1 mm = 0.05 a.
A_MM = sc.a * 1e3  # 20


def mm_to_a(mm: float) -> float:
    return mm / A_MM


MESH_GRADES = {
    # h_crit / h_air / h_pml in mm
    "FEM-L": dict(h_crit_mm=0.8, h_air_mm=4.0, h_pml_mm=2.0),
    "FEM-M": dict(h_crit_mm=0.4, h_air_mm=2.5, h_pml_mm=1.2),
    "FEM-H": dict(h_crit_mm=0.2, h_air_mm=2.0, h_pml_mm=0.8),
}


def build_mesh(grade: str, device_mode: str = "full") -> Dict[str, Any]:
    import gmsh

    cfg = MESH_GRADES[grade]
    h_crit = mm_to_a(cfg["h_crit_mm"])
    h_air = mm_to_a(cfg["h_air_mm"])
    h_pml = mm_to_a(cfg["h_pml_mm"])

    nx, ny = float(sc.nx_ports), float(sc.ny_ports)
    dpml = float(sc.dpml_ports)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"plasmeep_{grade}")

    # Domain rectangle [0,nx] x [0,ny]
    factory = gmsh.model.occ
    rect = factory.addRectangle(0, 0, 0, nx, ny)

    # Collect critical curves: bulb circles + approximate horn wall segments
    set_geometry_context(res=50, horn_walls="prism")
    rho = default_uniform_rho()
    B = np.zeros(3)
    pmm, _, _ = build_circulator_device(rho, B, res=50, device_mode=device_mode, wall_pec=True)

    crit_curves = []
    # Bulbs
    if device_mode == "full":
        locs = np.asarray(pmm.train_elem_locs, dtype=float)
        r_q = sc.r_bulb_outer
        r_p = 4.6 * sc.r_bulb_inner / 6.5
        for i in range(len(locs)):
            cx, cy = float(locs[i, 0]), float(locs[i, 1])
            # quartz outer + plasma core as size-control circles (embedded)
            c1 = factory.addCircle(cx, cy, 0, r_q)
            c2 = factory.addCircle(cx, cy, 0, r_p)
            crit_curves.extend([c1, c2])

    # Horn polylines as thin wire curves for size field
    for horn in sc.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], dtype=float)
            pts = []
            for p in poly:
                pts.append(factory.addPoint(p[0] + nx / 2, p[1] + ny / 2, 0, h_crit))
            for a, b in zip(pts[:-1], pts[1:]):
                crit_curves.append(factory.addLine(a, b))

    factory.synchronize()

    # Distance + threshold size fields
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", crit_curves)
    gmsh.model.mesh.field.setNumber(f_dist, "Sampling", 80)

    f_th = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_th, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMin", h_crit)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMax", h_air)
    gmsh.model.mesh.field.setNumber(f_th, "DistMin", 2 * h_crit)
    gmsh.model.mesh.field.setNumber(f_th, "DistMax", 0.15 * min(nx, ny))

    # PML box size
    f_box = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(f_box, "VIn", h_pml)
    gmsh.model.mesh.field.setNumber(f_box, "VOut", h_air)
    gmsh.model.mesh.field.setNumber(f_box, "XMin", 0)
    gmsh.model.mesh.field.setNumber(f_box, "XMax", nx)
    gmsh.model.mesh.field.setNumber(f_box, "YMin", 0)
    gmsh.model.mesh.field.setNumber(f_box, "YMax", ny)
    # Actually Box VIn is inside box — set thin PML frames via Min of threshold and box edges.
    # Simpler: use MathEval for PML strips
    f_pml = gmsh.model.mesh.field.add("MathEval")
    # default air size, smaller in PML bands
    expr = f"{h_air}"
    gmsh.model.mesh.field.setString(f_pml, "F", expr)

    f_min = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", [f_th])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_min)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay

    t0 = time.perf_counter()
    gmsh.model.mesh.generate(2)
    t_mesh = time.perf_counter() - t0

    # Stats
    nodes = gmsh.model.mesh.getNodes()
    n_nodes = len(nodes[0])
    etypes, elemTags, _ = gmsh.model.mesh.getElements(2)
    n_tri = 0
    for et, tags in zip(etypes, elemTags):
        if et == 2:  # triangle
            n_tri += len(tags)

    # Element size estimate from mesh
    coords = nodes[1].reshape(-1, 3)[:, :2]
    # sample edge lengths from triangles
    _, _, nodeTags = gmsh.model.mesh.getElements(2)
    # get triangles connectivity
    tris = None
    for et, tags, conn in zip(*gmsh.model.mesh.getElements(2)):
        if et == 2:
            tris = np.array(conn).reshape(-1, 3) - 1  # 0-based if tags are 1-based node nums
            # gmsh node tags may not be contiguous — map
            tag_to_idx = {int(t): i for i, t in enumerate(nodes[0])}
            tris = np.vectorize(lambda t: tag_to_idx[int(t)])(np.array(conn).reshape(-1, 3))
            break

    edge_lens = []
    if tris is not None and len(tris):
        for tri in tris[:: max(1, len(tris) // 5000)]:  # subsample
            p = coords[tri]
            for a, b in ((0, 1), (1, 2), (2, 0)):
                edge_lens.append(float(np.linalg.norm(p[a] - p[b])))
    edge_lens = np.array(edge_lens) if edge_lens else np.array([h_air])

    # Save mesh
    OUT.mkdir(parents=True, exist_ok=True)
    msh_path = OUT / f"mesh_{grade}_{device_mode}.msh"
    gmsh.write(str(msh_path))
    # also vtk via meshio if possible
    try:
        import meshio
        m = meshio.read(str(msh_path))
        meshio.write(str(OUT / f"mesh_{grade}_{device_mode}.xdmf"), m)
    except Exception as e:
        pass

    gmsh.finalize()

    # Air-region max: use 90th percentile of large edges as proxy
    stats = {
        "grade": grade,
        "device_mode": device_mode,
        "h_crit_mm": cfg["h_crit_mm"],
        "h_air_mm": cfg["h_air_mm"],
        "h_crit_a": h_crit,
        "h_air_a": h_air,
        "n_nodes": int(n_nodes),
        "n_triangles": int(n_tri),
        "scalar_unknowns_cg1": int(n_nodes),  # CG1 = nodes
        "min_edge_a": float(edge_lens.min()),
        "median_edge_a": float(np.median(edge_lens)),
        "max_edge_a": float(edge_lens.max()),
        "min_edge_mm": float(edge_lens.min() * A_MM),
        "median_edge_mm": float(np.median(edge_lens) * A_MM),
        "max_edge_mm": float(edge_lens.max() * A_MM),
        "mesh_time_s": t_mesh,
        "msh_path": str(msh_path),
        "vs_uniform_50ppc_8_4M": float(n_nodes) / 8.4e6,
        "memory_est_GiB_csr": float(n_nodes) * 15 * 16 / (1024**3),  # rough ~15 nnz/row complex
    }
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", choices=list(MESH_GRADES) + ["all"], default="all")
    ap.add_argument("--device-mode", default="full", choices=["full", "horns_only"])
    args = ap.parse_args()
    grades = list(MESH_GRADES) if args.grade == "all" else [args.grade]
    all_stats = []
    for g in grades:
        print(f"Building {g} ...", flush=True)
        try:
            st = build_mesh(g, device_mode=args.device_mode)
            print(json.dumps({k: st[k] for k in st if k != "msh_path"}, indent=2), flush=True)
            all_stats.append(st)
            (OUT / f"stats_{g}_{args.device_mode}.json").write_text(json.dumps(st, indent=2) + "\n")
        except Exception as e:
            err = {"grade": g, "error": str(e)}
            print(err, flush=True)
            all_stats.append(err)
    (OUT / f"mesh_stats_{args.device_mode}.json").write_text(json.dumps(all_stats, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
