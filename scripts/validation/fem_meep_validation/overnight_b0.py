#!/usr/bin/env python3
"""Quartz-only and B=0-plasma FEM on the canonical horn operator.

PEC is excised (no stiffness, no mass). PML is the Meep quadratic profile.
The mesh is the mirrored constrained triangulation, with concentric rings
through each 1 mm quartz wall.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from b0_campaign import build_symmetric_mesh, feed_rho, horn_rho  # noqa: E402
from fem_validated_solver import (  # noqa: E402
    LAST_OPERATOR,
    ElementSampler,
    apply_metal_walls_rho,
    assemble_anisotropic,
    assign_rho,
    eps_tensor_at_bias,
    guide_normal_flux,
    inject_mode_consistent,
    load_numerical_mode,
    meep_sigma_max,
    solve_system,
    wall_triangle_mask,
)
from sixport_common import set_geometry_context  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "overnight_b0"
LEVELS = [
    ("FEM-H", 0.070, 0.025, 0.045, 0.050, 0.016),
    ("FEM-VH", 0.050, 0.016, 0.030, 0.035, 0.010),
]


def _symmetrize_rho(points, tris, rho) -> None:
    """Copy upper-half materials onto exact mirror triangles."""
    from scipy.spatial import cKDTree

    y0 = float(sc.ny_ports) / 2.0
    cents = points[tris].mean(1)
    upper = np.flatnonzero(cents[:, 1] >= y0 - 1e-9)
    lower = np.flatnonzero(cents[:, 1] < y0 - 1e-9)
    if len(lower) == 0:
        return
    tree = cKDTree(cents[upper])
    query = np.column_stack([cents[lower, 0], 2.0 * y0 - cents[lower, 1]])
    dist, partner = tree.query(query, k=1)
    if float(dist.max()) > 1e-6:
        raise RuntimeError(f"mesh is not a mirror: max centroid mismatch {dist.max()}")
    src = upper[partner]
    for comp in rho:
        comp[lower] = comp[src]


def quartz_mask(points, tris, centers):
    cents = points[tris].mean(1)
    r_in = float(sc.r_bulb_inner)
    r_out = float(sc.r_bulb_outer)
    mask = np.zeros(len(tris), dtype=bool)
    for cx, cy in centers:
        d2 = (cents[:, 0] - cx) ** 2 + (cents[:, 1] - cy) ** 2
        mask |= (d2 <= r_out**2) & (d2 > r_in**2)
    return mask


def material_rho(points, tris, kind: str):
    """kind is quartz (air interiors) or plasma (uniform fp inside r_plasma)."""
    from fem_validated_solver import bulb_centers_mesh

    centers = bulb_centers_mesh()
    if kind == "quartz":
        T = len(tris)
        rho = [
            np.ones(T, dtype=np.complex128),
            np.zeros(T, dtype=np.complex128),
            np.zeros(T, dtype=np.complex128),
            np.ones(T, dtype=np.complex128),
        ]
        q = quartz_mask(points, tris, centers)
        rq = 1.0 / (3.8 + 0j)
        for a in rho:
            pass
        rho[0][q] = rq
        rho[3][q] = rq
        rho[1][q] = 0.0
        rho[2][q] = 0.0
    elif kind == "plasma":
        rho = list(assign_rho(points, tris, bias_a=0.0))
        q = quartz_mask(points, tris, centers)
    else:
        raise ValueError(kind)
    metal = wall_triangle_mask(points, tris)
    rho = list(apply_metal_walls_rho(*rho, metal))
    _symmetrize_rho(points, tris, rho)
    return rho, q, metal, centers


def wall_resolution(points, tris, sampler, centers):
    r_in = float(sc.r_bulb_inner)
    r_out = float(sc.r_bulb_outer)
    spans = []
    # A few bulbs: on-axis if present, and the first several upper-half centers.
    y0 = sc.ny_ports / 2.0
    chosen = [c for c in centers if c[1] >= y0 - 1e-8][:8]
    for cx, cy in chosen:
        for ang in np.linspace(0.15, np.pi - 0.15, 4):
            rr = np.linspace(r_in + 1e-4, r_out - 1e-4, 48)
            xy = np.column_stack([cx + rr * np.cos(ang), cy + rr * np.sin(ang)])
            xy = xy[xy[:, 1] >= y0 - 1e-6]
            if len(xy) < 4:
                continue
            ids = sampler.locate(xy)
            spans.append(int(len(set(int(i) for i in ids if i >= 0))))
    q = quartz_mask(points, tris, centers)
    if not np.any(q):
        return {"n_quartz_tris": 0}
    e = []
    qt = tris[q]
    for a, b in ((0, 1), (1, 2), (2, 0)):
        e.append(np.linalg.norm(points[qt[:, a]] - points[qt[:, b]], axis=1))
    edges = np.concatenate(e)
    return {
        "wall_thickness_a": float(r_out - r_in),
        "wall_thickness_mm": float((r_out - r_in) * sc.a * 1e3),
        "n_quartz_tris": int(q.sum()),
        "edge_min_a": float(edges.min()),
        "edge_median_a": float(np.median(edges)),
        "edge_max_a": float(edges.max()),
        "min_elements_across_wall": int(min(spans)) if spans else None,
        "median_elements_across_wall": float(np.median(spans)) if spans else None,
    }


def solve_level(name, spec, kind, src):
    h_bulk, h_wall, h_center, h_pml, h_q = spec
    print(f"\n=== {kind} {name} h_q={h_q} ===", flush=True)
    t0 = time.perf_counter()
    points, tris = build_symmetric_mesh(
        h_bulk,
        h_wall,
        h_center,
        h_pml,
        h_quartz=h_q,
        h_plasma=(h_q if kind == "plasma" else None),
    )
    t_mesh = time.perf_counter() - t0
    k0 = 2 * np.pi * float(sc.fs_a)
    pec = np.zeros(len(points), dtype=bool)
    t0 = time.perf_counter()
    rho, qmask, metal, centers = material_rho(points, tris, kind)
    t_mat = time.perf_counter() - t0
    sampler = ElementSampler(points, tris)
    qstat = wall_resolution(points, tris, sampler, centers)
    print("QUARTZ", json.dumps(qstat), flush=True)
    if qstat.get("min_elements_across_wall", 0) is None or qstat["min_elements_across_wall"] < 2:
        raise RuntimeError(f"quartz wall is not spanned by at least 2 elements: {qstat}")
    t0 = time.perf_counter()
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    t_asm = time.perf_counter() - t0
    op = dict(LAST_OPERATOR)
    print("OPERATOR", json.dumps(op), flush=True)
    if op.get("pec_model") != "excised_neumann_no_mass":
        raise RuntimeError(f"PEC formulation fell back: {op}")
    if op.get("pml_profile") != "meep_quadratic_R1e-15" or op["sigma_max"] < 20.0:
        raise RuntimeError(f"PML formulation fell back: {op}")
    if abs(op["sigma_max"] - meep_sigma_max(sc.dpml_ports)) > 1e-9:
        raise RuntimeError("sigma_max does not match the Meep formula")
    if op["pec_elements"] < 100 or int(metal.sum()) < 100:
        raise RuntimeError("PEC excision missed the prisms")
    eps = eps_tensor_at_bias(0.0)
    print(
        "EPS",
        json.dumps(
            {
                "fs_a": float(sc.fs_a),
                "fp_a": float(sc.fp_a),
                "gamma_a": float(sc.gamma_a),
                "eps_xx": [eps[0].real, eps[0].imag],
                "eps_xy": [eps[1].real, eps[1].imag],
                "quartz_eps": 3.8,
                "kind": kind,
            }
        ),
        flush=True,
    )
    if kind == "plasma":
        import meep as mp

        med = mp.Medium(
            epsilon=1.0,
            E_susceptibilities=[
                mp.DrudeSusceptibility(frequency=float(sc.fp_a), gamma=float(sc.gamma_a), sigma=1.0)
            ],
        )
        eps_meep = complex(med.epsilon(float(sc.fs_a))[0, 0])
        diff = complex(eps[0]) - eps_meep
        print(
            "PLASMA_EPS",
            json.dumps(
                {
                    "fs_Hz": float(sc.fs_Hz),
                    "fp_Hz": float(sc.fp_Hz),
                    "gamma_Hz": float(sc.gamma_Hz),
                    "fs_a": float(sc.fs_a),
                    "fp_a": float(sc.fp_a),
                    "gamma_a": float(sc.gamma_a),
                    "eps_FEM": [eps[0].real, eps[0].imag],
                    "eps_Meep": [eps_meep.real, eps_meep.imag],
                    "difference": [diff.real, diff.imag],
                }
            ),
            flush=True,
        )
        if abs(diff) > 1e-8:
            raise RuntimeError(f"FEM and Meep Drude epsilon differ by {diff}")
    b = inject_mode_consistent(points, tris, src, sampler)
    x, lu, t_fac, t_sol, resid = solve_system(A, b, pec)
    raw = np.array(
        [guide_normal_flux(points, x, None, None, p, sampler=sampler, rho=rho, omega=k0) for p in range(6)]
    )
    rho_f = feed_rho(points, tris)
    Aref = assemble_anisotropic(points, tris, *rho_f, k0, pec)
    xref, *_ = solve_system(Aref, b, pec)
    pinc = abs(guide_normal_flux(points, xref, None, None, 0, sampler=sampler, rho=rho_f, omega=k0))
    if sampler.misses:
        raise RuntimeError(f"monitor samples missed the mesh: {sampler.misses}")
    norm = raw / pinc
    rec = {
        "kind": kind,
        "level": name,
        "n_nodes": int(len(points)),
        "n_tris": int(len(tris)),
        "mesh_s": t_mesh,
        "material_s": t_mat,
        "assemble_s": t_asm,
        "factor_s": t_fac,
        "solve_s": t_sol,
        "residual": resid,
        "P_inc": float(pinc),
        "normalized": norm.tolist(),
        "quartz": qstat,
        "operator": op,
        "n_quartz_tris": int(qmask.sum()),
        "n_metal_tris": int(metal.sum()),
    }
    print(name, "DOFs", rec["n_nodes"], "P", [round(v, 6) for v in norm], flush=True)
    return rec, lu, A, points, tris, rho, sampler, k0, b, x


def main() -> int:
    kind = sys.argv[1] if len(sys.argv) > 1 else "quartz"
    levels = LEVELS
    if len(sys.argv) > 2 and sys.argv[2] == "fine":
        levels = [("FEM-QF", 0.050, 0.016, 0.030, 0.035, 0.006)]
    if len(sys.argv) > 2 and sys.argv[2] == "bulk":
        levels = [("FEM-VB", 0.040, 0.012, 0.022, 0.028, 0.010)]
    if len(sys.argv) > 2 and sys.argv[2] == "finer":
        levels = [("FEM-VC", 0.032, 0.010, 0.016, 0.022, 0.008)]
    if kind == "plasma" and len(sys.argv) <= 2:
        levels = [
            ("FEM-VH", 0.050, 0.016, 0.030, 0.035, 0.010),
            ("FEM-VB", 0.040, 0.012, 0.022, 0.028, 0.010),
        ]
    OUT.mkdir(parents=True, exist_ok=True)
    set_geometry_context(res=50, horn_walls="prism", grid_offset_cells=(0, 0), coord_rotation_deg=0.0)
    src = load_numerical_mode()
    rows = []
    for name, *spec in levels:
        rec, *_rest = solve_level(name, spec, kind, src)
        rows.append(rec)
        path = OUT / f"fem_{kind}.json"
        base = []
        if path.exists() and not (kind == "quartz" and levels is LEVELS):
            base = json.loads(path.read_text())
            base = [row for row in base if row.get("level") != name]
        path.write_text(json.dumps(base + rows, indent=2) + "\n")
    # Internal convergence
    if len(rows) == 2:
        a = np.array(rows[0]["normalized"])
        b = np.array(rows[1]["normalized"])
        dB = []
        for i in range(1, 6):
            if a[i] > 0 and b[i] > 0:
                dB.append(float(10 * np.log10(b[i] / a[i])))
        print("VH minus H dB", [round(v, 4) for v in dB], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
