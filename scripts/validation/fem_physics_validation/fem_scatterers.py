#!/usr/bin/env python3
"""Body-fitted scatterer sweeps against the independent T-matrix.

Same point-source convention as analytic_sweeps.py: H0 at (-4.5, 0), ratios
Hz/Hz_vacuum. No fit to Meep or to a previous FEM solve.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu
from scipy.special import hankel1

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, solve_clusters  # noqa: E402
from analytic_sweeps import R_CORE, R_GAP, R_SHELL, SRC, plasma, probes  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from plasma_interface import _tri  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
LX, LY, DPML = 14.0, 10.0, 1.2


def to_mesh(xy):
    return np.asarray(xy, float) + np.array([LX / 2.0, LY / 2.0])


def mesh_circles(centers, radii, h, h_edge, h_near=None, r_near=None):
    tr = _tri()
    verts = []
    segs = []
    index = {}

    def add_vertex(p):
        key = (round(float(p[0]), 8), round(float(p[1]), 8))
        if key in index:
            return index[key]
        index[key] = len(verts)
        verts.append([key[0], key[1]])
        return index[key]

    def add_chain(pts, closed):
        ids = [add_vertex(p) for p in pts]
        pairs = list(zip(ids, ids[1:]))
        if closed:
            pairs.append((ids[-1], ids[0]))
        for a, b in pairs:
            if a != b:
                segs.append([a, b])

    corners = [(0.0, 0.0), (LX, 0.0), (LX, LY), (0.0, LY)]
    for i in range(4):
        p0 = np.array(corners[i], float)
        p1 = np.array(corners[(i + 1) % 4], float)
        m = max(2, int(np.ceil(np.linalg.norm(p1 - p0) / h)))
        t = np.linspace(0.0, 1.0, m)
        add_chain((1 - t)[:, None] * p0 + t[:, None] * p1, closed=False)
    regions = []
    for c, rads in zip(centers, radii):
        mc = to_mesh(c)
        for radius in rads:
            m = max(48, int(np.ceil(2 * np.pi * radius / h_edge)))
            ang = np.linspace(0.0, 2 * np.pi, m, endpoint=False)
            ring = mc + radius * np.column_stack([np.cos(ang), np.sin(ang)])
            add_chain(ring, closed=True)
        regions.append([mc[0], mc[1], 1, 0.45 * h_edge * h_edge])
        if len(rads) == 3:
            mid = 0.5 * (rads[1] + rads[2])
            regions.append([mc[0] + mid, mc[1], 2, 0.45 * h_edge * h_edge])
        if h_near is not None and r_near is not None:
            m = max(64, int(np.ceil(2 * np.pi * r_near / h_near)))
            ang = np.linspace(0.0, 2 * np.pi, m, endpoint=False)
            ring = mc + r_near * np.column_stack([np.cos(ang), np.sin(ang)])
            add_chain(ring, closed=True)
            outer = float(rads[-1])
            air_r = 0.5 * (outer + r_near)
            regions.append([mc[0] + air_r, mc[1], 3, 0.45 * h_near * h_near])
            if len(rads) == 3:
                gap_r = 0.5 * (float(rads[0]) + float(rads[1]))
                regions.append([mc[0] + gap_r, mc[1], 4, 0.45 * h_near * h_near])
    if h_near is None:
        amax = 0.45 * h * h
        opts = f"pq20a{amax:.8f}A"
    else:
        # A numeric area after -a makes Triangle ignore the per-region areas.
        regions.append([0.25, 0.25, 5, 0.45 * h * h])
        opts = "pq20aA"
    mesh = tr.triangulate(
        {"vertices": np.asarray(verts, float), "segments": np.asarray(segs, np.int32), "regions": np.asarray(regions, float)},
        opts,
    )
    return np.asarray(mesh["vertices"], float), np.asarray(mesh["triangles"], int)


def rho_of(pts, tris, centers, kind, eps, coated):
    T = len(tris)
    rho = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    cents = pts[tris].mean(1) - np.array([LX / 2.0, LY / 2.0])
    plasma_area = 0.0
    quartz_area = 0.0
    areas = 0.5 * np.abs(
        (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1])
        - (pts[tris[:, 2], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
    )
    r_plasma = R_CORE if kind != "custom_radius" else None
    for i, c in enumerate(centers):
        d = np.linalg.norm(cents - c, axis=1)
        if coated:
            plasma = d <= R_CORE
            quartz = (d > R_GAP) & (d <= R_SHELL)
        else:
            plasma = d <= eps  # placeholder overwritten below
            quartz = np.zeros(T, dtype=bool)
        # The radius is passed separately. This branch is replaced by the caller.
        del plasma, quartz
    return rho, areas


def assign(pts, tris, centers, radius, material, coated, areas_out):
    T = len(tris)
    rho = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    origin = pts[tris].mean(1) - np.array([LX / 2.0, LY / 2.0])
    a = 0.5 * np.abs(
        (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1])
        - (pts[tris[:, 2], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
    )
    plasma_mask = np.zeros(T, dtype=bool)
    quartz_mask = np.zeros(T, dtype=bool)
    for c in centers:
        d = np.linalg.norm(origin - c, axis=1)
        plasma_mask |= d <= radius
        if coated:
            quartz_mask |= (d > R_GAP) & (d <= R_SHELL)
    rho[0][plasma_mask] = 1.0 / material
    rho[3][plasma_mask] = 1.0 / material
    rho[0][quartz_mask] = 1.0 / 3.8
    rho[3][quartz_mask] = 1.0 / 3.8
    areas_out["plasma"] = float(a[plasma_mask].sum())
    areas_out["quartz"] = float(a[quartz_mask].sum())
    areas_out["plasma_exact"] = float(len(centers) * np.pi * radius ** 2)
    areas_out["quartz_exact"] = float(len(centers) * np.pi * (R_SHELL ** 2 - R_GAP ** 2)) if coated else 0.0
    return rho


def solve_case(name, centers, radius, material, coated, levels):
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = LX, LY, DPML
    sc.fs_a = float(np.real(material)) * 0 + float(sc.fs_a)  # keep fs; k0 set below
    k0 = 2 * np.pi * float(sc.fs_a)
    # Material cases that are not the production plasma still use fs as the frequency.
    rows = []
    outer = R_SHELL if coated else radius
    probe = probes(outer)
    names = list(probe)
    xy = np.array(list(probe.values()), float)
    radii = [ (R_CORE, R_GAP, R_SHELL) if coated else (radius,) ]
    rad_list = [radii[0] for _ in centers]
    for spec in levels:
        h, h_edge = spec[0], spec[1]
        h_near = spec[2] if len(spec) > 2 else None
        r_near = spec[3] if len(spec) > 3 else None
        print(f"=== {name} h={h} h_near={h_near} r_near={r_near} ===", flush=True)
        pts, tris = mesh_circles(centers, rad_list, h, h_edge, h_near, r_near)
        h_stats = {}
        if h_near is not None:
            cents = pts[tris].mean(1) - np.array([LX / 2.0, LY / 2.0])
            edges = np.stack(
                [
                    np.linalg.norm(pts[tris[:, 1]] - pts[tris[:, 0]], axis=1),
                    np.linalg.norm(pts[tris[:, 2]] - pts[tris[:, 1]], axis=1),
                    np.linalg.norm(pts[tris[:, 0]] - pts[tris[:, 2]], axis=1),
                ],
                axis=1,
            ).mean(1)
            rad = np.linalg.norm(cents, axis=1)
            quartz = (rad > R_GAP) & (rad <= R_SHELL)
            h_stats = {
                "h_plasma_median": float(np.median(edges[rad <= radius])),
                "h_quartz_median": float(np.median(edges[quartz])) if np.any(quartz) else None,
                "h_near_median": float(np.median(edges[(rad > R_SHELL) & (rad < r_near)])),
                "h_far_median": float(np.median(edges[rad > r_near + 0.5])) if np.any(rad > r_near + 0.5) else None,
                "quartz_elements_across": float((R_SHELL - R_GAP) / np.median(edges[quartz])) if np.any(quartz) else None,
            }
        print(f"  nodes {len(pts)} tris {len(tris)}", flush=True)
        areas = {}
        rho_obj = assign(pts, tris, centers, radius, material, coated, areas)
        T = len(tris)
        rho_vac = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
        sampler = ElementSampler(pts, tris)
        src = to_mesh(SRC)
        b = np.zeros(len(pts), np.complex128)
        t_src = int(sampler.locate(src[None, :])[0])
        if t_src < 0:
            raise RuntimeError("source missed the mesh")
        w = sampler._bary(t_src, src)
        for k in range(3):
            b[int(tris[t_src, k])] += complex(w[k])
        pec = np.zeros(len(pts), dtype=bool)
        fields = {}
        residuals = {}
        for tag, rho in (("vac", rho_vac), ("obj", rho_obj)):
            A = assemble_anisotropic(pts, tris, *rho, k0, pec)
            x = splu(A.tocsc()).solve(b)
            residuals[tag] = float(np.linalg.norm(A @ x - b) / (np.linalg.norm(b) + 1e-30))
            hz, ex, ey = sampler.fields(x, *rho, k0, to_mesh(xy))
            fields[tag] = (hz, ex, ey, x)
        ratio = fields["obj"][0] / fields["vac"][0]
        bcoef = solve_clusters(centers, k0, radius, complex(material), SRC, 12, coated=([R_CORE, R_GAP, R_SHELL], [complex(material), 1.0, 3.8]) if coated else None)
        analytic = cluster_field(xy, centers, k0, outer, complex(material), SRC, bcoef, 12)
        vac_a = hankel1(0, k0 * np.linalg.norm(xy - SRC, axis=1))
        ar = analytic / vac_a
        per = {}
        for i, pname in enumerate(names):
            if not np.isfinite(ar[i]) or not np.isfinite(ratio[i]):
                per[pname] = None
                continue
            d = ratio[i] - ar[i]
            per[pname] = {
                "abs_err": float(abs(d)),
                "db": float(20 * np.log10(abs(ratio[i]) / abs(ar[i]))) if abs(ar[i]) > 1e-8 else None,
                "phase_deg": float(np.angle(ratio[i] / ar[i]) * 180 / np.pi),
                "fem": [float(ratio[i].real), float(ratio[i].imag)],
                "analytic": [float(ar[i].real), float(ar[i].imag)],
            }
        # Mirror residual for geometries symmetric in y: side_p60 versus side_m60 of the FEM ratio.
        mirror = None
        if "side_p60" in per and per["side_p60"] and per["side_m60"]:
            a = complex(*per["side_p60"]["fem"])
            c = complex(*per["side_m60"]["fem"])
            mirror = float(abs(a - c))
        ang_r = np.linspace(0.0, 2 * np.pi, 96, endpoint=False)
        ring = np.column_stack([np.cos(ang_r), np.sin(ang_r)])
        hz_obj = sampler.fields(fields["obj"][3], *rho_obj, k0, to_mesh(ring))[0]
        hz_vac = sampler.fields(fields["vac"][3], *rho_vac, k0, to_mesh(ring))[0]
        ratio_r = hz_obj / hz_vac
        b_ring = cluster_field(ring, centers, k0, outer, complex(material), SRC, bcoef, 12)
        vac_r = hankel1(0, k0 * np.linalg.norm(ring - SRC, axis=1))
        ar_r = b_ring / vac_r
        ring_l2 = float(np.linalg.norm(ratio_r - ar_r) / np.linalg.norm(ar_r))
        rec = {
            "case": name,
            "h": h,
            "h_edge": h_edge,
            "h_near": h_near,
            "r_near": r_near,
            "dofs": int(len(pts)),
            "n_tris": int(len(tris)),
            "areas": areas,
            "residual_obj": residuals["obj"],
            "probes": per,
            "mirror_side_abs": mirror,
            "quartz_cells_across_wall": (R_SHELL - R_GAP) / h_edge if coated else None,
            "mesh_stats": h_stats,
            "ring_l2": ring_l2,
        }
        rows.append(rec)
        fwd = per.get("forward")
        print("  forward", fwd, "resid", residuals["obj"], "area", areas, flush=True)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Keep the physical fs. material epsilon is separate.
    saved = float(sc.fs_a)
    eps = plasma(saved)
    pair = np.array([[0.0, 0.0], [1.0, 0.0]])
    ang = np.deg2rad(30.0)
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    three = np.array(json.loads((ROOT / "outputs/validation/fem_meep_validation/plasma_ring_campaign/fem_cluster3.json").read_text())[-1]["centers_a"])
    seven = np.array(json.loads((ROOT / "outputs/validation/fem_meep_validation/plasma_ring_campaign/fem_cluster7.json").read_text())[-1]["centers_a"])
    coarse = [(0.04, 0.012)]
    fine = [(0.04, 0.012), (0.02, 0.008)]
    # Which cases: environment CASE=name or all.
    catalog = {
        "bare_r012": (np.zeros((1, 2)), 0.12, eps, False, fine),
        "bare_prod": (np.zeros((1, 2)), R_CORE, eps, False, fine + [(0.012, 0.005)]),
        "bare_r040": (np.zeros((1, 2)), 0.40, eps, False, fine),
        "dielectric": (np.zeros((1, 2)), R_CORE, 3.8 + 0j, False, fine),
        "neg_1.5": (np.zeros((1, 2)), R_CORE, -1.5 + 0j, False, coarse),
        "neg_8": (np.zeros((1, 2)), R_CORE, -8.0 + 0.01j, False, coarse),
        "coated_1": (np.zeros((1, 2)), R_CORE, eps, True, fine),
        "bare_2": (pair, R_CORE, eps, False, fine),
        "bare_2_rot30": (pair @ rot.T, R_CORE, eps, False, coarse),
        "coated_2": (pair, R_CORE, eps, True, coarse),
        "bare_3": (three, R_CORE, eps, False, fine),
        "coated_3": (three, R_CORE, eps, True, fine),
        "bare_7": (seven, R_CORE, eps, False, fine),
        "coated_7": (seven, R_CORE, eps, True, fine),
    }
    chosen = sys.argv[1:] or ["bare_prod", "dielectric", "coated_1", "bare_2", "coated_7"]
    all_rows = []
    for name in chosen:
        centers, radius, material, coated, levels = catalog[name]
        all_rows.extend(solve_case(name, centers, radius, complex(material), coated, levels))
        (OUT / "scatter_partial.json").write_text(json.dumps(all_rows, indent=2) + "\n")
    sc.fs_a = saved
    print("SCATTER_DONE", len(all_rows), flush=True)


if __name__ == "__main__":
    main()
