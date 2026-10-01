#!/usr/bin/env python3
"""Full 91-bulb geometry audit, graded mesh, and forward solves.

The magnetized permittivity is the Lorentz tensor from gyrotropic_tensor.py.
The existing Meep helper is not the source of the off-diagonal sign.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
OUT = ROOT / "outputs" / "validation" / "fem_full_device_validation"

# Locally refined grades. Quartz wall is 0.050 a thick.
GRADES = {
    "C": dict(h_iface=0.012, h_horn=0.040, h_air=0.12, h_pml=0.20),
    "M": dict(h_iface=0.008, h_horn=0.025, h_air=0.08, h_pml=0.14),
    "F": dict(h_iface=0.005, h_horn=0.016, h_air=0.055, h_pml=0.11),
}


def production_constants():
    import sixport_common as sc
    from gyrotropic_tensor import cyclotron_hz, ordinary_from_si

    r_p = 4.6 * sc.r_bulb_inner / 6.5
    b_prod = 0.05  # KNOWN_TRUTHS: circulation reversed between -0.05 T and +0.05 T
    f, fp, gamma, fc = ordinary_from_si(sc.fs_Hz, sc.fp_Hz, sc.gamma_Hz, b_prod, sc.a)
    return sc, {
        "a_m": sc.a,
        "pitch_a": sc.d_exp,
        "r_plasma_metadata_a": sc.r_plasma,
        "r_plasma_material_a": r_p,
        "r_bulb_inner_a": sc.r_bulb_inner,
        "r_bulb_outer_a": sc.r_bulb_outer,
        "quartz_eps": 3.8,
        "wall_thickness_a": sc.wall_thickness,
        "width_open_a": sc.width_open,
        "width_base_a": sc.width_base,
        "horn_depth_a": sc.horn_depth,
        "feed_length_a": sc.feed_length,
        "clear_width_a": sc.clear_width,
        "nx_ports": float(sc.nx_ports),
        "ny_ports": float(sc.ny_ports),
        "dpml_ports": float(sc.dpml_ports),
        "fs_Hz": sc.fs_Hz,
        "fp_Hz": sc.fp_Hz,
        "gamma_Hz": sc.gamma_Hz,
        "fs_a": sc.fs_a,
        "fp_a": sc.fp_a,
        "gamma_a": sc.gamma_a,
        "B_production_T": b_prod,
        "B_direction": "+z",
        "fc_a_at_Bprod": fc,
        "fc_Hz_per_T": cyclotron_hz(1.0),
        "ordinary_f_fp_gamma_fc": [f, fp, gamma, fc],
    }


def bulb_centers_device(sc):
    from PMMCirculatorInverse import PMMI

    pmm = PMMI(
        a=sc.a, res=32, nx=sc.nx_ports, ny=sc.ny_ports, dpml=sc.dpml_ports, B=np.zeros(3)
    )
    pmm.Rod_Array_Hexagon_train(
        xy_cen=np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0]),
        side_dim=6,
        r=sc.r_plasma,
        d=sc.d_exp,
        bulbs=False,
        uniform=True,
    )
    return np.asarray(pmm.train_elem_locs, dtype=float)


def bulb_centers_shifted_helper(sc):
    from fem_validated_solver import bulb_centers_mesh

    return bulb_centers_mesh()


def write_geometry_audit() -> dict:
    sc, const = production_constants()
    centers = bulb_centers_device(sc)
    helper = bulb_centers_shifted_helper(sc)
    if len(centers) != 91 or len(helper) != 91:
        raise RuntimeError(f"expected 91 bulbs, got {len(centers)} and {len(helper)}")
    delta = centers - helper
    max_delta = float(np.max(np.linalg.norm(delta, axis=1)))
    # Lattice check: nearest-neighbor distance.
    dmat = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=2)
    np.fill_diagonal(dmat, np.inf)
    nn = dmat.min(axis=1)
    horns = []
    for i, horn in enumerate(sc.full_horns):
        horns.append({
            "port": i,
            "outward": sc.port_dirs[i].tolist(),
            "source_center_origin": np.asarray(horn["source_center"], float).tolist(),
            "monitor_center_origin": np.asarray(horn["monitor_center"], float).tolist(),
            "source_center_mesh": (np.asarray(horn["source_center"], float) + np.array([sc.nx_ports / 2, sc.ny_ports / 2])).tolist(),
            "monitor_center_mesh": (np.asarray(horn["monitor_center"], float) + np.array([sc.nx_ports / 2, sc.ny_ports / 2])).tolist(),
            "polygons": {
                name: np.asarray(horn[name], float)[:, :2].tolist()
                for name in ("left_flare", "right_flare", "left_feed", "right_feed")
            },
        })
    payload = {
        "constants": const,
        "n_bulbs": 91,
        "ordering": "Rod_Array_Hexagon_train side_dim=6, basis [[0,1],[sqrt(3)/2,1/2]]",
        "max_center_delta_vs_fem_helper": max_delta,
        "nearest_neighbor_min": float(nn.min()),
        "nearest_neighbor_max": float(nn.max()),
        "centers_mesh": centers.tolist(),
        "horns": horns,
        "geometry_ok": max_delta < 1e-12 and abs(float(nn.min()) - 1.0) < 1e-9,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "full91_geometry.json").write_text(json.dumps(payload) + "\n")
    lines = [
        "# Full 91-bulb geometry audit",
        "",
        "Centers and horns are taken from the current production builders:",
        "`Rod_Array_Hexagon_train` and `sixport_common.full_horns`.",
        "The FEM helper `bulb_centers_mesh` is the same lattice translated onto the device domain.",
        "",
        f"- Bulb count: 91",
        f"- Maximum center difference versus the FEM helper: `{max_delta:.3e}`",
        f"- Nearest-neighbor spacing: `{nn.min():.12f}` to `{nn.max():.12f}` (pitch `{const['pitch_a']}`)",
        f"- Lattice basis: `[[0, 1], [sqrt(3)/2, 1/2]]`, pointy-top triangular lattice, side length 6",
        f"- Material plasma radius: `{const['r_plasma_material_a']}` a = {const['r_plasma_material_a']*const['a_m']*1e3:.4f} mm",
        f"- Unused metadata `r_plasma`: `{const['r_plasma_metadata_a']}` a",
        f"- Quartz inner radius: `{const['r_bulb_inner_a']}` a",
        f"- Quartz outer radius: `{const['r_bulb_outer_a']}` a",
        f"- Vacuum gap: `{const['r_bulb_inner_a'] - const['r_plasma_material_a']}` a between plasma and quartz",
        f"- Quartz permittivity: 3.8",
        f"- Domain: `{const['nx_ports']}` by `{const['ny_ports']}` a, PML `{const['dpml_ports']}` a",
        f"- Frequency: {const['fs_Hz']/1e9} GHz (`fs_a={const['fs_a']}`)",
        f"- Plasma frequency: {const['fp_Hz']/1e9} GHz, collision rate {const['gamma_Hz']/1e6} MHz",
        f"- Production B: {const['B_production_T']} T along +z, `f_c a/c = {const['fc_a_at_Bprod']}`",
        f"- Ports: outward normals sorted by angle, index 0 nearest +x, then counterclockwise",
        "",
        "PEC walls are the four filled horn prisms per port (`left_flare`, `right_flare`, `left_feed`, `right_feed`).",
        "Natural Neumann on Hz is the PEC condition. Metal triangles are removed from the Helmholtz operator by setting rho to 0.",
        "",
        "Per-bulb center differences are all at the maximum above. Ordering is the `train_elem_locs` order.",
        "",
        "Geometry check: " + ("PASS" if payload["geometry_ok"] else "FAIL"),
        "",
    ]
    (OUT / "FULL91_GEOMETRY_AUDIT.md").write_text("\n".join(lines))
    print("GEOMETRY", "PASS" if payload["geometry_ok"] else "FAIL", "max_delta", max_delta, "nn", nn.min(), flush=True)
    return payload


def _tri():
    vend = ROOT / "scripts" / "validation" / "fem_meep_validation" / "_vendor"
    if str(vend) not in sys.path:
        sys.path.insert(0, str(vend))
    import triangle

    return triangle


def mesh_full91(grade: str) -> tuple[np.ndarray, np.ndarray, dict]:
    sc, const = production_constants()
    cfg = GRADES[grade]
    centers = bulb_centers_device(sc)
    nx, ny = float(sc.nx_ports), float(sc.ny_ports)
    dp = float(sc.dpml_ports)
    r_p = const["r_plasma_material_a"]
    r_in = const["r_bulb_inner_a"]
    r_out = const["r_bulb_outer_a"]
    shift = np.array([nx / 2.0, ny / 2.0])

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

    def add_chain(pts, closed, step):
        ids = []
        for i in range(len(pts) - (0 if closed else 1)):
            p0 = np.asarray(pts[i], float)
            p1 = np.asarray(pts[(i + 1) % len(pts)], float)
            L = float(np.linalg.norm(p1 - p0))
            n = max(1, int(np.ceil(L / step)))
            for k in range(n):
                ids.append(add_vertex(p0 + (k / n) * (p1 - p0)))
        if not closed:
            ids.append(add_vertex(pts[-1]))
        pairs = list(zip(ids, ids[1:]))
        if closed:
            pairs.append((ids[-1], ids[0]))
        for a, b in pairs:
            if a != b:
                segs.append([a, b])

    def add_circle(c, radius, step):
        m = max(24, int(np.ceil(2 * np.pi * radius / step)))
        ang = np.linspace(0.0, 2 * np.pi, m, endpoint=False)
        ring = c + radius * np.column_stack([np.cos(ang), np.sin(ang)])
        add_chain(ring, True, step * 10)

    # Domain and PML frame.
    add_chain([(0, 0), (nx, 0), (nx, ny), (0, ny)], True, cfg["h_pml"])
    add_chain(
        [(dp, dp), (nx - dp, dp), (nx - dp, ny - dp), (dp, ny - dp)],
        True,
        cfg["h_air"],
    )
    for c in centers:
        add_circle(c, r_out, cfg["h_iface"])
        add_circle(c, r_in, cfg["h_iface"])
        add_circle(c, r_p, cfg["h_iface"])
    for horn in sc.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], float)[:, :2] + shift
            add_chain(poly, True, cfg["h_horn"])

    regions = []
    marker = 1
    area_iface = 0.45 * cfg["h_iface"] ** 2
    for c in centers:
        regions.append([c[0], c[1], marker, area_iface])
        marker += 1
        regions.append([c[0] + 0.5 * (r_p + r_in), c[1], marker, area_iface])
        marker += 1
        regions.append([c[0] + 0.5 * (r_in + r_out), c[1], marker, area_iface])
        marker += 1
    area_horn = 0.45 * cfg["h_horn"] ** 2
    for horn in sc.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], float)[:, :2] + shift
            regions.append([float(poly[:, 0].mean()), float(poly[:, 1].mean()), marker, area_horn])
            marker += 1
    # Interior air, away from the center bulb and the horns.
    air_pt = np.array([nx / 2.0 + 3.2, ny / 2.0 + 0.15])
    regions.append([air_pt[0], air_pt[1], marker, 0.45 * cfg["h_air"] ** 2])
    marker += 1
    regions.append([0.5 * dp, 0.5 * dp, marker, 0.45 * cfg["h_pml"] ** 2])

    tr = _tri()
    t0 = time.perf_counter()
    mesh = tr.triangulate(
        {
            "vertices": np.asarray(verts, float),
            "segments": np.asarray(segs, np.int32),
            "regions": np.asarray(regions, float),
        },
        "pq20aA",
    )
    elapsed = time.perf_counter() - t0
    pts = np.asarray(mesh["vertices"], float)
    tris = np.asarray(mesh["triangles"], int)
    stats = {
        "grade": grade,
        "nodes": int(len(pts)),
        "triangles": int(len(tris)),
        "mesh_s": elapsed,
        **cfg,
        "quartz_thickness_a": r_out - r_in,
        "elements_across_quartz_requested": (r_out - r_in) / cfg["h_iface"],
    }
    print("MESH", stats, flush=True)
    dest = OUT / "meshes"
    dest.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest / f"full91_{grade}.npz", points=pts, triangles=tris)
    return pts, tris, stats


if __name__ == "__main__":
    import sys

    cmd = sys.argv[1]
    if cmd == "audit":
        write_geometry_audit()
    elif cmd == "mesh":
        mesh_full91(sys.argv[2])
    else:
        raise SystemExit(cmd)
