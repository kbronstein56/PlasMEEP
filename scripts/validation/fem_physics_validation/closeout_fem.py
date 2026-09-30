#!/usr/bin/env python3
"""Close the remaining cheap analytic and FEM audit rows.

Subcommands: freq, orient, pml, power, recip, linalg, ports, cluster.
Each writes its own JSON. Thresholds are not modified.
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
import analytic_sweeps as ans  # noqa: E402
import fem_scatterers as fsc  # noqa: E402
from analytic_maxwell import cluster_field, fresnel_ht, kz_of, solve_clusters, stack_response  # noqa: E402
from analytic_sweeps import R_CORE, R_GAP, R_SHELL, plasma, probes  # noqa: E402
from fem_scatterers import solve_case  # noqa: E402
from fem_validated_solver import (  # noqa: E402
    ElementSampler,
    PML_SIGMA_SCALE,
    assemble_anisotropic,
    pml_sx_sy,
)
import fem_validated_solver as fv  # noqa: E402
from planar_fem import assign_rho, dirichlet, rect_mesh  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
FS0 = float(sc.fs_a)


def fa(ghz: float) -> float:
    return FS0 * ghz / 3.85


def set_source(xy):
    ans.SRC[:] = xy
    fsc.SRC[:] = xy


def dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=2) + "\n")
    print("WROTE", name, flush=True)


def freq():
    saved_src = ans.SRC.copy()
    set_source([-4.5, 0.0])
    rows = []
    levels = [(0.04, 0.012), (0.02, 0.008)]
    try:
        for ghz in (3.20, 3.85, 4.50, 5.30, 6.00):
            sc.fs_a = fa(ghz)
            eps = plasma(sc.fs_a)
            lv = list(levels)
            if abs(ghz - 5.30) < 1e-6:
                lv.append((0.012, 0.005))
            rows.extend(solve_case(f"bare_{ghz:.2f}GHz", np.zeros((1, 2)), R_CORE, eps, False, lv))
        for ghz in (3.50, 3.85, 4.20, 5.30):
            sc.fs_a = fa(ghz)
            eps = plasma(sc.fs_a)
            rows.extend(solve_case(f"coated_{ghz:.2f}GHz", np.zeros((1, 2)), R_CORE, eps, True, levels))
    finally:
        sc.fs_a = FS0
        set_source(saved_src)
    dump("freq_fem.json", rows)


def pair_centers(angle_deg: float) -> np.ndarray:
    ang = np.deg2rad(angle_deg)
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    arm = rot @ np.array([0.5, 0.0])
    return np.vstack([-arm, arm])


def orient():
    saved_src = ans.SRC.copy()
    set_source([-4.5, 0.0])
    sc.fs_a = FS0
    eps = plasma(FS0)
    rows = []
    try:
        for ang in (0.0, 30.0, 45.0, 60.0, 90.0):
            rows.extend(solve_case(f"bare_pair_{ang:.0f}", pair_centers(ang), R_CORE, eps, False, [(0.04, 0.012), (0.02, 0.008)]))
        # Rotate the source with the pair. The forward probe in the rotated frame is the lab side probe only for 90 deg.
        # Compare the analytic ratio directly by a second source placement.
        set_source([0.0, -4.5])
        rows.extend(solve_case("bare_pair_90_source_rotated", pair_centers(90.0), R_CORE, eps, False, [(0.04, 0.012)]))
    finally:
        sc.fs_a = FS0
        set_source(saved_src)
    dump("orient_fem.json", rows)


def _disk_probe(dp, scale, lx, ly, h=0.04):
    fv.PML_SIGMA_SCALE = scale
    fsc.LX, fsc.LY, fsc.DPML = lx, ly, dp
    sc.fs_a = FS0
    eps = plasma(FS0)
    # Source stays 4.5 left of the origin. The mesh origin is the domain corner.
    src = np.array([-min(4.5, 0.35 * lx), 0.0])
    set_source(src)
    rows = solve_case(f"pml_disk_dp{dp}_s{scale}_L{lx}x{ly}", np.zeros((1, 2)), R_CORE, eps, False, [(h, 0.012)])
    fv.PML_SIGMA_SCALE = 1.0
    fsc.LX, fsc.LY, fsc.DPML = 14.0, 10.0, 1.2
    return rows


def oblique_wave(dp, scale, length, ly=1.6, h=0.02):
    """Oblique vacuum wave, Dirichlet on three sides, PML only on +x."""
    fv.PML_SIGMA_SCALE = scale
    k0 = 2 * np.pi * FS0
    ang = np.deg2rad(25.0)
    kx = k0 * np.cos(ang)
    ky = k0 * np.sin(ang)
    ly = float(ly)
    nx = length + 2 * dp
    ny = ly + 2 * dp
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = nx, ny, dp
    sc.fs_a = FS0
    x0, x1 = dp, dp + length
    y0, y1 = dp, dp + ly
    pts, tris, _, _ = rect_mesh(x0, x1 + dp, y0, y1, h)
    sy = pml_sx_sy(pts[tris].mean(1))[1]
    if np.max(np.abs(sy - 1)) > 1e-6:
        raise RuntimeError("unexpected y PML")

    def ufn(x, y):
        return np.exp(1j * (kx * x + ky * y))

    phys_right = x0 + length
    # Do not impose the unstretched plane wave on the top and bottom of the PML.
    bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (
        ((np.abs(pts[:, 1] - y0) < 1e-10) | (np.abs(pts[:, 1] - y1) < 1e-10)) & (pts[:, 0] <= phys_right + 1e-8)
    )
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    b = np.zeros(len(pts), np.complex128)
    A, b = dirichlet(A, b, bound, ufn(pts[:, 0], pts[:, 1]))
    uh = splu(A.tocsc()).solve(b)
    # Sample the physical region, before the PML.
    interior = (pts[:, 0] > x0 + 0.15) & (pts[:, 0] < x0 + length - 0.15) & (pts[:, 1] > y0 + 0.1) & (pts[:, 1] < y1 - 0.1)
    exact = ufn(pts[interior, 0], pts[interior, 1])
    err = uh[interior] - exact
    rel = float(np.linalg.norm(err) / np.linalg.norm(exact))
    # A probe near the PML face measures the reflected contamination.
    x_probe = x0 + length - 0.2
    i = int(np.argmin((pts[:, 0] - x_probe) ** 2 + (pts[:, 1] - (y0 + 0.5 * ly)) ** 2))
    ref = uh[i] / ufn(pts[i, 0], pts[i, 1]) - 1.0
    fv.PML_SIGMA_SCALE = 1.0
    return {
        "dpml": dp,
        "sigma_scale": scale,
        "length": length,
        "rel_L2": rel,
        "probe_complex_reflection": [float(ref.real), float(ref.imag)],
        "probe_abs": float(abs(ref)),
        "dofs": int(len(pts)),
    }


def pml():
    saved = ans.SRC.copy()
    rows = {"disk": [], "pair": [], "oblique": []}
    try:
        for dp, scale, lx, ly in (
            (1.2, 1.0, 14.0, 10.0),
            (0.8, 1.0, 14.0, 10.0),
            (1.6, 1.0, 14.0, 10.0),
            (1.2, 0.5, 14.0, 10.0),
            (1.2, 2.0, 14.0, 10.0),
            (1.2, 1.0, 11.0, 8.0),
            (1.2, 1.0, 18.0, 12.0),
        ):
            rows["disk"].extend(_disk_probe(dp, scale, lx, ly))
        fsc.LX, fsc.LY, fsc.DPML = 14.0, 10.0, 1.2
        set_source([-4.5, 0.0])
        eps = plasma(FS0)
        for dp, scale in ((1.2, 1.0), (0.8, 1.0), (1.6, 2.0)):
            fv.PML_SIGMA_SCALE = scale
            fsc.DPML = dp
            rows["pair"].extend(solve_case(f"pml_pair_dp{dp}_s{scale}", pair_centers(0.0), R_CORE, eps, False, [(0.04, 0.012)]))
        fv.PML_SIGMA_SCALE = 1.0
        fsc.DPML = 1.2
        for dp, scale, length in ((0.8, 1.0, 3.0), (1.2, 1.0, 3.0), (1.6, 1.0, 3.0), (1.2, 0.5, 3.0), (1.2, 2.0, 3.0), (1.2, 1.0, 4.2)):
            rec = oblique_wave(dp, scale, length)
            rows["oblique"].append(rec)
            print("oblique", rec, flush=True)
    finally:
        sc.fs_a = FS0
        fv.PML_SIGMA_SCALE = 1.0
        fsc.LX, fsc.LY, fsc.DPML = 14.0, 10.0, 1.2
        set_source(saved)
    dump("pml_domain.json", rows)


def contour_and_volume():
    """Analytic contour for a lossless pair, and FEM volume absorption for plasma."""
    k0 = 2 * np.pi * FS0
    ang = np.linspace(0, 2 * np.pi, 1441, endpoint=False)
    radius_c = 1.6
    xy = np.column_stack([radius_c * np.cos(ang), radius_c * np.sin(ang)])
    dl = radius_c * (ang[1] - ang[0])
    normal = xy / radius_c
    src = np.array([-4.5, 0.0])
    centers = pair_centers(0.0)

    def flux(field_fn):
        step = 1e-5
        hz = field_fn(xy)
        dx = (field_fn(xy + np.array([step, 0.0])) - field_fn(xy - np.array([step, 0.0]))) / (2 * step)
        dy = (field_fn(xy + np.array([0.0, step])) - field_fn(xy - np.array([0.0, step]))) / (2 * step)
        ex = (1j / k0) * dy
        ey = (1j / k0) * (-dx)
        sx = 0.5 * np.real(ey * np.conj(hz))
        sy = -0.5 * np.real(ex * np.conj(hz))
        return float(np.sum((sx * normal[:, 0] + sy * normal[:, 1]) * dl))

    def field_of(eps, coat, pts):
        b = solve_clusters(centers, k0, R_CORE, eps, src, 10, coated=coat)
        outer = R_CORE if coat is None else coat[0][-1]
        return cluster_field(pts, centers, k0, outer, eps, src, b, 10)

    dielectric = 3.8 + 0j
    p_die = flux(lambda pts: field_of(dielectric, None, pts))
    eps = plasma(FS0)
    coat = ([R_CORE, R_GAP, R_SHELL], [eps, 1.0 + 0j, 3.8 + 0j])
    # One coated cylinder, not the pair, for the volume comparison below.
    # Pair dielectric is the lossless cluster.

    # FEM volume absorption, one plasma cylinder.
    fsc.LX, fsc.LY, fsc.DPML = 14.0, 10.0, 1.2
    set_source(src)
    sc.fs_a = FS0
    # Reuse the scatterer assembly by a direct solve at h=0.03.
    from fem_scatterers import assign, mesh_circles, to_mesh

    centers1 = np.zeros((1, 2))
    pts, tris = mesh_circles(centers1, [(R_CORE,)], 0.03, 0.01)
    areas = {}
    rho = assign(pts, tris, centers1, R_CORE, eps, False, areas)
    sampler = ElementSampler(pts, tris)
    bvec = np.zeros(len(pts), np.complex128)
    src_m = to_mesh(src)
    t_src = int(sampler.locate(src_m[None, :])[0])
    w = sampler._bary(t_src, src_m)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    uh = splu(A.tocsc()).solve(bvec)
    cents = pts[tris].mean(1)
    twice = (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1]) - (pts[tris[:, 2], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
    area = 0.5 * np.abs(twice)
    # Element gradient.
    bx = np.stack([
        pts[tris[:, 1], 1] - pts[tris[:, 2], 1],
        pts[tris[:, 2], 1] - pts[tris[:, 0], 1],
        pts[tris[:, 0], 1] - pts[tris[:, 1], 1],
    ], axis=1) / twice[:, None]
    by = np.stack([
        pts[tris[:, 2], 0] - pts[tris[:, 1], 0],
        pts[tris[:, 0], 0] - pts[tris[:, 2], 0],
        pts[tris[:, 1], 0] - pts[tris[:, 0], 0],
    ], axis=1) / twice[:, None]
    dux = np.sum(bx * uh[tris], axis=1)
    duy = np.sum(by * uh[tris], axis=1)
    origin = cents - np.array([fsc.LX / 2.0, fsc.LY / 2.0])
    inside = np.linalg.norm(origin, axis=1) <= R_CORE
    # E = (i/k0) rho (dy Hz, -dx Hz), rho = 1/eps inside the plasma.
    rho_e = 1.0 / eps
    ex = (1j / k0) * rho_e * duy
    ey = (1j / k0) * rho_e * (-dux)
    e2 = np.abs(ex) ** 2 + np.abs(ey) ** 2
    # e^{-iωt}, passive Im(eps)>0: absorbed power (ω/2) ∫ Im(eps) |E|^2 dA.
    p_abs = float(0.5 * k0 * np.imag(eps) * np.sum(e2[inside] * area[inside]))
    # Outward flux on r=1.2 using the FEM field. Source is outside.
    ang_f = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    circ = np.column_stack([1.2 * np.cos(ang_f), 1.2 * np.sin(ang_f)])
    hz_c, ex_c, ey_c = sampler.fields(uh, *rho, k0, to_mesh(circ))
    nrm = circ / 1.2
    dl_f = 1.2 * (ang_f[1] - ang_f[0])
    sx = 0.5 * np.real(ey_c * np.conj(hz_c))
    sy = -0.5 * np.real(ex_c * np.conj(hz_c))
    p_out = float(np.sum((sx * nrm[:, 0] + sy * nrm[:, 1]) * dl_f))
    # Slab analytic absorption for the stack, as a second material case.
    slab = stack_response(k0, [(0.40, eps)])
    slab_abs = 1.0 - abs(slab["r"]) ** 2 - abs(slab["t_over_vacuum"]) ** 2
    fr = fresnel_ht(k0, 1.0, eps, 0.0)
    rec = {
        "lossless_pair_contour_power": p_die,
        "fem_plasma_volume_absorption": p_abs,
        "fem_plasma_outward_flux": p_out,
        "fem_balance_out_plus_abs": p_out + p_abs,
        "absorption_nonnegative": bool(p_abs >= 0.0),
        "plasma_slab_power_deficit": float(np.real(slab_abs)),
        "fresnel_power_sum": fr["power_sum"],
        "dofs": int(len(pts)),
        "coated_note": "pair contour is two dielectric cylinders, source outside",
    }
    # Coated single-cylinder analytic flux and a quartz/plasma/quartz check.
    def coated_field(pts):
        b = solve_clusters(np.zeros((1, 2)), k0, R_CORE, eps, src, 10, coated=coat)
        return cluster_field(pts, np.zeros((1, 2)), k0, R_SHELL, eps, src, b, 10)

    ang2 = np.linspace(0, 2 * np.pi, 1441, endpoint=False)
    circ2 = np.column_stack([1.2 * np.cos(ang2), 1.2 * np.sin(ang2)])
    rec["coated_analytic_flux"] = flux(coated_field) if False else None
    # flux() closes over radius 1.6 and the pair. Recompute for one coated cylinder on r=1.2.
    def flux_xy(field_fn, circ_xy, rad):
        step = 1e-5
        hz = field_fn(circ_xy)
        dx = (field_fn(circ_xy + np.array([step, 0.0])) - field_fn(circ_xy - np.array([step, 0.0]))) / (2 * step)
        dy = (field_fn(circ_xy + np.array([0.0, step])) - field_fn(circ_xy - np.array([0.0, step]))) / (2 * step)
        ex = (1j / k0) * dy
        ey = (1j / k0) * (-dx)
        sx = 0.5 * np.real(ey * np.conj(hz))
        sy = -0.5 * np.real(ex * np.conj(hz))
        nrm = circ_xy / rad
        return float(np.sum((sx * nrm[:, 0] + sy * nrm[:, 1]) * rad * (ang2[1] - ang2[0])))

    rec["coated_analytic_flux"] = flux_xy(coated_field, circ2, 1.2)
    rec["coated_absorption_proxy_sign"] = float(-rec["coated_analytic_flux"])
    print("power", rec, flush=True)
    dump("power_closeout.json", rec)


def guide_reciprocity():
    """Point-source exchange inside a parallel-plate guide with a quartz block."""
    k0 = 5.0
    width = 1.0
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    x0, y0, length = 2.0, 2.0, 2.4
    h = 0.02
    pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
    # Plates are natural Neumann (PEC for Hz). Ends are not Dirichlet; a small PML is not used.
    # The guide is closed by Dirichlet on the ends so the operator is the guide Green's function.
    bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)

    def eps_of(x, y):
        return 3.8 + 0j if (x0 + 1.0) <= x <= (x0 + 1.4) else 1.0 + 0j

    rho = assign_rho(pts[tris].mean(1), eps_of)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    A, _ = dirichlet(A, np.zeros(len(pts), np.complex128), bound, np.zeros(len(pts), np.complex128))
    lu = splu(A.tocsc())
    targets = [(x0 + 0.4, y0 + 0.25), (x0 + 0.4, y0 + 0.75), (x0 + 2.0, y0 + 0.25), (x0 + 2.0, y0 + 0.75)]
    idx = [int(np.argmin((pts[:, 0] - x) ** 2 + (pts[:, 1] - y) ** 2)) for x, y in targets]
    rows = []
    cols = []
    for a in idx:
        ea = np.zeros(len(pts), np.complex128)
        ea[a] = 1.0
        cols.append(lu.solve(ea))
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            diff = cols[i][idx[j]] - cols[j][idx[i]]
            rows.append({
                "a": targets[i],
                "b": targets[j],
                "abs_diff": float(abs(diff)),
                "rel": float(abs(diff) / (abs(cols[i][idx[j]]) + 1e-30)),
                "phase_deg": float(np.angle((cols[i][idx[j]] + 1e-30) / (cols[j][idx[i]] + 1e-30)) * 180 / np.pi),
            })
    rec = {"max_rel": max(r["rel"] for r in rows), "pairs": rows, "dofs": int(len(pts)), "k0": k0}
    print("guide recip", rec["max_rel"], flush=True)
    dump("guide_reciprocity.json", rec)


def horns_reciprocity():
    from fem_horns_only import metal_mask, quads_mesh
    from fem_validated_solver import load_mesh

    poly = ROOT / "outputs/validation/fem_meep_validation/horn_localization/meep_horn_polygons.json"
    doc = json.loads(poly.read_text())
    points, tris = load_mesh("FEM-L")
    quads = quads_mesh(doc)
    mask = metal_mask(points, tris, quads)
    T = len(tris)
    rho = [np.ones(T, np.complex128), np.zeros(T, np.complex128), np.zeros(T, np.complex128), np.ones(T, np.complex128)]
    for arr in rho:
        arr[mask] = 0.0
    k0 = 2 * np.pi * FS0
    A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    lu = splu(A.tocsc())
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    targets = []
    for horn in doc["horns"]:
        targets.append(tuple(np.asarray(horn["monitor_center_a"], float) + shift))
    idx = []
    for x, y in targets:
        i = int(np.argmin((points[:, 0] - x) ** 2 + (points[:, 1] - y) ** 2))
        idx.append(i)
    cols = []
    for a in idx:
        ea = np.zeros(len(points), np.complex128)
        ea[a] = 1.0
        cols.append(lu.solve(ea))
    n = len(idx)
    mat = np.zeros((n, n), np.complex128)
    for i in range(n):
        for j in range(n):
            mat[i, j] = cols[i][idx[j]]
    antisym = mat - mat.T
    rel = np.abs(antisym) / (np.abs(mat) + 1e-30)
    rec = {
        "grade": "FEM-L",
        "dofs": int(len(points)),
        "n_ports": n,
        "max_abs": float(np.max(np.abs(antisym))),
        "max_rel": float(np.max(rel)),
        "targets": targets,
    }
    print("horns recip", rec["max_rel"], "dofs", rec["dofs"], flush=True)
    dump("horns_reciprocity.json", rec)


def seven_reciprocity():
    sc.fs_a = FS0
    eps = plasma(FS0)
    k0 = 2 * np.pi * FS0
    centers = np.array(
        [[0.0, 0.0], [-0.86602540378, 0.5], [-0.86602540378, -0.5], [0.86602540378, 0.5], [0.86602540378, -0.5], [0.0, 1.0], [0.0, -1.0]]
    )
    from fem_scatterers import assign, mesh_circles, to_mesh

    fsc.LX, fsc.LY, fsc.DPML = 14.0, 10.0, 1.2
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 14.0, 10.0, 1.2
    pts, tris = mesh_circles(centers, [(R_CORE, R_GAP, R_SHELL)] * 7, 0.05, 0.016)
    areas = {}
    rho = assign(pts, tris, centers, R_CORE, eps, True, areas)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    lu = splu(A.tocsc())
    targets = [(-2.2, 0.3), (2.4, -0.2), (0.55, 0.15), (-0.4, 1.6), (1.3, -1.5), (0.0, -2.2)]
    idx = []
    for x, y in targets:
        m = to_mesh([x, y])
        idx.append(int(np.argmin((pts[:, 0] - m[0]) ** 2 + (pts[:, 1] - m[1]) ** 2)))
    cols = []
    for a in idx:
        ea = np.zeros(len(pts), np.complex128)
        ea[a] = 1.0
        cols.append(lu.solve(ea))
    rows = []
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            diff = cols[i][idx[j]] - cols[j][idx[i]]
            rows.append({"rel": float(abs(diff) / (abs(cols[i][idx[j]]) + 1e-30)), "abs": float(abs(diff))})
    rec = {"max_rel": max(r["rel"] for r in rows), "n_pairs": len(rows), "dofs": int(len(pts)), "pairs": rows}
    print("seven recip", rec["max_rel"], flush=True)
    dump("seven_reciprocity.json", rec)


def linalg():
    """Repeat factorizations and a 1-norm condition estimate."""
    from scipy.sparse.linalg import LinearOperator, onenormest

    sc.fs_a = FS0
    k0 = 2 * np.pi * FS0
    # Small guide.
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
    pts, tris, _, _ = rect_mesh(2.2, 3.8, 2.2, 3.2, 0.02)
    bound = (np.abs(pts[:, 0] - 2.2) < 1e-10) | (np.abs(pts[:, 0] - 3.8) < 1e-10)
    rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
    A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
    b = np.zeros(len(pts), np.complex128)
    A1, b1 = dirichlet(A, b, bound, np.zeros(len(pts), np.complex128))
    free = np.flatnonzero(~bound)
    b1[free[len(free) // 2]] = 1.0
    x1 = splu(A1.tocsc()).solve(b1)
    x2 = splu(A1.tocsc(), permc_spec="NATURAL").solve(b1)
    x3 = splu(A1.tocsc(), permc_spec="COLAMD").solve(b1)
    resid = float(np.linalg.norm(A1 @ x1 - b1) / np.linalg.norm(b1))

    def condest(mat):
        lu = splu(mat.tocsc())
        n = mat.shape[0]
        op = LinearOperator(
            (n, n),
            matvec=lambda v: lu.solve(np.asarray(v, np.complex128).ravel()),
            rmatvec=lambda v: lu.solve(np.asarray(v, np.complex128).ravel(), trans="H"),
            dtype=np.complex128,
        )
        return float(onenormest(mat) * onenormest(op))

    guide_cond = condest(A1)
    # Plasma disk at the coarse scatterer mesh.
    from fem_scatterers import assign, mesh_circles

    eps = plasma(FS0)
    fsc.LX, fsc.LY, fsc.DPML = 14.0, 10.0, 1.2
    sc.nx_ports, sc.ny_ports, sc.dpml_ports = 14.0, 10.0, 1.2
    dpts, dtris = mesh_circles(np.zeros((1, 2)), [(R_CORE,)], 0.05, 0.016)
    areas = {}
    drho = assign(dpts, dtris, np.zeros((1, 2)), R_CORE, eps, False, areas)
    Ad = assemble_anisotropic(dpts, dtris, *drho, k0, np.zeros(len(dpts), dtype=bool))
    disk_cond = condest(Ad)
    bd = np.zeros(len(dpts), np.complex128)
    bd[len(dpts) // 5] = 1.0
    xd1 = splu(Ad.tocsc()).solve(bd)
    xd2 = splu(Ad.tocsc()).solve(bd)
    rec = {
        "guide_dofs": int(len(pts)),
        "guide_residual": resid,
        "guide_reorder_rel": float(np.linalg.norm(x2 - x3) / np.linalg.norm(x1)),
        "guide_repeat_rel": float(np.linalg.norm(x1 - splu(A1.tocsc()).solve(b1)) / np.linalg.norm(x1)),
        "guide_cond_1": guide_cond,
        "disk_dofs": int(len(dpts)),
        "disk_cond_1": disk_cond,
        "disk_repeat_rel": float(np.linalg.norm(xd1 - xd2) / np.linalg.norm(xd1)),
        "disk_residual": float(np.linalg.norm(Ad @ xd1 - bd) / np.linalg.norm(bd)),
    }
    print("linalg", rec, flush=True)
    dump("linalg_closeout.json", rec)


def ports():
    """Guide-normal power versus station, mesh, and strip width."""
    width, k0 = 1.0, 5.0
    ky = np.pi / width
    beta = kz_of(k0, 1.0, ky)
    p_exact = 0.5 * np.real(beta) * (width / 2.0) / k0
    rows = []
    for h in (0.02, 0.01, 0.005):
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = 8.0, 8.0, 1.0
        x0, y0, length = 2.2, 2.2, 1.6

        def ufn(x, y, ky=ky, beta=beta, x0=x0, y0=y0):
            return np.cos(ky * (y - y0)) * np.exp(1j * beta * (x - x0))

        pts, tris, _, _ = rect_mesh(x0, x0 + length, y0, y0 + width, h)
        bound = (np.abs(pts[:, 0] - x0) < 1e-10) | (np.abs(pts[:, 0] - (x0 + length)) < 1e-10)
        rho = assign_rho(pts[tris].mean(1), lambda x, y: 1.0 + 0j)
        A = assemble_anisotropic(pts, tris, *rho, k0, np.zeros(len(pts), dtype=bool))
        b = np.zeros(len(pts), np.complex128)
        A, b = dirichlet(A, b, bound, ufn(pts[:, 0], pts[:, 1]))
        uh = splu(A.tocsc()).solve(b)
        cents = pts[tris].mean(1)
        twice = (pts[tris[:, 1], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 2], 1] - pts[tris[:, 0], 1]) - (pts[tris[:, 2], 0] - pts[tris[:, 0], 0]) * (pts[tris[:, 1], 1] - pts[tris[:, 0], 1])
        bx = np.stack([
            pts[tris[:, 1], 1] - pts[tris[:, 2], 1],
            pts[tris[:, 2], 1] - pts[tris[:, 0], 1],
            pts[tris[:, 0], 1] - pts[tris[:, 1], 1],
        ], axis=1) / twice[:, None]
        dux = np.sum(bx * uh[tris], axis=1)
        ey = (1j / k0) * (-dux)
        # Hz at the centroid from the nodal interpolant.
        hz = uh[tris].mean(1)
        sx = 0.5 * np.real(ey * np.conj(hz))
        stations = {}
        for frac in (0.25, 0.5, 0.75):
            xcut = x0 + frac * length
            on = np.abs(cents[:, 0] - xcut) <= 0.51 * h
            area = 0.5 * np.abs(twice)
            stations[str(frac)] = float(np.sum(sx[on] * area[on]) / np.sum(area[on]) * width)
        rels = {k: abs(v - p_exact) / abs(p_exact) for k, v in stations.items()}
        rows.append({"h": h, "dofs": int(len(pts)), "stations": stations, "power_rel": rels, "exact": p_exact})
        print("ports", rows[-1], flush=True)
    dump("port_stations.json", rows)


def cluster_finer():
    """One locally refined seven-coated mesh. Global h stays 0.02; the wall is finer."""
    saved = ans.SRC.copy()
    set_source([-4.5, 0.0])
    sc.fs_a = FS0
    eps = plasma(FS0)
    centers = np.array(
        [[0.0, 0.0], [-0.86602540378, 0.5], [-0.86602540378, -0.5], [0.86602540378, 0.5], [0.86602540378, -0.5], [0.0, 1.0], [0.0, -1.0]]
    )
    try:
        rows = solve_case("coated_7_local", centers, R_CORE, eps, True, [(0.02, 0.0045)])
    finally:
        sc.fs_a = FS0
        set_source(saved)
    dump("coated7_local.json", rows)


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "freq"
    {"freq": freq, "orient": orient, "pml": pml, "power": contour_and_volume, "recip": guide_reciprocity,
     "horns": horns_reciprocity, "seven": seven_reciprocity, "linalg": linalg, "ports": ports,
     "cluster": cluster_finer}[cmd]()


if __name__ == "__main__":
    main()
