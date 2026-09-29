#!/usr/bin/env python3
"""
PEC mass-term audit, excised-domain Neumann, PML strength, domain size, mesh seeds.

Does not change production defaults. Horn polygons come from meep_horn_polygons.json.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath
from scipy import sparse
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from fem_validated_solver import (  # noqa: E402
    fields_from_hz,
    load_mesh,
    load_numerical_mode,
    solve_system,
)
import sixport_common as sc  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "boundary_pml"
POLY = ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_localization" / "meep_horn_polygons.json"
MEEP_NORM = np.array([np.nan, 0.03763160465218182, 0.08553915041606872, 0.696684297884412, 0.08553915041606859, 0.037631604652181856])
K0 = 2 * np.pi * float(sc.fs_a)
# Meep quadratic PML: σ(u) = [-ln R /(2 d ∫u^2)] u^2, ∫_0^1 u^2 = 1/3.
MEEP_SIGMA_MAX = -np.log(1e-15) * 3.0 / (2.0 * float(sc.dpml_ports))


def quads_for_box(nx, ny):
    doc = json.loads(POLY.read_text())
    shift = np.array([nx / 2.0, ny / 2.0])
    return [np.asarray(w["vertices_a"], dtype=float) + shift for h in doc["horns"] for w in h["walls"]]


def metal_mask(points, tris, quads):
    cents = points[tris].mean(1)
    mask = np.zeros(len(tris), dtype=bool)
    for q in quads:
        mask |= MPath(np.asarray(q, dtype=float)).contains_points(cents, radius=1e-12)
    return mask


def stretches(cents, nx, ny, dp, sigma_max):
    x, y = cents[:, 0], cents[:, 1]
    sigx = np.zeros(len(cents))
    sigy = np.zeros(len(cents))
    m = x < dp
    sigx[m] = sigma_max * ((dp - x[m]) / dp) ** 2
    m = x > nx - dp
    sigx[m] = sigma_max * ((x[m] - (nx - dp)) / dp) ** 2
    m = y < dp
    sigy[m] = sigma_max * ((dp - y[m]) / dp) ** 2
    m = y > ny - dp
    sigy[m] = sigma_max * ((y[m] - (ny - dp)) / dp) ** 2
    omega = K0
    return 1.0 + 1j * sigx / omega, 1.0 + 1j * sigy / omega


def assemble(points, tris, k0, keep, mass_scale, nx, ny, dp, sigma_max, stiff_scale=None):
    """keep: elements to include. mass_scale multiplies k0^2. stiff_scale multiplies gradients."""
    n = len(points)
    use = np.flatnonzero(keep)
    tri = tris[use]
    ms = np.asarray(mass_scale, dtype=np.complex128)[use]
    if stiff_scale is None:
        ss = np.ones(len(use))
    else:
        ss = np.asarray(stiff_scale, dtype=np.complex128)[use]
    cents = points[tri].mean(1)
    sx, sy = stretches(cents, nx, ny, dp, sigma_max)
    x = points[tri, 0]
    y = points[tri, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    area = 0.5 * np.abs(twice)
    good = area > 1e-18
    bx = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1) / np.where(np.abs(twice) > 1e-18, twice, 1.0)[:, None]
    by = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1) / np.where(np.abs(twice) > 1e-18, twice, 1.0)[:, None]
    # Air ρ = I. Gradient weight is the PML anisotropic tensor.
    ax = sy / sx
    ay = sx / sy
    mass = (k0**2) * sx * sy * ms
    rows, cols, data = [], [], []
    for i in range(3):
        for j in range(3):
            term = ss * (ax * bx[:, j] * bx[:, i] + ay * by[:, j] * by[:, i])
            Ke = area * term
            Me = mass * area * (1.0 / 6.0 if i == j else 1.0 / 12.0)
            val = np.where(good, Ke - Me, 0.0)
            rows.append(tri[:, i])
            cols.append(tri[:, j])
            data.append(val)
    A = sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    row = np.abs(A).sum(axis=1).A1
    empty = row < 1e-14
    if np.any(empty):
        A = A.tolil()
        for i in np.flatnonzero(empty):
            A.rows[i] = [i]
            A.data[i] = [1.0 + 0j]
        A = A.tocsr()
    return A, empty


def classify(tris, metal):
    n = int(tris.max()) + 1
    n_m = np.zeros(n, dtype=np.int32)
    n_a = np.zeros(n, dtype=np.int32)
    for k in range(3):
        np.add.at(n_m, tris[:, k], metal.astype(np.int32))
        np.add.at(n_a, tris[:, k], (~metal).astype(np.int32))
    interior = (n_m > 0) & (n_a == 0)
    interface = (n_m > 0) & (n_a > 0)
    air = n_a > 0
    return interior, interface, air


def inject(points, center, tangent, offsets, amps, banned):
    b = np.zeros(len(points), dtype=np.complex128)
    for s, amp in zip(offsets, amps):
        xy = center + float(s) * tangent
        d2 = np.sum((points - xy) ** 2, axis=1)
        if banned is not None and np.any(banned):
            d2 = d2.copy()
            d2[banned] = np.inf
        b[int(np.argmin(d2))] += complex(amp)
    return b


def line_power(points, Hz, Ex, Ey, center, tangent, n_hat, span=None, n=31):
    if span is None:
        span = 0.96 * float(sc.clear_width)
    s = np.linspace(-span / 2, span / 2, n)
    ds = span / (n - 1)
    w = np.full(n, ds)
    w[0] *= 0.5
    w[-1] *= 0.5
    total = 0.0
    for si, wi in zip(s, w):
        xy = center + si * tangent
        i = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
        Sx = 0.5 * np.real(Ey[i] * np.conj(Hz[i]))
        Sy = -0.5 * np.real(Ex[i] * np.conj(Hz[i]))
        total += float(n_hat[0] * Sx + n_hat[1] * Sy) * float(wi)
    return total


def port_geometry(nx, ny):
    doc = json.loads(POLY.read_text())
    shift = np.array([nx / 2.0, ny / 2.0])
    ports = []
    for h in doc["horns"]:
        ports.append(
            {
                "center": np.asarray(h["monitor_center_a"], dtype=float) + shift,
                "source": np.asarray(h["source_center_a"], dtype=float) + shift,
                "n": np.asarray(h["outward"], dtype=float),
                "t": np.asarray(h["tangent"], dtype=float),
            }
        )
    return ports


def powers_from_solution(points, tris, x, rho_xx, rho_xy, rho_yx, rho_yy, nx, ny):
    Ex, Ey = fields_from_hz(points, tris, x, rho_xx, rho_xy, rho_yx, rho_yy, K0)
    out = []
    for p in port_geometry(nx, ny):
        out.append(line_power(points, x, Ex, Ey, p["center"], p["t"], p["n"]))
    return np.asarray(out, dtype=float)


def rho_for(tris, metal):
    T = len(tris)
    ones = np.ones(T, dtype=np.complex128)
    z = np.zeros(T, dtype=np.complex128)
    rxx, ryy = ones.copy(), ones.copy()
    rxx[metal] = 0.0
    ryy[metal] = 0.0
    return rxx, z.copy(), z.copy(), ryy


def solve_case(points, tris, metal, mode, nx, ny, dp, sigma_max, kind):
    """kind: rho0 | excised | nomass."""
    T = len(tris)
    if kind == "rho0":
        keep = np.ones(T, dtype=bool)
        mass_scale = np.ones(T)
        stiff = np.where(metal, 0.0, 1.0)
    elif kind == "excised":
        keep = ~metal
        mass_scale = np.ones(T)
        stiff = np.ones(T)
    elif kind == "nomass":
        keep = np.ones(T, dtype=bool)
        mass_scale = np.where(metal, 0.0, 1.0)
        stiff = np.where(metal, 0.0, 1.0)
    else:
        raise ValueError(kind)
    t0 = time.perf_counter()
    A, empty = assemble(points, tris, K0, keep, mass_scale, nx, ny, dp, sigma_max, stiff)
    b = inject(points, mode["center"], mode["tangent"], mode["offsets"], mode["amps"], empty)
    x, _lu, t_fac, t_sol, resid = solve_system(A, b, np.zeros(len(points), dtype=bool))
    # Incident straight feed on the same mesh / PML, excised slabs.
    feed_metal = feed_wall_mask(points, tris, nx, ny)
    Aref, empty_ref = assemble(points, tris, K0, ~feed_metal, np.ones(len(tris)), nx, ny, dp, sigma_max)
    bref = inject(points, mode["center"], mode["tangent"], mode["offsets"], mode["amps"], empty_ref)
    xref = solve_system(Aref, bref, np.zeros(len(points), dtype=bool))[0]
    raw = powers_from_solution(points, tris, x, *rho_for(tris, metal), nx, ny)
    pref = powers_from_solution(points, tris, xref, *rho_for(tris, feed_metal), nx, ny)
    pinc = abs(float(pref[0]))
    norm = raw / pinc
    return {
        "kind": kind,
        "sigma_max": float(sigma_max),
        "dp": float(dp),
        "nx": float(nx),
        "ny": float(ny),
        "n_nodes": int(len(points)),
        "factor_s": t_fac,
        "resid": resid,
        "P_inc": pinc,
        "raw": raw.tolist(),
        "normalized": norm.tolist(),
        "dB": [float(10 * np.log10(v)) if v > 0 else None for v in norm],
        "delta_dB_vs_meep": [
            None if (not np.isfinite(MEEP_NORM[i]) or norm[i] <= 0) else float(10 * np.log10(norm[i] / MEEP_NORM[i]))
            for i in range(6)
        ],
        "P2_minus_P6_dB": float(10 * np.log10(norm[1] / norm[5])) if norm[1] > 0 and norm[5] > 0 else None,
        "P3_minus_P5_dB": float(10 * np.log10(norm[2] / norm[4])) if norm[2] > 0 and norm[4] > 0 else None,
        "wall_s": time.perf_counter() - t0,
        "_x": x,
        "_interior_hz": None,
    }


def feed_wall_mask(points, tris, nx, ny):
    cents = points[tris].mean(1)
    shift = np.array([nx / 2.0, ny / 2.0])
    c = cents - shift
    half = 0.5 * float(sc.clear_width)
    th = float(sc.wall_thickness)
    # Slabs parallel to x, centered on y=0, long in x. P1 feed axis.
    return (np.abs(c[:, 1]) >= half) & (np.abs(c[:, 1]) <= half + th) & (np.abs(c[:, 0]) < 0.5 * nx)


def make_mesh(nx, ny, h, quads, jitter=0.0, seed=0):
    xs = np.arange(0.0, nx + 0.5 * h, h)
    ys = np.arange(0.0, ny + 0.5 * h, h)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    if jitter > 0:
        rng = np.random.default_rng(seed)
        # Keep the outer frame fixed so the PML box does not move.
        interior = (pts[:, 0] > h) & (pts[:, 0] < nx - h) & (pts[:, 1] > h) & (pts[:, 1] < ny - h)
        pts[interior] += rng.uniform(-jitter, jitter, size=(int(interior.sum()), 2))
    extra = np.vstack(quads)
    pts = np.vstack([pts, extra])
    _, keep = np.unique(np.round(pts, 8), axis=0, return_index=True)
    pts = pts[np.sort(keep)]
    tris = Delaunay(pts).simplices
    return pts, tris


def mode_for_box(src, nx, ny):
    shift = np.array([nx / 2.0, ny / 2.0]) - np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    return {
        "center": src["center_mesh"] + shift,
        "tangent": src["tangent"],
        "offsets": src["offsets"],
        "amps": src["amps"],
    }


def audit_matrix(points, tris, metal, nx, ny):
    interior, interface, air = classify(tris, metal)
    A_air, _ = assemble(points, tris, K0, ~metal, np.ones(len(tris)), nx, ny, sc.dpml_ports, 2.0)
    A_mass_only = mass_only(points, tris[metal], K0, len(points))
    d_air = np.array(A_air.diagonal())
    d_mass = np.array(A_mass_only.diagonal())
    rel = np.abs(d_mass[interface]) / np.maximum(np.abs(d_air[interface]), 1e-30)
    return {
        "n_interior_nodes": int(interior.sum()),
        "n_interface_nodes": int(interface.sum()),
        "n_air_nodes": int((air & ~interface).sum()),
        "interior_equation": "(-k0^2 M) Hz = 0, because stiffness is multiplied by rho=0 and the mass term is not",
        "mass_active_inside_pec": True,
        "interface_mass_perturbation_median": float(np.median(rel)) if rel.size else None,
        "interface_mass_perturbation_p90": float(np.percentile(rel, 90)) if rel.size else None,
        "note": "A nonzero metal mass on an interface node is absent from the excised-domain form.",
    }, interior, interface


def mass_only(points, tri, k0, n):
    if len(tri) == 0:
        return sparse.csr_matrix((n, n), dtype=np.complex128)
    x = points[tri, 0]
    y = points[tri, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    area = 0.5 * np.abs(twice)
    good = area > 1e-18
    rows, cols, data = [], [], []
    mass = k0**2
    for i in range(3):
        for j in range(3):
            Me = mass * area * (1.0 / 6.0 if i == j else 1.0 / 12.0)
            val = np.where(good, -Me, 0.0)
            rows.append(tri[:, i])
            cols.append(tri[:, j])
            data.append(val.astype(np.complex128))
    return sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()


def db_table(rec):
    return {k: rec[k] for k in ("kind", "sigma_max", "dp", "nx", "ny", "n_nodes", "factor_s", "P_inc", "normalized", "dB", "delta_dB_vs_meep", "P2_minus_P6_dB", "P3_minus_P5_dB", "wall_s")}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    src = load_numerical_mode()
    nx, ny = float(sc.nx_ports), float(sc.ny_ports)
    print("loading FEM-M", flush=True)
    points, tris = load_mesh("FEM-M")
    quads = quads_for_box(nx, ny)
    metal = metal_mask(points, tris, quads)
    info, interior, interface = audit_matrix(points, tris, metal, nx, ny)
    print("audit", json.dumps(info), flush=True)

    mode = mode_for_box(src, nx, ny)
    results = []
    for kind, sig in (
        ("rho0", 2.0),
        ("excised", 2.0),
        ("nomass", 2.0),
        ("rho0", MEEP_SIGMA_MAX),
        ("excised", MEEP_SIGMA_MAX),
    ):
        print(f"solve {kind} sigma={sig:.3f}", flush=True)
        rec = solve_case(points, tris, metal, mode, nx, ny, float(sc.dpml_ports), sig, kind)
        if kind == "rho0" and sig == 2.0:
            info["max_abs_Hz_interior"] = float(np.max(np.abs(rec["_x"][interior]))) if interior.any() else None
            info["max_abs_Hz_interface"] = float(np.max(np.abs(rec["_x"][interface]))) if interface.any() else None
            info["median_abs_Hz_air"] = float(np.median(np.abs(rec["_x"][~interior])))
        rec.pop("_x", None)
        rec.pop("_interior_hz", None)
        results.append(db_table(rec))
        print(json.dumps(results[-1]["delta_dB_vs_meep"]), flush=True)

    # Domain distance and PML thickness on a fresh mesh, excised + Meep-strength and baseline strength.
    h = 0.08
    domain_rows = []
    for extra, dp, sig, tag in (
        (0.0, 2.0, 2.0, "baseline_sigma2"),
        (0.125, 2.0, 2.0, "air_plus_25pct_of_0.5a"),
        (0.25, 2.0, 2.0, "air_plus_50pct_of_0.5a"),
        (2.0, 2.0, 2.0, "extra_2a_vacuum"),
        (4.0, 2.0, 2.0, "extra_4a_vacuum"),
        (2.0, 4.0, 2.0, "extra_2a_pml4_sigma2"),
        (2.0, 2.0, MEEP_SIGMA_MAX, "extra_2a_meep_sigma"),
        (2.0, 4.0, -np.log(1e-15) * 3.0 / (2.0 * 4.0), "extra_2a_pml4_meep_sigma"),
    ):
        # extra is added on each side, so the box grows by 2*extra. PML thickness dp sits inside the box.
        Lx, Ly = nx + 2 * extra, ny + 2 * extra
        if dp >= 0.5 * min(Lx, Ly):
            continue
        print("domain", tag, Lx, Ly, "dp", dp, flush=True)
        q = quads_for_box(Lx, Ly)
        pts, tri = make_mesh(Lx, Ly, h, q)
        met = metal_mask(pts, tri, q)
        md = mode_for_box(src, Lx, Ly)
        rec = solve_case(pts, tri, met, md, Lx, Ly, dp, sig, "excised")
        rec["tag"] = tag
        rec.pop("_x", None)
        rec.pop("_interior_hz", None)
        domain_rows.append(db_table(rec) | {"tag": tag})
        print(tag, rec["normalized"], rec["delta_dB_vs_meep"], flush=True)

    seed_rows = []
    for seed in range(4):
        print("seed", seed, flush=True)
        q = quads_for_box(nx, ny)
        pts, tri = make_mesh(nx, ny, h, q, jitter=0.35 * h, seed=seed)
        met = metal_mask(pts, tri, q)
        md = mode_for_box(src, nx, ny)
        rec = solve_case(pts, tri, met, md, nx, ny, 2.0, 2.0, "excised")
        seed_rows.append(
            {
                "seed": seed,
                "P2_minus_P6_dB": rec["P2_minus_P6_dB"],
                "P3_minus_P5_dB": rec["P3_minus_P5_dB"],
                "normalized": rec["normalized"],
                "delta_dB_vs_meep": rec["delta_dB_vs_meep"],
            }
        )
        print(seed_rows[-1], flush=True)

    payload = {
        "meep_sigma_max_at_d2": MEEP_SIGMA_MAX,
        "current_fem_sigma_max": 2.0,
        "audit": info,
        "fem_m_variants": results,
        "domain_sweep": domain_rows,
        "mesh_seeds": seed_rows,
    }
    (OUT / "pec_pml_fem.json").write_text(json.dumps(payload, indent=2) + "\n")
    print("wrote", OUT / "pec_pml_fem.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
