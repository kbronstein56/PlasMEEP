#!/usr/bin/env python3
"""Forward solves for the production 91-bulb device.

One factorization per mesh, frequency, and B. All six ports reuse that factorization.
The complex port matrix used for Onsager is the reaction matrix b_i^T A^{-1} b_j.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic, guide_normal_flux  # noqa: E402
from full91_geometry import OUT, production_constants  # noqa: E402
from gyrotropic_tensor import dissipation_matrix, matrix_of, ordinary_from_si, rho_xy, tensor_ordinary  # noqa: E402

MESH = OUT / "meshes"


def free_gib() -> float:
    avail = 0.0
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            avail = float(line.split()[1]) / (1024.0**2)
    return avail


def load_mesh(grade: str):
    data = np.load(MESH / f"full91_{grade}.npz")
    return data["points"], data["triangles"]


def material_masks(points, tris, sc_mod):
    from matplotlib.path import Path as MPath

    cents = points[tris].mean(1)
    r_p = 4.6 * sc_mod.r_bulb_inner / 6.5
    r_in, r_out = sc_mod.r_bulb_inner, sc_mod.r_bulb_outer
    from full91_geometry import bulb_centers_device

    centers = bulb_centers_device(sc_mod)
    plasma = np.zeros(len(tris), dtype=bool)
    quartz = np.zeros(len(tris), dtype=bool)
    for c in centers:
        d2 = (cents[:, 0] - c[0]) ** 2 + (cents[:, 1] - c[1]) ** 2
        plasma |= d2 <= r_p**2
        quartz |= (d2 <= r_out**2) & (d2 > r_in**2)
    shift = np.array([sc_mod.nx_ports / 2.0, sc_mod.ny_ports / 2.0])
    walls = np.zeros(len(tris), dtype=bool)
    for horn in sc_mod.full_horns:
        for name in ("left_flare", "right_flare", "left_feed", "right_feed"):
            poly = np.asarray(horn[name], float)[:, :2] + shift
            walls |= MPath(poly).contains_points(cents, radius=1e-12)
    return plasma, quartz, walls, centers


def rho_of(points, tris, f_hz: float, b_tesla: float, masks):
    plasma, quartz, walls, _ = masks
    f, fp, gamma, fc = ordinary_from_si(f_hz, sc.fp_Hz, sc.gamma_Hz, b_tesla, sc.a)
    exx, exy, eyx, eyy, _, _ = tensor_ordinary(f, fp, gamma, fc)
    rxx, rxy, ryx, ryy, _ = rho_xy(exx, exy, eyx, eyy)
    T = len(tris)
    rho = [
        np.ones(T, np.complex128),
        np.zeros(T, np.complex128),
        np.zeros(T, np.complex128),
        np.ones(T, np.complex128),
    ]
    rho[0][quartz] = 1.0 / 3.8
    rho[3][quartz] = 1.0 / 3.8
    rho[0][plasma] = rxx
    rho[1][plasma] = rxy
    rho[2][plasma] = ryx
    rho[3][plasma] = ryy
    for arr in rho:
        arr[walls] = 0.0
    eps = matrix_of(exx, exy, eyx, eyy)
    return rho, eps, fc, f


def port_lines(n_points: int = 41):
    """Uniform fundamental-mode samples across each feed. Hz PEC mode is flat."""
    lines = []
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    span = 0.92 * sc.clear_width
    offsets = np.linspace(-span / 2.0, span / 2.0, n_points)
    ds = span / (n_points - 1)
    weights = np.full(n_points, ds)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    for port in range(6):
        horn = sc.full_horns[port]
        n_hat = np.asarray(sc.port_dirs[port], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        src = np.asarray(horn["source_center"], float) + shift
        mon = np.asarray(horn["monitor_center"], float) + shift
        lines.append({
            "port": port,
            "outward": n_hat,
            "tangent": tang,
            "source_xy": np.vstack([src + float(s) * tang for s in offsets]),
            "monitor_xy": np.vstack([mon + float(s) * tang for s in offsets]),
            "weights": weights,
            "offsets": offsets,
        })
    return lines


def inject(points, tris, sampler, xy, weights):
    b = np.zeros(len(points), np.complex128)
    tri_id = sampler.locate(xy)
    for n in range(len(xy)):
        t = int(tri_id[n])
        if t < 0:
            i = int(np.argmin(np.sum((points - xy[n]) ** 2, axis=1)))
            b[i] += weights[n]
            continue
        w = sampler._bary(t, xy[n])
        nodes = tris[t]
        for k in range(3):
            b[int(nodes[k])] += float(w[k]) * float(weights[n])
    return b


def monitor_amplitude(sampler, uh, rho, k0, xy, weights):
    hz, _, _ = sampler.fields(uh, *rho, k0, xy)
    return complex(np.dot(weights, hz))


def volume_absorption(points, tris, uh, rho, k0, eps, plasma):
    diss = dissipation_matrix(eps)
    twice = (points[tris[:, 1], 0] - points[tris[:, 0], 0]) * (points[tris[:, 2], 1] - points[tris[:, 0], 1]) - (
        points[tris[:, 2], 0] - points[tris[:, 0], 0]
    ) * (points[tris[:, 1], 1] - points[tris[:, 0], 1])
    areas = 0.5 * np.abs(twice)
    bx = np.stack(
        [points[tris[:, 1], 1] - points[tris[:, 2], 1], points[tris[:, 2], 1] - points[tris[:, 0], 1], points[tris[:, 0], 1] - points[tris[:, 1], 1]],
        1,
    ) / twice[:, None]
    by = np.stack(
        [points[tris[:, 2], 0] - points[tris[:, 1], 0], points[tris[:, 0], 0] - points[tris[:, 2], 0], points[tris[:, 1], 0] - points[tris[:, 0], 0]],
        1,
    ) / twice[:, None]
    gx = np.sum(bx * uh[tris], 1)
    gy = np.sum(by * uh[tris], 1)
    ex = (1j / k0) * (rho[0] * gy + rho[1] * (-gx))
    ey = (1j / k0) * (rho[2] * gy + rho[3] * (-gx))
    # Dissipation uses the plasma tensor only on plasma elements.
    ev = np.stack([ex[plasma], ey[plasma]], 1)
    absorb = np.real(np.einsum("ti,ij,tj->t", ev.conj(), diss, ev))
    return float(np.sum(areas[plasma] * (k0 / 2.0) * absorb))


def rectangle_flux(sampler, uh, rho, k0):
    nx, ny, dp = float(sc.nx_ports), float(sc.ny_ports), float(sc.dpml_ports)
    n = 400
    xs = np.linspace(dp, nx - dp, n)
    ys = np.linspace(dp, ny - dp, n)
    segments = [
        (np.column_stack([xs, np.full(n, dp)]), np.array([0.0, -1.0]), xs),
        (np.column_stack([xs, np.full(n, ny - dp)]), np.array([0.0, 1.0]), xs),
        (np.column_stack([np.full(n, dp), ys]), np.array([-1.0, 0.0]), ys),
        (np.column_stack([np.full(n, nx - dp), ys]), np.array([1.0, 0.0]), ys),
    ]
    total = 0.0
    for xy, n_hat, coord in segments:
        hz, ex, ey = sampler.fields(uh, *rho, k0, xy)
        sx = 0.5 * np.real(ey * np.conj(hz))
        sy = -0.5 * np.real(ex * np.conj(hz))
        total += float(np.trapezoid((sx * n_hat[0] + sy * n_hat[1]), coord))
    return total


def solve_case(grade: str, b_tesla: float, f_hz: float, compare_scalar: bool = False) -> dict:
    points, tris = load_mesh(grade)
    sc.fs_a = f_hz * sc.a / 299792458.0
    k0 = 2 * np.pi * sc.fs_a
    masks = material_masks(points, tris, sc)
    rho, eps, fc, f_a = rho_of(points, tris, f_hz, b_tesla, masks)
    plasma, quartz, walls, _ = masks
    print(
        f"case grade={grade} B={b_tesla} f={f_hz/1e9:.4f} nodes={len(points)} "
        f"plasma={int(plasma.sum())} quartz={int(quartz.sum())} walls={int(walls.sum())} free_GiB={free_gib():.1f}",
        flush=True,
    )
    t0 = time.perf_counter()
    A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    t_asm = time.perf_counter() - t0
    scalar_resid = None
    if compare_scalar:
        rho_s, _, _, _ = rho_of(points, tris, f_hz, 0.0, masks)
        As = assemble_anisotropic(points, tris, *rho_s, k0, np.zeros(len(points), dtype=bool))
        diff = (A - As).tocsr()
        scalar_resid = float(np.linalg.norm(diff.data) / (np.linalg.norm(A.data) + 1e-30))
        del As, diff, rho_s
        print("scalar_operator_resid", scalar_resid, flush=True)
    nnz = int(A.nnz)
    # Complex CSR ~ 24 bytes/entry plus index overhead. LU fill is estimated from the
    # prior 5.94e6-DOF factorization at about 45 GiB, i.e. ~8e-6 GiB/DOF.
    est_lu = 8.0e-6 * len(points)
    print(f"nnz={nnz} matrix_GiB~{nnz*24/1024**3:.2f} est_LU_GiB~{est_lu:.1f} free_GiB={free_gib():.1f}", flush=True)
    if est_lu > 0.55 * free_gib():
        raise RuntimeError(f"estimated LU {est_lu:.1f} GiB is too close to {free_gib():.1f} GiB free")
    t1 = time.perf_counter()
    lu = splu(A.tocsc())
    t_fac = time.perf_counter() - t1
    print(f"factor_s={t_fac:.1f} free_GiB={free_gib():.1f}", flush=True)
    sampler = ElementSampler(points, tris)
    lines = port_lines()
    rhs = [inject(points, tris, sampler, line["source_xy"], line["weights"]) for line in lines]
    fields = []
    t_solve = 0.0
    for b in rhs:
        t2 = time.perf_counter()
        fields.append(lu.solve(b))
        t_solve += time.perf_counter() - t2
    # Reaction matrix. Column j is the response to port j.
    S = np.zeros((6, 6), np.complex128)
    monitor = np.zeros((6, 6), np.complex128)
    power = np.zeros((6, 6), float)
    absorb = []
    rect = []
    for j, uh in enumerate(fields):
        for i, b in enumerate(rhs):
            S[i, j] = np.dot(b, uh)
        for i, line in enumerate(lines):
            monitor[i, j] = monitor_amplitude(sampler, uh, rho, k0, line["monitor_xy"], line["weights"])
            # guide_normal_flux expects port index and uses sixport monitor geometry.
            power[i, j] = guide_normal_flux(
                points, uh, None, None, i, sampler=sampler, tris=tris, rho=rho, omega=k0
            )
        absorb.append(volume_absorption(points, tris, uh, rho, k0, eps, plasma))
        rect.append(rectangle_flux(sampler, uh, rho, k0))
        print(f"port {j} Pabs={absorb[-1]:.6e} Prect={rect[-1]:.6e} Pports={power[:, j].tolist()}", flush=True)
    def pack(z):
        return [[complex(v).real, complex(v).imag] for v in np.asarray(z).ravel()]

    out = {
        "grade": grade,
        "B_T": b_tesla,
        "f_Hz": f_hz,
        "fs_a": sc.fs_a,
        "fc_a": fc,
        "nodes": int(len(points)),
        "triangles": int(len(tris)),
        "nnz": nnz,
        "plasma_triangles": int(plasma.sum()),
        "quartz_triangles": int(quartz.sum()),
        "wall_triangles": int(walls.sum()),
        "scalar_operator_resid": scalar_resid,
        "factor_s": t_fac,
        "assemble_s": t_asm,
        "solve6_s": t_solve,
        "est_lu_gib": est_lu,
        "reaction_reim": pack(S),
        "monitor_reim": pack(monitor),
        "port_power": power.tolist(),
        "absorption": absorb,
        "pml_rectangle_outward": rect,
        "eps_xx": [complex(eps[0, 0]).real, complex(eps[0, 0]).imag],
        "eps_xy": [complex(eps[0, 1]).real, complex(eps[0, 1]).imag],
    }
    tag = f"full91_{grade}_B{b_tesla:+.4f}_f{f_hz/1e9:.4f}".replace("+", "p").replace("-", "m")
    (OUT / f"{tag}.json").write_text(json.dumps(out) + "\n")
    print("WROTE", tag, flush=True)
    return out


def main():
    grade = sys.argv[1]
    b_tesla = float(sys.argv[2])
    f_hz = float(sys.argv[3]) if len(sys.argv) > 3 else sc.fs_Hz
    compare = "--scalar" in sys.argv
    solve_case(grade, b_tesla, f_hz, compare_scalar=compare)


if __name__ == "__main__":
    main()
