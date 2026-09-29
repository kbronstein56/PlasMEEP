#!/usr/bin/env python3
"""Homogeneous-air phase constant at fs. FEM structured mesh and Meep res 50.

A point Hz source radiates a cylindrical wave. The far-field phase slope is k_num.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_continuum"


def fit_k(s, z):
    amp = np.abs(z)
    use = amp > 0.2 * np.median(amp[amp > 0])
    ph = np.unwrap(np.angle(z))
    # Drop the first and last 10% so the fit is away from the source and the PML.
    n = len(s)
    core = np.zeros(n, dtype=bool)
    core[int(0.25 * n) : int(0.85 * n)] = True
    use = use & core & np.isfinite(ph)
    coef = np.polyfit(s[use], ph[use], 1)
    pred = np.polyval(coef, s[use])
    rms = float(np.sqrt(np.mean((ph[use] - pred) ** 2)))
    return float(coef[0]), float(coef[1]), int(use.sum()), rms


def fem_lines():
    from scipy.spatial import Delaunay

    import sixport_common as sc
    from fem_validated_solver import ElementSampler, assemble_anisotropic, solve_system

    nx, ny, h = 28.0, 20.0, 0.04
    sc.nx_ports = nx
    sc.ny_ports = ny
    sc.dpml_ports = 2.0
    xs = np.arange(0.0, nx + 1e-9, h)
    ys = np.arange(0.0, ny + 1e-9, h)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    tris = Delaunay(pts).simplices
    k0 = 2 * np.pi * float(sc.fs_a)
    T = len(tris)
    rho = (
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    )
    pec = np.zeros(len(pts), dtype=bool)
    A = assemble_anisotropic(pts, tris, *rho, k0, pec)
    # Meep source at (-8, 0) in a cell [-14, 14] x [-10, 10] maps to mesh (6, 10).
    src_xy = np.array([6.0, 10.0])
    sampler = ElementSampler(pts, tris)
    b = np.zeros(len(pts), dtype=np.complex128)
    t = int(sampler.locate(src_xy[None, :])[0])
    w = sampler._bary(t, src_xy)
    for k in range(3):
        b[int(tris[t, k])] += float(w[k])
    x, *_ = solve_system(A, b, pec)

    def grab(xy):
        hz, _, _ = sampler.fields(x, *rho, k0, xy)
        return hz

    # Axial: mesh y=10, x from 10 to 24 (4a to 18a from the source at x=6).
    s = np.linspace(4.0, 16.0, 400)
    axial = np.column_stack([6.0 + s, np.full_like(s, 10.0)])
    # 45-degree ray from the source.
    diag = np.column_stack([6.0 + s / np.sqrt(2.0), 10.0 + s / np.sqrt(2.0)])
    # Keep the diagonal inside the air, short of the PML (PML starts at x=26, y=18).
    keep = (diag[:, 0] < 25.0) & (diag[:, 1] < 17.0)
    return {
        "k0": float(k0),
        "axial": (s, grab(axial)),
        "diagonal": (s[keep], grab(diag[keep])),
    }


def meep_lines():
    import meep as mp
    import sixport_common as sc

    fs = float(sc.fs_a)
    cell = mp.Vector3(28, 20)
    sim = mp.Simulation(
        cell_size=cell,
        resolution=50,
        boundary_layers=[mp.PML(2.0)],
        sources=[
            mp.Source(
                mp.GaussianSource(fs, fwidth=0.20 * fs),
                component=mp.Hz,
                center=mp.Vector3(-8, 0),
            )
        ],
        default_material=mp.Medium(epsilon=1),
        eps_averaging=False,
    )
    # Box covering the axial segment and the 45-degree ray. Coordinates are cell-centered.
    vol = mp.Volume(center=mp.Vector3(3.0, 3.0), size=mp.Vector3(18.0, 14.0))
    dft = sim.add_dft_fields([mp.Hz], fs, 0, 1, where=vol)
    sim.run(until_after_sources=20)
    hz = np.array(sim.get_dft_array(dft, mp.Hz, 0))
    xs, ys, *_ = sim.get_array_metadata(dft_cell=dft)
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    # hz layout follows Meep's (x, y) or (y, x). Match by shape.
    if hz.shape == (len(xs), len(ys)):
        field = hz
    elif hz.shape == (len(ys), len(xs)):
        field = hz.T
    else:
        raise RuntimeError(f"unexpected DFT shape {hz.shape} vs x {len(xs)} y {len(ys)}")

    def interp(pts):
        ix = np.interp(pts[:, 0], xs, np.arange(len(xs)))
        iy = np.interp(pts[:, 1], ys, np.arange(len(ys)))
        i0 = np.clip(np.floor(ix).astype(int), 0, len(xs) - 2)
        j0 = np.clip(np.floor(iy).astype(int), 0, len(ys) - 2)
        tx = np.clip(ix - i0, 0, 1)
        ty = np.clip(iy - j0, 0, 1)
        return (
            (1 - tx) * (1 - ty) * field[i0, j0]
            + tx * (1 - ty) * field[i0 + 1, j0]
            + (1 - tx) * ty * field[i0, j0 + 1]
            + tx * ty * field[i0 + 1, j0 + 1]
        )

    s = np.linspace(4.0, 16.0, 400)
    axial_xy = np.column_stack([-8.0 + s, np.zeros_like(s)])
    diag_xy = np.column_stack([-8.0 + s / np.sqrt(2.0), 0.0 + s / np.sqrt(2.0)])
    keep = (diag_xy[:, 0] < 11.0) & (diag_xy[:, 1] < 7.0) & (diag_xy[:, 0] > -13) & (diag_xy[:, 1] > -9)
    return {
        "axial": (s, interp(axial_xy)),
        "diagonal": (s[keep], interp(diag_xy[keep])),
        "x_range": [float(xs[0]), float(xs[-1])],
        "y_range": [float(ys[0]), float(ys[-1])],
    }


def pack(tag, s, z, k0):
    k, phi0, n, rms = fit_k(s, z)
    return {
        "k_num": k,
        "k0": k0,
        "k_over_k0": k / k0,
        "relative_error": k / k0 - 1.0,
        "phi0": phi0,
        "n_fit": n,
        "phase_fit_rms_rad": rms,
        "note": tag,
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "fem"
    if mode == "fem":
        data = fem_lines()
        k0 = data["k0"]
        out = {"solver": "FEM", "h_a": 0.04, "k0": k0}
        for name in ("axial", "diagonal"):
            s, z = data[name]
            out[name] = pack(name, s, z, k0)
            np.savez_compressed(OUT / f"freespace_fem_{name}.npz", s=s, Hz=z)
        (OUT / "freespace_fem.json").write_text(json.dumps(out, indent=2) + "\n")
        print(json.dumps(out, indent=2), flush=True)
        return 0
    if mode == "meep":
        import sixport_common as sc

        k0 = 2 * np.pi * float(sc.fs_a)
        data = meep_lines()
        out = {"solver": "Meep", "resolution": 50, "points_per_cm": 25.0, "k0": k0, "grid": data}
        for name in ("axial", "diagonal"):
            s, z = data[name]
            out[name] = pack(name, s, z, k0)
            np.savez_compressed(OUT / f"freespace_meep_{name}.npz", s=s, Hz=z)
        # grid ranges are lists already; drop non-json from a nested copy
        out["grid"] = {"x_range": data["x_range"], "y_range": data["y_range"]}
        (OUT / "freespace_meep.json").write_text(json.dumps(out, indent=2) + "\n")
        print(json.dumps(out, indent=2), flush=True)
        return 0
    raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    raise SystemExit(main())
