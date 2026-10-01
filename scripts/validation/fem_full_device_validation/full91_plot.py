#!/usr/bin/env python3
"""Downsampled |Hz| plots for one source port. Does not keep the nodal field."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from full91_solve import OUT, inject, load_mesh, material_masks, port_lines, rho_of  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from scipy.sparse.linalg import splu
import sixport_common as sc  # noqa: E402


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grade = "C"
    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    lines = port_lines()
    sampler = ElementSampler(points, tris)
    b = inject(points, tris, sampler, lines[0]["source_xy"], lines[0]["weights"])
    nx, ny = int(sc.nx_ports), int(sc.ny_ports)
    xs = np.linspace(2.0, nx - 2.0, 180)
    ys = np.linspace(2.0, ny - 2.0, 168)
    xx, yy = np.meshgrid(xs, ys)
    samples = np.column_stack([xx.ravel(), yy.ravel()])
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), sharey=True)
    for ax, b_tesla, title in zip(axes, (0.0, 0.05, -0.05), ("B = 0", "B = +0.05 T", "B = -0.05 T")):
        sc.fs_a = sc.fs_Hz * sc.a / 299792458.0
        k0 = 2 * np.pi * sc.fs_a
        rho, _, _, _ = rho_of(points, tris, sc.fs_Hz, b_tesla, masks)
        A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
        uh = splu(A.tocsc()).solve(b)
        hz, _, _ = sampler.fields(uh, *rho, k0, samples)
        img = np.log10(np.maximum(np.abs(hz).reshape(xx.shape), 1e-8))
        im = ax.imshow(
            img, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="viridis", vmin=-3.5, vmax=-1.0
        )
        ax.set_title(title)
        ax.set_aspect("equal")
        print("plot", title, flush=True)
    fig.colorbar(im, ax=axes, fraction=0.02, label="log10 |Hz|")
    fig.savefig(OUT / "full91_hz_port0.png", dpi=140)
    print("WROTE plot", flush=True)


if __name__ == "__main__":
    main()
