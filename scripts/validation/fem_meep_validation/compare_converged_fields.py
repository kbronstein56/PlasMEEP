#!/usr/bin/env python3
"""Sample the converged polygon FEM-H+ field onto the Meep Hz grid."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from b0_campaign import build_symmetric_mesh, horn_rho  # noqa: E402
from fem_validated_solver import (  # noqa: E402
    ElementSampler,
    assemble_anisotropic,
    inject_mode_consistent,
    load_numerical_mode,
    solve_system,
)

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "b0_campaign"
GRID = ROOT / "outputs" / "validation" / "fem_meep_validation" / "boundary_pml" / "meep_hz_grid.npz"


def main() -> int:
    points, tris = build_symmetric_mesh(0.05, 0.016, 0.030, 0.035)
    k0 = 2 * np.pi * float(sc.fs_a)
    pec = np.zeros(len(points), dtype=bool)
    rho, _mask = horn_rho(points, tris)
    sampler = ElementSampler(points, tris)
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    src = load_numerical_mode()
    b = inject_mode_consistent(points, tris, src, sampler)
    x, *_rest = solve_system(A, b, pec)
    hz = np.load(GRID)
    stride = 8
    xs, ys, MH = hz["x"][::stride], hz["y"][::stride], hz["Hz"][::stride, ::stride]
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    xy = np.column_stack([xx.ravel(), yy.ravel()])
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    fem_s, _ex, _ey = sampler.fields(x, *rho, k0, xy + shift)
    fem = fem_s.reshape(xx.shape)
    from b0_campaign import horn_prism_quads

    inside = np.zeros(xx.shape, dtype=bool)
    for q in horn_prism_quads(sc.nx_ports, sc.ny_ports):
        inside |= MPath(q - shift).contains_points(xy).reshape(xx.shape)
    sel = (np.abs(xx - 11.0) < 0.6) & (np.abs(yy) < 0.6) & ~inside
    den = np.vdot(fem[sel], fem[sel])
    fem = fem * (np.vdot(fem[sel], MH[sel]) / den)
    air = ~inside
    err = np.abs(fem - MH)
    rel = err / np.maximum(np.abs(MH), 1e-6)
    r = np.hypot(xx, yy)
    bins = []
    for r0, r1 in ((0, 3), (3, 6), (6, 9), (9, 11), (11, 13)):
        m = (r >= r0) & (r < r1) & air
        bins.append(
            {
                "r": [r0, r1],
                "median_rel": float(np.median(rel[m])),
                "p90_rel": float(np.percentile(rel[m], 90)),
            }
        )
    print("RADIAL", json.dumps(bins), flush=True)
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.2))
    extent = [xs[0], xs[-1], ys[0], ys[-1]]
    show_f = np.where(air, np.abs(fem), np.nan)
    show_m = np.where(air, np.abs(MH), np.nan)
    show_e = np.where(air, err, np.nan)
    for a, data, title in (
        (ax[0], show_m, "|Hz| Meep"),
        (ax[1], show_f, "|Hz| FEM-H+"),
        (ax[2], show_e, "|difference|"),
    ):
        im = a.imshow(data.T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
        a.set_title(title)
        fig.colorbar(im, ax=a, fraction=0.046)
    fig.tight_layout()
    fig.savefig(OUT / "field_diff_FEM-H+.png", dpi=120)
    (OUT / "field_diff_FEM-H+.json").write_text(json.dumps({"radial": bins, "n_nodes": int(len(points))}, indent=2) + "\n")
    print("wrote field_diff_FEM-H+.png", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
