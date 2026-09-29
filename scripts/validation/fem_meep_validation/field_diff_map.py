#!/usr/bin/env python3
"""Field-difference maps and one finer six-horn check with a Meep-matched PML."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, Rectangle
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from boundary_diagnosis import (  # noqa: E402
    MEEP_NORM,
    MEEP_SIGMA_MAX,
    make_mesh,
    metal_mask,
    mode_for_box,
    quads_for_box,
    solve_case,
)
from fem_validated_solver import load_mesh, load_numerical_mode  # noqa: E402
import sixport_common as sc  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "boundary_pml"
POLY = ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_localization" / "meep_horn_polygons.json"


def sample_fem(points, hz, xy):
    tree = cKDTree(points)
    _, idx = tree.query(xy, k=1)
    return hz[idx]


def align(ref, other):
    den = np.vdot(other, other)
    return other * (np.vdot(other, ref) / den)


def plot_pair(tag, meep_x, meep_y, meep_hz, fem_hz, quads):
    # meep_hz, fem_hz on the same flattened non-metal samples. Rebuild a grid for imshow
    # from the saved npz coordinates passed in as 1d axes and 2d arrays.
    fig, ax = plt.subplots(2, 3, figsize=(12.5, 8.2))
    panels = {
        (0, 0): (np.abs(meep_hz), "|Hz| Meep"),
        (0, 1): (np.abs(fem_hz), "|Hz| FEM"),
        (0, 2): (np.abs(fem_hz - meep_hz), "|Hz FEM − Hz Meep|"),
        (1, 0): (np.degrees(np.angle(fem_hz * np.conj(meep_hz))), "phase FEM−Meep (deg)"),
        (1, 1): (np.abs(fem_hz - meep_hz) / np.maximum(np.abs(meep_hz), 1e-6), "relative |error|"),
    }
    extent = [meep_x[0], meep_x[-1], meep_y[0], meep_y[-1]]
    for (i, j), (data, title) in panels.items():
        im = ax[i, j].imshow(data.T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
        ax[i, j].set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax[i, j], fraction=0.046)
    ax[1, 2].set_aspect("equal")
    ax[1, 2].set_title("geometry")
    for a in ax.ravel():
        for q in quads:
            a.add_patch(Polygon(q, closed=True, fill=False, ec="w", lw=0.4))
        a.add_patch(Rectangle((-13, -12), 26, 24, fill=False, ec="r", lw=0.8, ls="--"))
        a.plot([10.055, 10.055], [-0.96, 0.96], color="cyan", lw=1.0)
        a.plot([11.255, 11.255], [-0.96, 0.96], color="orange", lw=1.0)
        a.set_xlim(extent[0], extent[1])
        a.set_ylim(extent[2], extent[3])
    ax[1, 2].set_xlabel("red dashed = PML inner boundary\ncyan monitor, orange source")
    fig.suptitle(tag)
    fig.tight_layout()
    fig.savefig(OUT / f"field_diff_{tag}.png", dpi=120)
    plt.close(fig)


def grid_compare(points, x, meep, stride=8):
    hz = np.load(OUT / "meep_hz_grid.npz")
    mx, my, mHz = hz["x"], hz["y"], hz["Hz"]
    # Confirm orientation: shape (1500, 1400) matches (len x, len y).
    xs = mx[::stride]
    ys = my[::stride]
    MH = mHz[::stride, ::stride]
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    xy = np.column_stack([xx.ravel(), yy.ravel()])
    # FEM mesh is corner-origin. Meep grid is centered.
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    fem = sample_fem(points, x, xy + shift).reshape(xx.shape)
    doc = json.loads(POLY.read_text())
    quads = [np.asarray(w["vertices_a"], dtype=float) for h in doc["horns"] for w in h["walls"]]
    inside = np.zeros(xx.shape, dtype=bool)
    from matplotlib.path import Path as MPath

    for q in quads:
        inside |= MPath(q).contains_points(xy).reshape(xx.shape)
    # Scale on the source-monitor segment, outside metal.
    sel = (np.abs(xx - 10.5) < 0.8) & (np.abs(yy) < 0.7) & ~inside
    if sel.sum() < 10:
        sel = ~inside
    scale_src = fem[sel]
    scale_ref = MH[sel]
    den = np.vdot(scale_src, scale_src)
    fem = fem * (np.vdot(scale_src, scale_ref) / den)
    fem = np.where(inside, np.nan, fem)
    MH = np.where(inside, np.nan, MH)
    err = fem - MH
    rel = np.abs(err) / np.maximum(np.abs(MH), 1e-6)
    # Where is the error large, in radial bins from the origin.
    r = np.hypot(xx, yy)
    bins = []
    for r0, r1 in ((0, 3), (3, 6), (6, 9), (9, 11), (11, 13)):
        m = (r >= r0) & (r < r1) & np.isfinite(rel)
        bins.append({"r": [r0, r1], "median_rel": float(np.median(rel[m])) if m.any() else None, "p90_rel": float(np.percentile(rel[m], 90)) if m.any() else None})
    return xs, ys, MH, fem, quads, bins


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    src = load_numerical_mode()
    nx, ny = float(sc.nx_ports), float(sc.ny_ports)
    print("FEM-M rho0 current PML", flush=True)
    points, tris = load_mesh("FEM-M")
    metal = metal_mask(points, tris, quads_for_box(nx, ny))
    mode = mode_for_box(src, nx, ny)
    rec = solve_case(points, tris, metal, mode, nx, ny, 2.0, 2.0, "rho0")
    xs, ys, MH, fem, quads, bins = grid_compare(points, rec["_x"], None)
    plot_pair("rho0_sigma2", xs, ys, MH, fem, quads)
    print("bins rho0", json.dumps(bins), flush=True)

    print("finer mesh, +2a vacuum, Meep sigma, excised", flush=True)
    Lx, Ly = nx + 4.0, ny + 4.0
    q = quads_for_box(Lx, Ly)
    pts, tri = make_mesh(Lx, Ly, 0.055, q)
    met = metal_mask(pts, tri, q)
    md = mode_for_box(src, Lx, Ly)
    rec2 = solve_case(pts, tri, met, md, Lx, Ly, 2.0, MEEP_SIGMA_MAX, "excised")
    # Map this solution into the centered frame. points are in the enlarged box.
    # Meep xy + shift_large = mesh coordinate.
    hz = np.load(OUT / "meep_hz_grid.npz")
    stride = 8
    xs = hz["x"][::stride]
    ys = hz["y"][::stride]
    MH = hz["Hz"][::stride, ::stride]
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    xy = np.column_stack([xx.ravel(), yy.ravel()])
    shift = np.array([Lx / 2.0, Ly / 2.0])
    fem2 = sample_fem(pts, rec2["_x"], xy + shift).reshape(xx.shape)
    from matplotlib.path import Path as MPath

    inside = np.zeros(xx.shape, dtype=bool)
    for quad in quads:
        inside |= MPath(quad).contains_points(xy).reshape(xx.shape)
    sel = (np.abs(xx - 10.5) < 0.8) & (np.abs(yy) < 0.7) & ~inside
    den = np.vdot(fem2[sel], fem2[sel])
    fem2 = fem2 * (np.vdot(fem2[sel], MH[sel]) / den)
    fem2 = np.where(inside, np.nan, fem2)
    MHn = np.where(inside, np.nan, MH)
    plot_pair("excised_meep_sigma_box_plus_2a", xs, ys, MHn, fem2, quads)
    rel = np.abs(fem2 - MHn) / np.maximum(np.abs(MHn), 1e-6)
    r = np.hypot(xx, yy)
    bins2 = []
    for r0, r1 in ((0, 3), (3, 6), (6, 9), (9, 11), (11, 13)):
        m = (r >= r0) & (r < r1) & np.isfinite(rel)
        bins2.append({"r": [r0, r1], "median_rel": float(np.median(rel[m])) if m.any() else None})
    print("bins fixed", json.dumps(bins2), flush=True)
    out = {
        "h": 0.055,
        "box": [Lx, Ly],
        "sigma_max": MEEP_SIGMA_MAX,
        "kind": "excised",
        "normalized": rec2["normalized"],
        "dB": rec2["dB"],
        "delta_dB_vs_meep": rec2["delta_dB_vs_meep"],
        "P2_minus_P6_dB": rec2["P2_minus_P6_dB"],
        "P3_minus_P5_dB": rec2["P3_minus_P5_dB"],
        "n_nodes": rec2["n_nodes"],
        "factor_s": rec2["factor_s"],
        "meep": MEEP_NORM.tolist(),
        "radial_rel_current_rho0": bins,
        "radial_rel_excised_meep_pml": bins2,
    }
    (OUT / "confirmed_sixhorn.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out["delta_dB_vs_meep"]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
