#!/usr/bin/env python3
"""Plot |Hz| for source 0: uniform vs optimized design (compact PNG, not giant fields)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from circulator_core import B_PLUS, F_HZ, load_device, solve_six  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"


def plot_hz(dev, x, path, title):
    # sample on a coarse Cartesian grid for a small PNG
    pts = dev.points
    xmin, xmax = pts[:, 0].min(), pts[:, 0].max()
    ymin, ymax = pts[:, 1].min(), pts[:, 1].max()
    nx, ny = 240, 220
    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)
    XX, YY = np.meshgrid(xs, ys)
    xy = np.column_stack([XX.ravel(), YY.ravel()])
    # vacuum rho for display interpolation of Hz only
    rho = [
        np.ones(len(dev.tris), np.complex128),
        np.zeros(len(dev.tris), np.complex128),
        np.zeros(len(dev.tris), np.complex128),
        np.ones(len(dev.tris), np.complex128),
    ]
    hz, _, _ = dev.sampler.fields(x, *rho, 1.0, xy)
    amp = np.abs(hz).reshape(XX.shape)
    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=120)
    im = ax.imshow(
        amp,
        origin="lower",
        extent=[xmin, xmax, ymin, ymax],
        cmap="magma",
        aspect="equal",
    )
    ax.set_title(title)
    ax.set_xlabel("x / a")
    ax.set_ylabel("y / a")
    fig.colorbar(im, ax=ax, fraction=0.046, label="|Hz|")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print("WROTE", path, flush=True)


def main():
    opt = json.loads(Path(sys.argv[1]).read_text())
    grade = sys.argv[2] if len(sys.argv) > 2 else "C"
    dev = load_device(grade)
    P_ref = float(opt["P_ref"])
    s1 = np.ones(len(dev.rod_masks))
    sopt = np.asarray(opt["best"]["s"], float)
    b0 = solve_six(dev, s1, F_HZ, B_PLUS, P_ref=P_ref)
    b1 = solve_six(dev, sopt, F_HZ, B_PLUS, P_ref=P_ref)
    plot_hz(dev, b0.fields[0], OUT / f"field_uniform_{grade}_src0.png", f"Uniform s=1 |Hz| (src0, {grade})")
    plot_hz(dev, b1.fields[0], OUT / f"field_opt_{grade}_src0.png", f"Optimized |Hz| (src0, {grade})")
    # also save orbit q bar plot
    q = np.asarray(opt["best"]["q"], float)
    fig, ax = plt.subplots(figsize=(7, 3.2), dpi=120)
    ax.bar(np.arange(len(q)), q, color="#335577")
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("C6 orbit index")
    ax.set_ylabel("density scale q")
    ax.set_title("Optimized sixfold density scales")
    fig.tight_layout()
    fig.savefig(OUT / f"orbits_q_{Path(sys.argv[1]).stem}.png")
    plt.close(fig)
    print("WROTE orbits plot", flush=True)


if __name__ == "__main__":
    main()
