#!/usr/bin/env python3
"""Before/after formulation comparison figure from modal ablation JSON."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ABL = ROOT / "outputs" / "validation" / "modal_ports" / "modal_ablation_summary.json"
GRID = ROOT / "outputs" / "validation" / "modal_ports"
OUT = ROOT / "outputs" / "validation" / "modal_ports"


def main() -> None:
    data = json.loads(ABL.read_text(encoding="utf-8"))
    p12 = data["ablation_P1P2"]
    forms = list(p12.keys())
    labels = ["TE1+flux", "Num+flux", "TE1+modal", "Num+modal"]
    p12_vals = [p12["te1_hz_line"], p12["num_mode_hz_line"], p12["te1_hz_modal"], p12["num_mode_modal"]]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

    ax = axes[0]
    colors = ["#E45756", "#54A24B", "#F58518", "#4C78A8"]
    bars = ax.bar(labels, p12_vals, color=colors)
    ax.axhline(0.2, color="k", ls="--", lw=1, label="0.2 dB gate")
    ax.axhline(0.625, color="gray", ls=":", lw=1, label="TE1 baseline")
    ax.set_ylabel("P1↔P2 |P12−P21| (dB)")
    ax.set_title("Launch × receiver ablation (horns_only, res32, rt=5)")
    for b, v in zip(bars, p12_vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.3f}", ha="center", fontsize=9)
    ax.legend(fontsize=8)
    ax.set_ylim(0, max(p12_vals) * 1.25 + 0.1)

    ax = axes[1]
    grid = data["grid_offset_P1P2_num_mode_hz_line"]
    grid_offsets = list(grid.keys())
    grid_vals = list(grid.values())
    te1_range = (0.05, 0.79)
    ax.bar(grid_offsets, grid_vals, color="#54A24B", label="num_mode+flux")
    ax.axhspan(te1_range[0], te1_range[1], color="#E45756", alpha=0.2, label="TE1+flux prior range")
    ax.set_ylabel("P1↔P2 |P12−P21| (dB)")
    ax.set_title("Grid-offset robustness (num_mode_hz_line)")
    ax.legend(fontsize=8)

    fig.tight_layout()
    out = OUT / "modal_port_ablation.png"
    fig.savefig(out, dpi=160)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
