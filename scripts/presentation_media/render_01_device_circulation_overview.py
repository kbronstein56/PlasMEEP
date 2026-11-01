#!/usr/bin/env python3
"""01 — SCHEMATIC device / circulation overview (no EM field data)."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Wedge

sys.path.insert(0, str(Path(__file__).resolve().parent))
from media_common import (  # noqa: E402
    C_ACCENT,
    C_BG,
    C_GOOD,
    C_INK,
    C_MUTED,
    C_SCHEMATIC,
    OUT_ROOT,
    clear_frames,
    encode_mp4,
    frames_to_gif,
    save_frame,
    stamp_kind,
    write_poster,
)

OUT = OUT_ROOT / "01_device_circulation_overview"
FRAMES = OUT / "frames"


def port_xy(i: int, r: float = 3.4):
    ang = i * 60.0 * math.pi / 180.0
    return r * math.cos(ang), r * math.sin(ang)


def draw_base(ax, phase: float):
    ax.set_xlim(-5.2, 5.2)
    ax.set_ylim(-4.4, 4.8)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor(C_BG)
    stamp_kind(ax, "SCHEMATIC", "not a field solve")

    # cavity
    ax.add_patch(Circle((0, 0), 2.55, fc="#e8eef8", ec=C_INK, lw=2.2, zorder=1))
    # plasma bulb ring (schematic)
    for k in range(12):
        ang = (k + phase * 0.15) * 2 * math.pi / 12
        ax.add_patch(
            Circle(
                (1.55 * math.cos(ang), 1.55 * math.sin(ang)),
                0.18,
                fc="#9ec5ff",
                ec="#2458a8",
                lw=1.0,
                zorder=2,
            )
        )
    ax.add_patch(Circle((0, 0), 0.22, fc="#2458a8", ec="none", zorder=3))
    ax.text(0, 0.55, "91 plasma rods\n(C6-tied design)", ha="center", va="center",
            fontsize=14, color=C_MUTED, zorder=4)

    # B field badge
    ax.add_patch(
        FancyBboxPatch(
            (-4.9, 3.7), 2.6, 0.85, boxstyle="round,pad=0.08",
            fc="#fff3d6", ec=C_SCHEMATIC, lw=1.5, zorder=5,
        )
    )
    ax.text(-3.6, 4.12, r"$B_z = +B$  (B_PLUS)", ha="center", va="center",
            fontsize=16, fontweight="bold", color=C_SCHEMATIC, zorder=6)
    ax.annotate(
        "",
        xy=(-3.6, 3.55),
        xytext=(-3.6, 3.9),
        arrowprops=dict(arrowstyle="-|>", color=C_SCHEMATIC, lw=2),
        zorder=6,
    )
    ax.text(-3.6, 3.35, "out of page", ha="center", fontsize=12, color=C_SCHEMATIC)

    # ports
    for i in range(6):
        x, y = port_xy(i, 3.55)
        hx, hy = port_xy(i, 2.7)
        ax.plot([hx, x], [hy, y], color="#334155", lw=10, solid_capstyle="butt", zorder=3)
        ax.plot([hx, x], [hy, y], color="#cbd5e1", lw=6, solid_capstyle="butt", zorder=4)
        ax.text(
            1.18 * x, 1.18 * y, f"P{i+1}",
            ha="center", va="center", fontsize=20, fontweight="bold", color=C_INK, zorder=7,
        )
        # 40 mm guide callout on P1
        if i == 0:
            ax.text(x + 0.15, y - 0.55, "40 mm guide", fontsize=13, color=C_MUTED, ha="left")

    ax.text(
        0, -4.15,
        "Convention (code): for source j, desired = (j+1) mod 6, reverse = (j−1) mod 6  ·  B_PLUS",
        ha="center", fontsize=14, color=C_MUTED,
    )


def draw_circulation(ax, t: float):
    # animate circulating arrows P1→P2→… for B+
    for i in range(6):
        a0 = (i + t) * 60 * math.pi / 180
        a1 = (i + 0.55 + t) * 60 * math.pi / 180
        r = 2.95
        x0, y0 = r * math.cos(a0), r * math.sin(a0)
        x1, y1 = r * math.cos(a1), r * math.sin(a1)
        ax.add_patch(
            FancyArrowPatch(
                (x0, y0), (x1, y1),
                arrowstyle="-|>", mutation_scale=18,
                color=C_GOOD, lw=2.4, zorder=8,
            )
        )
    ax.text(
        0, -3.55,
        r"Circulation sense for $B_z=+B$:  P1→P2→P3→P4→P5→P6→P1",
        ha="center", fontsize=18, fontweight="bold", color=C_GOOD,
    )


def draw_roles(ax, t: float):
    # highlight one source example: P1 → desired P2, reverse P6
    src = 0
    des = 1
    rev = 5
    for i, (lab, col) in zip([src, des, rev], [("SRC", C_ACCENT), ("DES", C_GOOD), ("REV", "#c2410c")]):
        x, y = port_xy(i, 3.55)
        ax.add_patch(Circle((x, y), 0.42, fc=col, alpha=0.25 + 0.1 * math.sin(2 * math.pi * t), ec=col, lw=2))
        ax.text(x, y + 0.72, lab, ha="center", fontsize=14, fontweight="bold", color=col)
    ax.annotate(
        "",
        xy=port_xy(des, 3.1),
        xytext=port_xy(src, 3.1),
        arrowprops=dict(arrowstyle="-|>", color=C_GOOD, lw=3, connectionstyle="arc3,rad=0.25"),
    )
    ax.annotate(
        "",
        xy=port_xy(rev, 3.1),
        xytext=port_xy(src, 3.1),
        arrowprops=dict(arrowstyle="-|>", color="#c2410c", lw=2.2, connectionstyle="arc3,rad=-0.25", linestyle="--"),
    )
    ax.text(
        0, -3.55,
        "Example excitation at P1: reward desired P2, penalize reverse P6 (mode-preserving TE0 FoM)",
        ha="center", fontsize=16, color=C_INK,
    )


def main():
    clear_frames(FRAMES)
    paths = []
    n = 90  # 3 s @ 30 fps
    for i in range(n):
        t = i / n
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.90])
        fig.suptitle(
            "PlasMEEP six-port plasma circulator — geometry & circulation (SCHEMATIC)",
            fontsize=28, fontweight="bold", color=C_INK, y=0.97,
        )
        draw_base(ax, phase=t)
        if t < 0.45:
            draw_circulation(ax, t * 2)
        else:
            draw_roles(ax, (t - 0.45) / 0.55)
        p = FRAMES / f"frame_{i:04d}.png"
        save_frame(fig, p)
        paths.append(p)

    mp4 = OUT / "01_device_circulation_overview.mp4"
    encode_mp4(str(FRAMES / "frame_*.png"), mp4)
    write_poster(paths[len(paths) // 2], OUT / "01_device_circulation_overview_poster.png")
    frames_to_gif(paths, OUT / "01_device_circulation_overview.gif", fps=12)
    print("wrote", mp4)


if __name__ == "__main__":
    main()
