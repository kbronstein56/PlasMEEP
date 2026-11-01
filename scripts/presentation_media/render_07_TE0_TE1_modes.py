#!/usr/bin/env python3
"""07 — TE0 / TE1 mode explainer for a 40 mm parallel-plate / H-plane guide.

ANALYTIC transverse profiles only — not a device field export.
Animates Re[Hz exp(-iωt)] for the analytic modal shapes (correct oscillation),
never |H| as if it were oscillating.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from media_common import (  # noqa: E402
    C_ACCENT,
    C_GOOD,
    C_INK,
    C_MUTED,
    OUT_ROOT,
    clear_frames,
    encode_mp4,
    frames_to_gif,
    save_frame,
    stamp_kind,
    write_poster,
)

OUT = OUT_ROOT / "07_TE0_TE1_modes"
FRAMES = OUT / "frames"

W_MM = 40.0  # guide width


def main():
    clear_frames(FRAMES)
    y = np.linspace(0, W_MM, 400)
    # Parallel-plate / H-plane guide: TE_n ~ cos(n π y / a) for Hz (n=0,1,...)
    # TE0: uniform Hz; TE1: one half-wave.
    te0 = np.ones_like(y)
    te1 = np.cos(np.pi * y / W_MM)

    paths = []
    nframes = 120  # 4 s
    for i in range(nframes):
        phase = 2 * np.pi * i / nframes
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(
            "Waveguide modes in the 40 mm PlasMEEP port guides (ANALYTIC)",
            fontsize=26, fontweight="bold", color=C_INK, y=0.97,
        )
        ax0 = fig.add_axes([0.06, 0.42, 0.42, 0.45])
        ax1 = fig.add_axes([0.56, 0.42, 0.42, 0.45])
        axn = fig.add_axes([0.06, 0.08, 0.92, 0.25])
        stamp_kind(ax0, "ANALYTIC", "Re[Hz e^{-iωt}]")
        stamp_kind(ax1, "ANALYTIC", "Re[Hz e^{-iωt}]")
        stamp_kind(axn, "ANALYTIC")

        # 2D schematic strip: x=propagation, y=transverse
        x = np.linspace(0, 1, 80)
        X, Y = np.meshgrid(x, y)
        # complex amplitude with e^{-iωt} → Re = amp(y) cos(ωt) assuming real mode profile
        Z0 = np.outer(te0, np.ones_like(x)) * np.cos(phase)
        Z1 = np.outer(te1, np.ones_like(x)) * np.cos(phase)
        vmax = 1.15
        im0 = ax0.imshow(
            Z0, origin="lower", aspect="auto", cmap="RdBu_r",
            vmin=-vmax, vmax=vmax, extent=[0, 1, 0, W_MM],
        )
        im1 = ax1.imshow(
            Z1, origin="lower", aspect="auto", cmap="RdBu_r",
            vmin=-vmax, vmax=vmax, extent=[0, 1, 0, W_MM],
        )
        ax0.set_title("TE0 (fundamental) — uniform $H_z$", fontsize=18)
        ax1.set_title("TE1 — one transverse half-wave", fontsize=18)
        for ax in (ax0, ax1):
            ax.set_xlabel("propagation (schematic)")
            ax.set_ylabel("transverse y (mm)")
            ax.axhline(0, color="k", lw=2)
            ax.axhline(W_MM, color="k", lw=2)
        fig.colorbar(im0, ax=ax0, fraction=0.046, pad=0.03).set_label(r"Re $H_z$ (arb.)")
        fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.03).set_label(r"Re $H_z$ (arb.)")

        # line profiles
        axn.plot(y, te0 * np.cos(phase), color=C_GOOD, lw=3, label="TE0 Re[Hz]")
        axn.plot(y, te1 * np.cos(phase), color=C_ACCENT, lw=3, label="TE1 Re[Hz]")
        axn.plot(y, np.abs(te0), color=C_GOOD, ls=":", lw=2, label="|TE0| (static envelope)")
        axn.plot(y, np.abs(te1), color=C_ACCENT, ls=":", lw=2, label="|TE1| (static envelope)")
        axn.set_xlim(0, W_MM)
        axn.set_ylim(-1.3, 1.3)
        axn.set_xlabel("y (mm) across 40 mm guide")
        axn.set_ylabel("amplitude")
        axn.legend(loc="upper right", ncol=2, frameon=False, fontsize=12)
        axn.set_title(
            "Do not animate |H| as if it oscillates — envelopes are static; Re[Hz e^{-iωt}] oscillates",
            fontsize=14, color=C_MUTED,
        )
        axn.grid(True, alpha=0.25)

        p = FRAMES / f"frame_{i:04d}.png"
        save_frame(fig, p)
        paths.append(p)

    mp4 = OUT / "07_TE0_TE1_modes.mp4"
    encode_mp4(str(FRAMES / "frame_*.png"), mp4)
    write_poster(paths[0], OUT / "07_TE0_TE1_modes_poster.png")
    frames_to_gif(paths, OUT / "07_TE0_TE1_modes.gif", fps=12)
    print("wrote", mp4)


if __name__ == "__main__":
    main()
