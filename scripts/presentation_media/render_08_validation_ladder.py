#!/usr/bin/env python3
"""08 — Validation ladder movie from EXISTING overview figure / documented statuses.

No invented claims. Reuses the presentation ladder content from
make_project_overview_figures.py / VALIDATION_FIGURES_SUMMARY.md.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from media_common import (  # noqa: E402
    C_BAD,
    C_GOOD,
    C_INK,
    C_MUTED,
    OUT_ROOT,
    PROD_PRES_FIGS,
    clear_frames,
    encode_mp4,
    frames_to_gif,
    save_frame,
    stamp_kind,
    write_poster,
)

OUT = OUT_ROOT / "08_validation_ladder"
FRAMES = OUT / "frames"

C_PASS = "#2e7d32"
C_FAIL = "#c62828"
C_OPEN = "#9e9e9e"

# Exact rows from make_project_overview_figures.py (do not invent new statuses)
ROWS = [
    ("Homogeneous Faraday cell (Meep vs Lorentz)", "PASS", "|κ| err 0.024%"),
    ("FEM canonical B=0 benchmarks (Mie, slabs, guides)", "PASS", "≤0.01–0.05 dB"),
    ("FEM gyrotropic benchmarks (Voigt, T-matrix, Onsager)", "PASS", "≤2.6e-4 dB"),
    ("Meep↔FEM: horns only, B=0", "PASS", "≤0.14 dB @50 ppc, converging"),
    ("Meep↔FEM: quartz six-port, B=0", "PASS", "≤0.074 dB"),
    ("Meep↔FEM: 91 plasma bulbs, B=0", "FAIL", "diverges with ppc"),
    ("Meep↔FEM: full device, B≠0", "OPEN", "never run"),
    ("FEM full device B≠0 (mesh, Onsager, balance)", "PASS", "≤0.10 dB, 1e-13"),
    ("FEM 12×12 multimode DtN S (TE0+TE1)", "PASS", "Onsager 8e-14; bal ≤3.8%"),
    ("Adjoint gradient vs FD (multimode objective)", "PASS", "fwd≡adj 1e-13"),
    ("TE1-aware multimode optimization", "OPEN", "not started"),
]


def main():
    clear_frames(FRAMES)
    color = {"PASS": C_PASS, "FAIL": C_FAIL, "OPEN": C_OPEN}
    paths = []
    n = len(ROWS)
    # reveal rows one-by-one, then hold
    for reveal in range(1, n + 1):
        for _hold in range(12):
            fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
            ax = fig.add_axes([0.05, 0.08, 0.90, 0.82])
            fig.suptitle(
                "Project validation ladder (documented statuses — no new claims)",
                fontsize=26, fontweight="bold", color=C_INK, y=0.965,
            )
            stamp_kind(ax, "DATA", "from presentation_figures / overview script")
            ax.set_xlim(0, 1)
            ax.set_ylim(-0.2, n)
            ax.axis("off")
            for i, (label, status, note) in enumerate(ROWS[:reveal]):
                y = n - 1 - i
                ax.add_patch(
                    FancyBboxPatch(
                        (0.0, y + 0.12), 0.11, 0.76, boxstyle="round,pad=0.02",
                        fc=color[status], ec="none",
                    )
                )
                ax.text(0.055, y + 0.5, status, ha="center", va="center", color="white",
                        fontsize=14, fontweight="bold")
                ax.text(0.14, y + 0.5, label, ha="left", va="center", fontsize=16)
                ax.text(0.99, y + 0.5, note, ha="right", va="center", fontsize=14, color="#444444")
            ax.axhline(n - 3, color="#cccccc", lw=1)
            ax.axhline(n - 7, color="#cccccc", lw=1)
            ax.text(
                0.5, -0.05,
                "Unresolved: plasma B=0 Meep↔FEM divergence; no B≠0 full-device Meep↔FEM comparison saved.",
                ha="center", fontsize=14, color=C_MUTED, transform=ax.get_xaxis_transform(),
            )
            p = FRAMES / f"frame_{len(paths):04d}.png"
            save_frame(fig, p)
            paths.append(p)

    # append existing ladder PNG as final hold if present
    src = PROD_PRES_FIGS / "project_validation_ladder.png"
    if src.exists():
        shutil.copy(src, OUT / "project_validation_ladder_source.png")

    mp4 = OUT / "08_validation_ladder.mp4"
    encode_mp4(str(FRAMES / "frame_*.png"), mp4)
    write_poster(paths[-1], OUT / "08_validation_ladder_poster.png")
    frames_to_gif(paths, OUT / "08_validation_ladder.gif", fps=10)
    print("wrote", mp4)


if __name__ == "__main__":
    main()
