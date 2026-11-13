#!/usr/bin/env python3
"""Slow Stage-3 vs Stage-4 broadband spectral sweep video + GIF + poster."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from media_common import (  # noqa: E402
    C_ACCENT,
    C_BAD,
    C_GOOD,
    C_INK,
    C_MUTED,
    C_WARN,
    FPS,
    OUT_ROOT,
    PROD_CAMPAIGN,
    clear_frames,
    encode_mp4,
    frames_to_gif,
    load_json,
    save_frame,
    stamp_kind,
    write_poster,
)

OUT = OUT_ROOT / "19_stage3_vs_stage4_broadband_sweep"
DURATION_S = 14.0


def xy(rows, key):
    f = np.array([r["f_GHz"] for r in rows], float)
    y = np.array([r.get(key) if r.get(key) is not None else np.nan for r in rows], float)
    return f, y


def main():
    s3 = load_json(PROD_CAMPAIGN / "bb_stage3_dense_25MHz.json")
    s4 = load_json(PROD_CAMPAIGN / "bb_stage4_dense_25MHz.json")
    f3, a3 = xy(s3["rows"], "J_A2")
    f4, a4 = xy(s4["rows"], "J_A2")
    _, iso3 = xy(s3["rows"], "isolation_dB")
    _, iso4 = xy(s4["rows"], "isolation_dB")

    n = int(round(DURATION_S * FPS))
    clear_frames(OUT / "frames")
    paths = []
    fmin, fmax = 4.0, 6.0
    for i in range(n):
        u = min(1.0, i / (n - 30)) if n > 30 else i / max(n - 1, 1)
        f_cur = fmin + (fmax - fmin) * u
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(
            "Stage-3 → Stage-4 broadband spectral comparison",
            fontsize=26,
            fontweight="bold",
            color=C_INK,
            y=0.97,
        )
        ax0 = fig.add_axes([0.07, 0.55, 0.88, 0.35])
        ax1 = fig.add_axes([0.07, 0.10, 0.88, 0.35])
        for ax in (ax0, ax1):
            stamp_kind(ax, "DATA", "DENSE AUDIT")
        ax0.plot(f3, a3, color=C_MUTED, lw=2.2, label="Stage-3")
        ax0.plot(f4, a4, color=C_ACCENT, lw=2.6, label="Stage-4")
        for f, c in (
            (5.08125, C_BAD),
            (4.725, C_WARN),
            (5.8125, C_GOOD),
            (4.775, C_GOOD),
        ):
            ax0.axvline(f, color=c, ls="--", lw=1.4, alpha=0.8)
        ax0.axvline(f_cur, color="#0f766e", lw=2.5)
        ax0.axhline(0, color="#94a3b8", lw=1)
        ax0.set_xlim(fmin, fmax)
        ax0.set_ylabel("A2")
        ax0.legend(loc="lower right", fontsize=14, ncol=2)

        ax1.plot(f3, iso3, color=C_MUTED, lw=2.2, label="Stage-3")
        ax1.plot(f4, iso4, color=C_ACCENT, lw=2.6, label="Stage-4")
        ax1.axvline(f_cur, color="#0f766e", lw=2.5)
        ax1.axhline(0, color="#94a3b8", lw=1)
        ax1.set_xlim(fmin, fmax)
        ax1.set_xlabel("Frequency (GHz)")
        ax1.set_ylabel("Isolation (dB)")
        ax1.legend(loc="lower right", fontsize=14)
        fig.text(
            0.5,
            0.02,
            f"cursor @ {f_cur:.3f} GHz  ·  markers: 5.08125/4.725 new holes · 5.8125/4.775 targets",
            ha="center",
            fontsize=13,
            color=C_MUTED,
        )
        fp = OUT / "frames" / f"frame_{i:04d}.png"
        save_frame(fig, fp)
        paths.append(fp)

    stem = "19_stage3_vs_stage4_broadband_sweep"
    mp4 = OUT / f"{stem}.mp4"
    encode_mp4(str(OUT / "frames" / "frame_*.png"), mp4, fps=FPS)
    write_poster(paths[len(paths) // 2], OUT / f"{stem}_poster.png")
    frames_to_gif(paths, OUT / f"{stem}_preview.gif", fps=6)
    print("wrote", mp4, flush=True)


if __name__ == "__main__":
    main()
