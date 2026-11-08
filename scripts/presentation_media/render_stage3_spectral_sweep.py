#!/usr/bin/env python3
"""Slow Stage-2 vs Stage-3 broadband spectral sweep video + GIF + poster."""
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

OUT = OUT_ROOT / "13_stage2_vs_stage3_broadband_sweep"
DURATION_S = 14.0


def xy(rows, key):
    f = np.array([r["f_GHz"] for r in rows], float)
    y = np.array([r.get(key) if r.get(key) is not None else np.nan for r in rows], float)
    return f, y


def main():
    s2 = load_json(PROD_CAMPAIGN / "bb_stage2_audit_50MHz.json")
    dense = PROD_CAMPAIGN / "bb_stage3_dense_25MHz.json"
    s3 = load_json(dense if dense.exists() else PROD_CAMPAIGN / "bb_stage3_audit_50MHz.json")
    f2, a2 = xy(s2["rows"], "J_A2")
    f3, a3 = xy(s3["rows"], "J_A2")
    _, iso2 = xy(s2["rows"], "isolation_dB")
    _, iso3 = xy(s3["rows"], "isolation_dB")

    n = int(round(DURATION_S * FPS))
    clear_frames(OUT / "frames")
    paths = []
    fmin, fmax = 4.0, 6.0
    for i in range(n):
        # Cursor sweeps slowly left→right then brief hold at end
        u = min(1.0, i / (n - 30)) if n > 30 else i / max(n - 1, 1)
        f_cur = fmin + (fmax - fmin) * u
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(
            "Stage-2 → Stage-3 broadband spectral comparison",
            fontsize=26,
            fontweight="bold",
            color=C_INK,
            y=0.97,
        )
        ax0 = fig.add_axes([0.07, 0.55, 0.88, 0.35])
        ax1 = fig.add_axes([0.07, 0.10, 0.88, 0.35])
        for ax in (ax0, ax1):
            stamp_kind(ax, "DATA", "AUDIT")
        ax0.plot(f2, a2, color=C_MUTED, lw=2.2, label="Stage-2")
        ax0.plot(f3, a3, color=C_ACCENT, lw=2.6, label="Stage-3")
        for f, c, lab in (
            (5.75, C_BAD, "5.75 hole"),
            (5.825, C_BAD, "5.825 hole"),
            (4.3, C_WARN, "4.3 residual"),
            (5.775, "#7c3aed", "5.775 residual"),
        ):
            ax0.axvline(f, color=c, ls="--", lw=1.4, alpha=0.8)
        ax0.axvline(f_cur, color=C_GOOD, lw=2.5)
        ax0.axhline(0, color="#94a3b8", lw=1)
        ax0.set_xlim(fmin, fmax)
        ax0.set_ylabel("A2")
        ax0.legend(loc="lower right", fontsize=13, ncol=2)

        ax1.plot(f2, iso2, color=C_MUTED, lw=2.2, label="Stage-2")
        ax1.plot(f3, iso3, color=C_ACCENT, lw=2.6, label="Stage-3")
        ax1.axvline(f_cur, color=C_GOOD, lw=2.5)
        ax1.axhline(0, color="#94a3b8", lw=1)
        ax1.set_xlim(fmin, fmax)
        ax1.set_ylabel("Isolation (dB)")
        ax1.set_xlabel("Frequency (GHz)")
        ax1.text(
            0.99,
            0.95,
            f"cursor {f_cur:.3f} GHz",
            transform=ax1.transAxes,
            ha="right",
            va="top",
            fontsize=14,
            color=C_GOOD,
        )
        fp = OUT / "frames" / f"frame_{i:04d}.png"
        save_frame(fig, fp)
        paths.append(fp)

    mp4 = OUT / "13_stage2_vs_stage3_broadband_sweep.mp4"
    encode_mp4(str(OUT / "frames" / "frame_*.png"), mp4, fps=FPS)
    write_poster(paths[n // 2], OUT / "13_stage2_vs_stage3_broadband_sweep_poster.png")
    frames_to_gif(paths, OUT / "13_stage2_vs_stage3_broadband_sweep_preview.gif")
    print("wrote", mp4, f"({DURATION_S}s)", flush=True)


if __name__ == "__main__":
    main()
