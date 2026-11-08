#!/usr/bin/env python3
"""Slow Stage-1/2/3 optimization evolution video from REAL histories."""
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
    C_WARN,
    FPS,
    OUT_ROOT,
    PROD_CAMPAIGN,
    VIZ_ROOT,
    accepted_stage2_records,
    clear_frames,
    encode_mp4,
    frames_to_gif,
    load_json,
    save_frame,
    stamp_kind,
    write_poster,
)

sys.path[:0] = [
    str(VIZ_ROOT / "scripts/circulator_inverse_design"),
    str(VIZ_ROOT / "scripts/validation/fem_full_device_validation"),
    str(VIZ_ROOT / "scripts/validation/fem_meep_validation"),
    str(VIZ_ROOT / "scripts/validation"),
    str(VIZ_ROOT / "scripts"),
]

OUT = OUT_ROOT / "18_stage3_optimization_evolution"
DURATION_S = 14.0


def load_geometry():
    import sixport_common as sc
    from full91_geometry import bulb_centers_device
    from sixfold_orbits import build_orbits, expand_tied

    centers = np.asarray(bulb_centers_device(sc), float)
    orbits, origin = build_orbits(centers)
    return centers, orbits, origin, expand_tied


def accepted(h):
    return accepted_stage2_records(h)


def main():
    centers, orbits, _, expand_tied = load_geometry()
    s2 = load_json(PROD_CAMPAIGN / "opt_broadband_stage2_A2_C.json")
    s3 = load_json(PROD_CAMPAIGN / "opt_broadband_stage3_A2_C.json")
    h2 = accepted(load_json(PROD_CAMPAIGN / "bb_stage2_history.json")["history"])
    h3 = accepted(load_json(PROD_CAMPAIGN / "bb_stage3_history.json")["history"])

    # Timeline: Stage-2 checkpoints then Stage-3
    timeline = [("S2", r) for r in h2] + [("S3", r) for r in h3]
    n = int(round(DURATION_S * FPS))
    clear_frames(OUT / "frames")
    paths = []

    def j_curve(stage_recs, upto_idx):
        xs, ys, best = [], [], None
        for r in stage_recs[: upto_idx + 1]:
            j = float(r["J_BB"])
            best = j if best is None else max(best, j)
            xs.append(int(r["iteration"]))
            ys.append(best)
        return xs, ys

    for i in range(n):
        # Hold briefly on first/last
        u = np.clip((i - 15) / max(n - 30, 1), 0, 1)
        idx = int(round(u * (len(timeline) - 1)))
        stage, rec = timeline[idx]
        q = np.asarray(rec["q"], float)
        s = expand_tied(q, orbits)

        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(
            "Optimization evolution · Stage-2 → Stage-3 (real checkpoints)",
            fontsize=24,
            fontweight="bold",
            color=C_INK,
            y=0.97,
        )
        ax0 = fig.add_axes([0.04, 0.12, 0.42, 0.78])
        ax1 = fig.add_axes([0.52, 0.55, 0.44, 0.35])
        ax2 = fig.add_axes([0.52, 0.12, 0.44, 0.35])
        stamp_kind(ax0, "DATA", stage)
        stamp_kind(ax1, "DATA")
        stamp_kind(ax2, "DATA")

        xy = centers - centers.mean(0)
        sca = ax0.scatter(
            xy[:, 0],
            xy[:, 1],
            c=s,
            s=60,
            cmap="viridis",
            vmin=0.4,
            vmax=1.45,
            edgecolors="k",
            linewidths=0.25,
        )
        ax0.set_aspect("equal")
        ax0.set_title(
            f"{'Stage-2' if stage == 'S2' else 'Stage-3'}  nit={rec['iteration']}  "
            f"J_BB={float(rec['J_BB']):.3f}",
            fontsize=16,
        )
        fig.colorbar(sca, ax=ax0, fraction=0.046, pad=0.04)

        # Best-so-far curves
        if stage == "S2":
            i2 = idx
            x2, y2 = j_curve(h2, i2)
            ax1.plot(x2, y2, color=C_MUTED, lw=2.5)
            ax1.set_title("Stage-2 best-so-far J_BB", fontsize=15)
        else:
            x2, y2 = j_curve(h2, len(h2) - 1)
            x3, y3 = j_curve(h3, idx - len(h2))
            ax1.plot(x2, y2, color=C_MUTED, lw=2.0, label="Stage-2")
            ax1.plot(x3, y3, color=C_ACCENT, lw=2.6, label="Stage-3")
            ax1.legend(fontsize=12)
            ax1.set_title("Stage-3 best-so-far J_BB", fontsize=15)
        ax1.set_ylabel("J_BB")

        # Final comparison bars (always)
        ax2.bar(
            [0, 1],
            [float(s2["J_BB_best"]), float(s3["J_BB_best"])],
            color=[C_MUTED, C_GOOD],
        )
        ax2.set_xticks([0, 1])
        ax2.set_xticklabels(["Stage-2 final", "Stage-3 final"])
        ax2.set_ylabel("J_BB")
        ax2.set_title("Final soft-min objectives", fontsize=15)

        fp = OUT / "frames" / f"frame_{i:04d}.png"
        save_frame(fig, fp)
        paths.append(fp)

    mp4 = OUT / "18_stage3_optimization_evolution.mp4"
    encode_mp4(str(OUT / "frames" / "frame_*.png"), mp4, fps=FPS)
    write_poster(paths[n // 2], OUT / "18_stage3_optimization_evolution_poster.png")
    frames_to_gif(paths, OUT / "18_stage3_optimization_evolution_preview.gif")
    print("wrote", mp4, flush=True)


if __name__ == "__main__":
    main()
