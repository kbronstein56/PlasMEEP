#!/usr/bin/env python3
"""Slow Stage-1…5 optimization evolution video from REAL histories."""
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

OUT = OUT_ROOT / "31_stage1_to_stage5_optimization_evolution"
DURATION_S = 16.0


def load_geometry():
    import sixport_common as sc
    from full91_geometry import bulb_centers_device
    from sixfold_orbits import build_orbits, expand_tied

    centers = np.asarray(bulb_centers_device(sc), float)
    orbits, origin = build_orbits(centers)
    return centers, orbits, origin, expand_tied


def accepted(h):
    out = accepted_stage2_records(h)
    if out:
        return out
    # Stage-1 production history often lacks kind tags — keep records with q + J_BB
    fallback = []
    for rec in h:
        if "q" in rec and "J_BB" in rec:
            fallback.append(rec)
    return fallback


def main():
    centers, orbits, _, expand_tied = load_geometry()
    # Stage-1 production history if present; else start at Stage-2
    timeline = []
    for tag, hist_name in (
        ("S1", "bb_production_history.json"),
        ("S2", "bb_stage2_history.json"),
        ("S3", "bb_stage3_history.json"),
        ("S4", "bb_stage4_history.json"),
        ("S5", "bb_stage5_history.json"),
    ):
        p = PROD_CAMPAIGN / hist_name
        if not p.exists():
            continue
        raw = load_json(p)
        hist = raw["history"] if isinstance(raw, dict) and "history" in raw else raw
        for r in accepted(hist):
            timeline.append((tag, r))

    if len(timeline) < 2:
        raise SystemExit(f"insufficient history points: {len(timeline)}")

    n = int(round(DURATION_S * FPS))
    clear_frames(OUT / "frames")
    paths = []

    # Precompute best-so-far J per stage for curve panels
    by_stage = {}
    for tag, r in timeline:
        by_stage.setdefault(tag, []).append(r)

    for i in range(n):
        u = np.clip((i - 15) / max(n - 30, 1), 0, 1)
        idx = int(round(u * (len(timeline) - 1)))
        stage, rec = timeline[idx]
        q = np.asarray(rec["q"], float)
        s = expand_tied(q, orbits)

        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(
            "Optimization evolution · Stages 1–5 (real checkpoints)",
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

        xy0 = centers - centers.mean(0)
        sca = ax0.scatter(
            xy0[:, 0],
            xy0[:, 1],
            c=s,
            s=70,
            cmap="viridis",
            vmin=0.4,
            vmax=1.45,
            edgecolors="k",
            linewidths=0.3,
        )
        ax0.set_aspect("equal")
        ax0.set_title(f"{stage}  ·  iter {rec.get('iteration')}  ·  J_BB={rec['J_BB']:.3f}", fontsize=16)
        fig.colorbar(sca, ax=ax0, fraction=0.046, pad=0.04, label="s")

        colors = {"S1": C_MUTED, "S2": "#64748b", "S3": "#94a3b8", "S4": C_WARN, "S5": C_ACCENT}
        for tag, recs in by_stage.items():
            # only show up to current global index within this stage's portion
            stage_upto = 0
            for j, (t, _) in enumerate(timeline[: idx + 1]):
                if t == tag:
                    stage_upto += 1
            if stage_upto == 0:
                continue
            xs, ys, best = [], [], None
            for k, r in enumerate(recs[:stage_upto]):
                jbb = float(r["J_BB"])
                best = jbb if best is None else max(best, jbb)
                xs.append(k)
                ys.append(best)
            ax1.plot(xs, ys, color=colors.get(tag, C_MUTED), lw=2.2, label=tag)
        ax1.set_ylabel("best J_BB")
        ax1.set_xlabel("accepted index within stage")
        ax1.legend(fontsize=11, ncol=5, loc="lower right")
        ax1.set_title("Accepted-objective progress by stage", fontsize=14)

        # worst A2 trail for current stage
        cur_recs = by_stage.get(stage, [])
        stage_upto = sum(1 for t, _ in timeline[: idx + 1] if t == stage)
        xs = list(range(stage_upto))
        wa = [float(r.get("worst_A2", np.nan)) for r in cur_recs[:stage_upto]]
        ax2.plot(xs, wa, color=colors.get(stage, C_ACCENT), lw=2.4)
        ax2.axhline(0, color="#94a3b8", lw=1)
        ax2.set_ylabel("worst A2 (obj grid)")
        ax2.set_xlabel(f"accepted index within {stage}")
        ax2.set_title(
            f"obj-grid worst A2={rec.get('worst_A2', float('nan')):.3f} @ "
            f"{rec.get('worst_f_GHz', float('nan'))} GHz",
            fontsize=13,
        )

        fp = OUT / "frames" / f"frame_{i:04d}.png"
        save_frame(fig, fp)
        paths.append(fp)

    stem = "31_stage1_to_stage5_optimization_evolution"
    mp4 = OUT / f"{stem}.mp4"
    encode_mp4(str(OUT / "frames" / "frame_*.png"), mp4, fps=FPS)
    write_poster(paths[len(paths) // 2], OUT / f"{stem}_poster.png")
    frames_to_gif(paths, OUT / f"{stem}_preview.gif", fps=6)
    (OUT / f"{stem}_PROVENANCE.txt").write_text(
        f"timeline_points={len(timeline)}\n"
        f"stages={sorted(by_stage)}\n"
        f"animation=device q + J_BB curves from real histories\n"
    )
    print("wrote", mp4, flush=True)


if __name__ == "__main__":
    main()
