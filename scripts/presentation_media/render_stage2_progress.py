#!/usr/bin/env python3
"""Render Stage-2 A2 progress from REAL accepted checkpoints (read-only campaign JSON).

Used for:
  - 04_inverse_design_evolution.mp4
  - 09_stage2_progress_latest.mp4 / .png
"""
from __future__ import annotations

import argparse
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


def load_geometry():
    import sixport_common as sc
    from full91_geometry import bulb_centers_device
    from sixfold_orbits import build_orbits, expand_tied

    centers = np.asarray(bulb_centers_device(sc), float)
    orbits, origin = build_orbits(centers)
    return centers, orbits, origin, expand_tied


def panel_frame(rec, centers, orbits, expand_tied, j_curve, title_suffix=""):
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.32, wspace=0.28,
                          left=0.06, right=0.98, top=0.90, bottom=0.08)
    fig.suptitle(
        f"Stage-2 A2 broadband inverse design — REAL checkpoints{title_suffix}",
        fontsize=26, fontweight="bold", color=C_INK, y=0.97,
    )

    # Design map
    ax0 = fig.add_subplot(gs[:, 0])
    stamp_kind(ax0, "DATA", "Stage-2 q → C6 rods")
    q = np.asarray(rec["q"], float)
    s = expand_tied(q, orbits)
    xy = centers - centers.mean(axis=0)
    sca = ax0.scatter(
        xy[:, 0], xy[:, 1], c=s, s=55, cmap="viridis",
        vmin=0.4, vmax=1.45, edgecolors="k", linewidths=0.25,
    )
    ax0.set_aspect("equal")
    ax0.set_title(f"Design density  nit={rec['iteration']}  ({rec['kind']})", fontsize=16)
    ax0.set_xlabel("x (a)")
    ax0.set_ylabel("y (a)")
    cb = fig.colorbar(sca, ax=ax0, fraction=0.046, pad=0.04)
    cb.set_label("rod scale s")

    # J_BB history
    ax1 = fig.add_subplot(gs[0, 1:])
    stamp_kind(ax1, "DATA", "bb_stage2_history.json")
    nit = int(rec["iteration"])
    xs, ys = j_curve
    ax1.plot(xs, ys, "-o", color=C_ACCENT, ms=5, label="best J_BB so far")
    ax1.axvline(nit, color=C_WARN, ls="--", lw=1.5)
    ax1.scatter([nit], [float(rec["J_BB"])], s=90, color=C_BAD, zorder=5, label="this eval J_BB")
    ax1.set_xlabel("accepted L-BFGS iteration")
    ax1.set_ylabel(r"$J_{BB}$ (higher better)")
    ax1.set_title(f"J_BB={rec['J_BB']:.3f}   best_so_far tracked on accepted iters", fontsize=16)
    ax1.legend(loc="lower right", frameon=False)
    ax1.grid(True, alpha=0.25)

    # Spectrum of J_A2
    ax2 = fig.add_subplot(gs[1, 1])
    freqs = np.asarray(rec["freqs_GHz"], float)
    jpf = np.asarray(rec["J_per_freq"], float)
    ax2.plot(freqs, jpf, "-o", color=C_ACCENT, ms=4)
    wf = float(rec["worst_f_GHz"])
    ax2.axvline(wf, color=C_BAD, ls="--", lw=1.6)
    ax2.axvline(4.1, color=C_GOOD, ls=":", lw=1.6, label="4.1 GHz")
    ax2.set_xlabel("f (GHz)")
    ax2.set_ylabel(r"$J_{A2}(f)$")
    ax2.set_title(f"Per-freq A2  worst@{wf:.3f} GHz", fontsize=15)
    ax2.legend(loc="lower left", frameon=False, fontsize=12)
    ax2.grid(True, alpha=0.25)

    # Metrics callouts
    ax3 = fig.add_subplot(gs[1, 2])
    ax3.axis("off")
    stamp_kind(ax3, "DATA")
    w = rec.get("worst", {})
    a41 = rec.get("at_4p1", {})
    lines = [
        f"elapsed: {rec.get('elapsed_s', 0)/3600:.2f} h",
        f"worst A2: {rec.get('worst_A2', float('nan')):.3f}",
        f"worst f:  {rec.get('worst_f_GHz', float('nan')):.3f} GHz",
        "",
        "At 4.1 GHz:",
        f"  T_des = {a41.get('T_des', float('nan')):.4f}",
        f"  IL    = {a41.get('IL_dB', float('nan')):.2f} dB",
        f"  iso   = {a41.get('isolation_dB', float('nan')):.2f} dB",
        "",
        "Worst-freq metrics:",
        f"  T_des = {w.get('T_des', float('nan')):.4f}",
        f"  IL    = {w.get('IL_dB', float('nan')):.2f} dB",
        f"  TE1   = {w.get('TE1', float('nan')):.3f}",
    ]
    ax3.text(0.02, 0.98, "\n".join(lines), transform=ax3.transAxes, va="top",
             fontsize=15, family="monospace", color=C_INK,
             bbox=dict(boxstyle="round,pad=0.5", fc="#f0f4ff", ec=C_ACCENT))
    return fig


def build_j_curve(accepted):
    xs, ys, best = [], [], None
    for rec in accepted:
        j = float(rec["J_BB"])
        best = j if best is None else max(best, j)
        xs.append(int(rec["iteration"]))
        ys.append(best)
    return np.asarray(xs), np.asarray(ys)


def render_evolution(hold: int = 4):
    out = OUT_ROOT / "04_inverse_design_evolution"
    frames = out / "frames"
    clear_frames(frames)
    hist = load_json(PROD_CAMPAIGN / "bb_stage2_history.json")["history"]
    accepted = accepted_stage2_records(hist)
    centers, orbits, _, expand_tied = load_geometry()
    j_curve = build_j_curve(accepted)
    paths = []
    fi = 0
    # subsample if many: keep all 31 (start+30)
    for rec in accepted:
        fig = panel_frame(rec, centers, orbits, expand_tied, j_curve, title_suffix="")
        for _ in range(hold):
            p = frames / f"frame_{fi:04d}.png"
            # re-save same content hold times without redraw: copy via save once then duplicate
            if _ == 0:
                save_frame(fig, p)
                first = p
            else:
                import shutil
                shutil.copy(first, p)
            paths.append(p)
            fi += 1
        plt.close("all")
    mp4 = out / "04_inverse_design_evolution.mp4"
    encode_mp4(str(frames / "frame_*.png"), mp4)
    write_poster(paths[-1], out / "04_inverse_design_evolution_poster.png")
    frames_to_gif(paths, out / "04_inverse_design_evolution.gif", fps=8)
    print("wrote", mp4, "n_checkpoints", len(accepted))


def render_latest():
    out = OUT_ROOT / "09_stage2_progress_latest"
    frames = out / "frames"
    clear_frames(frames)
    blob = load_json(PROD_CAMPAIGN / "opt_broadband_stage2_A2_C.json")
    hist = load_json(PROD_CAMPAIGN / "bb_stage2_history.json")["history"]
    accepted = accepted_stage2_records(hist)
    # Prefer last accepted_iter; fall back to rec with q_best
    rec = accepted[-1]
    # Ensure q is best if available
    if "q_best" in blob:
        rec = dict(rec)
        rec["q"] = blob["q_best"]
        rec["J_BB"] = blob.get("J_BB_best", rec["J_BB"])
        # use rec_best spectrum if present
        rb = blob.get("rec_best", {})
        if "J_per_freq" in rb and "freqs_GHz" in rb:
            rec["J_per_freq"] = rb["J_per_freq"]
            rec["freqs_GHz"] = rb["freqs_GHz"]
            rec["worst_f_GHz"] = rb.get("worst_f_GHz", rec.get("worst_f_GHz"))
            rec["worst_A2"] = rb.get("worst_J", rec.get("worst_A2"))
            if "worst_metrics" in rb:
                rec["worst"] = rb["worst_metrics"]
    centers, orbits, _, expand_tied = load_geometry()
    j_curve = build_j_curve(accepted)
    fig = panel_frame(
        rec, centers, orbits, expand_tied, j_curve,
        title_suffix=f"  |  FINAL best J_BB={blob.get('J_BB_best', float('nan')):.3f}  (30/30)",
    )
    paths = []
    for i in range(60):  # 2 s hold
        p = frames / f"frame_{i:04d}.png"
        if i == 0:
            save_frame(fig, p)
            first = p
        else:
            import shutil
            shutil.copy(first, p)
        paths.append(p)
    mp4 = out / "09_stage2_progress_latest.mp4"
    encode_mp4(str(frames / "frame_*.png"), mp4)
    write_poster(paths[0], out / "09_stage2_progress_latest.png")
    write_poster(paths[0], out / "09_stage2_progress_latest_poster.png")
    frames_to_gif(paths[:1] * 8, out / "09_stage2_progress_latest.gif", fps=4)
    print("wrote", mp4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["evolution", "latest", "both"], default="both")
    ap.add_argument("--hold", type=int, default=4, help="frames to hold each checkpoint")
    args = ap.parse_args()
    if args.mode in ("evolution", "both"):
        render_evolution(hold=args.hold)
    if args.mode in ("latest", "both"):
        render_latest()


if __name__ == "__main__":
    main()
